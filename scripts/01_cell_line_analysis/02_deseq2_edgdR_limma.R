# =========================================================
# Within-cell-line cfRNA vs cell_debris DE analysis
# Methods: limma-voom, DESeq2, edgeR
# Purpose:
#   For each cell line, compare cfRNA vs cell_debris
#   Output both up-regulated (cfRNA higher) and down-regulated genes
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
  library(edgeR)
  library(DESeq2)
})

# =========================================================
# 1. Path settings
# =========================================================

count_file <- "./expression_matrix/all_mlRNA_counts.txt"
sample_info_file <- "./information_table/cellcult_Sample_Information.txt"

outdir <- "./01_mlRNA_cfRNA_vs_cell_debris"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Parameter settings
# =========================================================

cell_lines <- c(
  "HEK293T",
  "HTR-8/SVneo",
  "Hep3B2.1-7",
  "HepG2",
  "K562"
)

logFC_cutoff <- 1
padj_cutoff <- 0.05

# =========================================================
# Safe filename handling function
# =========================================================

safe_name <- function(x) {
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[:*?\"<>|]", "_", x)
  return(x)
}

# =========================================================
# 3. Read expression matrix and sample information
# =========================================================

counts <- read.table(
  count_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

counts <- as.matrix(counts)

info <- read_tsv(sample_info_file, show_col_types = FALSE)

cat("Expression matrix dimensions:\n")
print(dim(counts))

# =========================================================
# 4. Sample information processing
# =========================================================

info2 <- info %>%
  filter(QC %in% c(1, 2)) %>%
  mutate(
    Component2 = case_when(
      Component %in% c("cell", "debris", "cell_debris") ~ "cell_debris",
      Component %in% c("cfRNA", "medium", "supernatant") ~ "cfRNA",
      TRUE ~ Component
    )
  ) %>%
  filter(
    Cell %in% cell_lines,
    Component2 %in% c("cfRNA", "cell_debris")
  )

valid_samples <- intersect(colnames(counts), info2$Sample)

if (length(valid_samples) == 0) {
  stop("No matching samples between expression matrix column names and sample information Sample column, please check sample names.")
}

counts <- counts[, valid_samples, drop = FALSE]

info2 <- info2 %>%
  filter(Sample %in% valid_samples) %>%
  distinct(Sample, .keep_all = TRUE)

info2 <- info2[match(colnames(counts), info2$Sample), ]

stopifnot(all(colnames(counts) == info2$Sample))

cat("\nNumber of retained samples: ", ncol(counts), "\n")
cat("\nComponent2 distribution:\n")
print(table(info2$Component2))
cat("\nCell distribution:\n")
print(table(info2$Cell))
cat("\nCell x Component2 distribution:\n")
print(table(info2$Cell, info2$Component2))

# =========================================================
# 5. Differential analysis function: limma-voom
# =========================================================

run_limma_voom <- function(count_sub, meta_sub, contrast_name, out_prefix) {

  group <- factor(meta_sub$DE_group, levels = c("cell_debris", "cfRNA"))

  dge <- DGEList(counts = count_sub)
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  if (nrow(dge) < 10) {
    warning("Too few genes after filtering, skipping limma: ", out_prefix)
    return(NULL)
  }

  dge <- calcNormFactors(dge)

  design <- model.matrix(~ group)
  colnames(design) <- make.names(colnames(design))

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  res <- topTable(
    fit,
    coef = "groupcfRNA",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  res$Gene <- rownames(res)
  res$method <- "limma_voom"
  res$comparison <- contrast_name

  res <- res %>%
    relocate(Gene, method, comparison)

  # All results
  write.csv(res, paste0(out_prefix, "_limma_voom_all.csv"), row.names = FALSE)

  # Up-regulated: cfRNA > cell_debris
  res_up <- res %>% filter(adj.P.Val < padj_cutoff, logFC > logFC_cutoff)
  write.csv(res_up, paste0(out_prefix, "_limma_voom_up.csv"), row.names = FALSE)

  # Down-regulated: cfRNA < cell_debris
  res_down <- res %>% filter(adj.P.Val < padj_cutoff, logFC < -logFC_cutoff)
  write.csv(res_down, paste0(out_prefix, "_limma_voom_down.csv"), row.names = FALSE)

  return(res)
}

# =========================================================
# 6. Differential analysis function: edgeR
# =========================================================

run_edgeR <- function(count_sub, meta_sub, contrast_name, out_prefix) {

  group <- factor(meta_sub$DE_group, levels = c("cell_debris", "cfRNA"))

  dge <- DGEList(counts = count_sub, group = group)
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  if (nrow(dge) < 10) {
    warning("Too few genes after filtering, skipping edgeR: ", out_prefix)
    return(NULL)
  }

  dge <- calcNormFactors(dge)

  design <- model.matrix(~ group)
  dge <- estimateDisp(dge, design)

  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)

  res <- topTags(qlf, n = Inf)$table

  res$Gene <- rownames(res)
  res$method <- "edgeR"
  res$comparison <- contrast_name

  res <- res %>%
    relocate(Gene, method, comparison)

  write.csv(res, paste0(out_prefix, "_edgeR_all.csv"), row.names = FALSE)

  res_up <- res %>% filter(FDR < padj_cutoff, logFC > logFC_cutoff)
  write.csv(res_up, paste0(out_prefix, "_edgeR_up.csv"), row.names = FALSE)

  res_down <- res %>% filter(FDR < padj_cutoff, logFC < -logFC_cutoff)
  write.csv(res_down, paste0(out_prefix, "_edgeR_down.csv"), row.names = FALSE)

  return(res)
}

# =========================================================
# 7. Differential analysis function: DESeq2
# =========================================================

run_DESeq2 <- function(count_sub, meta_sub, contrast_name, out_prefix) {

  meta_sub$DE_group <- factor(meta_sub$DE_group, levels = c("cell_debris", "cfRNA"))

  count_sub <- round(count_sub)
  count_sub <- count_sub[rowSums(count_sub) > 0, , drop = FALSE]

  if (nrow(count_sub) < 10) {
    warning("Too few genes after filtering, skipping DESeq2: ", out_prefix)
    return(NULL)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = count_sub,
    colData = meta_sub,
    design = ~ DE_group
  )

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]

  if (nrow(dds) < 10) {
    warning("Too few genes after filtering, skipping DESeq2: ", out_prefix)
    return(NULL)
  }

  dds <- DESeq(dds, quiet = TRUE)

  res <- results(
    dds,
    contrast = c("DE_group", "cfRNA", "cell_debris"),
    alpha = padj_cutoff
  )

  res <- as.data.frame(res)
  res$Gene <- rownames(res)
  res$method <- "DESeq2"
  res$comparison <- contrast_name

  res <- res %>%
    relocate(Gene, method, comparison) %>%
    arrange(padj)

  write.csv(res, paste0(out_prefix, "_DESeq2_all.csv"), row.names = FALSE)

  res_up <- res %>% filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > logFC_cutoff)
  write.csv(res_up, paste0(out_prefix, "_DESeq2_up.csv"), row.names = FALSE)

  res_down <- res %>% filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -logFC_cutoff)
  write.csv(res_down, paste0(out_prefix, "_DESeq2_down.csv"), row.names = FALSE)

  return(res)
}

# =========================================================
# 8. cfRNA vs cell_debris comparison for a single cell line
# =========================================================

run_cfRNA_vs_debris <- function(cell_line) {

  cat("\n====================================================\n")
  cat("Cell line:", cell_line, "\n")
  cat("Comparison: cfRNA vs cell_debris\n")
  cat("====================================================\n")

  cell_safe <- safe_name(cell_line)

  meta_sub <- info2 %>%
    filter(Cell == cell_line) %>%
    mutate(DE_group = Component2)

  count_sub <- counts[, meta_sub$Sample, drop = FALSE]

  group_table <- table(meta_sub$DE_group)
  print(group_table)

  if (!all(c("cfRNA", "cell_debris") %in% names(group_table))) {
    warning("Missing cfRNA or cell_debris, skipping: ", cell_line)
    return(NULL)
  }

  if (any(group_table[c("cfRNA", "cell_debris")] < 2)) {
    warning("cfRNA or cell_debris sample count less than 2, results may be unstable: ", cell_line)
  }

  contrast_name <- paste0(cell_line, "_cfRNA_vs_cell_debris")

  out_dir <- file.path(outdir, cell_safe)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix <- file.path(out_dir, safe_name(contrast_name))

  # Save sample grouping
  write.csv(
    meta_sub %>% select(Sample, Cell, Component, Component2, DE_group),
    paste0(out_prefix, "_sample_group.csv"),
    row.names = FALSE
  )

  res_limma <- run_limma_voom(
    count_sub = count_sub,
    meta_sub = meta_sub,
    contrast_name = contrast_name,
    out_prefix = out_prefix
  )

  res_edgeR <- run_edgeR(
    count_sub = count_sub,
    meta_sub = meta_sub,
    contrast_name = contrast_name,
    out_prefix = out_prefix
  )

  res_DESeq2 <- run_DESeq2(
    count_sub = count_sub,
    meta_sub = meta_sub,
    contrast_name = contrast_name,
    out_prefix = out_prefix
  )

  return(list(
    limma = res_limma,
    edgeR = res_edgeR,
    DESeq2 = res_DESeq2
  ))
}

# =========================================================
# 9. Loop over all cell lines
# =========================================================

all_results <- list()

for (cell in cell_lines) {
  res <- run_cfRNA_vs_debris(cell)
  all_results[[cell]] <- res
}

# =========================================================
# 10. Summarize number of significant DE genes
# =========================================================

summary_files_up <- list.files(
  outdir,
  pattern = "_up.csv$",
  recursive = TRUE,
  full.names = TRUE
)

summary_files_down <- list.files(
  outdir,
  pattern = "_down.csv$",
  recursive = TRUE,
  full.names = TRUE
)

summary_up <- data.frame(
  file      = summary_files_up,
  direction = "up",
  n_genes   = sapply(summary_files_up, function(f) nrow(read.csv(f)))
)

summary_down <- data.frame(
  file      = summary_files_down,
  direction = "down",
  n_genes   = sapply(summary_files_down, function(f) nrow(read.csv(f)))
)

summary_df <- bind_rows(summary_up, summary_down) %>%
  arrange(file, direction)

write.csv(
  summary_df,
  file.path(outdir, "summary_DE_gene_number.csv"),
  row.names = FALSE
)

cat("\n✅ All DE analyses finished.\n")
cat("Results saved to:\n")
cat(outdir, "\n")
cat("\nSummary of significant DE gene numbers:\n")
print(summary_df)