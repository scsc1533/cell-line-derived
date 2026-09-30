suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ============================================================
# RBP cell-type-specific expression screening with Wilcoxon test
#
# Samples retained:
#   QC == 1
#   Component == "cell"
#
# RBPs retained:
#   RBP_status == "RBP"
#   Composite_Score > 20
#   matched by SYMBOL
#   expressed in at least 20% retained samples with TPM > 1
#
# Comparisons:
#   1. Each individual cell type vs all other cell types
#   2. Hep3B2.1-7 + HepG2 combined as Hep_liver_like vs all others
#
# Notes:
#   Wilcoxon p-values are used as reference only because sample sizes are small.
#   log2FC is calculated from raw TPM means:
#     log2((mean_target_TPM + pseudocount) / (mean_other_TPM + pseudocount))
# ============================================================

# -------------------------
# 1. Input and output paths
# -------------------------

rbp_file <- "./High_confidence_human_RBPs0827.txt"
expr_file <- "./expression_matrix/all_mlRNA_TPM.txt"
meta_file <- "./cellcult_Sample_Information.txt"

outdir <- "./03_RBP_cell_type_specific_wilcox_cell_samples"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 2. Parameters
# -------------------------

component_keep <- "cell"
composite_cutoff <- 20
min_tpm <- 1
min_sample_fraction <- 0.2
pseudocount <- 0.1

log2fc_cutoff <- 1
mean_target_tpm_cutoff <- 1
target_positive_fraction_cutoff <- 0.5
fdr_cutoff <- 0.05

top_n_per_group <- 100
combined_group_name <- "Hep_liver_like"
combined_group_cells <- c("Hep3B2.1-7", "HepG2")

# -------------------------
# 3. Helper functions
# -------------------------

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x)
}

median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  median(x)
}

wilcox_greater <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 1 || length(y) < 1) {
    return(NA_real_)
  }
  if (length(unique(c(x, y))) < 2) {
    return(1)
  }

  tryCatch(
    wilcox.test(x, y, alternative = "greater", exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

write_tsv <- function(x, file) {
  write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_matrix <- function(mat, file, id_col = "SYMBOL") {
  write.table(
    data.frame(setNames(list(rownames(mat)), id_col), mat, check.names = FALSE),
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

run_one_vs_rest <- function(group_name, target_samples, meta_use, tpm_mat, log2_mat) {
  other_samples <- setdiff(colnames(tpm_mat), target_samples)

  target_tpm <- tpm_mat[, target_samples, drop = FALSE]
  other_tpm <- tpm_mat[, other_samples, drop = FALSE]
  target_log2 <- log2_mat[, target_samples, drop = FALSE]
  other_log2 <- log2_mat[, other_samples, drop = FALSE]

  mean_target_tpm <- apply(target_tpm, 1, mean_or_na)
  mean_other_tpm <- apply(other_tpm, 1, mean_or_na)
  median_target_tpm <- apply(target_tpm, 1, median_or_na)
  median_other_tpm <- apply(other_tpm, 1, median_or_na)
  mean_target_log2 <- apply(target_log2, 1, mean_or_na)
  mean_other_log2 <- apply(other_log2, 1, mean_or_na)
  target_positive_fraction <- rowMeans(target_tpm > min_tpm, na.rm = TRUE)
  other_positive_fraction <- rowMeans(other_tpm > min_tpm, na.rm = TRUE)

  wilcox_p <- sapply(rownames(tpm_mat), function(gene_i) {
    wilcox_greater(
      as.numeric(tpm_mat[gene_i, target_samples]),
      as.numeric(tpm_mat[gene_i, other_samples])
    )
  })

  data.frame(
    SYMBOL = rownames(tpm_mat),
    target_group = group_name,
    target_cells = paste(sort(unique(meta_use[target_samples, "Cell"])), collapse = ";"),
    other_cells = paste(sort(unique(meta_use[other_samples, "Cell"])), collapse = ";"),
    target_sample_n = length(target_samples),
    other_sample_n = length(other_samples),
    target_samples = paste(target_samples, collapse = ";"),
    other_samples = paste(other_samples, collapse = ";"),
    mean_target_TPM = as.numeric(mean_target_tpm),
    mean_other_TPM = as.numeric(mean_other_tpm),
    median_target_TPM = as.numeric(median_target_tpm),
    median_other_TPM = as.numeric(median_other_tpm),
    mean_target_log2TPM = as.numeric(mean_target_log2),
    mean_other_log2TPM = as.numeric(mean_other_log2),
    log2FC_target_vs_other_TPM = log2((mean_target_tpm + pseudocount) / (mean_other_tpm + pseudocount)),
    TPM_difference_target_minus_other = as.numeric(mean_target_tpm - mean_other_tpm),
    target_positive_fraction = as.numeric(target_positive_fraction),
    other_positive_fraction = as.numeric(other_positive_fraction),
    wilcox_p_greater = as.numeric(wilcox_p),
    stringsAsFactors = FALSE
  )
}

# -------------------------
# 4. Read files
# -------------------------

rbp <- fread(rbp_file, data.table = FALSE, check.names = FALSE)
expr <- fread(expr_file, data.table = FALSE, check.names = FALSE)
meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)

required_rbp_cols <- c("SYMBOL", "RBP_status", "Composite_Score")
required_expr_cols <- c("sym_id")
required_meta_cols <- c("Cell", "Sample", "Component", "day", "QC")

missing_rbp_cols <- setdiff(required_rbp_cols, colnames(rbp))
missing_expr_cols <- setdiff(required_expr_cols, colnames(expr))
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_rbp_cols) > 0) {
  stop("Missing columns in RBP table: ", paste(missing_rbp_cols, collapse = ", "))
}
if (length(missing_expr_cols) > 0) {
  stop("Missing columns in expression matrix: ", paste(missing_expr_cols, collapse = ", "))
}
if (length(missing_meta_cols) > 0) {
  stop("Missing columns in metadata: ", paste(missing_meta_cols, collapse = ", "))
}

rbp$Composite_Score <- as.numeric(rbp$Composite_Score)
if ("RBP2GO_Score" %in% colnames(rbp)) {
  rbp$RBP2GO_Score <- as.numeric(rbp$RBP2GO_Score)
} else {
  rbp$RBP2GO_Score <- NA_real_
}
meta$QC <- as.numeric(meta$QC)
meta$day <- as.numeric(meta$day)

# -------------------------
# 5. Select candidate RBPs and QC=1 cell samples
# -------------------------

rbp_candidate <- rbp %>%
  filter(
    RBP_status == "RBP",
    Composite_Score > composite_cutoff,
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  arrange(desc(Composite_Score), desc(RBP2GO_Score)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

rbp_symbols <- unique(rbp_candidate$SYMBOL)
expr_samples <- setdiff(colnames(expr), "sym_id")

meta_cell <- meta %>%
  filter(
    QC == 1,
    Component == component_keep,
    Sample %in% expr_samples
  ) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  arrange(Cell, day, Sample)

cell_samples <- meta_cell$Sample
cell_types <- unique(meta_cell$Cell)

if (length(cell_samples) < 2) {
  stop("Fewer than 2 QC=1 cell samples were found in both metadata and expression matrix.")
}
if (length(cell_types) < 2) {
  stop("Fewer than 2 cell types were found among QC=1 cell samples.")
}
if (!all(combined_group_cells %in% cell_types)) {
  stop(
    "Combined group cells not found among selected samples: ",
    paste(setdiff(combined_group_cells, cell_types), collapse = ", ")
  )
}

# -------------------------
# 6. Build RBP expression matrix
# -------------------------

expr_mat <- expr %>%
  select(sym_id, all_of(cell_samples))

expr_mat[, cell_samples] <- lapply(expr_mat[, cell_samples, drop = FALSE], function(x) {
  as.numeric(as.character(x))
})

expr_mat <- expr_mat %>%
  group_by(sym_id) %>%
  summarise(across(all_of(cell_samples), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

expr_rbp <- expr_mat %>%
  filter(sym_id %in% rbp_symbols)

if (nrow(expr_rbp) < 2) {
  stop("Fewer than 2 RBP SYMBOLs matched the expression matrix.")
}

rbp_tpm <- as.data.frame(expr_rbp)
rownames(rbp_tpm) <- rbp_tpm$sym_id
rbp_tpm$sym_id <- NULL
rbp_tpm <- as.matrix(rbp_tpm)
mode(rbp_tpm) <- "numeric"
rbp_tpm <- rbp_tpm[rownames(rbp_tpm) %in% rbp_symbols, cell_samples, drop = FALSE]

min_sample_n <- ceiling(ncol(rbp_tpm) * min_sample_fraction)
keep_expr <- rowSums(rbp_tpm > min_tpm, na.rm = TRUE) >= min_sample_n
rbp_tpm_filt <- rbp_tpm[keep_expr, , drop = FALSE]

if (nrow(rbp_tpm_filt) < 2) {
  stop("Fewer than 2 expressed RBPs remained after TPM filtering.")
}

rbp_log2 <- log2(rbp_tpm_filt + 1)
finite_gene <- apply(rbp_log2, 1, function(x) all(is.finite(x)))
finite_sample <- apply(rbp_log2, 2, function(x) all(is.finite(x)))
rbp_log2 <- rbp_log2[finite_gene, finite_sample, drop = FALSE]
rbp_tpm_filt <- rbp_tpm_filt[rownames(rbp_log2), colnames(rbp_log2), drop = FALSE]

meta_cell <- meta_cell %>%
  filter(Sample %in% colnames(rbp_tpm_filt)) %>%
  distinct(Sample, .keep_all = TRUE)
rownames(meta_cell) <- meta_cell$Sample
meta_cell <- meta_cell[colnames(rbp_tpm_filt), , drop = FALSE]

if (!identical(rownames(meta_cell), colnames(rbp_tpm_filt))) {
  stop("Sample metadata order does not match expression matrix columns.")
}

cell_types <- unique(meta_cell$Cell)

# -------------------------
# 7. Run one-vs-rest tests
# -------------------------

comparison_list <- list()

for (cell_i in cell_types) {
  target_samples <- rownames(meta_cell)[meta_cell$Cell == cell_i]
  comparison_list[[cell_i]] <- run_one_vs_rest(
    group_name = cell_i,
    target_samples = target_samples,
    meta_use = meta_cell,
    tpm_mat = rbp_tpm_filt,
    log2_mat = rbp_log2
  )
}

combined_samples <- rownames(meta_cell)[meta_cell$Cell %in% combined_group_cells]
comparison_list[[combined_group_name]] <- run_one_vs_rest(
  group_name = combined_group_name,
  target_samples = combined_samples,
  meta_use = meta_cell,
  tpm_mat = rbp_tpm_filt,
  log2_mat = rbp_log2
)

all_results <- bind_rows(comparison_list) %>%
  group_by(target_group) %>%
  mutate(FDR_by_target_group = p.adjust(wilcox_p_greater, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global = p.adjust(wilcox_p_greater, method = "BH")) %>%
  left_join(rbp_candidate, by = "SYMBOL") %>%
  arrange(
    target_group,
    desc(log2FC_target_vs_other_TPM),
    desc(TPM_difference_target_minus_other),
    FDR_by_target_group,
    desc(mean_target_TPM)
  )

candidate_results <- all_results %>%
  filter(
    mean_target_TPM >= mean_target_tpm_cutoff,
    log2FC_target_vs_other_TPM >= log2fc_cutoff,
    target_positive_fraction >= target_positive_fraction_cutoff,
    FDR_by_target_group <= fdr_cutoff
  ) %>%
  arrange(
    target_group,
    FDR_by_target_group,
    desc(log2FC_target_vs_other_TPM),
    desc(mean_target_TPM)
  )

effect_size_candidates <- all_results %>%
  filter(
    mean_target_TPM >= mean_target_tpm_cutoff,
    log2FC_target_vs_other_TPM >= log2fc_cutoff,
    target_positive_fraction >= target_positive_fraction_cutoff
  ) %>%
  arrange(
    target_group,
    desc(log2FC_target_vs_other_TPM),
    FDR_by_target_group,
    desc(mean_target_TPM)
  )

top_by_group <- all_results %>%
  group_by(target_group) %>%
  arrange(
    desc(log2FC_target_vs_other_TPM),
    FDR_by_target_group,
    desc(mean_target_TPM),
    .by_group = TRUE
  ) %>%
  slice_head(n = top_n_per_group) %>%
  ungroup()

# -------------------------
# 8. Output
# -------------------------

write_tsv(
  all_results,
  file.path(outdir, "RBP_cell_type_one_vs_rest_wilcox_all_results.tsv")
)

write_tsv(
  candidate_results,
  file.path(outdir, "RBP_cell_type_specific_candidates_wilcox_FDR005.tsv")
)

write_tsv(
  effect_size_candidates,
  file.path(outdir, "RBP_cell_type_specific_candidates_effect_size_wilcox_reference.tsv")
)

write_tsv(
  top_by_group,
  file.path(outdir, paste0("Top", top_n_per_group, "_RBP_per_group_by_log2FC_wilcox_reference.tsv"))
)

write_tsv(
  all_results %>% filter(target_group == combined_group_name),
  file.path(outdir, "RBP_Hep_liver_like_vs_rest_wilcox_all_results.tsv")
)

write_tsv(
  candidate_results %>% filter(target_group == combined_group_name),
  file.path(outdir, "RBP_Hep_liver_like_specific_candidates_wilcox_FDR005.tsv")
)

write_tsv(
  effect_size_candidates %>% filter(target_group == combined_group_name),
  file.path(outdir, "RBP_Hep_liver_like_specific_candidates_effect_size_wilcox_reference.tsv")
)

write_matrix(
  rbp_tpm_filt,
  file.path(outdir, "RBP_TPM_QC1_cell_samples.tsv")
)

write_matrix(
  rbp_log2,
  file.path(outdir, "RBP_log2TPM_QC1_cell_samples.tsv")
)

write_tsv(
  meta_cell,
  file.path(outdir, "cell_sample_metadata_QC1.tsv")
)

sample_count <- meta_cell %>%
  count(Cell, name = "sample_n") %>%
  arrange(Cell)

write_tsv(
  sample_count,
  file.path(outdir, "cell_sample_count_by_cell_type.tsv")
)

summary_df <- data.frame(
  item = c(
    "Candidate RBP unique SYMBOLs",
    "QC=1 cell samples matched expression matrix",
    "Cell types among selected cell samples",
    "RBP SYMBOLs matched expression matrix",
    "RBP retained after TPM filtering",
    "Comparison groups",
    "All result rows",
    "Candidate rows after effect-size thresholds",
    "Candidate rows after effect-size and FDR thresholds",
    "Combined group name",
    "Combined group cells",
    "log2FC cutoff",
    "mean target TPM cutoff",
    "target positive fraction cutoff",
    "FDR cutoff",
    "Wilcoxon alternative"
  ),
  value = c(
    length(rbp_symbols),
    ncol(rbp_tpm_filt),
    length(cell_types),
    nrow(rbp_tpm),
    nrow(rbp_tpm_filt),
    length(comparison_list),
    nrow(all_results),
    nrow(effect_size_candidates),
    nrow(candidate_results),
    combined_group_name,
    paste(combined_group_cells, collapse = ";"),
    log2fc_cutoff,
    mean_target_tpm_cutoff,
    target_positive_fraction_cutoff,
    fdr_cutoff,
    "greater"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_df,
  file.path(outdir, "RBP_cell_type_specific_wilcox_summary.tsv")
)

cat("Candidate RBP SYMBOLs:", length(rbp_symbols), "\n")
cat("QC=1 cell samples matched expression matrix:", ncol(rbp_tpm_filt), "\n")
cat("Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("Sample count by cell type:\n")
print(sample_count)
cat("RBP retained after TPM filtering:", nrow(rbp_tpm_filt), "\n")
cat("Comparison groups:", paste(names(comparison_list), collapse = ", "), "\n")
cat("All result rows:", nrow(all_results), "\n")
cat("Candidate rows after effect-size thresholds:", nrow(effect_size_candidates), "\n")
cat("Candidate rows after effect-size and FDR thresholds:", nrow(candidate_results), "\n")

cat("\nTop candidates per group by log2FC:\n")
print(
  top_by_group %>%
    group_by(target_group) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    select(
      SYMBOL,
      target_group,
      target_sample_n,
      other_sample_n,
      mean_target_TPM,
      mean_other_TPM,
      log2FC_target_vs_other_TPM,
      wilcox_p_greater,
      FDR_by_target_group
    )
)

cat("\nDone.\n")
cat("Output directory:\n", outdir, "\n")

