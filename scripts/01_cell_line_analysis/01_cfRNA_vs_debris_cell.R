# ============================================================================
# Script: cfRNA vs debris_cell differential expression analysis (stratified by cell type)
#         Add gene annotation (RNA tissue specificity, nTPM, mRNA class)
#         Additionally add cfRNA vs cell separate comparison
# ============================================================================

# Load necessary packages
library(data.table)
library(dplyr)
library(tidyr)

# ------------------------------ Path settings ------------------------------------
expr_file <- "./expression_matrix/all_mRNA_TPM.txt"
info_file <- "./information_table/cellcult_Sample_Information.txt"
anno_file <- "./01_generate_annotation/gene_annotation_filtered.tsv"
out_dir   <- "./02_tpm_Component comparison/"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------ Data reading ------------------------------------
expr <- fread(expr_file, header = TRUE, sep = "\t", data.table = FALSE)
rownames(expr) <- expr[, 1]
expr <- expr[, -1]

info <- fread(info_file, header = TRUE, sep = "\t", data.table = FALSE)

if (file.exists(anno_file)) {
  gene_anno <- fread(anno_file, header = TRUE, sep = "\t", data.table = FALSE, 
                     stringsAsFactors = FALSE, quote = "")
  colnames(gene_anno) <- make.names(colnames(gene_anno))
  cat("Gene annotation file loaded, total", nrow(gene_anno), "records\n")
  required_cols <- c("Gene", "RNA.tissue.specificity", "RNA.tissue.specific.nTPM", "mRNA.class")
  if (!all(required_cols %in% colnames(gene_anno))) {
    warning("Annotation file missing expected columns, actual column names: ", paste(colnames(gene_anno), collapse = ", "))
    for (col in required_cols) {
      if (!col %in% colnames(gene_anno)) gene_anno[[col]] <- NA_character_
    }
  }
  gene_anno <- gene_anno[, required_cols, drop = FALSE]
  gene_anno <- gene_anno[!duplicated(gene_anno$Gene), ]
} else {
  warning("Gene annotation file does not exist: ", anno_file, ", annotation information will not be added")
  gene_anno <- data.frame(Gene = character(), RNA.tissue.specificity = character(),
                          RNA.tissue.specific.nTPM = character(), mRNA.class = character(),
                          stringsAsFactors = FALSE)
}

# ------------------------------ Sample filtering (modified)-----------------------------
info_qc1 <- info %>% filter(QC == 1)
# Keep three components: cfRNA, cell, debris
info_qc1 <- info_qc1 %>% filter(Component %in% c("cfRNA", "cell", "debris"))

# Construct grouping variables for two comparisons
info_qc1 <- info_qc1 %>%
  mutate(
    # Comparison 1: cfRNA vs debris_cell (debris+cell)
    Group_debris_cell = ifelse(Component %in% c("debris", "cell"), "debris_cell", Component),
    # Comparison 2: cfRNA vs cell (keep only cfRNA and cell, exclude debris)
    Group_cell_only = ifelse(Component == "debris", NA_character_, Component)
  )

# Intersection samples
samples_keep <- intersect(info_qc1$Sample, colnames(expr))
if (length(samples_keep) == 0) stop("No common samples!")

expr_sub <- expr[, samples_keep, drop = FALSE]
info_sub <- info_qc1 %>% filter(Sample %in% samples_keep)

# ------------------------------ Gene filtering (based on overall samples)----------------------
prop_expr <- rowSums(expr_sub > 0.2) / ncol(expr_sub)
genes_keep <- prop_expr >= 0.1
expr_filt <- expr_sub[genes_keep, , drop = FALSE]

cat("\nOverall gene filtering: original gene count", nrow(expr_sub), 
    ", retained after filtering", nrow(expr_filt), "genes (TPM>0.1 in at least 20% of samples)\n")
if (nrow(expr_filt) == 0) stop("No genes passed filtering!")

# ------------------------------ Loop analysis by cell type (modified)--------------------
cell_types <- unique(info_sub$Cell)

# Define configurations for two comparisons
comparisons <- list(
  "cfRNA_vs_debris_cell" = list(
    group_col = "Group_debris_cell",
    ref = "cfRNA",
    target = "debris_cell",
    target_label = "debris_cell"
  ),
  "cfRNA_vs_cell" = list(
    group_col = "Group_cell_only",
    ref = "cfRNA",
    target = "cell",
    target_label = "cell"
  )
)

# Store combined results for two comparisons
combined_list <- list()

for (comp_name in names(comparisons)) {
  comp_config <- comparisons[[comp_name]]
  cat("\n====== Starting comparison: ", comp_name, " ======\n")
  
  results_list <- list()  # Store results by cell type for current comparison
  
  for (cell in cell_types) {
    cat("\nProcessing cell type:", cell, "\n")
    
    # Extract sample information for current cell type
    info_cell <- info_sub %>% filter(Cell == cell)
    
    # Use grouping column of current comparison, and remove samples with NA in that column (e.g., debris in cfRNA vs cell comparison)
    info_cell_comp <- info_cell %>% filter(!is.na(.data[[comp_config$group_col]]))
    
    # Check whether both groups exist
    groups_present <- info_cell_comp[[comp_config$group_col]]
    if (!(comp_config$ref %in% groups_present) || !(comp_config$target %in% groups_present)) {
      cat("  Skipped: missing ", comp_config$ref, " or ", comp_config$target, " group\n")
      next
    }
    
    # Expression matrix for current cell type used in this comparison
    samples_comp <- info_cell_comp$Sample
    expr_cell <- expr_filt[, samples_comp, drop = FALSE]
    
    # Separate two groups of samples
    ref_samples <- info_cell_comp %>% filter(.data[[comp_config$group_col]] == comp_config$ref) %>% pull(Sample)
    target_samples <- info_cell_comp %>% filter(.data[[comp_config$group_col]] == comp_config$target) %>% pull(Sample)
    
    expr_ref <- expr_cell[, ref_samples, drop = FALSE]
    expr_target <- expr_cell[, target_samples, drop = FALSE]
    
    # Initialize result vectors
    gene_names <- rownames(expr_cell)
    n_genes <- length(gene_names)
    
    pvals   <- numeric(n_genes)
    stats   <- numeric(n_genes)
    mean_ref <- numeric(n_genes)
    mean_target <- numeric(n_genes)
    median_ref <- numeric(n_genes)
    median_target <- numeric(n_genes)
    
    for (i in seq_along(gene_names)) {
      x <- as.numeric(expr_ref[i, ])
      y <- as.numeric(expr_target[i, ])
      
      mean_ref[i]    <- mean(x, na.rm = TRUE)
      mean_target[i] <- mean(y, na.rm = TRUE)
      median_ref[i]  <- median(x, na.rm = TRUE)
      median_target[i] <- median(y, na.rm = TRUE)
      
      if (length(x) >= 2 && length(y) >= 2) {
        wt <- wilcox.test(x, y, exact = FALSE)
        pvals[i] <- wt$p.value
        stats[i] <- wt$statistic
      } else {
        pvals[i] <- NA
        stats[i] <- NA
      }
    }
    
    # Calculate log2FC: target / ref
    log2fc <- log2( (mean_target + 1e-5) / (mean_ref + 1e-5) )
    
    padj <- rep(NA, n_genes)
    non_na <- !is.na(pvals)
    if (sum(non_na) > 0) {
      padj[non_na] <- p.adjust(pvals[non_na], method = "BH")
    }
    
    # Direction
    direction <- ifelse(log2fc > 0, paste0("up_in_", comp_config$target_label), 
                        paste0("up_in_", comp_config$ref))
    direction[is.na(log2fc)] <- NA
    
    # Build result data frame
    res_cell <- data.frame(
      gene = gene_names,
      mean_ref = mean_ref,
      mean_target = mean_target,
      median_ref = median_ref,
      median_target = median_target,
      log2FC = log2fc,
      wilcox_statistic = stats,
      p_value = pvals,
      p_adjusted = padj,
      direction = direction,
      stringsAsFactors = FALSE
    ) %>% arrange(p_value)
    
    # Rename columns to make them more readable
    colnames(res_cell)[colnames(res_cell) == "mean_ref"] <- paste0("mean_", comp_config$ref)
    colnames(res_cell)[colnames(res_cell) == "mean_target"] <- paste0("mean_", comp_config$target_label)
    colnames(res_cell)[colnames(res_cell) == "median_ref"] <- paste0("median_", comp_config$ref)
    colnames(res_cell)[colnames(res_cell) == "median_target"] <- paste0("median_", comp_config$target_label)
    
    # Add gene annotation
    if (exists("gene_anno") && nrow(gene_anno) > 0) {
      res_cell <- res_cell %>%
        left_join(gene_anno, by = c("gene" = "Gene")) %>%
        mutate(
          RNA.tissue.specificity = if_else(is.na(RNA.tissue.specificity), "none", RNA.tissue.specificity),
          RNA.tissue.specific.nTPM = if_else(is.na(RNA.tissue.specific.nTPM), "none", RNA.tissue.specific.nTPM),
          mRNA.class = if_else(is.na(mRNA.class), "none", mRNA.class)
        )
    } else {
      res_cell$RNA.tissue.specificity <- "none"
      res_cell$RNA.tissue.specific.nTPM <- "none"
      res_cell$mRNA.class <- "none"
    }
    
    res_cell <- res_cell %>%
      select(gene, RNA.tissue.specificity, RNA.tissue.specific.nTPM, mRNA.class, everything())
    
    results_list[[cell]] <- res_cell
    
    # Output single cell type result
    safe_cell_name <- gsub("/", "_", cell)
    out_file <- file.path(out_dir, paste0(safe_cell_name, "_", comp_name, "_Wilcoxon_results.txt"))
    fwrite(res_cell, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("  Results saved to:", out_file, "\n")
  }
  
  # Merge all cell type results for current comparison
  if (length(results_list) > 0) {
    combined_res <- bind_rows(results_list, .id = "Cell_type")
    combined_file <- file.path(out_dir, paste0("all_celltypes_", comp_name, "_Wilcoxon_results.txt"))
    fwrite(combined_res, file = combined_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\n", comp_name, "combined results saved to:", combined_file, "\n")
    combined_list[[comp_name]] <- combined_res
  } else {
    cat("\n", comp_name, "no results generated\n")
  }
}

cat("\nAll analyses completed!\n")