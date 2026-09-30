# ============================================================
# Script: PCA plots for multiple RNA types
# Purpose:
#   Read expression matrices of lncRNA, mRNA, mlncRNA and miRNA
#   plus sample information, keep QC == 1 samples, perform PCA
#   on log2-transformed expression values, and save 2D PCA
#   scatter plots (PC1 vs PC2) colored by Cell and shaped by
#   Component for each RNA type.
# ============================================================

# Load necessary packages
library(ggplot2)
library(dplyr)

# Set file paths
lncRNA_file  <- "./all_lncRNA_TPM.txt"
mRNA_file    <- "./all_mRNA_TPM.txt"
mlncRNA_file <- "./all_mlRNA_TPM.txt"
miRNA_file   <- "./all_miRNA_rpm.txt"
sample_info_file <- "./cellcult_Sample_Information.txt"
output_dir   <- "./01_figure1/05_PCA/"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read sample information table, keep only QC == 1 samples
sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_info_filtered <- sample_info %>%
  filter(QC == 1) %>%
  mutate(
    Cell      = factor(Cell, levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")),
    Component = factor(Component, levels = c("cfRNA", "debris", "cell"))
  )

# ================== PCA plotting function ==================
plot_pca <- function(expr_file, title, sample_info_filtered) {
  expr_mat <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  expr_mat <- as.matrix(expr_mat)
  valid_samples <- intersect(colnames(expr_mat), sample_info_filtered$Sample)
  if (length(valid_samples) == 0) stop("No valid samples found.")
  expr_mat_filtered <- expr_mat[, valid_samples, drop = FALSE]
  expr_log <- log2(expr_mat_filtered + 1)
  pca_result <- prcomp(t(expr_log), center = TRUE, scale. = FALSE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$Sample <- rownames(pca_df)
  pca_df <- pca_df %>% left_join(sample_info_filtered, by = "Sample")
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell, shape = Component)) +
    geom_point(size = 3, stroke = 1.2) +
    scale_shape_manual(values = c("cfRNA"  = 16,  # solid circle
                                  "debris" = 3,   # hollow circle with cross
                                  "cell"   = 1),  # hollow circle
                       breaks  = c("cfRNA", "debris", "cell")) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = title,
         color = "Cell Type",
         shape = "Component") +
    theme_minimal() +
    theme(legend.position = "right")
  return(p)
}

# ================== Plot and save ==================
# lncRNA
p_lncRNA <- plot_pca(lncRNA_file, "PCA of lncRNA Expression (log2(TPM+1))", sample_info_filtered)
print(p_lncRNA)
ggsave(file.path(output_dir, "PCA_lncRNA.pdf"), plot = p_lncRNA, width = 10, height = 8)

# mRNA
p_mRNA <- plot_pca(mRNA_file, "PCA of mRNA Expression (log2(TPM+1))", sample_info_filtered)
print(p_mRNA)
ggsave(file.path(output_dir, "PCA_mRNA.pdf"), plot = p_mRNA, width = 10, height = 8)

# mlncRNA
if (file.exists(mlncRNA_file)) {
  p_mlncRNA <- plot_pca(mlncRNA_file, "PCA of mlncRNA Expression (log2(TPM+1))", sample_info_filtered)
  print(p_mlncRNA)
  ggsave(file.path(output_dir, "PCA_mlncRNA.pdf"), plot = p_mlncRNA, width = 10, height = 8)
} else {
  warning("mlncRNA file not found: ", mlncRNA_file)
}

# miRNA
if (file.exists(miRNA_file)) {
  p_miRNA <- plot_pca(miRNA_file, "PCA of miRNA Expression (log2(RPM+1))", sample_info_filtered)
  print(p_miRNA)
  ggsave(file.path(output_dir, "PCA_miRNA.pdf"), plot = p_miRNA, width = 10, height = 8)
} else {
  warning("miRNA file not found: ", miRNA_file)
}

cat("All PCA plots saved to:", output_dir, "\n")