# 加载必要的包
library(ggplot2)
library(dplyr)

# 设置文件路径
lncRNA_file  <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_lncRNA_TPM.txt"
mRNA_file    <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mRNA_TPM.txt"
mlncRNA_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mlRNA_TPM.txt"
miRNA_file   <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_miRNA_rpm.txt"
sample_info_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
output_dir   <- "/data/work/01_2603cell_culture/07_figure/01_figure1/05_PCA/"

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取样本信息表，只保留 QC == 1 的样本
sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_info_filtered <- sample_info %>%
  filter(QC == 1) %>%
  mutate(
    Cell      = factor(Cell, levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")),
    Component = factor(Component, levels = c("cfRNA", "debris", "cell"))
  )

# ================== PCA 绘图函数 ==================
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
    scale_shape_manual(values = c("cfRNA"  = 16,  # 实心圆
                                  "debris" = 3,   # 空心圆加十字
                                  "cell"   = 1),  # 空心圆
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

# ================== 绘制并保存 ==================
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
