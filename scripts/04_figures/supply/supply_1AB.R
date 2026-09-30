# 加载必要的包
library(ggplot2)
library(ggrepel)
library(dplyr)

# 设置文件路径
lncRNA_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_lncRNA_TPM.txt"
mRNA_file   <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mRNA_TPM.txt"
sample_info_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
output_dir  <- "/data/work/01_2603cell_culture/02_qc_batch/02_batch_PCA/"

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取样本信息表
sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 筛选QC为1或2的样本
sample_info_filtered <- sample_info %>% filter(QC %in% c(1,2))

# ================== 原有PCA图生成（保持不变） ==================
# 定义函数：绘制PCA图（无凸包）
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
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell, shape = taskID, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    scale_shape_manual(values = c("W250219001" = 16, "W250319023" = 1)) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = title,
         color = "Cell Type",
         shape = "Task ID") +
    theme_minimal() +
    theme(legend.position = "right")
  return(p)
}

# 绘制并保存原有PCA图
p_lncRNA <- plot_pca(lncRNA_file, "PCA of lncRNA Expression (log2(TPM+1))", sample_info_filtered)
p_mRNA   <- plot_pca(mRNA_file,   "PCA of mRNA Expression (log2(TPM+1))",   sample_info_filtered)
print(p_lncRNA)
print(p_mRNA)
ggsave(file.path(output_dir, "PCA_lncRNA.pdf"), plot = p_lncRNA, width = 10, height = 8)
ggsave(file.path(output_dir, "PCA_mRNA.pdf"),   plot = p_mRNA,   width = 10, height = 8)

# ================== 新增：只包含blank和water的PCA图（带凸包） ==================
# 筛选blank和water样本
sample_info_blank_water <- sample_info_filtered %>% filter(Cell %in% c("blank", "water"))

# 检查是否有数据
if (nrow(sample_info_blank_water) == 0) {
  warning("No blank or water samples found after QC filtering. Skip drawing hull plots.")
} else {
  # 定义带凸包的PCA绘图函数
  plot_pca_with_hull <- function(expr_file, title, sample_info) {
    expr_mat <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    expr_mat <- as.matrix(expr_mat)
    valid_samples <- intersect(colnames(expr_mat), sample_info$Sample)
    if (length(valid_samples) == 0) stop("No valid samples for blank/water.")
    expr_mat_filtered <- expr_mat[, valid_samples, drop = FALSE]
    expr_log <- log2(expr_mat_filtered + 1)
    pca_result <- prcomp(t(expr_log), center = TRUE, scale. = FALSE)
    pca_df <- as.data.frame(pca_result$x[, 1:2])
    pca_df$Sample <- rownames(pca_df)
    pca_df <- pca_df %>% left_join(sample_info, by = "Sample")
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
    
    # 计算凸包：按Cell和taskID分组，提取每个组在PC1-PC2平面上的凸包顶点
    hull_data <- pca_df %>%
      group_by(Cell, taskID) %>%
      # 至少需要3个点才能构成凸包
      filter(n() >= 3) %>%
      slice(chull(PC1, PC2))
    
    # 绘图：先画凸包（透明填充），再画点、标签等
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell, shape = taskID, label = Sample)) +
      # 凸包层：按Cell和taskID分组，颜色继承自Cell（与点一致），填充透明
      geom_polygon(data = hull_data, aes(group = interaction(Cell, taskID), color = Cell),
                   fill = NA, linewidth = 0.8, alpha = 0.5) +
      geom_point(size = 3) +
      geom_text_repel(size = 3, max.overlaps = 20) +
      scale_shape_manual(values = c("W250219001" = 16, "W250319023" = 1)) +
      labs(x = paste0("PC1 (", var_explained[1], "%)"),
           y = paste0("PC2 (", var_explained[2], "%)"),
           title = title,
           color = "Cell Type",
           shape = "Task ID") +
      theme_minimal() +
      theme(legend.position = "right")
    
    return(p)
  }
  
  # 绘制blank+water的凸包PCA图
  p_lncRNA_hull <- plot_pca_with_hull(lncRNA_file, "PCA of lncRNA Expression (blank & water)", sample_info_blank_water)
  p_mRNA_hull   <- plot_pca_with_hull(mRNA_file,   "PCA of mRNA Expression (blank & water)",   sample_info_blank_water)
  
  # 显示并保存
  print(p_lncRNA_hull)
  print(p_mRNA_hull)
  ggsave(file.path(output_dir, "PCA_lncRNA_blank_water_with_hull.pdf"), plot = p_lncRNA_hull, width = 10, height = 8)
  ggsave(file.path(output_dir, "PCA_mRNA_blank_water_with_hull.pdf"),   plot = p_mRNA_hull,   width = 10, height = 8)
  
  cat("Blank & water PCA plots with convex hulls have been saved to:", output_dir, "\n")
}


# ================== 新增：mlncRNA的PCA图 ==================
mlncRNA_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mlRNA_TPM.txt"

# 检查文件是否存在，若存在则生成图
if (file.exists(mlncRNA_file)) {
  # 全样本PCA图
  p_mlncRNA <- plot_pca(mlncRNA_file, "PCA of mlncRNA Expression (log2(TPM+1))", sample_info_filtered)
  print(p_mlncRNA)
  ggsave(file.path(output_dir, "PCA_mlncRNA.pdf"), plot = p_mlncRNA, width = 10, height = 8)
  
  # blank和water组的凸包PCA图（复用之前筛选的sample_info_blank_water）
  if (nrow(sample_info_blank_water) > 0) {
    p_mlncRNA_hull <- plot_pca_with_hull(mlncRNA_file, "PCA of mlncRNA Expression (blank & water)", sample_info_blank_water)
    print(p_mlncRNA_hull)
    ggsave(file.path(output_dir, "PCA_mlncRNA_blank_water_with_hull.pdf"), plot = p_mlncRNA_hull, width = 10, height = 8)
    cat("mlncRNA PCA plots with convex hulls have been saved to:", output_dir, "\n")
  } else {
    warning("No blank or water samples for mlncRNA hull plots.")
  }
} else {
  warning("mlncRNA file not found: ", mlncRNA_file)
}