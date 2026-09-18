library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# ================== 读取数据 ==================
exp <- read.table("/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_miRNA_rpm.txt",
                  row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)

info <- read_tsv("/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt")

# ================== QC筛选 ==================
info_filtered <- info %>%
  filter(QC %in% c(1, 2) & !Component %in% c("blank", "water")) %>%
  distinct(Sample, .keep_all = TRUE)

# ================== 以样本为单位保留表达矩阵 ==================
keep_samples <- intersect(info_filtered$Sample, colnames(exp))

info_filtered <- info_filtered %>%
  filter(Sample %in% keep_samples)

exp_mat <- exp[, info_filtered$Sample, drop = FALSE]
exp_log <- log2(as.matrix(exp_mat) + 1)

# ================== Spearman 相关 ==================
cor_mat <- cor(exp_log,
               method = "spearman",
               use = "pairwise.complete.obs")

cat("相关性范围:",
    round(min(cor_mat, na.rm = TRUE), 2), "到",
    round(max(cor_mat, na.rm = TRUE), 2), "\n")

# ================== 排序规则 ==================
cell_order <- c("Hep3B2.1-7", "HepG2", "HTR-8/SVneo", "K562", "HEK293T")
comp_order <- c("cfRNA", "debris", "cell")

annotation_col <- info_filtered %>%
  select(Sample, Cell, Component) %>%
  mutate(
    Cell = factor(Cell, levels = cell_order),
    Component = factor(Component, levels = comp_order)
  ) %>%
  arrange(Cell, Component, Sample) %>%
  column_to_rownames("Sample")

ordered_samples <- rownames(annotation_col)
ordered_samples <- ordered_samples[ordered_samples %in% colnames(cor_mat)]

cor_mat_ordered <- cor_mat[ordered_samples, ordered_samples]
annotation_col <- annotation_col[ordered_samples, , drop = FALSE]
annotation_row <- annotation_col

# ================== 细胞系分隔线位置 ==================
cell_sizes <- table(annotation_col$Cell)
cell_sizes <- cell_sizes[cell_sizes > 0]

gap_positions <- cumsum(as.numeric(cell_sizes))
gap_positions <- gap_positions[-length(gap_positions)]

# ================== 配色 ==================
col <- colorRampPalette(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                          "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3",
                          "#2166AC", "#053061"))(200)

breaks <- seq(-1, 1, length.out = 201)

ann_colors <- list(
  Cell = setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(unique(annotation_col$Cell))),
    unique(annotation_col$Cell)
  ),
  Component = setNames(
    colorRampPalette(brewer.pal(8, "Set1"))(length(unique(annotation_col$Component))),
    unique(annotation_col$Component)
  )
)

# ================== 总图：不做层次聚类 ==================
pheatmap(cor_mat_ordered,
         color = col,
         breaks = breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row = gap_positions,
         gaps_col = gap_positions,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 6,
         fontsize_row = 7,
         fontsize_col = 7,
         angle_col = 45,
         border_color = NA,
         filename = "All_Cell_Component_Correlation_miRNA_sample_level_ordered_no_cluster.pdf",
         width = 16,
         height = 16)

# ================== 分图：每个细胞系内部层次聚类 ==================
for (cl in cell_order) {
  samples_cl <- rownames(annotation_col)[annotation_col$Cell == cl]
  samples_cl <- samples_cl[samples_cl %in% colnames(cor_mat)]

  if (length(samples_cl) < 2) {
    next
  }

  cor_cl <- cor_mat[samples_cl, samples_cl, drop = FALSE]
  annotation_cl <- annotation_col[samples_cl, , drop = FALSE]

  hc_cl <- hclust(as.dist(1 - cor_cl), method = "average")

  out_name <- paste0("Correlation_miRNA_", gsub("[^A-Za-z0-9]", "_", cl), "_sample_level_hclust.pdf")

  pheatmap(cor_cl,
           color = col,
           breaks = breaks,
           cluster_rows = hc_cl,
           cluster_cols = hc_cl,
           annotation_row = annotation_cl,
           annotation_col = annotation_cl,
           annotation_colors = ann_colors,
           display_numbers = TRUE,
           number_color = "white",
           fontsize_number = 7,
           fontsize_row = 8,
           fontsize_col = 8,
           angle_col = 45,
           border_color = NA,
           treeheight_row = 70,
           treeheight_col = 70,
           filename = out_name,
           width = 10,
           height = 10)
}

cat("\nDone!\n")
cat("总图输出: All_Cell_Component_Correlation_miRNA_sample_level_ordered_no_cluster.pdf\n")
cat("分图输出: 每个细胞系一个 Correlation_miRNA_*_sample_level_hclust.pdf\n")