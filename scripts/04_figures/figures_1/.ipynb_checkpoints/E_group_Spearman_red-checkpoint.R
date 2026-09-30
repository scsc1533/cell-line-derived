#!/usr/bin/env Rscript

# ============================================================
# Figure 10
# 4-mer motif proportion correlation between
# Cell line × Component groups
#
# 图形要求：
# 1. 与原图一样：每个分组为 Cell line | Component
# 2. 绘制完整相关性矩阵
# 3. 使用 Spearman 相关系数
# 4. 保留行列聚类树
# 5. 保留行列分组名称
# 6. 每个方格显示相关系数
# 7. 只使用白色—红色渐变，不使用蓝色
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# ============================================================
# 1. 输入和输出路径
# ============================================================

qc_file <- paste0(
  "/data/work/1/cellcult_Sample_Information.txt"
)

# 已复制到工作目录的4-mer文件
motif_file <- "/data/work/sc/sample_mRNA_4motif_fre.txt"

outdir <- paste0(
  "/data/work/sc/redraw_figures_0717/",
  "10_fig10_group_Spearman_red"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(qc_file)) {
  stop("找不到样本信息文件：", qc_file)
}

if (!file.exists(motif_file)) {
  stop("找不到4-mer频率文件：", motif_file)
}

cat("样本信息文件：", qc_file, "\n")
cat("4-mer频率文件：", motif_file, "\n")
cat("输出目录：", outdir, "\n")

# ============================================================
# 2. 读取样本信息
# ============================================================

qc_data <- read.delim(
  qc_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_qc_columns <- c(
  "Sample",
  "Cell",
  "Component",
  "QC"
)

missing_qc_columns <- setdiff(
  required_qc_columns,
  colnames(qc_data)
)

if (length(missing_qc_columns) > 0) {
  stop(
    "样本信息表缺少以下列：",
    paste(missing_qc_columns, collapse = ", ")
  )
}

qc_use <- qc_data %>%
  transmute(
    Sample = trimws(as.character(Sample)),
    Cell = trimws(as.character(Cell)),
    Component = trimws(as.character(Component)),
    QC = suppressWarnings(
      as.numeric(as.character(QC))
    )
  ) %>%
  filter(
    QC == 1,
    Component %in% c(
      "cfRNA",
      "debris",
      "cell"
    ),
    !is.na(Sample),
    Sample != "",
    !is.na(Cell),
    Cell != ""
  ) %>%
  distinct(
    Sample,
    .keep_all = TRUE
  )

if (nrow(qc_use) == 0) {
  stop("QC==1且属于目标组分的样本数为0。")
}

cat("\nQC筛选后的样本数：", nrow(qc_use), "\n")

cat("\n各 Cell × Component 样本数：\n")
print(
  table(
    qc_use$Cell,
    qc_use$Component,
    useNA = "ifany"
  )
)

# ============================================================
# 3. 读取4-mer频率矩阵
#
# 文件格式：
# 第一列：Motif
# 后续各列：样本
# ============================================================

motif_data <- read.delim(
  motif_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(motif_data) < 3) {
  stop(
    "4-mer文件列数过少，",
    "应至少包括1列Motif和2列样本。"
  )
}

motif_colname <- colnames(motif_data)[1]
sample_columns <- colnames(motif_data)[-1]

cat("\nMotif列名：", motif_colname, "\n")
cat("4-mer文件中的样本列数：", length(sample_columns), "\n")

# ============================================================
# 4. 匹配QC样本和4-mer样本
# ============================================================

valid_samples <- intersect(
  qc_use$Sample,
  sample_columns
)

if (length(valid_samples) == 0) {
  stop(
    "4-mer文件中的样本列与QC样本没有匹配。\n",
    "请检查Sample名称是否一致。"
  )
}

# 保持QC表中的顺序
valid_samples <- qc_use$Sample[
  qc_use$Sample %in% valid_samples
]

cat("匹配到的样本数：", length(valid_samples), "\n")

unmatched_qc <- setdiff(
  qc_use$Sample,
  sample_columns
)

if (length(unmatched_qc) > 0) {
  cat(
    "QC表中未在4-mer文件中找到的样本数：",
    length(unmatched_qc),
    "\n"
  )
}

# ============================================================
# 5. 转换成长表
# ============================================================

motif_filtered <- motif_data %>%
  select(
    all_of(motif_colname),
    all_of(valid_samples)
  ) %>%
  rename(
    Motif = all_of(motif_colname)
  )

motif_long <- motif_filtered %>%
  pivot_longer(
    cols = all_of(valid_samples),
    names_to = "Sample",
    values_to = "Proportion"
  ) %>%
  mutate(
    Sample = trimws(as.character(Sample)),
    Motif = as.character(Motif),
    Proportion = suppressWarnings(
      as.numeric(as.character(Proportion))
    )
  ) %>%
  left_join(
    qc_use %>%
      select(
        Sample,
        Cell,
        Component
      ),
    by = "Sample"
  ) %>%
  filter(
    !is.na(Cell),
    !is.na(Component),
    is.finite(Proportion)
  )

if (nrow(motif_long) == 0) {
  stop("转换成长表并合并样本信息后没有有效数据。")
}

write.table(
  motif_long,
  file = file.path(
    outdir,
    "Fig10_motif_original_long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 6. 对每个 Cell × Component 分组计算每个Motif的平均比例
#
# 例如：
# HepG2 | cfRNA
# HepG2 | debris
# HepG2 | cell
# ============================================================

group_summary <- motif_long %>%
  group_by(
    Cell,
    Component,
    Motif
  ) %>%
  summarise(
    MeanProportion = mean(
      Proportion,
      na.rm = TRUE
    ),
    SD = sd(
      Proportion,
      na.rm = TRUE
    ),
    sample_n = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(
      Cell,
      Component,
      sep = " | "
    )
  )

write.table(
  group_summary,
  file = file.path(
    outdir,
    "Fig10_motif_mean_by_Cell_Component.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\n最终 Cell × Component 分组数：",
  n_distinct(group_summary$Group),
  "\n"
)

cat("\n最终分组名称：\n")
print(
  sort(
    unique(group_summary$Group)
  )
)

# ============================================================
# 7. 构建 Motif × Group 矩阵
# ============================================================

group_matrix_df <- group_summary %>%
  select(
    Motif,
    Group,
    MeanProportion
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = MeanProportion
  )

group_matrix <- group_matrix_df %>%
  column_to_rownames(
    var = "Motif"
  ) %>%
  as.matrix()

mode(group_matrix) <- "numeric"

# 删除全部为NA的Motif
keep_motif <- rowSums(
  is.finite(group_matrix)
) > 0

group_matrix <- group_matrix[
  keep_motif,
  ,
  drop = FALSE
]

# 删除有效Motif过少或数值无变化的分组
group_qc <- tibble(
  Group = colnames(group_matrix),

  valid_motif_n = apply(
    group_matrix,
    2,
    function(x) {
      sum(is.finite(x))
    }
  ),

  motif_sd = apply(
    group_matrix,
    2,
    function(x) {
      sd(
        x[is.finite(x)],
        na.rm = TRUE
      )
    }
  )
)

write.table(
  group_qc,
  file = file.path(
    outdir,
    "Fig10_group_matrix_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

valid_groups <- group_qc %>%
  filter(
    valid_motif_n >= 10,
    is.finite(motif_sd),
    motif_sd > 0
  ) %>%
  pull(Group)

if (length(valid_groups) < 2) {
  stop("通过质量检查的分组少于2个，无法计算相关性。")
}

group_matrix <- group_matrix[
  ,
  valid_groups,
  drop = FALSE
]

cat("用于相关性分析的Motif数：", nrow(group_matrix), "\n")
cat("用于相关性分析的分组数：", ncol(group_matrix), "\n")

write.table(
  group_matrix,
  file = file.path(
    outdir,
    "Fig10_motif_by_group_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ============================================================
# 8. 计算分组之间的Spearman相关性
#
# 行列均为：
# Cell line | Component
# ============================================================

cor_rho <- cor(
  group_matrix,
  use = "pairwise.complete.obs",
  method = "spearman"
)

# 消除极小浮点误差
cor_rho <- (
  cor_rho +
    t(cor_rho)
) / 2

diag(cor_rho) <- 1

if (anyNA(cor_rho)) {
  stop(
    "Spearman相关性矩阵中存在NA，",
    "请检查部分分组的有效Motif数量。"
  )
}

cor_rho[cor_rho > 1] <- 1
cor_rho[cor_rho < -1] <- -1

write.csv(
  cor_rho,
  file = file.path(
    outdir,
    "Fig10_group_Spearman_correlation_matrix.csv"
  ),
  row.names = TRUE
)

cor_long <- as.data.frame(
  as.table(cor_rho),
  stringsAsFactors = FALSE
) %>%
  rename(
    Group1 = Var1,
    Group2 = Var2,
    Spearman_rho = Freq
  )

write.table(
  cor_long,
  file = file.path(
    outdir,
    "Fig10_group_Spearman_correlation_long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nSpearman相关性范围：",
  min(cor_rho, na.rm = TRUE),
  "至",
  max(cor_rho, na.rm = TRUE),
  "\n"
)

# ============================================================
# 9. 设置单一红色渐变
#
# 低相关：白色或极浅红
# 高相关：深红色
# 不再使用蓝色
# ============================================================

red_palette <- colorRampPalette(
  c(
    "#FFFFFF",
    "#FFF5F0",
    "#FEE0D2",
    "#FCBBA1",
    "#FC9272",
    "#FB6A4A",
    "#EF3B2C",
    "#CB181D",
    "#A50F15",
    "#67000D"
  )
)(100)

# 与参考图相近，将0.2设为图例最低值
# 低于0.2的值按0.2颜色显示，但格子数字仍显示真实相关系数
color_min <- 0.2
color_max <- 1.0

cor_plot <- cor_rho

cor_plot[cor_plot < color_min] <- color_min
cor_plot[cor_plot > color_max] <- color_max

heatmap_breaks <- seq(
  color_min,
  color_max,
  length.out = length(red_palette) + 1
)

legend_breaks <- seq(
  color_min,
  color_max,
  by = 0.1
)

legend_labels <- sprintf(
  "%.1f",
  legend_breaks
)

# 格子中显示真实的Spearman相关系数，而不是截断后的值
number_matrix <- matrix(
  sprintf("%.2f", cor_rho),
  nrow = nrow(cor_rho),
  ncol = ncol(cor_rho),
  dimnames = dimnames(cor_rho)
)

# ============================================================
# 10. 行列使用同一个聚类顺序
#
# 与原图一样使用average linkage
# ============================================================

group_hclust <- hclust(
  dist(cor_plot),
  method = "average"
)

n_group <- ncol(cor_plot)

heatmap_width <- max(
  8,
  n_group * 0.58
)

heatmap_height <- max(
  8,
  n_group * 0.58
)

# 避免生成过大的图片
heatmap_width <- min(
  heatmap_width,
  14
)

heatmap_height <- min(
  heatmap_height,
  14
)

# ============================================================
# 11. 绘制PDF
# ============================================================

pdf_file <- file.path(
  outdir,
  "Fig10_group_Spearman_heatmap_red.pdf"
)

pdf(
  pdf_file,
  width = heatmap_width,
  height = heatmap_height
)

pheatmap(
  cor_plot,

  color = red_palette,
  breaks = heatmap_breaks,

  # 行列使用同一个聚类树
  cluster_rows = group_hclust,
  cluster_cols = group_hclust,

  clustering_method = "average",

  # 保留格子中的真实相关系数
  display_numbers = number_matrix,
  number_color = "grey30",
  fontsize_number = 7,

  # 保留完整矩阵及方格边框
  border_color = "grey75",

  # 保留分组名称
  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,

  angle_col = 90,

  legend_breaks = legend_breaks,
  legend_labels = legend_labels,

  main = "Correlation of motif proportions between groups (Spearman)"
)

dev.off()

# ============================================================
# 12. 绘制PNG
# ============================================================

png_file <- file.path(
  outdir,
  "Fig10_group_Spearman_heatmap_red.png"
)

png(
  png_file,
  width = heatmap_width,
  height = heatmap_height,
  units = "in",
  res = 300,
  bg = "white"
)

pheatmap(
  cor_plot,

  color = red_palette,
  breaks = heatmap_breaks,

  cluster_rows = group_hclust,
  cluster_cols = group_hclust,

  clustering_method = "average",

  display_numbers = number_matrix,
  number_color = "grey30",
  fontsize_number = 7,

  border_color = "grey75",

  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,

  angle_col = 90,

  legend_breaks = legend_breaks,
  legend_labels = legend_labels,

  main = "Correlation of motif proportions between groups (Spearman)"
)

dev.off()

# ============================================================
# 13. 完成信息
# ============================================================

cat("\n图10绘制完成。\n")
cat("PDF：", pdf_file, "\n")
cat("PNG：", png_file, "\n")

cat(
  "相关性矩阵：",
  file.path(
    outdir,
    "Fig10_group_Spearman_correlation_matrix.csv"
  ),
  "\n"
)

cat(
  "分组均值表：",
  file.path(
    outdir,
    "Fig10_motif_mean_by_Cell_Component.tsv"
  ),
  "\n"
)