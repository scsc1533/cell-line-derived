#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ============================================================
# 1. 输入和输出路径
# ============================================================

ROOT <- "/data/users/shichu/shichu_02c8549606a945a79d8e5e9c1bf9cb1c/online"
PROJ <- file.path(ROOT, "01_2603cell_culture")

GENE_FILE <- file.path(
  PROJ,
  "01_raw_date/length_motif/genenum_all.txt"
)

INFO_FILE <- file.path(
  PROJ,
  "01_raw_date/information_table/cellcult_Sample_Information.txt"
)

OUTDIR <- "/data/work/01_2603cell_culture/07_figure/01_figure1/03_detected_RNA_types"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

QC_KEEP <- c(1, 2)

# ============================================================
# 2. 配色和顺序
# ============================================================

rna_colors <- c(
  "mRNA" = "#E64B35",
  "lncRNA" = "#4DBBD5",
  "miRNA" = "#00A087",
  "tRNA" = "#3C5488",
  "piRNA" = "#F39B7F",
  "other RNA" = "#8491B4"
)

rna_order <- c(
  "mRNA",
  "lncRNA",
  "miRNA",
  "tRNA",
  "piRNA",
  "other RNA"
)

component_order <- c(
  "Cell",
  "Debris",
  "Supernatant"
)

# 柱子从下到上的顺序
stack_order_bottom <- c(
  "other RNA",
  "piRNA",
  "tRNA",
  "miRNA",
  "lncRNA",
  "mRNA"
)

# 图例顺序
legend_order <- c(
  "mRNA",
  "lncRNA",
  "miRNA",
  "tRNA",
  "piRNA",
  "other RNA"
)

# ============================================================
# 3. 图形参数
# ============================================================

# 柱子更窄
bar_width <- 0.30

# 连接带透明度
band_alpha <- 0.16

# 三根柱子更靠近
x_pos_map <- c(
  "Cell" = 1.00,
  "Debris" = 1.40,
  "Supernatant" = 1.80
)

# 只显示 >= 5% 的标签
label_cutoff <- 0.05

# ============================================================
# 4. 辅助函数
# ============================================================

find_column <- function(x, patterns) {
  for (pattern in patterns) {
    hit <- grep(pattern, x, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

map_component <- function(x) {
  x2 <- tolower(trimws(as.character(x)))

  case_when(
    x2 == "cell" ~ "Cell",
    grepl("debris", x2) ~ "Debris",
    grepl("cfrna|supernatant", x2) ~ "Supernatant",
    TRUE ~ NA_character_
  )
}

# ============================================================
# 5. 检查输入文件
# ============================================================

if (!file.exists(GENE_FILE)) {
  stop("找不到：", GENE_FILE)
}

if (!file.exists(INFO_FILE)) {
  stop("找不到：", INFO_FILE)
}

# ============================================================
# 6. 读取数据
# ============================================================

gene <- fread(
  GENE_FILE,
  data.table = FALSE,
  check.names = FALSE
)

info <- fread(
  INFO_FILE,
  data.table = FALSE,
  check.names = FALSE
)

# ============================================================
# 7. 整理样本信息
# ============================================================

sample_info_col <- find_column(names(info), c("^Sample$"))
component_col <- find_column(names(info), c("^Component$"))
qc_col <- find_column(names(info), c("^QC$"))

if (is.na(sample_info_col) || is.na(component_col)) {
  stop("样本信息表中缺少 Sample 或 Component 列。")
}

info_use <- info %>%
  transmute(
    Sample = trimws(as.character(.data[[sample_info_col]])),
    Component = map_component(.data[[component_col]]),
    QC = if (!is.na(qc_col)) {
      suppressWarnings(as.numeric(as.character(.data[[qc_col]])))
    } else {
      1
    }
  ) %>%
  filter(
    QC %in% QC_KEEP,
    Component %in% component_order,
    !is.na(Sample),
    Sample != ""
  ) %>%
  distinct(Sample, .keep_all = TRUE)

cat("样本信息筛选后的样本数：", nrow(info_use), "\n")

# ============================================================
# 8. 识别 genenum_all.txt 中的列
# ============================================================

gene_sample_col <- find_column(names(gene), c("^Sample$"))
mrna_col <- find_column(names(gene), c("^mRNA\\(TPM>0\\)$"))
lncrna_col <- find_column(names(gene), c("^lncRNA\\(TPM>0\\)$"))
mirna_col <- find_column(names(gene), c("^miRNA\\(RPM>0\\)$"))
trna_col <- find_column(names(gene), c("^tRNA\\(RPM>0\\)$"))
pirna_col <- find_column(names(gene), c("^piRNA\\(RPM>0\\)$"))
other_col <- find_column(names(gene), c("^otherRNA$", "^other RNA$"))

required <- c(
  gene_sample_col,
  mrna_col,
  lncrna_col,
  mirna_col,
  trna_col,
  pirna_col
)

if (any(is.na(required))) {
  stop(
    "genenum_all.txt 缺少必要列。当前列名：\n",
    paste(names(gene), collapse = "\n")
  )
}

# 处理 other RNA
if (!is.na(other_col)) {
  other_value <- suppressWarnings(
    as.numeric(as.character(gene[[other_col]]))
  )
} else {
  circ_cols <- grep(
    "^circRNA.*\\(.*>0\\)",
    names(gene),
    ignore.case = TRUE,
    value = TRUE
  )

  microbe_cols <- grep(
    "^Microbe.*\\(.*>0\\)",
    names(gene),
    ignore.case = TRUE,
    value = TRUE
  )

  other_cols <- unique(c(circ_cols, microbe_cols))

  if (length(other_cols) == 0) {
    other_value <- rep(0, nrow(gene))
  } else {
    other_matrix <- as.data.frame(
      lapply(gene[, other_cols, drop = FALSE], function(x) {
        suppressWarnings(as.numeric(as.character(x)))
      })
    )

    other_value <- rowSums(
      other_matrix,
      na.rm = TRUE
    )
  }
}

# ============================================================
# 9. 构建样本级检出数宽表
# ============================================================

gene_use <- data.frame(
  Sample = trimws(as.character(gene[[gene_sample_col]])),
  mRNA = suppressWarnings(as.numeric(as.character(gene[[mrna_col]]))),
  lncRNA = suppressWarnings(as.numeric(as.character(gene[[lncrna_col]]))),
  miRNA = suppressWarnings(as.numeric(as.character(gene[[mirna_col]]))),
  tRNA = suppressWarnings(as.numeric(as.character(gene[[trna_col]]))),
  piRNA = suppressWarnings(as.numeric(as.character(gene[[pirna_col]]))),
  `other RNA` = other_value,
  check.names = FALSE
) %>%
  inner_join(info_use, by = "Sample")

if (nrow(gene_use) == 0) {
  stop("genenum_all.txt 与样本信息表合并后没有有效数据。")
}

cat("最终匹配到的样本数：", n_distinct(gene_use$Sample), "\n")

cat("\n各组分样本数：\n")
print(
  gene_use %>%
    distinct(Sample, Component) %>%
    count(Component)
)

# ============================================================
# 10. 样本级长表（保留原始检出数）
# ============================================================

gene_long <- gene_use %>%
  pivot_longer(
    cols = all_of(rna_order),
    names_to = "RNA_type",
    values_to = "Detected_number"
  ) %>%
  mutate(
    Detected_number = replace_na(
      suppressWarnings(as.numeric(Detected_number)),
      0
    )
  ) %>%
  mutate(
    Component = factor(Component, levels = component_order),
    RNA_type = factor(RNA_type, levels = rna_order)
  )

# ============================================================
# 11. 检查每个样本检出总数
# ============================================================

sample_sum_check <- gene_long %>%
  group_by(Sample, Component) %>%
  summarise(
    total_detected = sum(Detected_number, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n每个样本检出总数范围：\n")
print(range(sample_sum_check$total_detected, na.rm = TRUE))

# ============================================================
# 12. 按组分计算每类RNA检出数的平均值
# ============================================================

summary_table <- gene_long %>%
  group_by(Component, RNA_type) %>%
  summarise(
    mean_count = mean(Detected_number, na.rm = TRUE),
    sd_count = sd(Detected_number, na.rm = TRUE),
    n = sum(!is.na(Detected_number)),
    .groups = "drop"
  )

# ============================================================
# 13. 检查每个组分平均检出总数
# ============================================================

count_sum_qc <- summary_table %>%
  group_by(Component) %>%
  summarise(
    sum_mean_count = sum(mean_count, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n各组分平均检出总数：\n")
print(count_sum_qc)

# ============================================================
# 14. 保存样本级数据和汇总数据
# ============================================================

write.table(
  gene_long,
  file.path(OUTDIR, "Fig3_sample_level_detected_RNA_ratio.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_table,
  file.path(OUTDIR, "Fig3_summary_mean_ratio.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  count_sum_qc,
  file.path(OUTDIR, "Fig3_ratio_sum_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 15. 构建作图数据（按平均检出数堆叠）
# ============================================================

full_grid_fig3 <- expand_grid(
  Component = component_order,
  RNA_type = stack_order_bottom
)

summary_fig3 <- summary_table %>%
  transmute(
    Component = as.character(Component),
    RNA_type = as.character(RNA_type),
    mean_count = as.numeric(mean_count)
  )

plot_data_fig3 <- full_grid_fig3 %>%
  left_join(
    summary_fig3,
    by = c("Component", "RNA_type")
  ) %>%
  mutate(
    mean_count = replace_na(mean_count, 0),
    Component = factor(Component, levels = component_order),
    RNA_type = factor(RNA_type, levels = stack_order_bottom),
    stack_id = match(as.character(RNA_type), stack_order_bottom)
  ) %>%
  group_by(Component) %>%
  arrange(stack_id, .by_group = TRUE) %>%
  mutate(
    ymin = cumsum(mean_count) - mean_count,
    ymax = cumsum(mean_count),
    xpos = unname(x_pos_map[as.character(Component)])
  ) %>%
  ungroup()

# ============================================================
# 16. 检查柱高
# ============================================================

fig3_height_check <- plot_data_fig3 %>%
  group_by(Component) %>%
  summarise(
    sum_mean_count = sum(mean_count, na.rm = TRUE),
    bar_top = max(ymax, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n图3每个组分柱高检查（bar_top 应等于 sum_mean_count）：\n")
print(fig3_height_check)

if (any(abs(fig3_height_check$bar_top - fig3_height_check$sum_mean_count) > 1e-8)) {
  stop("图3作图数据柱高不正确。")
}

write.table(
  fig3_height_check,
  file.path(OUTDIR, "Fig3_plot_height_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 17. 构建柱子坐标
# ============================================================

bar_data_fig3 <- plot_data_fig3 %>%
  mutate(
    xmin = xpos - bar_width / 2,
    xmax = xpos + bar_width / 2
  )

# ============================================================
# 17.1 构建标签数据（占比 >=5% 的区段显示计数）
# ============================================================

label_data_fig3 <- bar_data_fig3 %>%
  group_by(Component) %>%
  mutate(
    bar_total = sum(mean_count, na.rm = TRUE),
    seg_frac = ifelse(bar_total > 0, mean_count / bar_total, 0)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(seg_frac),
    seg_frac >= label_cutoff
  ) %>%
  mutate(
    y_label = (ymin + ymax) / 2,
    count_label = as.character(round(mean_count, 0)),
    text_color = case_when(
      as.character(RNA_type) %in% c(
        "mRNA", "miRNA", "tRNA", "other RNA"
      ) ~ "white",
      TRUE ~ "black"
    )
  )

cat("\n将显示的计数标签：\n")
print(
  label_data_fig3 %>%
    select(Component, RNA_type, mean_count, count_label)
)

cat("\n标签数量：", nrow(label_data_fig3), "\n")

# ============================================================
# 18. 构建直线透明连接带
# ============================================================

make_straight_band_fig3 <- function(
  left_row,
  right_row,
  x_left,
  x_right,
  rna_name,
  band_name
) {
  data.frame(
    x = c(x_left, x_right, x_right, x_left),
    y = c(
      left_row$ymax[[1]],
      right_row$ymax[[1]],
      right_row$ymin[[1]],
      left_row$ymin[[1]]
    ),
    RNA_type = rna_name,
    band_id = band_name,
    stringsAsFactors = FALSE
  )
}

band_list_fig3 <- list()
band_index_fig3 <- 1

for (component_index in seq_len(length(component_order) - 1)) {

  left_component <- component_order[component_index]
  right_component <- component_order[component_index + 1]

  left_center <- unname(x_pos_map[left_component])
  right_center <- unname(x_pos_map[right_component])

  x_left <- left_center + bar_width / 2
  x_right <- right_center - bar_width / 2

  for (rna_name in stack_order_bottom) {

    left_row <- plot_data_fig3 %>%
      filter(Component == left_component, RNA_type == rna_name)

    right_row <- plot_data_fig3 %>%
      filter(Component == right_component, RNA_type == rna_name)

    band_list_fig3[[band_index_fig3]] <- make_straight_band_fig3(
      left_row = left_row,
      right_row = right_row,
      x_left = x_left,
      x_right = x_right,
      rna_name = rna_name,
      band_name = paste(left_component, right_component, rna_name, sep = "__")
    )

    band_index_fig3 <- band_index_fig3 + 1
  }
}

band_data_fig3 <- bind_rows(band_list_fig3) %>%
  mutate(
    RNA_type = factor(RNA_type, levels = stack_order_bottom)
  )

cat("\n图3直线连接带数据行数，应为48：", nrow(band_data_fig3), "\n")
stopifnot(nrow(band_data_fig3) == 48)

# ============================================================
# 19. 保存最终作图数据
# ============================================================

write.table(
  plot_data_fig3,
  file.path(OUTDIR, "Fig3_plot_data_strict_100_percent_compact.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 20. 绘图
# ============================================================

p_fig3 <- ggplot() +

  # 连接带
  geom_polygon(
    data = band_data_fig3,
    aes(
      x = x,
      y = y,
      group = band_id,
      fill = RNA_type
    ),
    alpha = band_alpha,
    color = NA
  ) +

  # 柱子
  geom_rect(
    data = bar_data_fig3,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = RNA_type
    ),
    color = "white",
    linewidth = 0.28
  ) +

  # 计数标签（仅占比 >=5% 的区段）
  geom_text(
    data = label_data_fig3,
    aes(
      x = xpos,
      y = y_label,
      label = count_label,
      color = text_color
    ),
    inherit.aes = FALSE,
    size = 2.7,
    fontface = "bold",
    show.legend = FALSE
  ) +

  scale_color_identity() +

  scale_fill_manual(
    values = rna_colors,
    breaks = legend_order,
    drop = FALSE
  ) +

  scale_x_continuous(
    breaks = unname(x_pos_map),
    labels = names(x_pos_map),
    limits = c(0.82, 2.02),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +

  labs(
    title = "Composition of detected RNA types across components",
    x = NULL,
    y = "Mean number of detected RNAs",
    fill = NULL
  ) +

  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "right",
    plot.margin = margin(
      t = 8, r = 8, b = 8, l = 8
    )
  )

# ============================================================
# 21. 保存图片
# ============================================================

pdf_file <- file.path(
  OUTDIR,
  "Fig3_detected_RNA_type_narrow_close_straight_bands_label5pct.pdf"
)

png_file <- file.path(
  OUTDIR,
  "Fig3_detected_RNA_type_narrow_close_straight_bands_label5pct.png"
)

ggsave(
  filename = pdf_file,
  plot = p_fig3,
  width = 5.4,
  height = 5.6,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = p_fig3,
  width = 5.4,
  height = 5.6,
  dpi = 300,
  bg = "white"
)

cat("\n图3绘制完成。\n")
cat("PDF：", pdf_file, "\n")
cat("PNG：", png_file, "\n")
