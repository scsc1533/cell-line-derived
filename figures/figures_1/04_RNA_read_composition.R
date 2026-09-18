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

COUNT_FILE <- file.path(
  PROJ,
  "01_raw_date/length_motif/RNA_totalCounts.txt"
)

INFO_FILE <- file.path(
  PROJ,
  "01_raw_date/information_table/cellcult_Sample_Information.txt"
)

OUTDIR <- "/data/work/sc/redraw_figures_0717/04_fig4_RNA_read_composition"
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

# 柱子从下到上的堆叠顺序
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

# 百分比标签阈值：只显示 >= 5%
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

if (!file.exists(COUNT_FILE)) {
  stop("找不到 RNA reads 文件：", COUNT_FILE)
}

if (!file.exists(INFO_FILE)) {
  stop("找不到样本信息文件：", INFO_FILE)
}

# ============================================================
# 6. 读取数据
# ============================================================

counts <- fread(
  COUNT_FILE,
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
component_col   <- find_column(names(info), c("^Component$"))
qc_col          <- find_column(names(info), c("^QC$"))

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
# 8. 识别 RNA_totalCounts.txt 中的列
# ============================================================

sample_col <- find_column(names(counts), c("^Sample$"))
mrna_col   <- find_column(names(counts), c("^mRNA$"))
lncrna_col <- find_column(names(counts), c("^lncRNA$"))
mirna_col  <- find_column(names(counts), c("^miRNA$"))
trna_col   <- find_column(names(counts), c("^tRNA$"))
pirna_col  <- find_column(names(counts), c("^piRNA$"))
other_col  <- find_column(names(counts), c("^otherRNA$", "^other RNA$"))

required_columns <- c(
  Sample   = sample_col,
  mRNA     = mrna_col,
  lncRNA   = lncrna_col,
  miRNA    = mirna_col,
  tRNA     = trna_col,
  piRNA    = pirna_col,
  otherRNA = other_col
)

if (any(is.na(required_columns))) {
  stop(
    "RNA_totalCounts.txt 缺少以下列：",
    paste(names(required_columns)[is.na(required_columns)], collapse = ", "),
    "\n当前列名：\n",
    paste(names(counts), collapse = "\n")
  )
}

# ============================================================
# 9. 构建样本级 RNA reads 宽表
# ============================================================

count_use <- data.frame(
  Sample = trimws(as.character(counts[[sample_col]])),
  mRNA = suppressWarnings(as.numeric(as.character(counts[[mrna_col]]))),
  lncRNA = suppressWarnings(as.numeric(as.character(counts[[lncrna_col]]))),
  miRNA = suppressWarnings(as.numeric(as.character(counts[[mirna_col]]))),
  tRNA = suppressWarnings(as.numeric(as.character(counts[[trna_col]]))),
  piRNA = suppressWarnings(as.numeric(as.character(counts[[pirna_col]]))),
  `other RNA` = suppressWarnings(as.numeric(as.character(counts[[other_col]]))),
  check.names = FALSE
) %>%
  inner_join(info_use, by = "Sample")

if (nrow(count_use) == 0) {
  stop("RNA_totalCounts.txt 与样本信息表合并后没有有效数据。")
}

cat("最终匹配到的样本数：", n_distinct(count_use$Sample), "\n")

cat("\n各组分样本数：\n")
print(
  count_use %>%
    distinct(Sample, Component) %>%
    count(Component)
)

# ============================================================
# 10. 每个样本内计算 reads 比例
# ============================================================

sample_proportion <- count_use %>%
  pivot_longer(
    cols = all_of(rna_order),
    names_to = "RNA_type",
    values_to = "Reads"
  ) %>%
  mutate(
    Reads = replace_na(Reads, 0)
  ) %>%
  group_by(Sample, Component) %>%
  mutate(
    Total_reads = sum(Reads, na.rm = TRUE),
    Proportion = ifelse(
      Total_reads > 0,
      Reads / Total_reads,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    Component = factor(Component, levels = component_order),
    RNA_type = factor(RNA_type, levels = rna_order)
  )

# ============================================================
# 11. 检查每个样本比例之和
# ============================================================

sample_sum_check <- sample_proportion %>%
  group_by(Sample, Component) %>%
  summarise(
    proportion_sum = sum(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n每个样本比例总和范围：\n")
print(range(sample_sum_check$proportion_sum, na.rm = TRUE))

# ============================================================
# 12. 按组分计算各 RNA 类型平均比例
# ============================================================

summary_table <- sample_proportion %>%
  group_by(Component, RNA_type) %>%
  summarise(
    mean_proportion = mean(Proportion, na.rm = TRUE),
    sd_proportion = sd(Proportion, na.rm = TRUE),
    n = sum(!is.na(Proportion)),
    .groups = "drop"
  )

# ============================================================
# 13. 保存样本级数据
# ============================================================

write.table(
  sample_proportion,
  file.path(OUTDIR, "Fig4_sample_level_read_proportions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_table,
  file.path(OUTDIR, "Fig4_summary_mean_read_proportions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 14. 构建严格 100% 的作图数据
# ============================================================

full_grid_fig4 <- expand_grid(
  Component = component_order,
  RNA_type = stack_order_bottom
)

summary_fig4 <- summary_table %>%
  transmute(
    Component = as.character(Component),
    RNA_type = as.character(RNA_type),
    mean_proportion = as.numeric(mean_proportion)
  )

plot_data_fig4 <- full_grid_fig4 %>%
  left_join(
    summary_fig4,
    by = c("Component", "RNA_type")
  ) %>%
  mutate(
    mean_proportion = replace_na(mean_proportion, 0),
    Component = factor(Component, levels = component_order),
    RNA_type = factor(RNA_type, levels = stack_order_bottom),
    stack_id = match(as.character(RNA_type), stack_order_bottom)
  ) %>%
  group_by(Component) %>%
  mutate(
    component_total = sum(mean_proportion, na.rm = TRUE),
    proportion_100 = ifelse(
      component_total > 0,
      mean_proportion / component_total,
      0
    )
  ) %>%
  arrange(Component, stack_id) %>%
  mutate(
    ymin = cumsum(proportion_100) - proportion_100,
    ymax = cumsum(proportion_100),
    xpos = unname(x_pos_map[as.character(Component)])
  ) %>%
  ungroup()

# ============================================================
# 15. 检查柱高是否为 100%
# ============================================================

fig4_height_check <- plot_data_fig4 %>%
  group_by(Component) %>%
  summarise(
    proportion_sum = sum(proportion_100, na.rm = TRUE),
    bar_top = max(ymax, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n图4每个组分柱高检查，应全部为 1：\n")
print(fig4_height_check)

if (
  any(abs(fig4_height_check$proportion_sum - 1) > 1e-8) ||
  any(abs(fig4_height_check$bar_top - 1) > 1e-8)
) {
  stop("图4作图数据没有正确归一化到 100%。")
}

write.table(
  fig4_height_check,
  file.path(OUTDIR, "Fig4_plot_height_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 16. 构建柱子坐标
# ============================================================

bar_data_fig4 <- plot_data_fig4 %>%
  mutate(
    xmin = xpos - bar_width / 2,
    xmax = xpos + bar_width / 2
  )

# ============================================================
# 16.1 构建百分比标签数据（>=5% 才显示）
# ============================================================

label_data_fig4 <- bar_data_fig4 %>%
  filter(
    !is.na(proportion_100),
    proportion_100 >= label_cutoff
  ) %>%
  mutate(
    y_label = (ymin + ymax) / 2,
    percent_label = sprintf("%.0f%%", proportion_100 * 100),
    text_color = case_when(
      as.character(RNA_type) %in% c(
        "mRNA", "miRNA", "tRNA", "other RNA"
      ) ~ "white",
      TRUE ~ "black"
    )
  )

cat("\n将显示的百分比标签：\n")
print(
  label_data_fig4 %>%
    select(Component, RNA_type, proportion_100, percent_label)
)

cat("\n标签数量：", nrow(label_data_fig4), "\n")

# ============================================================
# 17. 构建直线透明连接带
# ============================================================

make_straight_band_fig4 <- function(
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

band_list_fig4 <- list()
band_index_fig4 <- 1

for (component_index in seq_len(length(component_order) - 1)) {

  left_component <- component_order[component_index]
  right_component <- component_order[component_index + 1]

  left_center <- unname(x_pos_map[left_component])
  right_center <- unname(x_pos_map[right_component])

  x_left <- left_center + bar_width / 2
  x_right <- right_center - bar_width / 2

  for (rna_name in stack_order_bottom) {

    left_row <- plot_data_fig4 %>%
      filter(Component == left_component, RNA_type == rna_name)

    right_row <- plot_data_fig4 %>%
      filter(Component == right_component, RNA_type == rna_name)

    band_list_fig4[[band_index_fig4]] <- make_straight_band_fig4(
      left_row = left_row,
      right_row = right_row,
      x_left = x_left,
      x_right = x_right,
      rna_name = rna_name,
      band_name = paste(left_component, right_component, rna_name, sep = "__")
    )

    band_index_fig4 <- band_index_fig4 + 1
  }
}

band_data_fig4 <- bind_rows(band_list_fig4) %>%
  mutate(
    RNA_type = factor(RNA_type, levels = stack_order_bottom)
  )

cat("\n图4直线连接带数据行数，应为 48：", nrow(band_data_fig4), "\n")
stopifnot(nrow(band_data_fig4) == 48)

# ============================================================
# 18. 保存最终作图数据
# ============================================================

write.table(
  plot_data_fig4,
  file.path(OUTDIR, "Fig4_plot_data_strict_100_percent_compact.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 19. 绘图
# ============================================================

p <- ggplot() +

  # 连接带
  geom_polygon(
    data = band_data_fig4,
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
    data = bar_data_fig4,
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

  # 百分比标签（仅 >=5%）
  geom_text(
    data = label_data_fig4,
    aes(
      x = xpos,
      y = y_label,
      label = percent_label,
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
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +

  labs(
    title = "RNA read composition across components",
    x = NULL,
    y = "Mean proportion of RNA reads",
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
# 20. 保存图片
# ============================================================

pdf_file <- file.path(
  OUTDIR,
  "Fig4_RNA_read_composition_narrow_close_straight_bands_label5pct.pdf"
)

png_file <- file.path(
  OUTDIR,
  "Fig4_RNA_read_composition_narrow_close_straight_bands_label5pct.png"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 5.4,
  height = 5.6,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 5.4,
  height = 5.6,
  dpi = 300,
  bg = "white"
)

cat("\n图4绘制完成。\n")
cat("PDF：", pdf_file, "\n")
cat("PNG：", png_file, "\n")
