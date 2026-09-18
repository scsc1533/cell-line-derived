#!/usr/bin/env Rscript

# ============================================================
# Figure 9
# Comparison of freqTC among Cell, Debris and Supernatant
# All cell lines are combined
# One point represents one sample
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggsignif)
})

# ============================================================
# 1. 输入和输出路径
# ============================================================

qc_file <- paste0(
  "/data/work/1/cellcult_Sample_Information.txt"
)

motif_file <- paste0(
  "/data/work/sc/sample_motif_ratios_summary.txt"
)

outdir <- "/data/work/sc/redraw_figures_0717/09_fig9_freqTC"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(qc_file)) {
  stop("找不到样本信息文件：", qc_file)
}

if (!file.exists(motif_file)) {
  stop("找不到 freqTC 输入文件：", motif_file)
}

# ============================================================
# 2. 配色
# ============================================================

component_colors <- c(
  "Cell" = "#4DBBD5",
  "Debris" = "#00A087",
  "Supernatant" = "#E64B35"
)

component_order <- c(
  "Cell",
  "Debris",
  "Supernatant"
)

# ============================================================
# 3. 读取数据
# ============================================================

qc_data <- read.delim(
  qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

motif_data <- read.delim(
  motif_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("样本信息表维度：", nrow(qc_data), "×", ncol(qc_data), "\n")
cat("freqTC数据表维度：", nrow(motif_data), "×", ncol(motif_data), "\n")

required_qc_columns <- c(
  "Sample",
  "QC",
  "Component"
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

required_motif_columns <- c(
  "sample_name",
  "freqTC"
)

missing_motif_columns <- setdiff(
  required_motif_columns,
  colnames(motif_data)
)

if (length(missing_motif_columns) > 0) {
  stop(
    "freqTC数据表缺少以下列：",
    paste(missing_motif_columns, collapse = ", ")
  )
}

# ============================================================
# 4. 样本筛选与组分名称整理
# ============================================================

valid_qc <- qc_data %>%
  mutate(
    Sample = trimws(as.character(Sample)),
    QC = suppressWarnings(
      as.numeric(as.character(QC))
    ),
    Component = trimws(as.character(Component))
  ) %>%
  filter(
    QC == 1,
    Component %in% c(
      "cell",
      "debris",
      "cfRNA"
    )
  ) %>%
  transmute(
    Sample,
    Component = case_when(
      Component == "cell" ~ "Cell",
      Component == "debris" ~ "Debris",
      Component == "cfRNA" ~ "Supernatant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Sample),
    Sample != "",
    !is.na(Component)
  ) %>%
  distinct(
    Sample,
    .keep_all = TRUE
  )

# ============================================================
# 5. 整理freqTC数据
# ============================================================

motif_use <- motif_data %>%
  transmute(
    Sample = trimws(
      as.character(sample_name)
    ),
    freqTC = suppressWarnings(
      as.numeric(as.character(freqTC))
    )
  ) %>%
  filter(
    !is.na(Sample),
    Sample != "",
    is.finite(freqTC)
  )

# ============================================================
# 6. 合并数据
#
# 如果一个样本在motif文件中出现多次，
# 对这个样本的freqTC取平均。
# 最终一个样本只对应一个点。
# ============================================================

plot_df <- valid_qc %>%
  inner_join(
    motif_use,
    by = "Sample"
  ) %>%
  group_by(
    Sample,
    Component
  ) %>%
  summarise(
    freqTC = mean(
      freqTC,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Component = factor(
      Component,
      levels = component_order
    )
  ) %>%
  filter(
    !is.na(Component),
    is.finite(freqTC)
  )

if (nrow(plot_df) == 0) {
  stop(
    "样本信息表与freqTC数据合并后没有有效数据，",
    "请检查Sample与sample_name是否一致。"
  )
}

if (n_distinct(plot_df$Component) < 2) {
  stop("有效组分少于两个，无法进行组间比较。")
}

cat("\n最终用于图9的数据概况：\n")
print(
  table(
    plot_df$Component,
    useNA = "ifany"
  )
)

cat("\n每个样本最多出现次数：\n")
print(
  plot_df %>%
    count(Sample) %>%
    summarise(max_rows_per_sample = max(n))
)

# ============================================================
# 7. 保存样本级绘图数据
# ============================================================

write.table(
  plot_df,
  file = file.path(
    outdir,
    "Fig9_freqTC_sample_level_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_count <- plot_df %>%
  count(
    Component,
    name = "sample_n"
  )

write.table(
  sample_count,
  file = file.path(
    outdir,
    "Fig9_sample_count_by_component.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 8. 两两Wilcoxon检验
#
# 使用BH方法校正三次两两比较。
# ============================================================

stat_res <- plot_df %>%
  pairwise_wilcox_test(
    freqTC ~ Component,
    p.adjust.method = "BH",
    exact = FALSE
  ) %>%
  add_significance(
    "p.adj"
  )

# 将P值格式化成图中显示的文本
format_adjusted_p <- function(p) {

  if (is.na(p)) {
    return("BH-adjusted P = NA")
  }

  if (p < 0.0001) {
    return("BH-adjusted P < 0.0001")
  }

  if (p < 0.001) {
    return(
      paste0(
        "BH-adjusted P = ",
        formatC(
          p,
          format = "e",
          digits = 2
        )
      )
    )
  }

  paste0(
    "BH-adjusted P = ",
    formatC(
      p,
      format = "f",
      digits = 3
    )
  )
}

stat_res <- stat_res %>%
  mutate(
    p_label = vapply(
      p.adj,
      format_adjusted_p,
      character(1)
    )
  )

write.csv(
  stat_res,
  file = file.path(
    outdir,
    "Fig9_pairwise_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)

cat("\n两两比较结果：\n")
print(
  stat_res %>%
    select(
      group1,
      group2,
      n1,
      n2,
      statistic,
      p,
      p.adj,
      p.adj.signif
    )
)

# ============================================================
# 9. 设置P值横线位置
# ============================================================

plot_min <- min(
  plot_df$freqTC,
  na.rm = TRUE
)

plot_max <- max(
  plot_df$freqTC,
  na.rm = TRUE
)

plot_span <- plot_max - plot_min

if (!is.finite(plot_span) || plot_span == 0) {
  plot_span <- max(
    abs(plot_max),
    0.1
  )
}

annotation_step <- plot_span * 0.15

stat_res <- stat_res %>%
  filter(
    !is.na(p.adj)
  ) %>%
  mutate(
    y_position = plot_max +
      annotation_step *
      seq_len(n())
  )

comparisons_list <- purrr::map2(
  stat_res$group1,
  stat_res$group2,
  c
)

# ============================================================
# 10. 绘图：自定义三个箱线图的实际位置
# ============================================================

# 真正控制三个箱子之间的距离
x_map <- c(
  "Cell" = 1.00,
  "Debris" = 1.55,
  "Supernatant" = 2.10
)

plot_df_positioned <- plot_df %>%
  mutate(
    x_pos = unname(
      x_map[as.character(Component)]
    )
  )

# 将P值横线的组名转换成自定义横坐标
stat_plot <- stat_res %>%
  mutate(
    xmin = unname(
      x_map[as.character(group1)]
    ),
    xmax = unname(
      x_map[as.character(group2)]
    ),
    x_text = (xmin + xmax) / 2
  ) %>%
  filter(
    !is.na(xmin),
    !is.na(xmax),
    !is.na(y_position)
  )

# P值括号两端小竖线的长度
bracket_tip <- annotation_step * 0.10

set.seed(123)

p <- ggplot(
  plot_df_positioned,
  aes(
    x = x_pos,
    y = freqTC,
    group = Component,
    color = Component
  )
) +

  # 透明箱线图
  geom_boxplot(
    width = 0.22,
    fill = NA,
    outlier.shape = NA,
    linewidth = 0.75
  ) +

  # 一个点代表一个样本
  geom_point(
    position = position_jitter(
      width = 0.035,
      height = 0
    ),
    size = 1.7,
    alpha = 0.62
  ) +

  scale_color_manual(
    values = component_colors,
    drop = FALSE
  ) +

  # 使用自定义位置显示组名
  scale_x_continuous(
    breaks = unname(x_map),
    labels = names(x_map),
    limits = c(0.78, 2.32),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +

  labs(
    title = NULL,
    x = NULL,
    y = "freqTC"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid.major = element_line(
      color = "#E8E8E8",
      linewidth = 0.45
    ),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      color = "grey35",
      fill = NA,
      linewidth = 0.7
    ),

    axis.title.y = element_text(
      size = 10.5,
      color = "black",
      margin = margin(r = 6)
    ),

    axis.text.x = element_text(
      size = 8.8,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),

    axis.ticks = element_line(
      color = "grey35",
      linewidth = 0.5
    ),

    legend.position = "none",

    plot.margin = margin(
      t = 15,
      r = 12,
      b = 7,
      l = 12
    )
  ) +

  coord_cartesian(
    ylim = c(
      plot_min,
      if (nrow(stat_plot) > 0) {
        max(stat_plot$y_position) +
          annotation_step * 0.55
      } else {
        plot_max + annotation_step
      }
    ),
    clip = "off"
  )

# ============================================================
# 11. 手动添加P值括号
# 因为横坐标已经改成自定义数值位置
# ============================================================

if (nrow(stat_plot) > 0) {

  p <- p +

    # 括号横线
    geom_segment(
      data = stat_plot,
      aes(
        x = xmin,
        xend = xmax,
        y = y_position,
        yend = y_position
      ),
      inherit.aes = FALSE,
      color = "grey25",
      linewidth = 0.45
    ) +

    # 左侧小竖线
    geom_segment(
      data = stat_plot,
      aes(
        x = xmin,
        xend = xmin,
        y = y_position,
        yend = y_position - bracket_tip
      ),
      inherit.aes = FALSE,
      color = "grey25",
      linewidth = 0.45
    ) +

    # 右侧小竖线
    geom_segment(
      data = stat_plot,
      aes(
        x = xmax,
        xend = xmax,
        y = y_position,
        yend = y_position - bracket_tip
      ),
      inherit.aes = FALSE,
      color = "grey25",
      linewidth = 0.45
    ) +

    # P值文字
    geom_text(
      data = stat_plot,
      aes(
        x = x_text,
        y = y_position +
          annotation_step * 0.06,
        label = p_label
      ),
      inherit.aes = FALSE,
      color = "grey25",
      size = 3.0,
      vjust = 0
    )
}

# ============================================================
# 12. 保存
# ============================================================

pdf_file <- file.path(
  outdir,
  "Fig9_freqTC_transparent_boxplot_compact.pdf"
)

png_file <- file.path(
  outdir,
  "Fig9_freqTC_transparent_boxplot_compact.png"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 4.3,
  height = 5.2,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 4.3,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

cat("\n紧凑版透明箱线图绘制完成。\n")
cat("PDF：", pdf_file, "\n")
cat("PNG：", png_file, "\n")