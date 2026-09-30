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
# 1. Input and output paths
# ============================================================

qc_file <- paste0(
  "./cellcult_Sample_Information.txt"
)

motif_file <- paste0(
  "./sample_motif_ratios_summary.txt"
)

outdir <- "./redraw_figures_0717/09_fig9_freqTC"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(qc_file)) {
  stop("Sample information file not found: ", qc_file)
}

if (!file.exists(motif_file)) {
  stop("freqTC input file not found: ", motif_file)
}

# ============================================================
# 2. Colors
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
# 3. Read data
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

cat("Sample information table dimensions: ", nrow(qc_data), " x ", ncol(qc_data), "\n")
cat("freqTC data table dimensions: ", nrow(motif_data), " x ", ncol(motif_data), "\n")

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
    "Sample information table is missing the following columns: ",
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
    "freqTC data table is missing the following columns: ",
    paste(missing_motif_columns, collapse = ", ")
  )
}

# ============================================================
# 4. Sample filtering and component name cleaning
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
# 5. Clean freqTC data
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
# 6. Merge data
#
# If a sample appears multiple times in the motif file,
# average freqTC for that sample.
# Ultimately one sample corresponds to one point.
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
    "No valid data after merging sample information table with freqTC data. ",
    "Please check whether Sample and sample_name are consistent."
  )
}

if (n_distinct(plot_df$Component) < 2) {
  stop("Fewer than two valid components, cannot perform between-group comparison.")
}

cat("\nFinal data overview used for Figure 9:\n")
print(
  table(
    plot_df$Component,
    useNA = "ifany"
  )
)

cat("\nMaximum number of occurrences per sample:\n")
print(
  plot_df %>%
    count(Sample) %>%
    summarise(max_rows_per_sample = max(n))
)

# ============================================================
# 7. Save sample-level plotting data
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
# 8. Pairwise Wilcoxon tests
#
# Use BH method to adjust the three pairwise comparisons.
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

# Format P values into the text displayed in the plot
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

cat("\nPairwise comparison results:\n")
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
# 9. Set P-value bracket positions
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
# 10. Plotting: customize the actual positions of the three boxplots
# ============================================================

# Actually control the distance between the three boxes
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

# Convert P-value bracket group names to custom x coordinates
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

# Length of the small vertical ticks at both ends of the P-value bracket
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

  # Transparent boxplot
  geom_boxplot(
    width = 0.22,
    fill = NA,
    outlier.shape = NA,
    linewidth = 0.75
  ) +

  # One point represents one sample
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

  # Use custom positions to display group names
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
# 11. Manually add P-value brackets
# Because the x-axis has been changed to custom numeric positions
# ============================================================

if (nrow(stat_plot) > 0) {

  p <- p +

    # Bracket horizontal line
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

    # Left small vertical tick
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

    # Right small vertical tick
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

    # P-value text
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
# 12. Save
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

cat("\nCompact transparent boxplot completed.\n")
cat("PDF: ", pdf_file, "\n")
cat("PNG: ", png_file, "\n")