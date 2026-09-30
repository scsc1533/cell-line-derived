#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

###############################################################################
# Butterfly plot for gene-level D_target direction
#
# Input:
#   gene_level_distances.tsv.gz
#   gene_level_group_tests.tsv.gz
#
# Output:
#   1. All genes butterfly plot
#   2. Significant genes (wilcox_p_D_target < 0.05) butterfly plot
#
# Direction:
#   gdm_minus_healthy_D_target < 0
#       -> GDM closer to target reference
#
#   gdm_minus_healthy_D_target > 0
#       -> GDM farther from target reference
###############################################################################


# ======================== 1. Paths ========================

validation_dir <- "/data/work/01_2603cell_culture/07_figure/04_figure4/02_HBV_validation_SH/all"

distance_file <- file.path(
  validation_dir,
  "gene_level_distances.tsv.gz"
)

gene_test_file <- file.path(
  validation_dir,
  "gene_level_group_tests.tsv.gz"
)

output_dir <- file.path(
  validation_dir,
  "../../05_butterfly_plot/HBV_SH/"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================== 2. Parameters ========================

comparison_order <- c(
  "K562_vs_rest",
  "HTR_8_SVneo_vs_rest",
  "HepG2_Hep3B2.1_7_vs_rest"
)


comparison_labels <- c(
  "K562_vs_rest" = "K562",
  "HTR_8_SVneo_vs_rest" = "HTR-8/SVneo",
  "HepG2_Hep3B2.1_7_vs_rest" = "Liver (HepG2+Hep3B2.1-7)"
)


plot_colors <- c(
  "Closer to target reference" = "#4DBBD5FF",
  "Farther from target reference" = "#E64B35FF"
)


p_cutoff <- 0.05



# ======================== 3. Read data ========================

message("[1/5] Reading distance matrix ...")

gene_distance_dt <- fread(
  distance_file
)

required_cols <- c(
  "comparison",
  "Gene",
  "D_target"
)

missing_cols <- setdiff(
  required_cols,
  names(gene_distance_dt)
)

if (length(missing_cols) > 0) {
  stop(
    "gene_level_distances missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ======================== 4. Calculate gene direction ========================

message("[2/5] Calculating gene direction ...")


gene_direction_dt <- gene_distance_dt[
  ,
  .(
    healthy_median_D_target =
      median(D_target[group == 0], na.rm = TRUE),

    gdm_median_D_target =
      median(D_target[group == 1], na.rm = TRUE),

    n_healthy =
      sum(group == 0 & !is.na(D_target)),

    n_gdm =
      sum(group == 1 & !is.na(D_target))
  ),
  by = .(
    comparison,
    Gene
  )
]


gene_direction_dt[,
  gdm_minus_healthy_D_target :=
    gdm_median_D_target -
    healthy_median_D_target
]


gene_direction_dt[,
  Direction :=
    fifelse(
      gdm_minus_healthy_D_target < 0,
      "Closer to target reference",
      "Farther from target reference"
    )
]


gene_direction_dt[
  ,
  comparison := factor(
    comparison,
    levels = comparison_order
  )
]


fwrite(
  gene_direction_dt,
  file.path(
    output_dir,
    "gene_direction_all.tsv"
  ),
  sep = "\t"
)



# ======================== 5. Add significance information ========================

message("[3/5] Adding gene significance ...")


gene_test_dt <- fread(
  gene_test_file
)


sig_dt <- gene_test_dt[
  comparison %in% comparison_order,
  .(
    comparison,
    Gene,
    wilcox_p_D_target,
    wilcox_fdr_D_target
  )
]


gene_direction_sig_dt <- merge(
  gene_direction_dt,
  sig_dt,
  by = c(
    "comparison",
    "Gene"
  ),
  all.x = TRUE
)


fwrite(
  gene_direction_sig_dt,
  file.path(
    output_dir,
    "gene_direction_with_statistics.tsv"
  ),
  sep = "\t"
)



# ======================== 6. Butterfly plot function ========================

draw_butterfly <- function(
    input_dt,
    title_text,
    output_prefix
) {


  plot_dt <- copy(input_dt)


  plot_dt <- plot_dt[
    !is.na(Direction)
  ]


  count_dt <- plot_dt[
    ,
    .(
      Count = .N
    ),
    by = .(
      comparison,
      Direction
    )
  ]


  # 保证三个comparison和两个方向都存在
  complete_dt <- CJ(
    comparison = factor(
      comparison_order,
      levels = comparison_order
    ),
    Direction = names(plot_colors)
  )


  count_dt <- merge(
    complete_dt,
    count_dt,
    by = c(
      "comparison",
      "Direction"
    ),
    all.x = TRUE
  )


  count_dt[
    is.na(Count),
    Count := 0
  ]


  # 左侧显示 closer
  # 右侧显示 farther

  count_dt[
    Direction ==
      "Closer to target reference",
    Count_plot := -Count
  ]


  count_dt[
    Direction ==
      "Farther from target reference",
    Count_plot := Count
  ]


  count_dt[
    ,
    comparison_label :=
      comparison_labels[
        as.character(comparison)
      ]
  ]


  count_dt[
    ,
    comparison_label :=
      factor(
        comparison_label,
        levels =
          comparison_labels[comparison_order]
      )
  ]



  p <- ggplot(
    count_dt,
    aes(
      x = Count_plot,
      y = comparison_label,
      fill = Direction
    )
  ) +

    geom_col(
      width = 0.65
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.5
    ) +

    geom_text(
      aes(
        label = abs(Count_plot)
      ),
      size = 4,
      fontface = "bold",
      hjust =
        ifelse(
          count_dt$Count_plot > 0,
          -0.2,
          1.2
        )
    ) +

    scale_fill_manual(
      values = plot_colors
    ) +

    scale_x_continuous(
      labels = abs
    ) +

    labs(
      x = "Number of genes",
      y = NULL,
      fill = NULL,
      title = title_text
    ) +

    theme_bw(
      base_size = 12
    ) +

    theme(
      panel.grid.major.y =
        element_blank(),

      panel.grid.minor =
        element_blank(),

      legend.position =
        "bottom",

      axis.text =
        element_text(
          color = "black"
        ),

      plot.title =
        element_text(
          size = 14,
          face = "bold",
          hjust = 0.5
        )
    )


  ggsave(
    file.path(
      output_dir,
      paste0(
        output_prefix,
        ".pdf"
      )
    ),
    p,
    width = 9,
    height = 4.5
  )


  ggsave(
    file.path(
      output_dir,
      paste0(
        output_prefix,
        ".png"
      )
    ),
    p,
    width = 9,
    height = 4.5,
    dpi = 300
  )


  fwrite(
    count_dt,
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_count.tsv"
      )
    ),
    sep = "\t"
  )


  return(p)

}



# ======================== 7. Plot all genes ========================

message("[4/5] Plotting all genes ...")


draw_butterfly(
  gene_direction_dt,
  title_text =
    "Direction of plasma length similarity to cell-line reference (All genes)",
  output_prefix =
    "butterfly_all_genes"
)



# ======================== 8. Plot significant genes ========================

message("[5/5] Plotting significant genes ...")


gene_direction_sig <- gene_direction_sig_dt[
  !is.na(wilcox_p_D_target) &
    wilcox_p_D_target < p_cutoff
]


draw_butterfly(
  gene_direction_sig,
  title_text =
    "Direction of plasma length similarity to cell-line reference (Significant genes)",
  output_prefix =
    "butterfly_significant_genes"
)



message("========================================")
message("Butterfly plots completed.")
message(
  "Output directory: ",
  output_dir
)
message("========================================")