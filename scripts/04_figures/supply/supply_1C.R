# ============================================================
# Script: Detected RNA counts boxplots by component with IQR,
#         blank ceiling, and pairwise Wilcoxon tests
# Purpose:
#   Read genenum_all.txt and sample information, keep the six
#   RNA types (mRNA, lncRNA, mlncRNA, miRNA, tRNA, piRNA), order
#   components (cfRNA / debris / cell / blank / water), compute
#   per-group quartiles, fences, and blank M+3SD ceiling, draw
#   boxplots with jitter, blank ceiling line, IQR fences and Q1/Q3
#   labels, add significance marks for selected comparisons, and
#   export pairwise Wilcoxon test tables (BH-adjusted). mRNA,
#   lncRNA and mlncRNA plots are also combined into one figure.
# ============================================================

# Load necessary packages
library(tidyverse)
library(ggpubr)
library(gridExtra)

# Define group order
group_order <- c("cfRNA", "debris", "cell", "blank", "water")

# ------------------------------
# 1. Read and tidy data
# ------------------------------

genenum <- read_tsv("./genenum_all.txt",
                    col_select = c(Sample, 
                                   `mRNA(TPM>0)`, 
                                   `lncRNA(TPM>0)`,
                                   `miRNA(RPM>0)`,
                                   `tRNA(RPM>0)`,
                                   `piRNA(RPM>0)`)) %>%
  rename(mRNA   = `mRNA(TPM>0)`,
         lncRNA = `lncRNA(TPM>0)`,
         miRNA  = `miRNA(RPM>0)`,
         tRNA   = `tRNA(RPM>0)`,
         piRNA  = `piRNA(RPM>0)`)

info <- read_tsv("./cellcult_Sample_Information.txt") %>%
  select(Sample, Component)

data <- left_join(genenum, info, by = "Sample") %>%
  mutate(mlncRNA = mRNA + lncRNA,
         Component = factor(Component, levels = group_order))

print(unique(data$Component))

# ------------------------------
# 2. Define function: pairwise significance tests between groups and export table
# ------------------------------
pairwise_test <- function(df, y_var, y_name, output_prefix) {
  groups <- levels(df$Component)
  groups <- groups[groups %in% unique(df$Component)]
  combs <- combn(groups, 2, simplify = FALSE)
  
  results <- map_dfr(combs, function(g) {
    x <- df %>% filter(Component == g[1]) %>% pull({{ y_var }})
    y <- df %>% filter(Component == g[2]) %>% pull({{ y_var }})
    test <- wilcox.test(x, y, exact = FALSE)
    tibble(
      group1 = g[1],
      group2 = g[2],
      statistic = test$statistic,
      p_value = test$p.value,
      method = test$method
    )
  }) %>%
    mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  filename <- paste0(output_prefix, "_pairwise_test.txt")
  write_tsv(results, filename)
  message("Saved pairwise test results to ", filename)
  return(results)
}

# ------------------------------
# 3. Define plotting function (add blank background ceiling and experimental group IQR boundaries, with numeric labels)
# ------------------------------
plot_box_with_iqr <- function(df, y_var, y_label, comparisons = NULL, filename = NULL) {
  
  # Get group levels
  group_levels <- levels(df$Component)
  n_groups <- length(group_levels)
  
  # Compute statistics per group
  iqr_stats <- df %>%
    group_by(Component) %>%
    summarise(
      q1 = quantile({{ y_var }}, 0.25, na.rm = TRUE),
      q3 = quantile({{ y_var }}, 0.75, na.rm = TRUE),
      median = median({{ y_var }}, na.rm = TRUE),
      sd = sd({{ y_var }}, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      IQR = q3 - q1,
      lower_bound = q1 - 1.5 * IQR,
      upper_bound = q3 + 1.5 * IQR,
      background_ceiling = median + 3 * sd,
      # Assign x position for each group (1, 2, 3...)
      x_pos = match(Component, group_levels)
    )
  
  # Base boxplot
  p <- ggplot(df, aes(x = Component, y = {{ y_var }})) +
    geom_boxplot(fill = "skyblue", outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1)
  
  # Annotate blank group background ceiling (red dashed line, only covering blank group width)
  blank_stats <- iqr_stats %>% filter(Component == "blank")
  if (nrow(blank_stats) > 0) {
    x_blank <- blank_stats$x_pos[1]
    ceiling_val <- blank_stats$background_ceiling[1]
    
    # Short line: only covers blank group width (0.3 on each side)
    p <- p +
      annotate("segment", 
               x = x_blank - 0.3, xend = x_blank + 0.3,
               y = ceiling_val, yend = ceiling_val,
               linetype = "dashed", color = "red", linewidth = 1) +
      # Numeric label
      annotate("text",
               x = x_blank + 0.35,  # Slightly to the right
               y = ceiling_val,
               label = paste0("Ceiling=", round(ceiling_val, 0)),
               hjust = 0, vjust = 0.5, size = 3, color = "red", fontface = "bold")
  }
  
  # Annotate upper and lower boundaries of cfRNA, debris, cell groups (only covering each group width)
  exp_groups <- c("cfRNA", "debris", "cell", "blank")
  exp_stats <- iqr_stats %>% filter(Component %in% exp_groups)
  
  for (i in seq_len(nrow(exp_stats))) {
    group_name <- exp_stats$Component[i]
    x_pos <- exp_stats$x_pos[i]
    lb <- exp_stats$lower_bound[i]
    ub <- exp_stats$upper_bound[i]
    q1_val <- exp_stats$q1[i]
    q3_val <- exp_stats$q3[i]
    
    # Lower fence (green dashed line, only this group width)
    p <- p + 
      annotate("segment",
               x = x_pos - 0.3, xend = x_pos + 0.3,
               y = lb, yend = lb,
               linetype = "dashed", color = "darkgreen", linewidth = 0.8, alpha = 0.8) +
      # Lower fence numeric label (left side)
      annotate("text",
               x = x_pos - 0.35,
               y = lb,
               label = paste0("L=", round(lb, 0)),
               hjust = 1, vjust = 0.5, size = 2.8, color = "darkgreen")
    
    # Upper fence (purple dashed line, only this group width)
    p <- p + 
      annotate("segment",
               x = x_pos - 0.3, xend = x_pos + 0.3,
               y = ub, yend = ub,
               linetype = "dashed", color = "purple", linewidth = 0.8, alpha = 0.8) +
      # Upper fence numeric label (right side)
      annotate("text",
               x = x_pos + 0.35,
               y = ub,
               label = paste0("U=", round(ub, 0)),
               hjust = 0, vjust = 0.5, size = 2.8, color = "purple")
  }
  
  # Add Q1 and Q3 numeric labels (on the boxplot)
  p <- p +
    geom_text(data = iqr_stats, aes(x = Component, y = q1, label = round(q1, 0)),
              vjust = 1.5, size = 3, color = "darkred", fontface = "bold") +
    geom_text(data = iqr_stats, aes(x = Component, y = q3, label = round(q3, 0)),
              vjust = -0.8, size = 3, color = "darkred", fontface = "bold")
  
  # Add legend description
  p <- p +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = "Red: Blank Ceiling (M+3SD)\nGreen: Lower Fence (Q1-1.5IQR)\nPurple: Upper Fence (Q3+1.5IQR)",
      hjust = 1, vjust = 1, size = 2.5, color = "gray30",
      xjust = 1, yjust = 1
    ) +
    labs(title = paste(y_label, "detection number by group"),
         x = "Group", y = y_label) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 100, 10, 60)  # Leave space for left/right labels
    )
  
  # If a comparison list is provided, add significance marks
  if (!is.null(comparisons) && length(comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons,
                                 method = "wilcox.test",
                                 label = "p.signif",
                                 tip.length = 0.02,
                                 bracket.size = 0.3,
                                 size = 4,
                                 step.increase = 0.05)
  }
  
  # If a filename is specified, save figures (PNG + PDF)
  if (!is.null(filename)) {
    # Save PNG (original behavior)
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
    # ---- Modification start: also save PDF ----
    pdf_filename <- sub("\\.png$", ".pdf", filename)
    ggsave(pdf_filename, p, width = 12, height = 8)
    # ---- Modification end ----
  }
  
  return(p)
}

# ------------------------------
# 4. Build the required comparison lists
# ------------------------------
actual_groups <- levels(data$Component)[levels(data$Component) %in% unique(data$Component)]

# Comparisons of water with other groups
if ("water" %in% actual_groups) {
  water_comparisons <- map(setdiff(actual_groups, "water"), ~ c("water", .x))
} else {
  water_comparisons <- list()
}

# Pairwise comparisons among cell, debris, cfRNA
cell_debris_cfrna <- intersect(c("cell", "debris", "cfRNA"), actual_groups)
if (length(cell_debris_cfrna) >= 2) {
  cell_debris_comparisons <- combn(cell_debris_cfrna, 2, simplify = FALSE)
} else {
  cell_debris_comparisons <- list()
}

# Merge all comparisons that need to be marked (for mRNA, lncRNA, mlncRNA)
all_comparisons <- c(water_comparisons, cell_debris_comparisons)

# For other RNA types, only use comparisons of water with other groups
other_comparisons <- water_comparisons

# ------------------------------
# 5. Perform analysis and plotting for each RNA type
# ------------------------------

rna_list <- list(
  list(var = "mRNA",    label = "mRNA count"),
  list(var = "lncRNA",  label = "lncRNA count"),
  list(var = "mlncRNA", label = "mlncRNA count"),
  list(var = "miRNA",   label = "miRNA count"),
  list(var = "tRNA",    label = "tRNA count"),
  list(var = "piRNA",   label = "piRNA count")
)

plot_list <- list()

for (rna in rna_list) {
  var_name <- rna$var
  y_label  <- rna$label
  file_prefix <- var_name
  
  # Select comparison list based on RNA type
  if (var_name %in% c("mRNA", "lncRNA", "mlncRNA")) {
    comps <- all_comparisons
  } else {
    comps <- other_comparisons
  }
  
  # Generate figure and save individual files (function already saves both PNG and PDF)
  p <- plot_box_with_iqr(data, !!sym(var_name), y_label, 
                         comparisons = comps, 
                         filename = paste0(file_prefix, "_boxplot.png"))
  
  # If it is one of the first three, save to list for combining
  if (var_name %in% c("mRNA", "lncRNA", "mlncRNA")) {
    plot_list[[var_name]] <- p
  }
  
  # Perform pairwise significance tests and export table
  pairwise_test(data, !!sym(var_name), var_name, file_prefix)
}

# ------------------------------
# 6. Combine the three plots of mRNA, lncRNA, mlncRNA
# ------------------------------
g <- arrangeGrob(plot_list[["mRNA"]], 
                 plot_list[["lncRNA"]], 
                 plot_list[["mlncRNA"]], 
                 ncol = 3)

# Save combined figure (PNG)
ggsave("combined_mlRNA_boxplot.png", g, width = 28, height = 9, dpi = 300)
# ---- Modification start: also save PDF ----
ggsave("combined_mlRNA_boxplot.pdf", g, width = 28, height = 9)
# ---- Modification end ----