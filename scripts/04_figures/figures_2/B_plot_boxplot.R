# ============================================================
# Script: Boxplot for a specified cell line set and gene(s)
# Purpose:
#   Combine Hep3B2.1-7 and HepG2 into a target group "Hep3B2.1-7_HepG2",
#   read the mlRNA TPM matrix and QC-passed sample metadata, compare
#   the expression of the specified gene(s) (e.g. IGF2) between the
#   target group and all other cell lines, facet by Component
#   (cell+debris vs cfRNA), and add Wilcoxon significance marks.
#   Both per-gene and combined multi-gene plots are produced.
# ============================================================

###############################################################################
# Boxplots for a specified cell line set + specified genes (target vs other)
# Cell lines: Hep3B2.1-7 + HepG2 (merged as Hep3B2.1-7_HepG2)
# Gene: IGF2
###############################################################################

# ======================== 0. Load dependencies ========================
library(data.table)
library(ggplot2)

# ======================== 1. Configure paths ========================
qc_file <- "./information_table/cellcult_Sample_Information.txt"
expr_file <- "./all_mlRNA_TPM.txt"
out_dir <- "02_boxplot_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Colors (consistent with the original)
cell_colors <- c(
  "Hep3B2.1-7_HepG2" = "#4DBBD5FF",
  "Hep3B2.1-7"       = "#4DBBD5FF",
  "HepG2"            = "#00A087FF",
  "K562"             = "#E64B35FF",
  "HTR-8/SVneo"      = "#3C5488FF",
  "HEK293T"          = "#F39B7FFF",
  "Other"            = "#000000"
)

# ======================== 2. Specify target cell line set and target genes ========================
target_cells <- c("Hep3B2.1-7", "HepG2")   # Target cell line list
target_label <- "Hep3B2.1-7_HepG2"         # Used for output file names and titles
target_genes <- c("IGF2")                  # Target gene(s) (can be multiple)

# ======================== 3. Read QC metadata ========================
message("[1/4] Reading QC metadata ...")
qc <- fread(qc_file)
qc_pass <- qc[QC == 1]
message(sprintf("  QC-passed samples: %d", nrow(qc_pass)))

# Merge cell and debris into cell+debris
qc_pass[Component %in% c("cell", "debris"), Component := "cell+debris"]

# ======================== 4. Read expression matrix and convert to long format ========================
message("[2/4] Reading expression matrix ...")
expr <- fread(expr_file)
message(sprintf("  Expression matrix: %d genes x %d samples",
                nrow(expr), ncol(expr) - 1L))

sample_cols <- intersect(names(expr)[-1], qc_pass$Sample)
expr_long <- melt(expr,
                  id.vars       = "sym_id",
                  measure.vars  = sample_cols,
                  variable.name = "Sample",
                  value.name    = "TPM")
setnames(expr_long, "sym_id", "Gene")

expr_long <- merge(expr_long, qc_pass[, .(Sample, Cell, Component)],
                   by = "Sample", all.x = FALSE)
message(sprintf("  Long table rows: %d", nrow(expr_long)))

# ======================== 5. Filter expression data for target genes ========================
message("[3/4] Filtering expression data for target genes ...")

existing_genes <- intersect(target_genes, expr_long$Gene)
if (length(existing_genes) == 0) {
  stop("None of the target genes are found in the expression matrix.")
}
missing_genes <- setdiff(target_genes, existing_genes)
if (length(missing_genes) > 0) {
  warning(sprintf("Genes not found: %s", paste(missing_genes, collapse = ", ")))
}
target_genes <- existing_genes

plot_data <- expr_long[Gene %in% target_genes]

# Grouping: target cell line set vs other
plot_data[, Group := ifelse(Cell %in% target_cells, target_label, "Other")]

# Drop gene x Component x Group combinations with all-zero expression
plot_data[, has_expr := any(TPM > 0), by = .(Gene, Component, Group)]
plot_data <- plot_data[has_expr == TRUE]
plot_data[, has_expr := NULL]

# ======================== 6. Plotting ========================
message("[4/4] Generating plots ...")

# Color mapping
group_colors <- c(setNames(cell_colors[target_label], target_label),
                  "Other" = "#000000")

# Wilcoxon test
stat_dt <- plot_data[, {
  target_vals <- TPM[Group == target_label]
  other_vals  <- TPM[Group == "Other"]
  if (length(target_vals) >= 3 && length(other_vals) >= 3) {
    wt <- wilcox.test(target_vals, other_vals, exact = FALSE)
    .(p_value = wt$p.value, y_pos = max(TPM, na.rm = TRUE) * 1.15)
  } else {
    .(p_value = NA_real_, y_pos = NA_real_)
  }
}, by = .(Gene, Component)]

stat_dt[, sig := fcase(
  is.na(p_value),                    "",
  p_value < 0.0001,                  "****",
  p_value < 0.001,                   "***",
  p_value < 0.01,                    "**",
  p_value < 0.05,                    "*",
  default                            = "ns"
)]
stat_dt <- stat_dt[!is.na(p_value)]

# ---- 6.1 Plot per gene ----
for (gene in target_genes) {
  gene_data <- plot_data[Gene == gene]
  gene_stat <- stat_dt[Gene == gene]

  p <- ggplot(gene_data,
              aes(x     = Group,
                  y     = TPM + 0.01,
                  color = Group,
                  fill  = Group)) +

    geom_boxplot(outlier.shape = NA,
                 alpha         = 0,
                 linewidth     = 0.6,
                 width         = 0.6,
                 coef          = 1.5) +

    geom_jitter(width = 0.2,
                alpha = 0.7,
                size  = 1.5) +

    geom_text(data = gene_stat,
              aes(x = 1.5, y = y_pos, label = sig, color = NULL, fill = NULL),
              size      = 3.5,
              fontface  = "bold",
              inherit.aes = FALSE) +

    facet_wrap(~ Component, scales = "free_y") +

    scale_color_manual(values = group_colors, guide = "none") +
    scale_fill_manual(values = group_colors, guide = "none") +

    scale_y_log10(breaks = scales::log_breaks(n = 5),
                  labels = scales::label_number()) +

    labs(x        = NULL,
         y        = "TPM + 0.01 (log10)",
         title    = sprintf("%s - %s", target_label, gene)) +

    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = "grey80"),
      strip.text         = element_text(size = 12, face = "bold"),
      axis.text.x        = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 11),
      legend.position    = "none",
      plot.title         = element_text(size = 12, face = "bold", hjust = 0.5)
    )

  out_pdf <- file.path(out_dir,
                       sprintf("boxplot_%s_%s.pdf",
                               gsub("/", "_", target_label), gene))
  ggsave(out_pdf, plot = p, width = 6, height = 5, device = "pdf")
  message(sprintf("  Saved: %s", out_pdf))
}

# ---- 6.2 Combined plot (all target genes in one figure) ----
if (length(target_genes) > 1) {
  p_all <- ggplot(plot_data,
                  aes(x     = Group,
                      y     = TPM + 0.01,
                      color = Group,
                      fill  = Group)) +

    geom_boxplot(outlier.shape = NA,
                 alpha         = 0,
                 linewidth     = 0.4,
                 width         = 0.6,
                 coef          = 1.5) +

    geom_jitter(width = 0.2,
                alpha = 0.5,
                size  = 0.8) +

    geom_text(data = stat_dt,
              aes(x = 1.5, y = y_pos, label = sig, color = NULL, fill = NULL),
              size      = 2.5,
              fontface  = "bold",
              inherit.aes = FALSE) +

    facet_grid(Gene ~ Component, scales = "free_y") +

    scale_color_manual(values = group_colors, guide = "none") +
    scale_fill_manual(values = group_colors, guide = "none") +

    scale_y_log10(breaks = scales::log_breaks(n = 4),
                  labels = scales::label_number()) +

    labs(x        = NULL,
         y        = "TPM + 0.01 (log10)",
         title    = target_label) +

    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = "grey80"),
      strip.text         = element_text(size = 12, face = "bold"),
      axis.text.x        = element_text(size = 12),
      axis.text.y        = element_text(size = 11),
      legend.position    = "none",
      plot.title         = element_text(size = 12, face = "bold", hjust = 0.5)
    )

  n_genes <- length(target_genes)
  fig_height <- max(4, n_genes * 2.2)
  out_all <- file.path(out_dir,
                       sprintf("boxplot_all_%s.pdf", gsub("/", "_", target_label)))
  ggsave(out_all, plot = p_all, width = 4, height = fig_height, device = "pdf",
         limitsize = FALSE)
  message(sprintf("  Saved combined plot: %s", out_all))
}

message("\nDone. All outputs in: ", out_dir)