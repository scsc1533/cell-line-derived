suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# Offline plotting: pheatmap for A1, ggplot2 for A2 and B.
# Optional command-line arguments: enrichment_file expression_file sample_file outdir
enrichment_file <- "/data/work/01_2603cell_culture/03_analysis/05_RBP/07_RBP_motif_enrichment_new_back/RBP_level_motif_enrichment_Fisher_results.tsv"
expression_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mlRNA_TPM.txt"
sample_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
outdir <- "/data/work/01_2603cell_culture/03_analysis/05_RBP/08_plot"

# Select group x RBP pairs positive in at least one database x region.
fdr_column <- "FDR_by_database_global"
fdr_cutoff <- 0.05
or_cutoff <- 1
database_order <- c("CisBP_RNA", "ATtRACT")
region_order <- c("five_prime_UTR", "CDS", "three_prime_UTR")
group_order <- c("HTR-8-like", "liver-like", "K562-like")
cell_order <- c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")
group_palette <- c("HTR-8-like" = "#00B0F6", "liver-like" = "#F8766D", "K562-like" = "#7474B4")
bubble_width <- 8
bubble_row_height <- 0.24
expression_mode <- "cell_mean" # Alternative: "sample"
supplement_regions <- "three_prime_UTR" # Or use region_order for all regions.
log2_or_color_limit <- 3
neglog10_fdr_size_limit <- 6
draw_dual_support_line <- TRUE
png_dpi <- 300
heatmap_palette <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#EF3B2C", "#CB181D", "#67000D")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L && length(args) != 4L) {
  stop("Supply either no arguments, or enrichment_file expression_file sample_file outdir.")
}
if (length(args) == 4L) {
  enrichment_file <- args[1]; expression_file <- args[2]
  sample_file <- args[3]; outdir <- args[4]
}
if (!expression_mode %in% c("cell_mean", "sample")) stop("Invalid expression_mode")
if (!all(supplement_regions %in% region_order)) stop("Invalid supplement_regions")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) fwrite(x, file.path(outdir, name), sep = "\t", na = "NA")
read_table <- function(path, required) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing/empty file: ", path)
  x <- fread(path, sep = "\t", header = TRUE, check.names = FALSE, na.strings = c("", "NA"))
  if (anyDuplicated(names(x))) stop("Duplicate column names in ", path)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x
}
as_num_checked <- function(x, label) {
  y <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & is.na(y))) stop("Non-numeric values in ", label)
  y
}
normalize_symbol <- function(x) toupper(trimws(as.character(x)))
region_labels <- c(five_prime_UTR = "5'UTR", CDS = "CDS", three_prime_UTR = "3'UTR")
database_labels <- c(CisBP_RNA = "CisBP-RNA", ATtRACT = "ATtRACT")

# 1. Validate statistics and define positive rows without recomputing FDR.
required <- c("database", "feature_group", "region_type", "RBP_SYMBOL", "motif_n",
              "feature_gene_n", "background_gene_n", "feature_hit_gene_n",
              "background_hit_gene_n", "odds_ratio", "p_value", fdr_column)
enr <- read_table(enrichment_file, required)
for (nm in c("database", "feature_group", "region_type", "RBP_SYMBOL")) {
  set(enr, j = nm, value = trimws(as.character(enr[[nm]])))
  if (anyNA(enr[[nm]]) || any(enr[[nm]] == "")) stop("Missing identifiers in ", nm)
}
enr[, RBP_SYMBOL := normalize_symbol(RBP_SYMBOL)]
unknown <- setdiff(unique(enr$database), database_order)
if (length(unknown)) stop("Unrecognized databases: ", paste(unknown, collapse = ", "))
if (!all(enr$region_type %in% region_order)) stop("Unrecognized region_type")
key_cols <- c("database", "feature_group", "region_type", "RBP_SYMBOL")
if (anyDuplicated(enr, by = key_cols)) stop("Duplicate database/group/region/RBP rows")
count_cols <- c("motif_n", "feature_gene_n", "background_gene_n",
                "feature_hit_gene_n", "background_hit_gene_n")
for (nm in unique(c(count_cols, "odds_ratio", "p_value", fdr_column))) {
  set(enr, j = nm, value = as_num_checked(enr[[nm]], nm))
}
for (nm in count_cols) {
  if (any(!is.finite(enr[[nm]]) | enr[[nm]] < 0 | enr[[nm]] != floor(enr[[nm]]))) {
    stop("Invalid nonnegative integer counts in ", nm)
  }
}
if (any(enr$feature_hit_gene_n > enr$feature_gene_n) ||
    any(enr$background_hit_gene_n > enr$background_gene_n)) stop("Hit counts exceed denominators")
for (nm in unique(c("p_value", fdr_column))) {
  z <- enr[[nm]]
  if (any(!is.na(z) & (!is.finite(z) | z < 0 | z > 1))) stop("Invalid probability: ", nm)
}
if (any(enr$odds_ratio < 0, na.rm = TRUE)) stop("Negative odds ratios")
enr[, FDR_plot := get(fdr_column)]
if (any(enr$FDR_plot + 1e-12 < enr$p_value, na.rm = TRUE)) stop("FDR smaller than raw p-value")
enr[, test_available := !is.na(odds_ratio) & !is.na(FDR_plot) & !is.na(p_value) &
      feature_gene_n > 0 & background_gene_n > 0 & motif_n > 0]
enr[, positive := test_available & odds_ratio > or_cutoff & FDR_plot < fdr_cutoff]
enr[, source_row := TRUE]
selected <- enr[positive == TRUE, .(best_fdr = min(FDR_plot)), by = .(feature_group, RBP_SYMBOL)]
if (!nrow(selected)) stop("No positive group x RBP pairs at the specified FDR and OR cutoffs.")
group_order <- c(intersect(group_order, unique(selected$feature_group)),
                 sort(setdiff(unique(selected$feature_group), group_order)))
selected[, group_index := match(feature_group, group_order)]
setorder(selected, group_index, best_fdr, RBP_SYMBOL)
selected[, row_id := sprintf("row_%04d", .I)]
selected[, plot_y := .N - seq_len(.N) + 1L]
# Numeric coordinates prevent layer subsets from retraining discrete row order.
row_axis <- selected[order(plot_y)]
write_tsv(selected, "selected_positive_group_RBPs.tsv")
write_tsv(selected[, .(plot_row_top_to_bottom = .I, feature_group, RBP_SYMBOL, row_id, plot_y)],
          "A_shared_RBP_row_order.tsv")
write_tsv(enr[positive == TRUE], "positive_enrichment_rows.tsv")

# Complete the display grid, keeping absent results distinct from negative tests.
plot_grid <- selected[, .(database = rep(database_order, times = length(region_order)),
                         region_type = rep(region_order, each = length(database_order))),
                      by = .(feature_group, RBP_SYMBOL, row_id)]
bubble <- merge(plot_grid, enr, by = key_cols, all.x = TRUE, sort = FALSE)
bubble[is.na(test_available), test_available := FALSE]
bubble[is.na(positive), positive := FALSE]
bubble[, plot_status := fifelse(!test_available, "No usable result",
                               fifelse(positive, "Positive", "Not significant"))]
bubble[, log2_OR := log2(odds_ratio)]
bubble[, log2_OR_display := pmax(-log2_or_color_limit, pmin(log2_or_color_limit, log2_OR))]
bubble[, neglog10_FDR_display := pmin(neglog10_fdr_size_limit,
                                     -log10(pmax(FDR_plot, .Machine$double.xmin))) ]
support <- bubble[, .(positive_database_n = sum(positive),
                       usable_database_n = sum(test_available),
                       positive_databases = paste(database[positive], collapse = ";")),
                  by = .(feature_group, RBP_SYMBOL, row_id, region_type)]
support[, support_class := fifelse(positive_database_n == 2, "Both databases",
                                   fifelse(positive_database_n == 1, positive_databases,
                                           "No positive database"))]
write_tsv(support, "database_support_by_RBP_region.tsv")
write_tsv(bubble, "A_enrichment_plot_data.tsv")

# 2. Read expression using gene SYMBOL; only QC=1 and Component=cell.
meta <- read_table(sample_file, c("Cell", "Sample", "Component", "QC"))
meta[, `:=`(Sample = trimws(as.character(Sample)), Cell = trimws(as.character(Cell)),
            Component = trimws(as.character(Component)))]
meta[, QC := as_num_checked(QC, "QC")]
meta <- meta[!is.na(QC) & QC == 1 & tolower(Component) == "cell"]
if (!nrow(meta)) stop("No QC=1 cell samples")
if (anyNA(meta$Sample) || anyNA(meta$Cell) || anyDuplicated(meta$Sample)) stop("Invalid sample metadata")
expr <- read_table(expression_file, "sym_id")
meta[, in_expression := Sample %in% names(expr)]
write_tsv(meta, "cell_sample_matching_QC.tsv")
if (any(!meta$in_expression)) warning("Some QC=1 cell samples are absent from expression; see sample matching QC.")
meta <- meta[in_expression == TRUE]
if (!nrow(meta)) stop("No QC=1 cell samples matched the expression matrix")
cell_order <- c(intersect(cell_order, unique(meta$Cell)), sort(setdiff(unique(meta$Cell), cell_order)))
meta[, cell_index := match(Cell, cell_order)]
setorder(meta, cell_index, Sample)
expr[, RBP_SYMBOL := normalize_symbol(sym_id)]
expr <- expr[RBP_SYMBOL %in% selected$RBP_SYMBOL, c("RBP_SYMBOL", meta$Sample), with = FALSE]
if (anyDuplicated(expr$RBP_SYMBOL)) stop("Duplicate selected SYMBOLs in expression matrix; resolve upstream, do not silently sum.")
for (nm in meta$Sample) {
  values <- as_num_checked(expr[[nm]], nm)
  if (any(!is.finite(values) | values < 0)) stop("Missing/invalid TPM in ", nm)
  set(expr, j = nm, value = values)
}
missing_expr <- selected[!RBP_SYMBOL %in% expr$RBP_SYMBOL]
write_tsv(missing_expr, "selected_RBPs_missing_expression.tsv")
if (nrow(missing_expr)) warning("Some selected RBPs lack expression rows; heatmap cells will be gray.")
long <- melt(expr, id.vars = "RBP_SYMBOL", variable.name = "Sample", value.name = "TPM",
             variable.factor = FALSE)
long <- merge(long, meta[, .(Sample, Cell)], by = "Sample", sort = FALSE)
long[, log2_TPM_plus1 := log2(TPM + 1)]
means <- long[, .(mean_TPM = mean(TPM), mean_log2_TPM_plus1 = mean(log2_TPM_plus1),
                   sample_n = .N), by = .(RBP_SYMBOL, Cell)]
write_tsv(long, "cell_RBP_expression_per_sample.tsv")
write_tsv(means, "cell_RBP_expression_means.tsv")
if (expression_mode == "cell_mean") {
  display_cols <- cell_order
  cell_n <- meta[, .N, by = Cell]
  expression_labels <- setNames(paste0(cell_n$Cell, "\n(n=", cell_n$N, ")"), cell_n$Cell)
  heat_values <- means[, .(RBP_SYMBOL, display_column = Cell, value = mean_log2_TPM_plus1)]
  expression_legend <- "Mean log2(TPM + 1)"
} else {
  display_cols <- meta$Sample
  expression_labels <- setNames(paste(meta$Sample, meta$Cell, sep = "\n"), meta$Sample)
  heat_values <- long[, .(RBP_SYMBOL, display_column = Sample, value = log2_TPM_plus1)]
  expression_legend <- "log2(TPM + 1)"
}
heat <- selected[, .(display_column = display_cols), by = .(feature_group, RBP_SYMBOL, row_id)]
heat <- merge(heat, heat_values, by = c("RBP_SYMBOL", "display_column"), all.x = TRUE, sort = FALSE)
write_tsv(heat, "A_expression_plot_data.tsv")

# 3. Shared row ordering makes separately exported panels easy to align.
factor_rows <- function(x) {
  x <- copy(x)
  x[, feature_group := factor(feature_group, levels = group_order)]
  x[, plot_y := selected$plot_y[match(as.character(row_id), selected$row_id)]]
  if (anyNA(x$plot_y)) stop("Unmapped RBP row in plotting data")
  if ("region_type" %in% names(x)) x[, region_type := factor(region_type, levels = region_order)]
  x
}
heat <- factor_rows(heat)
bubble <- factor_rows(bubble)
support <- factor_rows(support)
heat[, display_column := factor(display_column, levels = display_cols)]
bubble[, db_x := match(database, database_order)]
common_theme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(fill = "#F0F2F3", color = NA),
        strip.text.y = element_text(angle = 0),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = "bottom", panel.spacing.y = grid::unit(0.18, "in"),
        plot.title = element_text(face = "bold", size = 13),
        plot.caption = element_text(hjust = 0, size = 9),
        plot.margin = margin(10, 12, 10, 10))
save_plot <- function(p, stem, width, height) {
  ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = width, height = height,
         device = "pdf", useDingbats = FALSE, limitsize = FALSE)
  ggsave(file.path(outdir, paste0(stem, ".png")), p, width = width, height = height,
         dpi = png_dpi, bg = "white", limitsize = FALSE)
}
fig_height <- max(4.8, 0.30 * nrow(selected) + 3.0 + 0.18 * length(group_order))
group_colors <- setNames(rep(c("#BA8C32", "#467FAD", "#689454"),
                             length.out = length(group_order)), group_order)
known_groups <- intersect(group_order, names(group_palette))
group_colors[known_groups] <- group_palette[known_groups]

# Explicit row/column indexing preserves the order shared with A2.
heat_matrix <- matrix(NA_real_, nrow = nrow(selected), ncol = length(display_cols),
                      dimnames = list(selected$row_id, display_cols))
heat_matrix[cbind(match(as.character(heat$row_id), selected$row_id),
                   match(as.character(heat$display_column), display_cols))] <- heat$value
row_annotation <- data.frame(
  Group = factor(selected$feature_group, levels = group_order), row.names = selected$row_id)
group_gaps <- which(diff(match(selected$feature_group, group_order)) != 0)
if (!length(group_gaps)) group_gaps <- NULL
heat_max <- max(c(1, heat_matrix[is.finite(heat_matrix)]))
heat_breaks <- seq(0, heat_max, length.out = 101)
heat_colors <- colorRampPalette(heatmap_palette)(100)
p_heat <- pheatmap::pheatmap(
  heat_matrix, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
  color = heat_colors, breaks = heat_breaks, na_col = "#D0D0D0",
  annotation_row = row_annotation, annotation_colors = list(Group = group_colors),
  annotation_names_row = FALSE, drop_levels = FALSE, gaps_row = group_gaps,
  labels_row = selected$RBP_SYMBOL, labels_col = unname(expression_labels[display_cols]),
  show_rownames = TRUE, show_colnames = TRUE, angle_col = "45",
  fontsize = 10, fontsize_row = 10, fontsize_col = 9, border_color = "white",
  main = paste("A1  RBP expression", expression_legend, sep = "\n"), silent = TRUE)
expr_width <- max(6.5, length(display_cols) * 0.62 + 3.2)
expr_height <- max(3.5, nrow(selected) * 0.26 + 2.5)
save_heatmap <- function(path, device) {
  if (device == "pdf") {
    grDevices::pdf(path, width = expr_width, height = expr_height, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = expr_width, height = expr_height, units = "in", res = png_dpi)
  }
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
}
save_heatmap(file.path(outdir, "A1_RBP_expression_heatmap.pdf"), "pdf")
save_heatmap(file.path(outdir, "A1_RBP_expression_heatmap.png"), "png")
write_tsv(data.table(Group = names(group_colors), color = unname(group_colors)), "group_annotation_colors.tsv")
write_tsv(data.table(display_column = display_cols, label = unname(expression_labels[display_cols])),
          "A1_expression_column_order.tsv")

p_bubble <- ggplot(bubble, aes(db_x, plot_y)) +
  geom_blank(aes(y = plot_y - 0.5)) +
  geom_blank(aes(y = plot_y + 0.5)) +
  facet_grid(feature_group ~ region_type, scales = "free_y", space = "free_y",
             labeller = labeller(region_type = as_labeller(region_labels)))
if (draw_dual_support_line) {
  p_bubble <- p_bubble + geom_segment(data = support[positive_database_n == 2],
                                     aes(x = 1, xend = 2, y = plot_y, yend = plot_y),
                                     inherit.aes = FALSE, color = "#505050", linewidth = 0.45)
}
p_bubble <- p_bubble +
  geom_point(data = bubble[test_available == TRUE],
             aes(fill = log2_OR_display, size = neglog10_FDR_display, color = plot_status),
             shape = 21, stroke = 0.8) +
  geom_point(data = bubble[test_available == FALSE], shape = 4, size = 2.6, color = "#888888") +
  scale_y_continuous(breaks = row_axis$plot_y, labels = row_axis$RBP_SYMBOL,
                     expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = 1:2, labels = unname(database_labels[database_order]), limits = c(0.5, 2.5)) +
  scale_fill_gradient2(low = "#3978AC", mid = "#F7F7F7", high = "#BE354C", midpoint = 0,
                       limits = c(-log2_or_color_limit, log2_or_color_limit), name = "log2(OR)") +
  scale_size_continuous(range = c(1.8, 7), limits = c(0, neglog10_fdr_size_limit),
                        breaks = seq(0, neglog10_fdr_size_limit, length.out = 4),
                        name = "-log10(FDR)") +
  scale_color_manual(values = c(Positive = "#202020", `Not significant` = "#BDBDBD"),
                      name = "Test result", drop = FALSE) +
  labs(title = "A2  Regional enrichment and database support",
       subtitle = paste0("Positive: OR > ", or_cutoff, "; FDR < ", fdr_cutoff),
       caption = paste0("Black outline: positive; gray: not significant; x: no usable result.\n",
                        "Line: both databases positive. Color capped at +/-", log2_or_color_limit,
                        "; size capped at ", neglog10_fdr_size_limit, ".")) +
  common_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
                       panel.spacing.x = grid::unit(0.10, "in"),
                       panel.spacing.y = grid::unit(0.10, "in"),
                       legend.text = element_text(size = 8), legend.title = element_text(size = 9),
                       legend.key.width = grid::unit(0.32, "cm"),
                       legend.spacing.x = grid::unit(0.1, "cm"),
                       plot.caption = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3, fill = "white")),
         size = guide_legend(nrow = 1, override.aes = list(fill = "#888888", color = "#555555")))
# Check the rendered facet scales, not the physical order of merged table rows.
bubble_build <- ggplot_build(p_bubble)
stopifnot(identical(rownames(heat_matrix), selected$row_id))
for (panel_i in seq_len(nrow(bubble_build$layout$layout))) {
  group_i <- as.character(bubble_build$layout$layout$feature_group[panel_i])
  displayed_y <- bubble_build$layout$panel_params[[panel_i]]$y$get_breaks()
  displayed_y <- sort(as.numeric(displayed_y[is.finite(displayed_y)]), decreasing = TRUE)
  expected_y <- selected[feature_group == group_i, plot_y]
  if (!isTRUE(all.equal(displayed_y, as.numeric(expected_y), check.attributes = FALSE))) {
    stop("A1/A2 row order mismatch in group: ", group_i)
  }
}
save_plot(p_bubble, "A2_RBP_region_database_bubble", bubble_width,
          max(3.5, bubble_row_height * nrow(selected) + 2.6 + 0.10 * length(group_order)))

# 4. Supplement B: same selected RBP rows, with both database estimates shown.
# Each database uses its own counts. Do not pool hits across databases.
supp <- bubble[as.character(region_type) %in% supplement_regions]
supp[, `:=`(feature_fraction = feature_hit_gene_n / feature_gene_n,
             background_fraction = background_hit_gene_n / background_gene_n)]
supp[, `:=`(feature_label = paste0(feature_hit_gene_n, "/", feature_gene_n),
             background_label = paste0(background_hit_gene_n, "/", background_gene_n))]
supp[, database := factor(database, levels = database_order)]
write_tsv(supp, "B_hit_proportion_plot_data.tsv")
for (reg in supplement_regions) {
  sub <- supp[as.character(region_type) == reg]
  valid <- sub[test_available == TRUE]
  p_supp <- ggplot(sub, aes(y = plot_y)) +
    geom_blank(aes(x = 0, y = plot_y - 0.5)) +
    geom_blank(aes(x = 0, y = plot_y + 0.5)) +
    geom_segment(data = valid, aes(x = background_fraction, xend = feature_fraction,
                                   yend = plot_y), color = "#B9BDC1", linewidth = 0.65) +
    geom_point(data = valid, aes(x = background_fraction), shape = 21,
               fill = "white", color = "#858C93", size = 2.7, stroke = 0.8) +
    geom_point(data = valid, aes(x = feature_fraction, fill = feature_group),
               shape = 21, color = "white", size = 3.2, stroke = 0.4) +
    geom_text(data = valid, aes(x = feature_fraction, label = feature_label, color = feature_group),
              nudge_y = 0.23, size = 2.6, show.legend = FALSE) +
    geom_text(data = valid, aes(x = background_fraction, label = background_label),
              nudge_y = -0.23, size = 2.6, color = "#727A81") +
    geom_text(data = sub[test_available == FALSE], aes(x = 0.5, label = "No usable result"),
              size = 3, color = "#888888") +
    facet_grid(feature_group ~ database, scales = "free_y", space = "free_y",
               labeller = labeller(database = as_labeller(database_labels))) +
    scale_y_continuous(breaks = row_axis$plot_y, labels = row_axis$RBP_SYMBOL,
                       expand = expansion(mult = 0)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"),
                       limits = c(-0.08, 1.08), expand = expansion(mult = 0)) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_color_manual(values = group_colors, guide = "none") +
    labs(title = paste0("B  Motif-hit gene proportions: ", region_labels[[reg]]),
         subtitle = "Colored point: feature genes; open gray point: background (feature set excluded).",
         x = "Genes with at least one RBP motif hit (%)",
         caption = "Labels: hit genes / tested genes. Rows follow A.\nResults from both databases shown, including nonsignificant tests.") +
    common_theme + theme(axis.title.x = element_text(),
                         panel.grid.major.x = element_line(color = "#EBEDEF", linewidth = 0.3))
  save_plot(p_supp, paste0("B_motif_hit_proportions_", reg), 10, max(fig_height, 0.45 * nrow(selected) + 3))
}
settings <- data.table(item = c("enrichment_file", "expression_file", "sample_file", "FDR_column",
                               "FDR_cutoff", "OR_cutoff", "expression_mode", "expression_transform",
                               "positive_group_RBP_pairs", "unique_RBPs", "QC_cell_samples", "groups",
                               "missing_result_cells", "missing_expression_RBP_pairs"),
                      value = as.character(c(enrichment_file, expression_file, sample_file, fdr_column,
                                              fdr_cutoff, or_cutoff, expression_mode,
                                              if (expression_mode == "cell_mean") "mean(log2(TPM+1)); no z-score" else "log2(TPM+1); no z-score",
                                              nrow(selected), uniqueN(selected$RBP_SYMBOL), nrow(meta),
                                              paste(group_order, collapse = ";"), sum(!bubble$test_available), nrow(missing_expr))))
write_tsv(settings, "plot_run_summary.tsv")
writeLines(c("These figures use observed input data; p-values/FDRs are not recalculated.",
             "An absent result is not proof of an absent database motif or a negative test.",
             "Rows are selected group x RBP pairs, positive in at least one tested region/database.",
             "The expression heatmap uses cell samples only. Enrichment gene sets retain their input definitions.",
             "Cross-database consistency is not necessarily independent experimental replication."),
           file.path(outdir, "figure_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("Selected group x RBP pairs:", nrow(selected), "\n")
cat("QC=1 cell samples:", nrow(meta), "\n")
cat("Figures and plotting tables:", outdir, "\n")
