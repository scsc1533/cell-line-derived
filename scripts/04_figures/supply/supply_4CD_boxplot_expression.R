# ============================================================
# Script: Boxplots of expression for All_three overlap genes
#         (Hep3B2.1-7 vs HepG2)
# Purpose:
#   From overlap_Hep3B2.1-7_HepG2.csv, take genes whose
#   Intersection == "All_three" and further keep genes with
#   count > 30 in all QC-passed liver cfRNA samples, then compare
#   Hep3B2.1-7 vs HepG2 expression levels (log2(TPM+1)) across
#   cfRNA / cell / debris / cell+debris components. Both overall
#   (all genes pooled) and per-gene boxplots are produced with
#   Wilcoxon significance marks, and the corresponding test
#   result tables are exported.
# ============================================================

###############################################################################
# All_three overlap gene expression boxplots
# Goal: for genes with Intersection == "All_three" in overlap_Hep3B2.1-7_HepG2.csv,
#       after count > 30 filtering, compare Hep3B2.1-7 vs HepG2 expression levels
#       in cfRNA / cell / debris / cell+debris (log2(TPM+1))
#       Includes overall combined plot and per-gene plots, exports Wilcoxon test result tables
###############################################################################

library(data.table)
library(ggplot2)

# ======================== File paths ========================
overlap_csv  <- "01_venn_output/overlap_Hep3B2.1-7_HepG2.csv"
info_file    <- "./cellcult_Sample_Information.txt"
expr_file    <- "./all_mlRNA_TPM.txt"

out_dir <- "07_boxplot_expression"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Colors
cell_colors <- c("Hep3B2.1-7" = "#E64B35FF", "HepG2" = "#4DBBD5FF")
all_comps    <- c("cfRNA", "cell", "debris", "cell+debris")
liver_cells  <- c("Hep3B2.1-7", "HepG2")

# ======================== General functions ========================

# Wilcoxon test (Hep3B2.1-7 vs HepG2)
calc_wilcoxon <- function(data, value_col) {
  data[, {
    vals <- get(value_col)
    h1 <- vals[Cell == "Hep3B2.1-7"]
    h2 <- vals[Cell == "HepG2"]
    if (length(h1) >= 3 && length(h2) >= 3) {
      wt <- wilcox.test(h1, h2, exact = FALSE)
      .(p_value = wt$p.value, y_pos = max(vals, na.rm = TRUE) * 1.12)
    } else {
      .(p_value = NA_real_, y_pos = NA_real_)
    }
  }, by = .(Component)]
}

# Significance markers
sig_label <- function(p) {
  fcase(
    is.na(p),          "",
    p < 0.0001,        "****",
    p < 0.001,         "***",
    p < 0.01,          "**",
    p < 0.05,          "*",
    default            = "ns"
  )
}

# Plotting function
draw_boxplot <- function(plot_data, stat_dt, title_text, y_label) {

  p <- ggplot(plot_data,
              aes(x     = Cell,
                  y     = log2TPM,
                  color = Cell,
                  fill  = Cell)) +

    geom_boxplot(outlier.shape = NA,
                 alpha         = 0,
                 linewidth     = 0.6,
                 width         = 0.6,
                 coef          = 1.5) +

    geom_jitter(width = 0.2,
                alpha = 0.7,
                size  = 1.5) +

    geom_text(data = stat_dt,
              aes(x = 1.5, y = y_pos, label = sig, color = NULL, fill = NULL),
              size      = 3.5,
              fontface  = "bold",
              inherit.aes = FALSE) +

    facet_wrap(~ Component, scales = "free_y", ncol = 4) +

    scale_color_manual(values = cell_colors, guide = "none") +
    scale_fill_manual(values  = cell_colors, guide = "none") +

    labs(x        = NULL,
         y        = y_label,
         title    = title_text) +

    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = "grey80"),
      strip.text         = element_text(size = 10, face = "bold"),
      axis.text.x        = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 8),
      legend.position    = "none",
      plot.title         = element_text(size = 11, face = "bold", hjust = 0.5)
    )

  p
}

# Run a full analysis: draw plots + export test table
run_analysis <- function(gene_agg, y_label, prefix) {

  message(sprintf("  Running Wilcoxon tests for %s ...", prefix))

  # Overall: pool all genes, take mean per sample x Component
  agg_overall <- gene_agg[, .(log2TPM = mean(log2TPM, na.rm = TRUE)),
                          by = .(sample, Cell, Component)]

  stat_overall <- calc_wilcoxon(agg_overall, "log2TPM")
  stat_overall[, sig := sig_label(p_value)]

  p_overall <- draw_boxplot(agg_overall, stat_overall,
                            sprintf("All_three genes (n=%d)", length(all_three_genes)),
                            y_label)
  overall_pdf <- file.path(out_dir, sprintf("boxplot_%s_overall.pdf", prefix))
  ggsave(overall_pdf, plot = p_overall, width = 16, height = 5, device = "pdf")
  message(sprintf("  Saved: %s", overall_pdf))

  # Per gene
  all_stats <- list()
  for (g in all_three_genes) {
    gene_data <- gene_agg[Gene == g]
    if (nrow(gene_data) == 0) next

    stat_gene <- calc_wilcoxon(gene_data, "log2TPM")
    stat_gene[, sig := sig_label(p_value)]
    stat_gene[, Gene := g]
    all_stats[[g]] <- stat_gene

    p_gene <- draw_boxplot(gene_data, stat_gene, g, y_label)
    safe_gene <- gsub("[/\\:*?\"<>| ]", "_", g)
    gene_pdf <- file.path(out_dir, sprintf("boxplot_%s_%s.pdf", prefix, safe_gene))
    ggsave(gene_pdf, plot = p_gene, width = 16, height = 5, device = "pdf")
    message(sprintf("  Saved: %s", gene_pdf))
  }

  # Table
  stat_overall[, Gene := "All_genes_pooled"]
  stat_overall[, sig := sig_label(p_value)]
  results_dt <- rbindlist(c(list(stat_overall), all_stats), use.names = TRUE, fill = TRUE)
  setcolorder(results_dt, c("Gene", "Component", "p_value", "sig", "y_pos"))
  setorder(results_dt, Gene, Component)
  results_dt[, p_formatted := formatC(p_value, format = "e", digits = 3)]

  out_csv <- file.path(out_dir, sprintf("wilcoxon_%s_results.csv", prefix))
  fwrite(results_dt[, .(Gene, Component, p_value = p_formatted, Significance = sig)], out_csv)
  message(sprintf("  Saved: %s", out_csv))
  message(sprintf("\n=== Wilcoxon Test Results (%s) ===", prefix))
  print(results_dt[, .(Gene, Component, p_value = p_formatted, sig)])
}

# ======================== 1. Read overlap, filter All_three genes ========================
message("[1/5] Reading overlap CSV and filtering All_three genes ...")
overlap_dt <- fread(overlap_csv)
all_three_genes <- overlap_dt[Intersection == "All_three", unique(Gene)]
message(sprintf("  All_three genes (before count filter): %d", length(all_three_genes)))
if (length(all_three_genes) == 0) stop("No All_three genes found.")

# ======================== 2. Read sample info + count filtering ========================
message("[2/5] Reading sample info ...")
qc_raw <- fread(info_file)
setnames(qc_raw, trimws(sub("\ufeff", "", names(qc_raw), fixed = TRUE)))
qc_pass <- qc_raw[QC == 1 & Cell %chin% liver_cells,
                  .(Sample, Cell, Component)]
message(sprintf("  Liver QC-passed samples: %d", nrow(qc_pass)))

# ------------------------------ count filtering ------------------------------
message("  Filtering genes by expression count > 30 in liver cfRNA samples ...")
counts_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mlRNA_counts.txt"

cfRNA_liver <- qc_pass[Component == "cfRNA", Sample]
message(sprintf("  Liver cfRNA QC-passed samples: %d", length(cfRNA_liver)))

counts_raw <- fread(counts_file)
gene_col <- intersect(c("gene_id", "sym_id", "Gene"), names(counts_raw))[1]
if (is.na(gene_col)) stop("Cannot find gene ID column in counts file.")
setnames(counts_raw, gene_col, "Gene")

count_cols <- intersect(names(counts_raw), cfRNA_liver)
message(sprintf("  Matched sample columns in counts file: %d / %d",
                length(count_cols), length(cfRNA_liver)))
if (length(count_cols) == 0) {
  stop("No cfRNA sample columns matched. Counts columns: ",
       paste(head(setdiff(names(counts_raw), "Gene"), 10), collapse = ", "))
}

counts_filt <- counts_raw[Gene %chin% all_three_genes,
                          .SD, .SDcols = c("Gene", count_cols)]
message(sprintf("  Genes found in counts file: %d / %d",
                nrow(counts_filt), length(all_three_genes)))

n_samp <- length(count_cols)
counts_filt[, keep := rowSums(.SD > 30, na.rm = TRUE) == n_samp, .SDcols = -"Gene"]
genes_pass_count <- counts_filt[keep == TRUE, unique(Gene)]
all_three_genes <- intersect(all_three_genes, genes_pass_count)
message(sprintf("  Genes passing count > 30 in all %d cfRNA samples: %d",
                n_samp, length(genes_pass_count)))
message(sprintf("  All_three genes (after count filter): %d", length(all_three_genes)))
if (length(all_three_genes) == 0) stop("No genes pass the count filter.")

# ======================== 3. Read expression matrix + long format conversion ========================
message("[3/5] Reading expression matrix (TPM) ...")

expr_raw <- fread(expr_file)
gene_col2 <- intersect(c("sym_id", "gene_id", "Gene"), names(expr_raw))[1]
if (is.na(gene_col2)) stop("Cannot find gene ID column in expression file.")
setnames(expr_raw, gene_col2, "Gene")

# Keep only All_three genes x liver samples
sample_cols <- intersect(names(expr_raw), qc_pass$Sample)
expr_filt <- expr_raw[Gene %chin% all_three_genes,
                      .SD, .SDcols = c("Gene", sample_cols)]
message(sprintf("  Expression matrix filtered: %d genes x %d samples",
                nrow(expr_filt), length(sample_cols)))

# Wide -> long
expr_long <- melt(expr_filt,
                  id.vars       = "Gene",
                  measure.vars  = sample_cols,
                  variable.name = "Sample",
                  value.name    = "TPM")

# Merge metadata
expr_long <- merge(expr_long,
                   qc_pass[, .(Sample, Cell, Component)],
                   by = "Sample", all.x = FALSE)
setnames(expr_long, "Sample", "sample")  # Unify to lowercase
expr_long <- expr_long[Gene %chin% all_three_genes]
message(sprintf("  Long table rows: %d", nrow(expr_long)))

# log2 transformation
expr_long[, log2TPM := log2(TPM + 1)]

# Create cell+debris pooled group
cd_pooled <- expr_long[Component %in% c("cell", "debris")]
cd_pooled[, Component := "cell+debris"]
expr_long <- rbind(expr_long, cd_pooled)

expr_long[, Cell      := factor(Cell,      levels = liver_cells)]
expr_long[, Component := factor(Component, levels = all_comps)]

# ======================== 4. Plotting + testing ========================
message("[4/5] Running Wilcoxon tests and plotting ...")

run_analysis(expr_long, "log2(TPM + 1)", "expression")

message("[5/5] Done. All outputs in: ", out_dir)