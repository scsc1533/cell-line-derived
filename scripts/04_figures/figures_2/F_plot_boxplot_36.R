# ============================================================
# Script: Boxplots of >36bp read proportion and freqTC for
#         All_three overlap genes (Hep3B2.1-7 vs HepG2)
# Purpose:
#   From overlap_Hep3B2.1-7_HepG2.csv, take genes whose
#   Intersection == "All_three" and further keep genes with
#   count > 30 in all QC-passed liver cfRNA samples, then for
#   each of the two metrics (proportion of reads > 36 bp, and
#   freqTC), compare Hep3B2.1-7 vs HepG2 across cfRNA / cell /
#   debris / cell+debris components, draw overall (all genes
#   pooled) and per-gene boxplots with Wilcoxon significance
#   marks, and export the corresponding test result tables.
# ============================================================

###############################################################################
# All_three overlap gene >36bp read proportion & freqTC boxplots
# Goal: for genes with Intersection == "All_three" in overlap_Hep3B2.1-7_HepG2.csv,
#       compare Hep3B2.1-7 vs HepG2 by boxplot
#       (1) >36bp read proportion  (2) freqTC
#       Includes overall combined plot and per-gene plots, exports Wilcoxon test result tables
###############################################################################

library(data.table)
library(ggplot2)

# ======================== File paths ========================
overlap_csv  <- "./overlap_Hep3B2.1-7_HepG2.csv"
info_file    <- "./cellcult_Sample_Information.txt"
length_file  <- "./all_samples_length_merged.mlncRNA.txt"
freqTC_file  <- "./freqTC.mRNA1.txt"

out_dir <- "./06_plot_boxplot_36bp"
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
draw_boxplot <- function(plot_data, stat_dt, title_text, value_col, y_label) {

  p <- ggplot(plot_data,
              aes(x     = Cell,
                  y     = .data[[value_col]],
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
run_analysis <- function(gene_agg, value_col, y_label, prefix) {

  message(sprintf("  Running Wilcoxon tests for %s ...", prefix))

  # Overall: pool all genes, take mean per sample x Component
  agg_overall <- gene_agg[, .(Value = mean(get(value_col), na.rm = TRUE)),
                          by = .(sample, Cell, Component)]
  setnames(agg_overall, "Value", value_col)

  stat_overall <- calc_wilcoxon(agg_overall, value_col)
  stat_overall[, sig := sig_label(p_value)]

  p_overall <- draw_boxplot(agg_overall, stat_overall,
                            sprintf("All_three genes (n=%d)", length(all_three_genes)),
                            value_col, y_label)
  overall_pdf <- file.path(out_dir, sprintf("boxplot_%s_overall.pdf", prefix))
  ggsave(overall_pdf, plot = p_overall, width = 12, height = 4, device = "pdf")
  message(sprintf("  Saved: %s", overall_pdf))

  # Per gene
  all_stats <- list()
  for (g in all_three_genes) {
    gene_data <- gene_agg[Gene == g]
    if (nrow(gene_data) == 0) next

    stat_gene <- calc_wilcoxon(gene_data, value_col)
    stat_gene[, sig := sig_label(p_value)]
    stat_gene[, Gene := g]
    all_stats[[g]] <- stat_gene

    p_gene <- draw_boxplot(gene_data, stat_gene, g, value_col, y_label)
    safe_gene <- gsub("[/\\:*?\"<>| ]", "_", g)
    gene_pdf <- file.path(out_dir, sprintf("boxplot_%s_%s.pdf", prefix, safe_gene))
    ggsave(gene_pdf, plot = p_gene, width = 12, height = 4, device = "pdf")
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
message("[1/6] Reading overlap CSV and filtering All_three genes ...")
overlap_dt <- fread(overlap_csv)
all_three_genes <- overlap_dt[Intersection == "All_three", unique(Gene)]
message(sprintf("  All_three genes (before count filter): %d", length(all_three_genes)))
if (length(all_three_genes) == 0) stop("No All_three genes found.")

# ======================== 2. Read sample info ========================
message("[2/6] Reading sample info ...")
qc_raw <- fread(info_file)
setnames(qc_raw, trimws(sub("\ufeff", "", names(qc_raw), fixed = TRUE)))
qc_pass <- qc_raw[QC == 1 & Cell %chin% liver_cells,
                  .(Sample, Cell, Component)]
message(sprintf("  Liver QC-passed samples: %d", nrow(qc_pass)))

# ======================== 2.5. Expression count filtering ========================
message("[2.5/6] Filtering genes by expression count > 30 in liver cfRNA samples ...")
counts_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mlRNA_counts.txt"

# All QC==1 liver cfRNA samples
cfRNA_liver <- qc_pass[Component == "cfRNA", Sample]
message(sprintf("  Liver cfRNA QC-passed samples: %d", length(cfRNA_liver)))

counts_raw <- fread(counts_file)

# Locate gene name column (may be gene_id or sym_id)
gene_col <- intersect(c("gene_id", "sym_id", "Gene"), names(counts_raw))[1]
if (is.na(gene_col)) stop("Cannot find gene ID column in counts file.")
setnames(counts_raw, gene_col, "Gene")

# Find liver cfRNA sample columns present in counts file
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

# Keep genes with count > 30 in all liver cfRNA samples
n_samp <- length(count_cols)
counts_filt[, keep := rowSums(.SD > 30, na.rm = TRUE) == n_samp, .SDcols = -"Gene"]
genes_pass_count <- counts_filt[keep == TRUE, unique(Gene)]
all_three_genes <- intersect(all_three_genes, genes_pass_count)
message(sprintf("  Genes passing count > 30 in all %d cfRNA samples: %d",
                n_samp, length(genes_pass_count)))
message(sprintf("  All_three genes (after count filter): %d", length(all_three_genes)))
if (length(all_three_genes) == 0) stop("No genes pass the count filter.")

# ======================== 3. Analysis (A): >36bp read proportion ========================
message("[3/6] ==== Analysis A: >36bp read proportion ====")

length_raw <- fread(length_file)
n_len_cols <- ncol(length_raw) - 4L
setnames(length_raw, seq_along(length_raw),
         c("sample", "Transcript", "Gene", "Type", as.character(seq_len(n_len_cols))))

length_filt <- length_raw[sample %in% qc_pass$Sample & Gene %chin% all_three_genes]
message(sprintf("  Filtered rows: %d (was %d)", nrow(length_filt), nrow(length_raw)))

length_filt <- merge(length_filt,
                     qc_pass[, .(Sample, Cell, Component)],
                     by.x = "sample", by.y = "Sample", all.x = FALSE)
length_filt[, c("Transcript", "Type") := NULL]

length_cols <- grep("^[0-9]+$", names(length_filt), value = TRUE)
length_long <- melt(length_filt,
                    id.vars       = c("sample", "Gene", "Cell", "Component"),
                    measure.vars  = length_cols,
                    variable.name = "Length",
                    value.name    = "Count")
length_long[, Length := as.integer(as.character(Length))]

gene_agg_36bp <- length_long[, .(
  TotalCount = sum(Count, na.rm = TRUE),
  LongCount  = sum(ifelse(Length > 36, Count, 0), na.rm = TRUE)
), by = .(sample, Gene, Cell, Component)]
gene_agg_36bp[, Prop36 := ifelse(TotalCount > 0, LongCount / TotalCount, 0)]

rm(length_raw, length_filt, length_long)
gc()

cd_pooled <- gene_agg_36bp[Component %in% c("cell", "debris")]
cd_pooled[, Component := "cell+debris"]
gene_agg_36bp <- rbind(gene_agg_36bp, cd_pooled)
gene_agg_36bp[, Cell      := factor(Cell,      levels = liver_cells)]
gene_agg_36bp[, Component := factor(Component, levels = all_comps)]

run_analysis(gene_agg_36bp, "Prop36", "Proportion of reads > 36 bp", "36bp")

# ======================== 4. Analysis (B): freqTC ========================
message("[4/6] ==== Analysis B: freqTC ====")

freqTC_raw <- fread(freqTC_file)
freqTC_long <- melt(freqTC_raw,
                    id.vars       = "Gene",
                    variable.name = "Sample",
                    value.name    = "freqTC")

# Merge sample info, filter liver samples + All_three genes
freqTC_merged <- merge(freqTC_long,
                       qc_pass[, .(Sample, Cell, Component)],
                       by = "Sample", all.x = FALSE)
setnames(freqTC_merged, "Sample", "sample")  # Unify to lowercase, consistent with by in run_analysis
freqTC_merged <- freqTC_merged[Gene %chin% all_three_genes]
message(sprintf("  freqTC rows after merge: %d", nrow(freqTC_merged)))

# Create cell+debris pooled group
cd_freq <- freqTC_merged[Component %in% c("cell", "debris")]
cd_freq[, Component := "cell+debris"]
freqTC_merged <- rbind(freqTC_merged, cd_freq)
freqTC_merged[, Cell      := factor(Cell,      levels = liver_cells)]
freqTC_merged[, Component := factor(Component, levels = all_comps)]

run_analysis(freqTC_merged, "freqTC", "freqTC", "freqTC")

# ======================== 5. Done ========================
message("[5/6] Done. All outputs in: ", out_dir)