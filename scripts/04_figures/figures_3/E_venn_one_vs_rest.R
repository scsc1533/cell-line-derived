# ============================================================
# Script: Venn analysis of one-vs-rest specific genes vs cfRNA
#         highly expressed genes
# Purpose:
#   From the one-vs-rest gene score table, take genes with
#   perm_fdr_wasserstein < 0.05 for selected comparisons
#   (HepG2 + Hep3B2.1-7, HTR-8/SVneo, K562), read the matching
#   overlap CSV to get cfRNA-high genes, and for each group
#   produce a two-set Venn diagram (One vs Rest vs cfRNA high)
#   together with a detailed overlap table.
# ============================================================

###############################################################################
# Venn analysis of one-vs-rest specific genes vs cfRNA highly expressed genes
# Goal: for genes with perm_fdr_wasserstein < 0.05 in the scores file,
#       overlap with cfRNA-high genes of the corresponding cell line, export Venn plots and tables
###############################################################################

library(data.table)
library(VennDiagram)

# ======================== File paths ========================
scores_file <- "./01_length_distribution_one_vs_rest/all_comparisons_gene_scores.tsv"
overlap_dir <- "./02_figure2/01_venn_output"

out_dir <- "02_venn_one_vs_rest"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======================== Config: comparison -> overlap CSV -> label ========================
comp_map <- data.table(
  comparison = c(
    "HepG2_Hep3B2.1_7_vs_rest",
    "HTR_8_SVneo_vs_rest",
    "K562_vs_rest"
  ),
  overlap_csv = c(
    "overlap_Hep3B2.1-7_HepG2.csv",
    "overlap_HTR-8_SVneo.csv",
    "overlap_K562.csv"
  ),
  label = c(
    "Hep3B2.1-7 & HepG2",
    "HTR-8/SVneo",
    "K562"
  )
)

# ======================== 1. Read scores file ========================
message("[1/3] Reading one-vs-rest scores ...")
scores <- fread(scores_file)
message(sprintf("  Total rows: %d", nrow(scores)))

# Filter: target comparisons + fdr < 0.05
scores_sig <- scores[comparison %chin% comp_map$comparison &
                     perm_fdr_wasserstein < 0.05]
message(sprintf("  Significant genes (fdr < 0.05): %d", nrow(scores_sig)))

# ======================== 2. Venn per group ========================
message("[2/3] Processing Venn overlaps ...")

# Venn plotting function (fixed circle size, no border, left One vs Rest / right cfRNA high)
draw_venn <- function(set1, set2, labels, title, out_pdf) {
  gene_list <- list(set1, set2)
  names(gene_list) <- labels

  fill_colors <- c("#E64B35BF", "#4DBBD5BF")
  futile.logger::flog.threshold(futile.logger::ERROR)

  pdf(out_pdf, width = 7, height = 7)

  venn.plot <- venn.diagram(
    x              = gene_list,
    filename       = NULL,
    category.names = labels,
    fill           = fill_colors,
    alpha          = 0.40,
    lty            = "blank",
    scaled         = FALSE,
    inverted       = FALSE,
    cex            = 1.8,
    cat.cex        = 1.4,
    cat.fontface   = "bold",
    main           = title,
    main.cex       = 1.4,
    margin         = 0.08
  )
  grid::grid.draw(venn.plot)

  dev.off()
}

# Loop over each group
for (i in seq_len(nrow(comp_map))) {

  comp    <- comp_map$comparison[i]
  csv_fn  <- comp_map$overlap_csv[i]
  label   <- comp_map$label[i]

  message(sprintf("\n--- %s ---", label))

  # Get one-vs-rest significant genes
  ovr_genes <- scores_sig[comparison == comp, unique(Gene)]
  message(sprintf("  One-vs-rest significant genes: %d", length(ovr_genes)))

  # Read overlap CSV, extract cfRNA-high genes
  csv_path <- file.path(overlap_dir, csv_fn)
  if (!file.exists(csv_path)) {
    message(sprintf("  Overlap CSV not found: %s, skipping.", csv_path))
    next
  }

  overlap_dt <- fread(csv_path)
  cfRNA_genes <- overlap_dt[cfRNA == TRUE, unique(Gene)]
  message(sprintf("  cfRNA-high genes: %d", length(cfRNA_genes)))

  # Skip groups where either set is empty
  if (length(ovr_genes) == 0 && length(cfRNA_genes) == 0) {
    message("  Both sets empty, skipping.")
    next
  }

  # Build overlap table
  all_genes <- unique(c(ovr_genes, cfRNA_genes))
  if (length(all_genes) == 0) next

  dt <- data.table(
    Gene       = all_genes,
    One_vs_rest = all_genes %chin% ovr_genes,
    cfRNA_high  = all_genes %chin% cfRNA_genes
  )

  dt[, Intersection := ""]
  dt[One_vs_rest == TRUE  & cfRNA_high == FALSE, Intersection := "One_vs_rest only"]
  dt[One_vs_rest == FALSE & cfRNA_high == TRUE,  Intersection := "cfRNA_high only"]
  dt[One_vs_rest == TRUE  & cfRNA_high == TRUE,  Intersection := "Overlap"]

  setorder(dt, Intersection, Gene)

  # Output overlap table
  safe_label <- gsub("[ /]", "_", label)
  table_out <- file.path(out_dir, sprintf("overlap_one_vs_rest_%s.csv", safe_label))
  fwrite(dt, table_out)
  message(sprintf("  Overlap table saved: %s (%d genes)", table_out, nrow(dt)))

  # Print intersection statistics
  intersect_counts <- dt[, .N, by = Intersection]
  for (j in seq_len(nrow(intersect_counts))) {
    message(sprintf("    %-18s : %d", intersect_counts$Intersection[j], intersect_counts$N[j]))
  }

  # Draw Venn diagram
  venn_out <- file.path(out_dir, sprintf("venn_one_vs_rest_%s.pdf", safe_label))
  draw_venn(
    set1    = ovr_genes,
    set2    = cfRNA_genes,
    labels  = c("One vs Rest\n(fdr<0.05)", "cfRNA high"),
    title   = label,
    out_pdf = venn_out
  )
  message(sprintf("  Venn plot saved: %s", venn_out))
}

message("\n[3/3] Done. All outputs in: ", out_dir)