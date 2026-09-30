# ============================================================
# Script: Venn diagrams of DE genes and HPA tissue-specific genes
# Purpose:
#   For each cell line, read cfRNA DEGs and cell_debris DEGs from
#   within-cell-line DE results and the HPA tissue-specific genes
#   for the corresponding tissue, build a three-set Venn diagram
#   (cell+debris, cfRNA, HPA) with and without edge lines, and
#   export detailed overlap tables per cell line.
# ============================================================

###############################################################################
# Venn diagram of DE genes and HPA tissue-specific genes
# Goal: for each cell line, plot a three-set Venn diagram of
#       cfRNA DEGs, cell_debris DEGs and HPA tissue-specific genes,
#       and export detailed overlap tables
###############################################################################

# ======================== 0. Load dependencies ========================
library(data.table)

# Try to load a Venn diagram package, prefer ggVennDiagram, fall back to VennDiagram
if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
  library(ggVennDiagram)
  use_ggvenn <- TRUE
} else if (requireNamespace("VennDiagram", quietly = TRUE)) {
  library(VennDiagram)
  use_ggvenn <- FALSE
} else {
  stop("Please install either 'ggVennDiagram' or 'VennDiagram' package.")
}

# ======================== 1. Configuration ========================
base_dir <- "./02_DE_cell_specific/mlRNA"
hpa_file <- "./02_database/HPA/tissue_enriched_genes.tsv"
out_dir  <- "01_venn_output"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Cell line configuration: R name -> folder name -> HPA tissue
cell_config <- list(
  "Hep3B2.1-7"       = list(folder = "Hep3B2.1-7",       tissue = "liver"),
  "HepG2"             = list(folder = "HepG2",             tissue = "liver"),
  "Hep3B2.1-7_HepG2" = list(folder = "Hep3B2.1-7_HepG2", tissue = "liver"),
  "K562"              = list(folder = "K562",              tissue = "bone marrow"),
  "HTR-8/SVneo"       = list(folder = "HTR-8_SVneo",       tissue = "placenta"),
  "HEK293T"           = list(folder = "HEK293T",           tissue = "adrenal gland")
)

# ======================== 2. Helper functions ========================

#' Read DE file, return genes with Combination == "DEL"
read_deg <- function(file_path) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(character(0))
  }
  dt <- fread(file_path)
      dt[Combination %in% c("DEL", "EL", "DL", "DE"), unique(Gene)]
}

#' Build overlap table
#' @param set_left  genes of left set (cell+debris)
#' @param set_right genes of right set (cfRNA)
#' @param set_hpa   HPA tissue-specific genes
#' @param set_labels labels of the three sets, order: left, right, hpa
build_overlap_table <- function(set_left, set_right, set_hpa, set_labels) {
  all_genes <- unique(c(set_left, set_right, set_hpa))
  if (length(all_genes) == 0) {
    return(data.table(Gene = character(0)))
  }

  dt <- data.table(
    Gene        = all_genes,
    cell_debris = all_genes %in% set_left,
    cfRNA       = all_genes %in% set_right,
    HPA         = all_genes %in% set_hpa
  )

  # Add intersection category labels
  dt[, Intersection := ""]
  dt[cell_debris == TRUE  & cfRNA == FALSE & HPA == FALSE, Intersection := set_labels[1]]
  dt[cell_debris == FALSE & cfRNA == TRUE  & HPA == FALSE, Intersection := set_labels[2]]
  dt[cell_debris == FALSE & cfRNA == FALSE & HPA == TRUE,  Intersection := set_labels[3]]
  dt[cell_debris == TRUE  & cfRNA == TRUE  & HPA == FALSE, Intersection := paste(set_labels[1], set_labels[2], sep = "&")]
  dt[cell_debris == TRUE  & cfRNA == FALSE & HPA == TRUE,  Intersection := paste(set_labels[1], set_labels[3], sep = "&")]
  dt[cell_debris == FALSE & cfRNA == TRUE  & HPA == TRUE,  Intersection := paste(set_labels[2], set_labels[3], sep = "&")]
  dt[cell_debris == TRUE  & cfRNA == TRUE  & HPA == TRUE,  Intersection := "All_three"]

  setorder(dt, Intersection, Gene)
  dt
}

#' Draw and save Venn diagram
draw_venn <- function(set1, set2, set3, labels, title, out_pdf, edge_size = 1.2) {
  gene_list <- list(set1, set2, set3)
  names(gene_list) <- labels

  pdf(out_pdf, width = 7, height = 7)

  if (use_ggvenn) {
    p <- ggVennDiagram(gene_list,
                       label_alpha = 0,
                       edge_size   = edge_size,
                       set_size    = 5) +
      scale_fill_gradient(low = "white", high = "#3C5488FF") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    print(p)
  } else {
    fill_colors <- c("#E64B35BF", "#4DBBD5BF", "#00A087BF")
    futile.logger::flog.threshold(futile.logger::ERROR)
    venn.plot <- venn.diagram(
      x              = gene_list,
      filename       = NULL,
      category.names = labels,
      fill           = fill_colors,
      alpha          = 0.40,
      lty            = if (edge_size == 0) "blank" else "solid",
      cex            = 1.8,
      cat.cex        = 1.4,
      cat.fontface   = "bold",
      main           = title,
      main.cex       = 1.4,
      margin         = 0.08
    )
    grid::grid.draw(venn.plot)
  }

  dev.off()
}

# ======================== 3. Read HPA data ========================
message("[1/3] Reading HPA tissue-enriched genes ...")
hpa <- fread(hpa_file)
message(sprintf("  HPA entries: %d", nrow(hpa)))

# ======================== 4. Process each cell line ========================
message("[2/3] Processing each cell line ...")

for (cell_name in names(cell_config)) {

  cfg     <- cell_config[[cell_name]]
  folder  <- cfg$folder
  tissue  <- cfg$tissue

  message(sprintf("\n--- %s (tissue: %s) ---", cell_name, tissue))

  # Read DE genes
  cf_file <- file.path(base_dir, "cfRNA", "single_cell_line", folder,
                       "gene_overlap_detailed_logFC2_adj0.01.csv")
  cd_file <- file.path(base_dir, "cell_debris", "single_cell_line", folder,
                       "gene_overlap_detailed_logFC2_adj0.01.csv")

  cf_genes <- read_deg(cf_file)
  cd_genes <- read_deg(cd_file)

  # Read HPA tissue-specific genes
  hpa_genes <- hpa[Tissues == tissue, unique(Gene)]

  message(sprintf("  cfRNA DEGs: %d, cell_debris DEGs: %d, HPA (%s): %d",
                  length(cf_genes), length(cd_genes), tissue, length(hpa_genes)))

  # Skip cell lines where all three sets are empty
  if (length(cf_genes) == 0 && length(cd_genes) == 0 && length(hpa_genes) == 0) {
    message("  All three sets empty, skipping.")
    next
  }

  # Three-set labels (left cell+debris, right cfRNA, bottom HPA)
  set_labels <- c("cell+debris", "cfRNA", paste0("HPA\n", tissue))

  # Output overlap table
  overlap_dt <- build_overlap_table(cd_genes, cf_genes, hpa_genes, set_labels)

  table_out <- file.path(out_dir, sprintf("overlap_%s.csv", gsub("/", "_", cell_name)))
  fwrite(overlap_dt, table_out)
  message(sprintf("  Overlap table saved: %s (%d genes)", table_out, nrow(overlap_dt)))

  # Print intersection statistics
  intersect_counts <- overlap_dt[, .N, by = Intersection]
  for (i in seq_len(nrow(intersect_counts))) {
    message(sprintf("    %-30s : %d", intersect_counts$Intersection[i], intersect_counts$N[i]))
  }

  # Draw Venn diagram (with edge)
  venn_out <- file.path(out_dir,
                        sprintf("venn_%s.pdf", gsub("/", "_", cell_name)))
  draw_venn(cd_genes, cf_genes, hpa_genes,
            labels    = set_labels,
            title     = cell_name,
            out_pdf   = venn_out,
            edge_size = 1.2)
  message(sprintf("  Venn plot saved: %s", venn_out))

  # Draw Venn diagram (no edge)
  venn_noedge_out <- file.path(out_dir,
                               sprintf("venn_noedge_%s.pdf", gsub("/", "_", cell_name)))
  draw_venn(cd_genes, cf_genes, hpa_genes,
            labels    = set_labels,
            title     = cell_name,
            out_pdf   = venn_noedge_out,
            edge_size = 0)
  message(sprintf("  Venn plot (no edge) saved: %s", venn_noedge_out))
}

message("\n[3/3] Done. All outputs in: ", out_dir)