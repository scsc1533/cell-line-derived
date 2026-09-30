# ============================================================
# Script: PCA of fragment length distribution
# Purpose:
#   Use fragment length proportion vectors of housekeeping genes
#   (RACK1, GAPDH, RPL23A) per sample as features, read the merged
#   length table and QC-passed sample information, compute
#   per-sample length proportion vectors by gene, and perform PCA
#   for each gene using all components and cfRNA-only samples,
#   coloring points by Cell and shaping by Component, with
#   confidence ellipses for cfRNA-only plots.
# ============================================================

###############################################################################
# Fragment length distribution PCA
# Perform PCA using length proportion vectors of RACK1 and GAPDH per sample as features
###############################################################################

library(ggplot2)
library(dplyr)
library(data.table)

# ======================== File paths ========================
length_file     <- "./all_samples_length_merged.mlncRNA.txt"
sample_info_file <- "./cellcult_Sample_Information.txt"
output_dir      <- "./04_plot_length_PCA"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ======================== Read sample information ========================
sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_info_filtered <- sample_info %>%
  filter(QC == 1) %>%
  mutate(
    Cell      = factor(Cell, levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")),
    Component = factor(Component, levels = c("cfRNA", "debris", "cell"))
  )

# ======================== Read length data ========================
message("[1/4] Reading length distribution data ...")
length_raw <- fread(length_file)

# Columns: first 4 columns are metadata, the rest are length positions 1-100
n_len_cols <- ncol(length_raw) - 4L
setnames(length_raw, seq_along(length_raw),
         c("sample", "Transcript", "Gene", "Type", as.character(seq_len(n_len_cols))))

message(sprintf("  Raw rows: %d, cols: %d", nrow(length_raw), ncol(length_raw)))

# ======================== Filtering and aggregation ========================
message("[2/4] Pre-filtering by samples and target genes ...")

# At the wide stage filter QC samples and RACK1/GAPDH simultaneously (before melt, greatly reduces memory)
target_genes <- c("RACK1", "GAPDH", "RPL23A")
length_filt <- length_raw[sample %in% sample_info_filtered$Sample & Gene %in% target_genes]
message(sprintf("  Filtered rows: %d (was %d)", nrow(length_filt), nrow(length_raw)))

# Merge metadata (use data.table join to avoid type loss caused by dplyr)
length_filt <- merge(length_filt,
                     sample_info_filtered[, c("Sample", "Cell", "Component")],
                     by.x = "sample", by.y = "Sample",
                     all.x = FALSE)

length_cols <- grep("^[0-9]+$", names(length_filt), value = TRUE)

# Convert to long format (data volume already greatly reduced at this point)
length_long <- melt(length_filt,
                    id.vars       = c("sample", "Cell", "Component", "Gene"),
                    measure.vars  = length_cols,
                    variable.name = "Length",
                    value.name    = "Count")
length_long[, Length := as.integer(as.character(Length))]

# Aggregate transcripts by sample x Gene and compute proportion
length_long <- length_long[, .(Count = sum(Count, na.rm = TRUE)),
                           by = .(sample, Cell, Component, Gene, Length)]
length_long[, Proportion := Count / sum(Count), by = .(sample, Gene)]

message(sprintf("  Aggregated rows: %d", nrow(length_long)))
rm(length_raw, length_filt)
gc()

# ======================== PCA plotting ========================
message("[3/4] Running PCA ...")

# PCA + plotting sub-function
run_pca <- function(gene_data, gene_name, subtitle, out_suffix, add_ellipse = FALSE) {

  if (nrow(gene_data) == 0) return(invisible(NULL))

  gene_wide <- dcast(gene_data, sample + Cell + Component ~ Length,
                     value.var = "Proportion", fill = 0)
  sample_cols <- setdiff(names(gene_wide), c("sample", "Cell", "Component"))
  feat_mat <- as.matrix(gene_wide[, ..sample_cols])
  rownames(feat_mat) <- gene_wide$sample

  pca_result <- prcomp(feat_mat, center = TRUE, scale. = FALSE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$sample <- rownames(pca_df)
  pca_df <- pca_df %>% left_join(
    gene_wide[, .(sample, Cell, Component)], by = "sample"
  )
  pca_df$Cell <- factor(pca_df$Cell,
                        levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T"))
  pca_df$Component <- factor(pca_df$Component, levels = c("cfRNA", "debris", "cell"))

  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  title_text <- paste0("PCA of ", gene_name, " Fragment Length Distribution")
  if (subtitle != "") title_text <- paste0(title_text, "\n(", subtitle, ")")

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell, shape = Component)) +
    geom_point(size = 3, stroke = 1.2)

  # Add confidence ellipse only for cfRNA plots
  if (add_ellipse) {
    p <- p + stat_ellipse(aes(group = Cell), level = 0.95,
                          linewidth = 0.6, linetype = "dashed", show.legend = FALSE)
  }

  p <- p +
    scale_shape_manual(values = c("cfRNA"  = 16,
                                  "debris" = 3,
                                  "cell"   = 1),
                       breaks  = c("cfRNA", "debris", "cell")) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = title_text,
         color = "Cell Type",
         shape = "Component") +
    theme_minimal() +
    theme(legend.position = "right")

  print(p)
  out_name <- paste0("PCA_length_", gene_name, out_suffix, ".pdf")
  ggsave(file.path(output_dir, out_name), plot = p, width = 10, height = 8)
  message(sprintf("  Saved: %s (PC1=%.1f%%, PC2=%.1f%%)", out_name, var_explained[1], var_explained[2]))
}

for (g in target_genes) {

  gene_data <- length_long[Gene == g]
  if (nrow(gene_data) == 0) {
    warning("Gene not found: ", g)
    next
  }

  # All samples
  run_pca(gene_data, g, "all components", "")

  # cfRNA only
  gene_cf <- gene_data[Component == "cfRNA"]
  run_pca(gene_cf, g, "cfRNA only", "_cfRNA", add_ellipse = TRUE)
}

message("[4/4] Done. Outputs: PCA_length_RACK1.pdf, PCA_length_GAPDH.pdf")