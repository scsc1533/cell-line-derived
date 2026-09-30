# ============================================================
# Script: Sample-level 4-mer correlation analysis
# Purpose:
#   Read sample-level terminal 4-mer frequency matrix and sample
#   information, keep QC == 1 cfRNA samples of the five cell
#   lines, compute sample-level Spearman (or Pearson) correlation
#   based on 4-mer composition, draw a correlation heatmap with
#   cell line annotation, perform a sample-level 4-mer PCA, build
#   a pairwise correlation long table (within vs between cell
#   lines), and compare correlation distributions between within-
#   cell-line and between-cell-line pairs.
# ============================================================

###############################################################################
# Part2: sample-level 4-mer correlation analysis
#
# Input:
#   sample_mRNA_4motif_fre.txt
#
# Rows:
#   Motif_4mer, 256 terminal 4-mer motifs
#
# Columns:
#   samples
#
# Sample inclusion:
#   Component == cfRNA
#   QC == 1
###############################################################################

library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ======================== 0. Path ========================

info_file <- "./cellcult_Sample_Information.txt"

motif4_file <- "./sample_mRNA_4motif_fre.txt"

outdir <- "./02_4mer_sample_correlation"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ======================== 1. Parameters ========================

cells <- c(
  "Hep3B2.1-7",
  "HepG2",
  "K562",
  "HTR-8/SVneo",
  "HEK293T"
)

cell_colors <- c(
  "Hep3B2.1-7" = "#E64B35FF",
  "HepG2"      = "#4DBBD5FF",
  "K562"       = "#00A087FF",
  "HTR-8/SVneo"= "#3C5488FF",
  "HEK293T"    = "#F39B7FFF"
)

cor_method <- "spearman"   # can be changed to "pearson"


theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(size = 12),
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    panel.border = element_rect(
      linewidth = 1,
      fill = NA
    )
  )


# ======================== 2. Read sample information ========================

message("[1/6] Reading sample information ...")

info <- fread(info_file)

setnames(
  info,
  trimws(
    sub(
      "\ufeff",
      "",
      names(info),
      fixed = TRUE
    )
  )
)

qc_cfRNA <- info[
  Component == "cfRNA" &
    QC == 1 &
    Cell %in% cells,
  .(
    Sample,
    Cell,
    Component,
    day,
    QC,
    taskID
  )
]

qc_cfRNA[, Cell := factor(Cell, levels = cells)]

fwrite(
  qc_cfRNA,
  file.path(outdir, "01_qc_pass_cfRNA_samples.tsv"),
  sep = "\t"
)

message("QC-passed cfRNA samples:")
print(qc_cfRNA[, .N, by = Cell])


# ======================== 3. Read 4-mer matrix ========================

message("[2/6] Reading 4-mer matrix ...")

motif4 <- fread(motif4_file)

motif_col <- intersect(
  c("Motif_4mer", "motif", "Motif", "X"),
  names(motif4)
)[1]

if (is.na(motif_col)) {
  stop("Cannot find motif column. Expected column name: Motif_4mer")
}

setnames(motif4, motif_col, "Motif_4mer")

sample_cols <- intersect(qc_cfRNA$Sample, names(motif4))

if (length(sample_cols) == 0) {
  stop("No QC cfRNA samples matched 4-mer matrix columns.")
}

missing_samples <- setdiff(qc_cfRNA$Sample, sample_cols)

if (length(missing_samples) > 0) {
  warning(
    "Some QC cfRNA samples are missing in 4-mer matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

motif4_sub <- motif4[
  ,
  c("Motif_4mer", sample_cols),
  with = FALSE
]

for (s in sample_cols) {
  motif4_sub[[s]] <- as.numeric(motif4_sub[[s]])
}

motif4_sub <- motif4_sub[
  complete.cases(motif4_sub[, ..sample_cols])
]

motif_sd <- apply(
  as.matrix(motif4_sub[, ..sample_cols]),
  1,
  sd,
  na.rm = TRUE
)

motif4_sub <- motif4_sub[motif_sd > 0]

message("Motifs used for correlation: ", nrow(motif4_sub))
message("Samples used for correlation: ", length(sample_cols))

fwrite(
  motif4_sub,
  file.path(outdir, "02_4mer_matrix_qc_cfRNA_samples.tsv"),
  sep = "\t"
)


# ======================== 4. Sample-level correlation ========================

message("[3/6] Calculating sample-level correlation ...")

mat_motif_sample <- as.matrix(motif4_sub[, ..sample_cols])

rownames(mat_motif_sample) <- motif4_sub$Motif_4mer

cor_mat <- cor(
  mat_motif_sample,
  method = cor_method,
  use = "pairwise.complete.obs"
)

cor_dt <- as.data.table(cor_mat, keep.rownames = "Sample")

fwrite(
  cor_dt,
  file.path(
    outdir,
    paste0("03_sample_4mer_", cor_method, "_correlation_matrix.tsv")
  ),
  sep = "\t"
)


# ======================== 5. Correlation heatmap ========================

message("[4/6] Drawing correlation heatmap ...")

anno <- qc_cfRNA[
  Sample %in% colnames(cor_mat),
  .(
    Sample,
    Cell
  )
]

anno <- as.data.frame(anno)
rownames(anno) <- anno$Sample
anno$Sample <- NULL

anno <- anno[colnames(cor_mat), , drop = FALSE]

annotation_colors <- list(
  Cell = cell_colors
)

heat_colors <- colorRampPalette(
  c("#2166AC", "white", "#B2182B")
)(100)

pdf(
  file.path(
    outdir,
    paste0("04_sample_4mer_", cor_method, "_correlation_heatmap.pdf")
  ),
  width = 8,
  height = 7
)

pheatmap(
  cor_mat,
  color = heat_colors,
  annotation_col = anno,
  annotation_row = anno,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  border_color = NA,
  main = paste0("Sample-level 4-mer ", cor_method, " correlation")
)

dev.off()


png(
  file.path(
    outdir,
    paste0("04_sample_4mer_", cor_method, "_correlation_heatmap.png")
  ),
  width = 2400,
  height = 2100,
  res = 300
)

pheatmap(
  cor_mat,
  color = heat_colors,
  annotation_col = anno,
  annotation_row = anno,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  border_color = NA,
  main = paste0("Sample-level 4-mer ", cor_method, " correlation")
)

dev.off()

# ======================== 5.1 Sample-level 4-mer PCA ========================

message("[PCA] Drawing sample-level 4-mer PCA ...")

# PCA input matrix: sample x motif
pca_input <- t(mat_motif_sample)

# Remove motifs with zero variance to avoid scale errors
motif_var <- apply(pca_input, 2, var, na.rm = TRUE)
pca_input <- pca_input[, motif_var > 0, drop = FALSE]

pca_res <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)

pca_df <- as.data.table(
  pca_res$x[, 1:2, drop = FALSE]
)

pca_df[, Sample := rownames(pca_res$x)]

pca_df <- merge(
  pca_df,
  qc_cfRNA[, .(Sample, Cell, day)],
  by = "Sample",
  all.x = TRUE
)

pca_df[, Cell := factor(Cell, levels = cells)]

var_explain <- round(
  100 * pca_res$sdev^2 / sum(pca_res$sdev^2),
  2
)

fwrite(
  pca_df,
  file.path(
    outdir,
    "05_sample_4mer_PCA_coordinates.tsv"
  ),
  sep = "\t"
)

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Cell
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  stat_ellipse(
    aes(group = Cell),
    level = 0.95,
    linewidth = 0.6,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = cell_colors
  ) +
  labs(
    x = paste0("PC1 (", var_explain[1], "%)"),
    y = paste0("PC2 (", var_explain[2], "%)"),
    title = "PCA of sample-level cfRNA 4-mer composition",
    color = "Cell line"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      color = "black"
    )
  )

ggsave(
  file.path(
    outdir,
    "05_sample_4mer_PCA.pdf"
  ),
  p_pca,
  width = 7,
  height = 5.5
)

ggsave(
  file.path(
    outdir,
    "05_sample_4mer_PCA.png"
  ),
  p_pca,
  width = 7,
  height = 5.5,
  dpi = 300
)

# Save loadings to inspect which 4-mers contribute most to PC1/PC2
loading_df <- as.data.table(
  pca_res$rotation[, 1:2, drop = FALSE],
  keep.rownames = "Motif_4mer"
)

fwrite(
  loading_df,
  file.path(
    outdir,
    "05_sample_4mer_PCA_loadings_PC1_PC2.tsv"
  ),
  sep = "\t"
)

# ======================== 6. Correlation long table ========================

message("[5/6] Building pairwise correlation table ...")

sample_meta <- qc_cfRNA[
  Sample %in% colnames(cor_mat),
  .(
    Sample,
    Cell,
    day
  )
]

pair_list <- list()
k <- 1

sample_names <- colnames(cor_mat)

for (i in seq_along(sample_names)) {
  
  for (j in seq_along(sample_names)) {
    
    if (i < j) {
      
      s1 <- sample_names[i]
      s2 <- sample_names[j]
      
      c1 <- sample_meta[Sample == s1, Cell]
      c2 <- sample_meta[Sample == s2, Cell]
      
      d1 <- sample_meta[Sample == s1, day]
      d2 <- sample_meta[Sample == s2, day]
      
      pair_type <- ifelse(
        as.character(c1) == as.character(c2),
        "Within cell line",
        "Between cell lines"
      )
      
      pair_label <- ifelse(
        as.character(c1) == as.character(c2),
        as.character(c1),
        paste(
          sort(c(as.character(c1), as.character(c2))),
          collapse = " vs "
        )
      )
      
      pair_list[[k]] <- data.table(
        sample1 = s1,
        sample2 = s2,
        cell1 = as.character(c1),
        cell2 = as.character(c2),
        day1 = d1,
        day2 = d2,
        pair_type = pair_type,
        pair_label = pair_label,
        correlation = cor_mat[s1, s2]
      )
      
      k <- k + 1
    }
  }
}

cor_long <- rbindlist(pair_list)

cor_long[, pair_type := factor(
  pair_type,
  levels = c(
    "Within cell line",
    "Between cell lines"
  )
)]

fwrite(
  cor_long,
  file.path(
    outdir,
    paste0("05_sample_4mer_", cor_method, "_correlation_long.tsv")
  ),
  sep = "\t"
)


# ======================== 7. Within vs between correlation boxplot ========================

message("[6/6] Drawing correlation distribution plots ...")

stat_within_between <- wilcox.test(
  correlation ~ pair_type,
  data = cor_long,
  exact = FALSE
)

stat_dt <- data.table(
  comparison = "Within cell line vs Between cell lines",
  method = "Wilcoxon rank-sum test",
  p_value = stat_within_between$p.value,
  p_value_format = formatC(
    stat_within_between$p.value,
    format = "e",
    digits = 3
  )
)

fwrite(
  stat_dt,
  file.path(
    outdir,
    paste0("06_sample_4mer_", cor_method, "_within_between_wilcoxon.tsv")
  ),
  sep = "\t"
)

p1 <- ggplot(
  cor_long,
  aes(
    x = pair_type,
    y = correlation,
    fill = pair_type
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.6
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.5,
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "Within cell line" = "#E64B35FF",
      "Between cell lines" = "#4DBBD5FF"
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = paste0("4-mer ", cor_method, " correlation"),
    title = "Sample-level 4-mer correlation",
    subtitle = paste0(
      "Wilcoxon P = ",
      stat_dt$p_value_format
    )
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      color = "black"
    )
  )

ggsave(
  file.path(
    outdir,
    paste0("06_sample_4mer_", cor_method, "_within_between_boxplot.pdf")
  ),
  p1,
  width = 5.5,
  height = 5
)

ggsave(
  file.path(
    outdir,
    paste0("06_sample_4mer_", cor_method, "_within_between_boxplot.png")
  ),
  p1,
  width = 5.5,
  height = 5,
  dpi = 300
)


# ======================== 8. Cell-pair correlation boxplot ========================

cor_long[, pair_label := factor(
  pair_label,
  levels = unique(
    cor_long[
      order(pair_type, pair_label),
      pair_label
    ]
  )
)]

p2 <- ggplot(
  cor_long,
  aes(
    x = pair_label,
    y = correlation,
    fill = pair_type
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "Within cell line" = "#E64B35FF",
      "Between cell lines" = "#4DBBD5FF"
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = paste0("4-mer ", cor_method, " correlation"),
    title = "Pairwise sample 4-mer correlation by cell-line pair"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      color = "black",
      size = 8
    )
  )

ggsave(
  file.path(
    outdir,
    paste0("07_sample_4mer_", cor_method, "_cellpair_boxplot.pdf")
  ),
  p2,
  width = 11,
  height = 5.5
)

ggsave(
  file.path(
    outdir,
    paste0("07_sample_4mer_", cor_method, "_cellpair_boxplot.png")
  ),
  p2,
  width = 11,
  height = 5.5,
  dpi = 300
)


# ======================== 9. Summary ========================

summary_dt <- data.table(
  item = c(
    "QC_pass_cfRNA_samples",
    "samples_matched_4mer_matrix",
    "motifs_used",
    "correlation_method",
    "within_between_wilcoxon_p"
  ),
  value = c(
    nrow(qc_cfRNA),
    length(sample_cols),
    nrow(motif4_sub),
    cor_method,
    stat_dt$p_value_format
  )
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    "00_sample_4mer_correlation_summary.tsv"
  ),
  sep = "\t"
)

message("Done.")
message("Output directory: ", outdir)