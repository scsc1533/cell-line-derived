#!/usr/bin/env Rscript

# ============================================================
# Figure 10
# 4-mer motif proportion correlation between
# Cell line x Component groups
#
# Figure requirements:
# 1. Same as the original figure: each group is Cell line | Component
# 2. Draw the full correlation matrix
# 3. Use Spearman correlation coefficient
# 4. Keep row and column clustering trees
# 5. Keep row and column group names
# 6. Display correlation coefficient in each cell
# 7. Use only a white-to-red gradient, no blue
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# ============================================================
# 1. Input and output paths
# ============================================================

qc_file <- paste0(
  "./cellcult_Sample_Information.txt"
)

# 4-mer file already copied to the working directory
motif_file <- "./sample_mRNA_4motif_fre.txt"

outdir <- paste0(
  "./redraw_figures_0717/",
  "10_fig10_group_Spearman_red"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(qc_file)) {
  stop("Sample information file not found: ", qc_file)
}

if (!file.exists(motif_file)) {
  stop("4-mer frequency file not found: ", motif_file)
}

cat("Sample information file: ", qc_file, "\n")
cat("4-mer frequency file: ", motif_file, "\n")
cat("Output directory: ", outdir, "\n")

# ============================================================
# 2. Read sample information
# ============================================================

qc_data <- read.delim(
  qc_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_qc_columns <- c(
  "Sample",
  "Cell",
  "Component",
  "QC"
)

missing_qc_columns <- setdiff(
  required_qc_columns,
  colnames(qc_data)
)

if (length(missing_qc_columns) > 0) {
  stop(
    "Sample information table is missing the following columns: ",
    paste(missing_qc_columns, collapse = ", ")
  )
}

qc_use <- qc_data %>%
  transmute(
    Sample = trimws(as.character(Sample)),
    Cell = trimws(as.character(Cell)),
    Component = trimws(as.character(Component)),
    QC = suppressWarnings(
      as.numeric(as.character(QC))
    )
  ) %>%
  filter(
    QC == 1,
    Component %in% c(
      "cfRNA",
      "debris",
      "cell"
    ),
    !is.na(Sample),
    Sample != "",
    !is.na(Cell),
    Cell != ""
  ) %>%
  distinct(
    Sample,
    .keep_all = TRUE
  )

if (nrow(qc_use) == 0) {
  stop("Number of samples with QC==1 and belonging to target components is 0.")
}

cat("\nNumber of samples after QC filtering: ", nrow(qc_use), "\n")

cat("\nNumber of samples per Cell x Component:\n")
print(
  table(
    qc_use$Cell,
    qc_use$Component,
    useNA = "ifany"
  )
)

# ============================================================
# 3. Read 4-mer frequency matrix
#
# File format:
# First column: Motif
# Subsequent columns: samples
# ============================================================

motif_data <- read.delim(
  motif_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(motif_data) < 3) {
  stop(
    "4-mer file has too few columns, ",
    "should contain at least 1 Motif column and 2 sample columns."
  )
}

motif_colname <- colnames(motif_data)[1]
sample_columns <- colnames(motif_data)[-1]

cat("\nMotif column name: ", motif_colname, "\n")
cat("Number of sample columns in 4-mer file: ", length(sample_columns), "\n")

# ============================================================
# 4. Match QC samples and 4-mer samples
# ============================================================

valid_samples <- intersect(
  qc_use$Sample,
  sample_columns
)

if (length(valid_samples) == 0) {
  stop(
    "Sample columns in 4-mer file do not match QC samples.\n",
    "Please check whether Sample names are consistent."
  )
}

# Keep the order in the QC table
valid_samples <- qc_use$Sample[
  qc_use$Sample %in% valid_samples
]

cat("Number of matched samples: ", length(valid_samples), "\n")

unmatched_qc <- setdiff(
  qc_use$Sample,
  sample_columns
)

if (length(unmatched_qc) > 0) {
  cat(
    "Number of QC samples not found in 4-mer file: ",
    length(unmatched_qc),
    "\n"
  )
}

# ============================================================
# 5. Convert to long format
# ============================================================

motif_filtered <- motif_data %>%
  select(
    all_of(motif_colname),
    all_of(valid_samples)
  ) %>%
  rename(
    Motif = all_of(motif_colname)
  )

motif_long <- motif_filtered %>%
  pivot_longer(
    cols = all_of(valid_samples),
    names_to = "Sample",
    values_to = "Proportion"
  ) %>%
  mutate(
    Sample = trimws(as.character(Sample)),
    Motif = as.character(Motif),
    Proportion = suppressWarnings(
      as.numeric(as.character(Proportion))
    )
  ) %>%
  left_join(
    qc_use %>%
      select(
        Sample,
        Cell,
        Component
      ),
    by = "Sample"
  ) %>%
  filter(
    !is.na(Cell),
    !is.na(Component),
    is.finite(Proportion)
  )

if (nrow(motif_long) == 0) {
  stop("No valid data after converting to long format and merging sample information.")
}

write.table(
  motif_long,
  file = file.path(
    outdir,
    "Fig10_motif_original_long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 6. Compute mean proportion of each Motif for each Cell x Component group
#
# For example:
# HepG2 | cfRNA
# HepG2 | debris
# HepG2 | cell
# ============================================================

group_summary <- motif_long %>%
  group_by(
    Cell,
    Component,
    Motif
  ) %>%
  summarise(
    MeanProportion = mean(
      Proportion,
      na.rm = TRUE
    ),
    SD = sd(
      Proportion,
      na.rm = TRUE
    ),
    sample_n = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(
      Cell,
      Component,
      sep = " | "
    )
  )

write.table(
  group_summary,
  file = file.path(
    outdir,
    "Fig10_motif_mean_by_Cell_Component.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFinal number of Cell x Component groups: ",
  n_distinct(group_summary$Group),
  "\n"
)

cat("\nFinal group names:\n")
print(
  sort(
    unique(group_summary$Group)
  )
)

# ============================================================
# 7. Build Motif x Group matrix
# ============================================================

group_matrix_df <- group_summary %>%
  select(
    Motif,
    Group,
    MeanProportion
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = MeanProportion
  )

group_matrix <- group_matrix_df %>%
  column_to_rownames(
    var = "Motif"
  ) %>%
  as.matrix()

mode(group_matrix) <- "numeric"

# Remove Motifs that are all NA
keep_motif <- rowSums(
  is.finite(group_matrix)
) > 0

group_matrix <- group_matrix[
  keep_motif,
  ,
  drop = FALSE
]

# Remove groups with too few valid Motifs or no variation in values
group_qc <- tibble(
  Group = colnames(group_matrix),

  valid_motif_n = apply(
    group_matrix,
    2,
    function(x) {
      sum(is.finite(x))
    }
  ),

  motif_sd = apply(
    group_matrix,
    2,
    function(x) {
      sd(
        x[is.finite(x)],
        na.rm = TRUE
      )
    }
  )
)

write.table(
  group_qc,
  file = file.path(
    outdir,
    "Fig10_group_matrix_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

valid_groups <- group_qc %>%
  filter(
    valid_motif_n >= 10,
    is.finite(motif_sd),
    motif_sd > 0
  ) %>%
  pull(Group)

if (length(valid_groups) < 2) {
  stop("Fewer than 2 groups passed quality check, cannot compute correlation.")
}

group_matrix <- group_matrix[
  ,
  valid_groups,
  drop = FALSE
]

cat("Number of Motifs used for correlation analysis: ", nrow(group_matrix), "\n")
cat("Number of groups used for correlation analysis: ", ncol(group_matrix), "\n")

write.table(
  group_matrix,
  file = file.path(
    outdir,
    "Fig10_motif_by_group_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ============================================================
# 8. Compute Spearman correlation between groups
#
# Both rows and columns are:
# Cell line | Component
# ============================================================

cor_rho <- cor(
  group_matrix,
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Eliminate tiny floating-point errors
cor_rho <- (
  cor_rho +
    t(cor_rho)
) / 2

diag(cor_rho) <- 1

if (anyNA(cor_rho)) {
  stop(
    "NA present in Spearman correlation matrix, ",
    "please check the number of valid Motifs for some groups."
  )
}

cor_rho[cor_rho > 1] <- 1
cor_rho[cor_rho < -1] <- -1

write.csv(
  cor_rho,
  file = file.path(
    outdir,
    "Fig10_group_Spearman_correlation_matrix.csv"
  ),
  row.names = TRUE
)

cor_long <- as.data.frame(
  as.table(cor_rho),
  stringsAsFactors = FALSE
) %>%
  rename(
    Group1 = Var1,
    Group2 = Var2,
    Spearman_rho = Freq
  )

write.table(
  cor_long,
  file = file.path(
    outdir,
    "Fig10_group_Spearman_correlation_long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nSpearman correlation range: ",
  min(cor_rho, na.rm = TRUE),
  " to ",
  max(cor_rho, na.rm = TRUE),
  "\n"
)

# ============================================================
# 9. Set single red gradient
#
# Low correlation: white or very light red
# High correlation: dark red
# No blue used
# ============================================================

red_palette <- colorRampPalette(
  c(
    "#FFFFFF",
    "#FFF5F0",
    "#FEE0D2",
    "#FCBBA1",
    "#FC9272",
    "#FB6A4A",
    "#EF3B2C",
    "#CB181D",
    "#A50F15",
    "#67000D"
  )
)(100)

# Similar to the reference figure, set 0.2 as the minimum of the legend
# Values below 0.2 are displayed with the 0.2 color, but cell numbers still show the true correlation
color_min <- 0.2
color_max <- 1.0

cor_plot <- cor_rho

cor_plot[cor_plot < color_min] <- color_min
cor_plot[cor_plot > color_max] <- color_max

heatmap_breaks <- seq(
  color_min,
  color_max,
  length.out = length(red_palette) + 1
)

legend_breaks <- seq(
  color_min,
  color_max,
  by = 0.1
)

legend_labels <- sprintf(
  "%.1f",
  legend_breaks
)

# Display the true Spearman correlation in cells, not the truncated value
number_matrix <- matrix(
  sprintf("%.2f", cor_rho),
  nrow = nrow(cor_rho),
  ncol = ncol(cor_rho),
  dimnames = dimnames(cor_rho)
)

# ============================================================
# 10. Use the same clustering order for rows and columns
#
# Use average linkage as in the original figure
# ============================================================

group_hclust <- hclust(
  dist(cor_plot),
  method = "average"
)

n_group <- ncol(cor_plot)

heatmap_width <- max(
  8,
  n_group * 0.58
)

heatmap_height <- max(
  8,
  n_group * 0.58
)

# Avoid generating an overly large image
heatmap_width <- min(
  heatmap_width,
  14
)

heatmap_height <- min(
  heatmap_height,
  14
)

# ============================================================
# 11. Draw PDF
# ============================================================

pdf_file <- file.path(
  outdir,
  "Fig10_group_Spearman_heatmap_red.pdf"
)

pdf(
  pdf_file,
  width = heatmap_width,
  height = heatmap_height
)

pheatmap(
  cor_plot,

  color = red_palette,
  breaks = heatmap_breaks,

  # Use the same clustering tree for rows and columns
  cluster_rows = group_hclust,
  cluster_cols = group_hclust,

  clustering_method = "average",

  # Keep the true correlation coefficients in cells
  display_numbers = number_matrix,
  number_color = "grey30",
  fontsize_number = 7,

  # Keep the full matrix and cell borders
  border_color = "grey75",

  # Keep group names
  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,

  angle_col = 90,

  legend_breaks = legend_breaks,
  legend_labels = legend_labels,

  main = "Correlation of motif proportions between groups (Spearman)"
)

dev.off()

# ============================================================
# 12. Draw PNG
# ============================================================

png_file <- file.path(
  outdir,
  "Fig10_group_Spearman_heatmap_red.png"
)

png(
  png_file,
  width = heatmap_width,
  height = heatmap_height,
  units = "in",
  res = 300,
  bg = "white"
)

pheatmap(
  cor_plot,

  color = red_palette,
  breaks = heatmap_breaks,

  cluster_rows = group_hclust,
  cluster_cols = group_hclust,

  clustering_method = "average",

  display_numbers = number_matrix,
  number_color = "grey30",
  fontsize_number = 7,

  border_color = "grey75",

  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,

  angle_col = 90,

  legend_breaks = legend_breaks,
  legend_labels = legend_labels,

  main = "Correlation of motif proportions between groups (Spearman)"
)

dev.off()

# ============================================================
# 13. Completion information
# ============================================================

cat("\nFigure 10 plotting completed.\n")
cat("PDF: ", pdf_file, "\n")
cat("PNG: ", png_file, "\n")

cat(
  "Correlation matrix: ",
  file.path(
    outdir,
    "Fig10_group_Spearman_correlation_matrix.csv"
  ),
  "\n"
)

cat(
  "Group mean table: ",
  file.path(
    outdir,
    "Fig10_motif_mean_by_Cell_Component.tsv"
  ),
  "\n"
)