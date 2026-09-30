#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# length_distribution_gdm_validation.R
#
# Goal:
#   Use cell-line significant gene sets to validate in plasma samples whether
#   length distributions are closer to the target reference distribution.
#
# Design:
#   1. Candidate gene sets:
#      - K562_vs_rest
#      - HTR_8_SVneo_vs_rest
#      - HepG2_Hep3B2.1_7_vs_rest
#      and require perm_fdr_wasserstein < 0.05
#   2. Reference distribution:
#      - Use the target reference length distribution from one-vs-rest output
#   3. Plasma grouping:
#      - group = 0: Healthy
#      - group = 1: GDM
#   4. Metrics:
#      - Per gene: D_target = Wasserstein(plasma, target_reference)
#      - Relative metric: DeltaD = mean(D_other_refs) - D_target
#      - Sample level: median/mean D_target, median/mean DeltaD
#   5. Test:
#      - GDM vs Healthy, using two-sided Wilcoxon test
#      - Direction output separately, no one-sided hypothesis preset
# =============================================================================

# ---------------------------- Configuration section ----------------------------
base_dir <- "./01_length_distribution_one_vs_rest"
gene_score_file <- file.path(base_dir, "all_comparisons_gene_scores.tsv")
reference_file <- file.path(base_dir, "all_cells_reference_length_distributions.tsv.gz")

plasma_length_file <- "./GDM_samples_length_merged.mlncRNA.txt"
plasma_group_file <- "./GDM_group_all.txt"

output_dir <- file.path(base_dir, "../../04_figure4/02_GDM_validation_SH/all")

selected_comparisons <- c(
  "K562_vs_rest",
  "HTR_8_SVneo_vs_rest",
  "HepG2_Hep3B2.1_7_vs_rest"
)

min_plasma_gene_count <- 10
min_genes_per_sample <- 10


# ---------------------------- Utility functions ----------------------------
parse_length_positions <- function(length_cols) {
  pos <- suppressWarnings(as.numeric(sub("^len_", "", sub("^X", "", length_cols))))
  if (anyNA(pos)) {
    pos <- seq_along(length_cols)
  }
  pos
}

read_tsv_header_fields <- function(file_path) {
  header_line <- readLines(file_path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (length(header_line) == 0) {
    stop(sprintf("File is empty, cannot read header: %s", file_path))
  }
  header_line <- sub("^\ufeff", "", header_line, useBytes = TRUE)
  trimws(strsplit(header_line, "\t", fixed = TRUE)[[1]])
}

detect_length_file_layout <- function(header_fields) {
  clean_names <- trimws(header_fields)
  sample_idx <- which(tolower(clean_names) == "sample")[1]
  gene_idx <- which(tolower(clean_names) == "gene")[1]
  transcript_idx <- which(tolower(clean_names) == "transcript")[1]
  type_idx <- which(tolower(clean_names) == "type")[1]

  numeric_like <- grepl("^X?[0-9]+$", clean_names)
  length_idx <- which(numeric_like)
  if (length(length_idx) == 0) {
    anchor_idx <- max(c(sample_idx, gene_idx, transcript_idx, type_idx), na.rm = TRUE)
    if (is.finite(anchor_idx) && anchor_idx < length(header_fields)) {
      length_idx <- seq.int(anchor_idx + 1L, length(header_fields))
    }
  }

  list(
    sample_idx = if (length(sample_idx) == 0 || is.na(sample_idx)) NA_integer_ else sample_idx,
    gene_idx = if (length(gene_idx) == 0 || is.na(gene_idx)) NA_integer_ else gene_idx,
    length_idx = length_idx,
    length_names = if (length(length_idx) > 0) clean_names[length_idx] else character(0)
  )
}

normalize_rows <- function(mat) {
  rs <- rowSums(mat, na.rm = TRUE)
  out <- mat
  valid <- rs > 0
  if (any(valid)) {
    out[valid, ] <- out[valid, , drop = FALSE] / rs[valid]
  }
  if (any(!valid)) {
    out[!valid, ] <- 0
  }
  out
}

wasserstein_1d <- function(p, q, positions) {
  if (length(p) != length(q) || length(p) != length(positions)) {
    stop("p, q, positions have inconsistent lengths.")
  }
  p <- as.numeric(p)
  q <- as.numeric(q)
  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)
  ord <- order(positions)
  p <- p[ord] / sum(p)
  q <- q[ord] / sum(q)
  pos <- positions[ord]
  if (length(pos) == 1) return(0)
  cdf_diff <- cumsum(p - q)
  sum(abs(cdf_diff[-length(cdf_diff)]) * diff(pos))
}

safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, alternative = "two.sided")$p.value)
}


# ---------------------------- Step 1. Candidate genes and reference distributions ----------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message(">> Reading candidate genes and reference distributions...")
gene_score_dt <- fread(gene_score_file)
candidate_dt <- gene_score_dt[
  comparison %in% selected_comparisons & perm_fdr_wasserstein < 0.05,
  .(
    comparison,
    target_cells,
    Gene,
    wasserstein_specificity_score,
    perm_fdr_wasserstein,
    rank_wasserstein
  )
]

if (nrow(candidate_dt) == 0) {
  stop("No candidate genes passed filtering, please check comparison or perm_fdr_wasserstein threshold.")
}

reference_dt <- fread(reference_file)
reference_dt <- reference_dt[
  comparison %in% selected_comparisons &
    reference_group == "target" &
    Gene %in% unique(candidate_dt$Gene)
]

ref_len_cols <- grep("^len_", names(reference_dt), value = TRUE)
if (length(ref_len_cols) == 0) {
  stop("No length columns with len_ prefix recognized in reference distribution file.")
}
ref_positions <- parse_length_positions(ref_len_cols)

fwrite(candidate_dt, file.path(output_dir, "candidate_gene_sets.tsv"), sep = "\t")


# ---------------------------- Step 2. Plasma sample length distributions ----------------------------
message(">> Reading plasma sample grouping and length distributions...")
group_dt <- fread(plasma_group_file)
setnames(group_dt, names(group_dt), trimws(names(group_dt)))
group_dt <- group_dt[group %in% c(0, 1)]
group_dt[, group_label := fifelse(group == 0, "Healthy", "GDM")]

if (!"sample" %in% names(group_dt)) {
  stop("Grouping file is missing the sample column.")
}

plasma_header_fields <- read_tsv_header_fields(plasma_length_file)
plasma_layout <- detect_length_file_layout(plasma_header_fields)
plasma_len_cols <- plasma_layout$length_names

if (is.na(plasma_layout$sample_idx) || is.na(plasma_layout$gene_idx) || length(plasma_len_cols) == 0) {
  stop("Plasma length file column layout cannot be recognized.")
}

plasma_positions <- parse_length_positions(plasma_len_cols)
if (length(plasma_positions) != length(ref_positions) || any(plasma_positions != ref_positions)) {
  stop("Plasma length columns are inconsistent with reference length columns.")
}

candidate_genes <- unique(candidate_dt$Gene)
plasma_samples <- unique(group_dt$sample)

plasma_select_idx <- c(plasma_layout$sample_idx, plasma_layout$gene_idx, plasma_layout$length_idx)
plasma_dt <- fread(
  plasma_length_file,
  sep = "\t",
  header = FALSE,
  skip = 1L,
  select = plasma_select_idx,
  showProgress = FALSE
)
setnames(plasma_dt, c("sample", "Gene", plasma_len_cols))

plasma_dt <- plasma_dt[sample %in% plasma_samples & Gene %in% candidate_genes]
if (nrow(plasma_dt) == 0) {
  stop("No matched grouping samples and candidate genes in plasma length file.")
}

plasma_gene_count_dt <- plasma_dt[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = .(sample, Gene),
  .SDcols = plasma_len_cols
]
plasma_gene_count_dt[, plasma_total_count := rowSums(.SD, na.rm = TRUE), .SDcols = plasma_len_cols]
plasma_gene_count_dt <- merge(plasma_gene_count_dt, group_dt[, .(sample, group, group_label)], by = "sample")
plasma_gene_count_dt <- plasma_gene_count_dt[plasma_total_count >= min_plasma_gene_count]

if (nrow(plasma_gene_count_dt) == 0) {
  stop("No genes in plasma samples satisfy the minimum length count threshold.")
}

plasma_prop_mat <- as.matrix(plasma_gene_count_dt[, ..plasma_len_cols])
storage.mode(plasma_prop_mat) <- "numeric"
plasma_prop_mat <- normalize_rows(plasma_prop_mat)
plasma_prop_dt <- copy(plasma_gene_count_dt)
plasma_prop_dt[, (plasma_len_cols) := as.data.table(plasma_prop_mat)]


# ---------------------------- Step 3. Compute D_target and DeltaD ----------------------------
message(">> Computing Wasserstein distances between plasma samples and reference distributions...")
distance_long_list <- list()

for (eval_cmp in selected_comparisons) {
  genes_eval <- candidate_dt[comparison == eval_cmp, unique(Gene)]
  plasma_eval <- plasma_prop_dt[Gene %in% genes_eval]
  if (nrow(plasma_eval) == 0) next

  for (ref_cmp in selected_comparisons) {
    ref_sub <- reference_dt[comparison == ref_cmp & Gene %in% genes_eval, c("Gene", ref_len_cols), with = FALSE]
    if (nrow(ref_sub) == 0) next

    merged_dt <- merge(
      plasma_eval[, c("sample", "group", "group_label", "Gene", "plasma_total_count", plasma_len_cols), with = FALSE],
      ref_sub,
      by = "Gene",
      suffixes = c("_plasma", "_ref"),
      all = FALSE
    )
    if (nrow(merged_dt) == 0) next

    plasma_cols_use <- plasma_len_cols
    ref_cols_use <- ref_len_cols
    dist_vals <- vapply(seq_len(nrow(merged_dt)), function(i) {
      p <- as.numeric(merged_dt[i, ..plasma_cols_use])
      q <- as.numeric(merged_dt[i, ..ref_cols_use])
      wasserstein_1d(p, q, positions = plasma_positions)
    }, numeric(1))

    distance_long_list[[paste(eval_cmp, ref_cmp, sep = "__")]] <- data.table(
      comparison = eval_cmp,
      reference_comparison = ref_cmp,
      sample = merged_dt$sample,
      group = merged_dt$group,
      group_label = merged_dt$group_label,
      Gene = merged_dt$Gene,
      plasma_total_count = merged_dt$plasma_total_count,
      D_ref = dist_vals
    )
  }
}

if (length(distance_long_list) == 0) {
  stop("No plasma-reference distance results generated.")
}

distance_long_dt <- rbindlist(distance_long_list, use.names = TRUE, fill = TRUE)

gene_distance_dt <- distance_long_dt[
  ,
  .(
    group = first(group),
    group_label = first(group_label),
    plasma_total_count = first(plasma_total_count),
    D_target = D_ref[reference_comparison == first(comparison)][1],
    mean_D_others = mean(D_ref[reference_comparison != first(comparison)], na.rm = TRUE),
    D_K562 = D_ref[reference_comparison == "K562_vs_rest"][1],
    D_HTR8 = D_ref[reference_comparison == "HTR_8_SVneo_vs_rest"][1],
    D_Liver = D_ref[reference_comparison == "HepG2_Hep3B2.1_7_vs_rest"][1]
  ),
  by = .(comparison, sample, Gene)
]
gene_distance_dt[, DeltaD := mean_D_others - D_target]


# ---------------------------- Step 4. Sample-level summary and group tests ----------------------------
message(">> Summarizing sample-level metrics and comparing GDM vs Healthy...")
sample_score_dt <- gene_distance_dt[
  ,
  .(
    group = first(group),
    group_label = first(group_label),
    n_genes_used = sum(!is.na(D_target)),
    median_D_target = median(D_target, na.rm = TRUE),
    mean_D_target = mean(D_target, na.rm = TRUE),
    median_DeltaD = median(DeltaD, na.rm = TRUE),
    mean_DeltaD = mean(DeltaD, na.rm = TRUE),
    median_D_K562 = median(D_K562, na.rm = TRUE),
    median_D_HTR8 = median(D_HTR8, na.rm = TRUE),
    median_D_Liver = median(D_Liver, na.rm = TRUE)
  ),
  by = .(comparison, sample)
]
sample_score_dt <- sample_score_dt[n_genes_used >= min_genes_per_sample]

metric_map <- c("median_D_target", "mean_D_target", "median_DeltaD", "mean_DeltaD")
sample_test_nested <- lapply(selected_comparisons, function(cmp) {
  sub_dt <- sample_score_dt[comparison == cmp]
  lapply(metric_map, function(metric_name) {
    healthy <- sub_dt[group == 0][[metric_name]]
    gdm <- sub_dt[group == 1][[metric_name]]
    data.table(
      comparison = cmp,
      metric = metric_name,
      n_healthy = sum(!is.na(healthy)),
      n_gdm = sum(!is.na(gdm)),
      healthy_median = median(healthy, na.rm = TRUE),
      gdm_median = median(gdm, na.rm = TRUE),
      gdm_minus_healthy = median(gdm, na.rm = TRUE) - median(healthy, na.rm = TRUE),
      wilcox_p = safe_wilcox(healthy, gdm)
    )
  })
})
sample_test_list <- rbindlist(lapply(sample_test_nested, rbindlist), use.names = TRUE, fill = TRUE)
sample_test_list[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]

gene_test_dt <- gene_distance_dt[
  ,
  .(
    n_healthy = sum(group == 0 & !is.na(D_target)),
    n_gdm = sum(group == 1 & !is.na(D_target)),
    healthy_median_D_target = median(D_target[group == 0], na.rm = TRUE),
    gdm_median_D_target = median(D_target[group == 1], na.rm = TRUE),
    gdm_minus_healthy_D_target = median(D_target[group == 1], na.rm = TRUE) - median(D_target[group == 0], na.rm = TRUE),
    wilcox_p_D_target = safe_wilcox(D_target[group == 0], D_target[group == 1]),
    healthy_median_DeltaD = median(DeltaD[group == 0], na.rm = TRUE),
    gdm_median_DeltaD = median(DeltaD[group == 1], na.rm = TRUE),
    gdm_minus_healthy_DeltaD = median(DeltaD[group == 1], na.rm = TRUE) - median(DeltaD[group == 0], na.rm = TRUE),
    wilcox_p_DeltaD = safe_wilcox(DeltaD[group == 0], DeltaD[group == 1])
  ),
  by = .(comparison, Gene)
]
gene_test_dt[, wilcox_fdr_D_target := p.adjust(wilcox_p_D_target, method = "BH"), by = comparison]
gene_test_dt[, wilcox_fdr_DeltaD := p.adjust(wilcox_p_DeltaD, method = "BH"), by = comparison]


# ---------------------------- Step 5. Output ----------------------------
summary_dt <- data.table(
  selected_comparisons = paste(selected_comparisons, collapse = ";"),
  n_candidate_genes_total = uniqueN(candidate_dt$Gene),
  n_candidate_gene_rows = nrow(candidate_dt),
  n_plasma_samples_total = uniqueN(group_dt$sample),
  n_healthy = uniqueN(group_dt[group == 0, sample]),
  n_gdm = uniqueN(group_dt[group == 1, sample]),
  min_plasma_gene_count = min_plasma_gene_count,
  min_genes_per_sample = min_genes_per_sample
)

fwrite(summary_dt, file.path(output_dir, "run_summary.tsv"), sep = "\t")
fwrite(candidate_dt, file.path(output_dir, "candidate_gene_sets.tsv"), sep = "\t")
fwrite(gene_distance_dt, file.path(output_dir, "gene_level_distances.tsv.gz"), sep = "\t")
fwrite(sample_score_dt, file.path(output_dir, "sample_level_scores.tsv"), sep = "\t")
fwrite(sample_test_list, file.path(output_dir, "sample_level_group_tests.tsv"), sep = "\t")
fwrite(gene_test_dt, file.path(output_dir, "gene_level_group_tests.tsv.gz"), sep = "\t")

message("========================================")
message("GDM vs Healthy plasma length distribution validation completed!")
message(sprintf("Output directory: %s", output_dir))
message("========================================")