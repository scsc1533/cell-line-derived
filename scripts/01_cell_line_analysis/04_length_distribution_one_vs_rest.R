#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# =============================================================================
# length_distribution_one_vs_rest.R
#
# Goal:
#   Step 1. In cell line supernatant cfRNA, construct a length distribution vector for each gene
#   Step 2. Use distribution-level one-vs-rest to screen genes with the most cell-line-specific length distributions
#
# Input:
#   1. Length matrix: all_samples_length_merged.mlncRNA.txt
#      - Row unit is sample + transcript
#      - The same Gene may correspond to multiple Transcripts, needs to be merged to gene level first
#   2. Count matrix: all_mlRNA_counts.txt
#      - Used for gene pre-filtering in Step 1
#   3. Sample information: cellcult_Sample_Information.txt
#      - Keep only samples with Component == cfRNA and QC == 1
#
# Step 1 rules:
#   - For each comparison, separately keep genes with count > 30 in all QC-passing samples of that comparison
#   - Merge transcript length counts to gene level
#   - Normalize each gene/sample length counts to relative abundance distribution
#
# Step 2 rules:
#   - Main metric: Wasserstein distance
#   - Auxiliary metric: JSD
#   - Main ranking: Wasserstein specificity score
#   - Significance: label permutation test
# =============================================================================

# ---------------------------- Configuration section ----------------------------
length_file <- "./all_samples_length_merged.mlncRNA.txt"
count_file <- "./all_mlRNA_counts.txt"
sample_info_file <- "./cellcult_Sample_Information.txt"
output_dir <- "./01_length_distribution_one_vs_rest"

component_filter <- "cfRNA"
qc_keep <- "1"
all_cell_types <- c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")

comparison_list <- list(
  list(name = "Hep3B2.1_7_vs_rest", target = c("Hep3B2.1-7"), rest = NULL),
  list(name = "HepG2_vs_rest", target = c("HepG2"), rest = NULL),
  list(name = "K562_vs_rest", target = c("K562"), rest = NULL),
  list(name = "HTR_8_SVneo_vs_rest", target = c("HTR-8/SVneo"), rest = NULL),
  list(name = "HEK293T_vs_rest", target = c("HEK293T"), rest = NULL),
  list(name = "HepG2_Hep3B2.1_7_vs_rest", target = c("HepG2", "Hep3B2.1-7"), rest = NULL),
  list(name = "HepG2_vs_Hep3B2.1_7", target = c("HepG2"), rest = c("Hep3B2.1-7"))
)

count_threshold <- 30
min_samples_per_group <- 2

n_perm <- 1000
random_seed <- 123
n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
eps <- 1e-12

write_step1_tables <- TRUE


# ---------------------------- Utility functions ----------------------------
sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

safe_filename <- function(x) {
  x <- sanitize_label(x)
  gsub("_+$", "", x)
}

parse_length_positions <- function(length_cols) {
  pos <- suppressWarnings(as.numeric(sub("^X", "", length_cols)))
  if (anyNA(pos)) {
    pos <- seq_along(length_cols)
  }
  pos
}

find_header_col <- function(colnames, target_name) {
  idx <- which(tolower(trimws(colnames)) == tolower(target_name))
  if (length(idx) == 0) {
    return(NA_character_)
  }
  colnames[idx[1]]
}

detect_length_cols <- function(colnames) {
  clean_names <- trimws(colnames)

  numeric_like <- grepl("^X?[0-9]+$", clean_names)
  if (any(numeric_like)) {
    return(colnames[numeric_like])
  }

  meta_idx <- which(tolower(clean_names) %in% c("sample", "transcript", "gene", "type"))
  if (length(meta_idx) > 0 && max(meta_idx) < length(colnames)) {
    return(colnames[(max(meta_idx) + 1):length(colnames)])
  }

  character(0)
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

jsd_distance <- function(p, q, eps = 1e-12) {
  if (length(p) != length(q)) stop("p and q have inconsistent lengths.")
  p <- as.numeric(p)
  q <- as.numeric(q)
  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)
  p <- p / sum(p)
  q <- q / sum(q)
  p <- pmax(p, eps)
  q <- pmax(q, eps)
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
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

mean_distance_to_centroid <- function(mat, centroid, metric = c("wasserstein", "jsd"), positions = NULL) {
  metric <- match.arg(metric)
  if (nrow(mat) < 2) return(NA_real_)
  vals <- apply(mat, 1, function(x) {
    if (metric == "wasserstein") {
      wasserstein_1d(x, centroid, positions = positions)
    } else {
      jsd_distance(x, centroid)
    }
  })
  mean(vals, na.rm = TRUE)
}

compute_group_metrics <- function(mat, labels, positions, eps = 1e-12) {
  target_mat <- mat[labels, , drop = FALSE]
  rest_mat <- mat[!labels, , drop = FALSE]

  if (nrow(target_mat) < 1 || nrow(rest_mat) < 1) {
    return(NULL)
  }

  target_centroid <- colMeans(target_mat)
  rest_centroid <- colMeans(rest_mat)

  between_wasserstein <- wasserstein_1d(target_centroid, rest_centroid, positions)
  between_jsd <- jsd_distance(target_centroid, rest_centroid)

  within_target_wasserstein <- mean_distance_to_centroid(
    target_mat, target_centroid, metric = "wasserstein", positions = positions
  )
  within_rest_wasserstein <- mean_distance_to_centroid(
    rest_mat, rest_centroid, metric = "wasserstein", positions = positions
  )
  within_target_jsd <- mean_distance_to_centroid(
    target_mat, target_centroid, metric = "jsd"
  )
  within_rest_jsd <- mean_distance_to_centroid(
    rest_mat, rest_centroid, metric = "jsd"
  )

  specificity_wasserstein <- between_wasserstein / (within_target_wasserstein + within_rest_wasserstein + eps)
  specificity_jsd <- between_jsd / (within_target_jsd + within_rest_jsd + eps)

  list(
    target_centroid = target_centroid,
    rest_centroid = rest_centroid,
    between_wasserstein = between_wasserstein,
    between_jsd = between_jsd,
    within_target_wasserstein = within_target_wasserstein,
    within_rest_wasserstein = within_rest_wasserstein,
    within_target_jsd = within_target_jsd,
    within_rest_jsd = within_rest_jsd,
    specificity_wasserstein = specificity_wasserstein,
    specificity_jsd = specificity_jsd
  )
}

run_permutation_test <- function(mat, labels, positions, n_perm, seed_offset = 0L, eps = 1e-12) {
  if (n_perm <= 0) {
    return(list(
      perm_p_wasserstein = NA_real_,
      perm_p_jsd = NA_real_,
      perm_mean_score_wasserstein = NA_real_,
      perm_mean_score_jsd = NA_real_
    ))
  }

  obs <- compute_group_metrics(mat, labels, positions, eps = eps)
  if (is.null(obs)) {
    return(list(
      perm_p_wasserstein = NA_real_,
      perm_p_jsd = NA_real_,
      perm_mean_score_wasserstein = NA_real_,
      perm_mean_score_jsd = NA_real_
    ))
  }

  set.seed(random_seed + seed_offset)
  perm_w <- rep(NA_real_, n_perm)
  perm_j <- rep(NA_real_, n_perm)

  for (i in seq_len(n_perm)) {
    perm_labels <- sample(labels, replace = FALSE)
    perm_stat <- compute_group_metrics(mat, perm_labels, positions, eps = eps)
    if (!is.null(perm_stat)) {
      perm_w[i] <- perm_stat$specificity_wasserstein
      perm_j[i] <- perm_stat$specificity_jsd
    }
  }

  valid_w <- !is.na(perm_w)
  valid_j <- !is.na(perm_j)

  list(
    perm_p_wasserstein = if (any(valid_w)) {
      (sum(perm_w[valid_w] >= obs$specificity_wasserstein) + 1) / (sum(valid_w) + 1)
    } else {
      NA_real_
    },
    perm_p_jsd = if (any(valid_j)) {
      (sum(perm_j[valid_j] >= obs$specificity_jsd) + 1) / (sum(valid_j) + 1)
    } else {
      NA_real_
    },
    perm_mean_score_wasserstein = if (any(valid_w)) mean(perm_w[valid_w]) else NA_real_,
    perm_mean_score_jsd = if (any(valid_j)) mean(perm_j[valid_j]) else NA_real_
  )
}

analyze_one_comparison <- function(comparison_name,
                                   target_cells,
                                   gene_levels,
                                   gene_split,
                                   prop_mat,
                                   cell_vec,
                                   sample_vec,
                                   positions,
                                   count_mat,
                                   sample_to_col,
                                   target_samples,
                                   rest_samples,
                                   n_perm,
                                   min_samples_per_group,
                                   n_cores,
                                   eps = 1e-12) {
  worker_fun <- function(i) {
    row_idx <- gene_split[[i]]
    gene_name <- gene_levels[i]

    sub_mat <- prop_mat[row_idx, , drop = FALSE]
    sub_cells <- cell_vec[row_idx]
    sub_samples <- sample_vec[row_idx]

    valid <- rowSums(sub_mat) > 0 & !is.na(sub_cells)
    sub_mat <- sub_mat[valid, , drop = FALSE]
    sub_cells <- sub_cells[valid]
    sub_samples <- sub_samples[valid]

    labels <- sub_cells %in% target_cells
    target_n <- sum(labels)
    rest_n <- sum(!labels)

    if (target_n < min_samples_per_group || rest_n < min_samples_per_group) {
      return(NULL)
    }

    stat <- compute_group_metrics(sub_mat, labels, positions, eps = eps)
    if (is.null(stat)) {
      return(NULL)
    }

    perm <- run_permutation_test(
      mat = sub_mat,
      labels = labels,
      positions = positions,
      n_perm = n_perm,
      seed_offset = i,
      eps = eps
    )

    count_row <- count_mat[gene_name, , drop = TRUE]
    target_cols <- sample_to_col[intersect(target_samples, names(sample_to_col))]
    rest_cols <- sample_to_col[intersect(rest_samples, names(sample_to_col))]

    list(
      result = data.table(
        comparison = comparison_name,
        target_cells = paste(target_cells, collapse = ";"),
        Gene = gene_name,
        n_valid_samples = nrow(sub_mat),
        n_target = target_n,
        n_rest = rest_n,
        target_samples = paste(sub_samples[labels], collapse = ";"),
        rest_samples = paste(sub_samples[!labels], collapse = ";"),
        target_mean_count = mean(count_row[target_cols], na.rm = TRUE),
        rest_mean_count = mean(count_row[rest_cols], na.rm = TRUE),
        target_detect_frac = mean(count_row[target_cols] > count_threshold, na.rm = TRUE),
        rest_detect_frac = mean(count_row[rest_cols] > count_threshold, na.rm = TRUE),
        between_wasserstein = stat$between_wasserstein,
        within_target_wasserstein = stat$within_target_wasserstein,
        within_rest_wasserstein = stat$within_rest_wasserstein,
        wasserstein_specificity_score = stat$specificity_wasserstein,
        between_jsd = stat$between_jsd,
        within_target_jsd = stat$within_target_jsd,
        within_rest_jsd = stat$within_rest_jsd,
        jsd_specificity_score = stat$specificity_jsd,
        perm_p_wasserstein = perm$perm_p_wasserstein,
        perm_p_jsd = perm$perm_p_jsd,
        perm_mean_score_wasserstein = perm$perm_mean_score_wasserstein,
        perm_mean_score_jsd = perm$perm_mean_score_jsd
      ),
      reference = rbind(
        data.table(comparison = comparison_name, target_cells = paste(target_cells, collapse = ";"), reference_group = "target", Gene = gene_name),
        data.table(comparison = comparison_name, target_cells = paste(target_cells, collapse = ";"), reference_group = "rest", Gene = gene_name)
      )
    ) -> out

    ref_mat <- rbind(stat$target_centroid, stat$rest_centroid)
    ref_dt <- as.data.table(ref_mat)
    setnames(ref_dt, names(ref_dt), paste0("len_", positions))
    out$reference <- cbind(out$reference, ref_dt)
    out
  }

  if (.Platform$OS.type != "windows" && n_cores > 1L) {
    res_list <- mclapply(seq_along(gene_levels), worker_fun, mc.cores = n_cores)
  } else {
    res_list <- lapply(seq_along(gene_levels), worker_fun)
  }

  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0) {
    return(list(results = data.table(), references = data.table()))
  }

  result_dt <- rbindlist(lapply(res_list, `[[`, "result"), use.names = TRUE, fill = TRUE)
  reference_dt <- rbindlist(lapply(res_list, `[[`, "reference"), use.names = TRUE, fill = TRUE)

  result_dt[, perm_fdr_wasserstein := p.adjust(perm_p_wasserstein, method = "BH")]
  result_dt[, perm_fdr_jsd := p.adjust(perm_p_jsd, method = "BH")]
  setorder(result_dt, -wasserstein_specificity_score, perm_p_wasserstein, -between_wasserstein, Gene)
  result_dt[, rank_wasserstein := seq_len(.N)]

  list(results = result_dt, references = reference_dt)
}


# ---------------------------- Step 0. Sample filtering ----------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_info <- fread(sample_info_file)
sample_info[, QC := as.character(QC)]
sample_info <- sample_info[Component == component_filter & QC == qc_keep]
sample_info <- sample_info[Cell %in% all_cell_types]

if (nrow(sample_info) == 0) {
  stop("No available samples after filtering, please check Component/QC/Cell settings in sample_info.")
}

message(">> Reading count matrix...")
count_dt <- fread(count_file)
count_gene_col <- names(count_dt)[1]

count_sample_cols <- intersect(sample_info$Sample, names(count_dt))
sample_info <- sample_info[Sample %in% count_sample_cols]
sample_info <- sample_info[match(count_sample_cols, Sample)]

if (length(count_sample_cols) == 0) {
  stop("No QC-passing cfRNA samples matched in count file.")
}

count_filtered <- count_dt[, c(count_gene_col, count_sample_cols), with = FALSE]
setnames(count_filtered, count_gene_col, "Gene")
count_mat <- as.matrix(count_filtered[, ..count_sample_cols])
rownames(count_mat) <- count_filtered$Gene
storage.mode(count_mat) <- "numeric"
sample_to_col <- setNames(seq_along(count_sample_cols), count_sample_cols)

filter_summary <- data.table(
  n_samples_qc_cfRNA = nrow(sample_info),
  n_target_cells = uniqueN(sample_info$Cell),
  n_genes_raw = nrow(count_dt),
  n_genes_in_count_matrix = nrow(count_filtered),
  count_threshold = count_threshold,
  component_filter = component_filter,
  qc_keep = qc_keep,
  n_perm = n_perm
)
fwrite(filter_summary, file.path(output_dir, "step1_filter_summary.tsv"), sep = "\t")
fwrite(sample_info, file.path(output_dir, "samples_used.tsv"), sep = "\t")


# ---------------------------- Step 1. Build gene length distribution ----------------------------
message(">> Reading length matrix and merging transcript -> gene ...")
length_header_fields <- read_tsv_header_fields(length_file)
length_layout <- detect_length_file_layout(length_header_fields)
length_cols <- length_layout$length_names

if (is.na(length_layout$sample_idx) || is.na(length_layout$gene_idx) || length(length_cols) == 0) {
  stop(sprintf(
    "Length file column names cannot be recognized. Detected column names: %s",
    paste(length_header_fields[seq_len(min(length(length_header_fields), 12))], collapse = ", ")
  ))
}

positions <- parse_length_positions(length_cols)

length_select_idx <- c(length_layout$sample_idx, length_layout$gene_idx, length_layout$length_idx)
length_dt <- fread(
  length_file,
  sep = "\t",
  header = FALSE,
  skip = 1L,
  select = length_select_idx,
  showProgress = FALSE
)
setnames(length_dt, c("sample", "Gene", length_cols))

length_dt <- length_dt[sample %in% sample_info$Sample]
if (nrow(length_dt) == 0) {
  stop("No records retained after filtering in length file, please check whether sample matches input.")
}

gene_len_count_dt <- length_dt[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = .(sample, Gene),
  .SDcols = length_cols
]

gene_len_count_dt <- merge(
  gene_len_count_dt,
  sample_info[, .(sample = Sample, Cell)],
  by = "sample",
  all.x = TRUE
)
setcolorder(gene_len_count_dt, c("sample", "Cell", "Gene", length_cols))

count_matrix_step1 <- as.matrix(gene_len_count_dt[, ..length_cols])
storage.mode(count_matrix_step1) <- "numeric"
prop_matrix_step1 <- normalize_rows(count_matrix_step1)

gene_len_prop_dt <- copy(gene_len_count_dt)
gene_len_prop_dt[, (length_cols) := as.data.table(prop_matrix_step1)]

if (write_step1_tables) {
  merged_summary <- data.table(
    n_rows_length_input_after_sample_filter = nrow(length_dt),
    n_gene_sample_rows_after_transcript_merge = nrow(gene_len_count_dt),
    n_unique_genes_after_transcript_merge = uniqueN(gene_len_count_dt$Gene),
    n_unique_samples_after_transcript_merge = uniqueN(gene_len_count_dt$sample)
  )
  fwrite(merged_summary, file.path(output_dir, "step1_merged_length_summary.tsv"), sep = "\t")
}

message(">> transcript -> gene merge completed, starting Step 1 + Step 2 per comparison ...")


# ---------------------------- Step 2. One-vs-Rest ----------------------------
prop_mat <- as.matrix(gene_len_prop_dt[, ..length_cols])
storage.mode(prop_mat) <- "numeric"
gene_vec <- gene_len_prop_dt$Gene
cell_vec <- gene_len_prop_dt$Cell
sample_vec <- gene_len_prop_dt$sample

all_results <- list()
all_references <- list()
run_summary <- list()

for (cmp in comparison_list) {
  comparison_name <- cmp$name
  target_cells <- cmp$target
  rest_cells <- cmp$rest

  if (is.null(rest_cells)) {
    rest_cells <- setdiff(all_cell_types, target_cells)
  }

  comparison_cells <- c(target_cells, rest_cells)
  comparison_samples <- sample_info[Cell %in% comparison_cells, Sample]
  target_samples <- sample_info[Cell %in% target_cells, Sample]
  rest_samples <- sample_info[Cell %in% rest_cells, Sample]

  message(sprintf(">> Analyzing: %s", comparison_name))

  if (length(target_samples) < min_samples_per_group || length(rest_samples) < min_samples_per_group) {
    warning(sprintf(
      "%s has insufficient valid samples, skipping. target=%d, rest=%d",
      comparison_name, length(target_samples), length(rest_samples)
    ))
    next
  }

  comparison_cols <- sample_to_col[comparison_samples]
  keep_genes_cmp <- rowSums(count_mat[, comparison_cols, drop = FALSE] > count_threshold, na.rm = TRUE) == length(comparison_cols)
  filtered_gene_ids <- rownames(count_mat)[keep_genes_cmp]

  if (length(filtered_gene_ids) == 0) {
    warning(sprintf("%s has no genes satisfying count > %d in all samples of this comparison.", comparison_name, count_threshold))
    next
  }

  row_keep <- sample_vec %in% comparison_samples & gene_vec %in% filtered_gene_ids & cell_vec %in% comparison_cells
  prop_mat_cmp <- prop_mat[row_keep, , drop = FALSE]
  gene_vec_cmp <- gene_vec[row_keep]
  cell_vec_cmp <- cell_vec[row_keep]
  sample_vec_cmp <- sample_vec[row_keep]

  gene_split_cmp <- split(seq_len(sum(row_keep)), gene_vec_cmp)
  gene_levels_cmp <- names(gene_split_cmp)

  out_dir_target <- file.path(output_dir, safe_filename(comparison_name))
  dir.create(out_dir_target, showWarnings = FALSE, recursive = TRUE)

  step1_cmp_summary <- data.table(
    comparison = comparison_name,
    target_cells = paste(target_cells, collapse = ";"),
    rest_cells = paste(rest_cells, collapse = ";"),
    n_target_samples = length(target_samples),
    n_rest_samples = length(rest_samples),
    n_comparison_samples = length(comparison_samples),
    n_genes_kept = length(filtered_gene_ids),
    count_threshold = count_threshold,
    gene_filter_rule = sprintf("all_%d_samples_count_gt_%d", length(comparison_samples), count_threshold)
  )
  fwrite(step1_cmp_summary, file.path(out_dir_target, "step1_filter_summary.tsv"), sep = "\t")
  fwrite(data.table(Gene = filtered_gene_ids), file.path(out_dir_target, "step1_filtered_genes.tsv"), sep = "\t")
  fwrite(
    sample_info[Cell %in% comparison_cells],
    file.path(out_dir_target, "samples_used.tsv"),
    sep = "\t"
  )

  if (write_step1_tables) {
    step1_count_out <- copy(gene_len_count_dt[sample %in% comparison_samples & Gene %in% filtered_gene_ids & Cell %in% comparison_cells])
    setnames(step1_count_out, length_cols, paste0("len_", length_cols))
    fwrite(step1_count_out, file.path(out_dir_target, "step1_gene_length_counts.tsv.gz"), sep = "\t")

    step1_prop_out <- copy(gene_len_prop_dt[sample %in% comparison_samples & Gene %in% filtered_gene_ids & Cell %in% comparison_cells])
    setnames(step1_prop_out, length_cols, paste0("len_", length_cols))
    fwrite(step1_prop_out, file.path(out_dir_target, "step1_gene_length_proportions.tsv.gz"), sep = "\t")
  }

  one_target <- analyze_one_comparison(
    comparison_name = comparison_name,
    target_cells = target_cells,
    gene_levels = gene_levels_cmp,
    gene_split = gene_split_cmp,
    prop_mat = prop_mat_cmp,
    cell_vec = cell_vec_cmp,
    sample_vec = sample_vec_cmp,
    positions = positions,
    count_mat = count_mat,
    sample_to_col = sample_to_col,
    target_samples = target_samples,
    rest_samples = rest_samples,
    n_perm = n_perm,
    min_samples_per_group = min_samples_per_group,
    n_cores = n_cores,
    eps = eps
  )

  res_dt <- one_target$results
  ref_dt <- one_target$references

  if (nrow(res_dt) == 0) {
    warning(sprintf("%s did not produce analyzable results.", comparison_name))
    next
  }

  fwrite(res_dt, file.path(out_dir_target, "comparison_gene_scores.tsv"), sep = "\t")
  fwrite(ref_dt, file.path(out_dir_target, "reference_length_distributions.tsv.gz"), sep = "\t")

  summary_dt <- data.table(
    comparison = comparison_name,
    target_cells = paste(target_cells, collapse = ";"),
    rest_cells = paste(rest_cells, collapse = ";"),
    n_target_samples = length(target_samples),
    n_rest_samples = length(rest_samples),
    n_genes_scored = nrow(res_dt),
    n_genes_perm_fdr_wasserstein_lt_0.05 = sum(res_dt$perm_fdr_wasserstein < 0.05, na.rm = TRUE),
    top_gene = res_dt$Gene[1],
    top_score = res_dt$wasserstein_specificity_score[1],
    top_perm_p = res_dt$perm_p_wasserstein[1]
  )
  fwrite(summary_dt, file.path(out_dir_target, "run_summary.tsv"), sep = "\t")

  all_results[[comparison_name]] <- res_dt
  all_references[[comparison_name]] <- ref_dt
  run_summary[[comparison_name]] <- summary_dt
}

if (length(all_results) == 0) {
  stop("All target cell lines produced no results, please check sample counts or input data.")
}

all_result_dt <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
all_reference_dt <- rbindlist(all_references, use.names = TRUE, fill = TRUE)
all_summary_dt <- rbindlist(run_summary, use.names = TRUE, fill = TRUE)

fwrite(all_result_dt, file.path(output_dir, "all_comparisons_gene_scores.tsv"), sep = "\t")
fwrite(all_reference_dt, file.path(output_dir, "all_cells_reference_length_distributions.tsv.gz"), sep = "\t")
fwrite(all_summary_dt, file.path(output_dir, "all_comparisons_run_summary.tsv"), sep = "\t")

message("========================================")
message("Length distribution one-vs-rest analysis completed!")
message(sprintf("Output directory: %s", output_dir))
message("========================================")