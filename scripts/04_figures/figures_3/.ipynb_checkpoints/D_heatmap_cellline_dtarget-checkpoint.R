#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
})

# =============================================================================
# plot_reference_emd_heatmaps.R
#
# 功能:
#   1. 从 one-vs-rest 结果中，分别为 K562 / Placenta / Liver 选取 top20 / top30 代表基因
#   2. 基于参考分布，计算每个样本、每个基因到 target reference 的 EMD 距离
#   3. 输出两类热图:
#      - mean EMD heatmap: 每个样本在三组基因上的平均 EMD 距离
#      - single-gene EMD heatmap: 每个样本在单基因层面的 EMD 距离
#   4. 分别针对:
#      - 仅 cfRNA 样本
#      - 所有通过 QC 的样本 (cfRNA / debris / cell)
#
# 依赖输入:
#   - all_comparisons_gene_scores.tsv
#   - all_cells_reference_length_distributions.tsv.gz
#   - all_samples_length_merged.mlncRNA.txt
#   - cellcult_Sample_Information.txt
#
# 输出:
#   - top20 / top30 的代表基因列表
#   - mean EMD matrix
#   - gene-level EMD matrix
#   - 热图 pdf / png
# =============================================================================


# ---------------------------- 配置区 ----------------------------

# one-vs-rest 结果目录
one_vs_rest_dir <- "/data/work/01_2603cell_culture/07_figure/03_figure3/01_length_distribution_one_vs_rest"

# 输入文件
gene_score_file <- file.path(one_vs_rest_dir, "all_comparisons_gene_scores.tsv")
reference_file  <- file.path(one_vs_rest_dir, "all_cells_reference_length_distributions.tsv.gz")

cell_length_file <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/all_samples_length_merged.mlncRNA.txt"
cell_info_file   <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"

# 输出目录
output_dir <- file.path(one_vs_rest_dir, "04_heatmap")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# feature_set 与 comparison 的映射
target_feature_map <- data.table(
  feature_set = c("K562", "Placenta", "Liver"),
  comparison  = c("K562_vs_rest", "HTR_8_SVneo_vs_rest", "HepG2_Hep3B2.1_7_vs_rest")
)

# 只保留显著基因
fdr_cutoff <- 0.05

# 每个基因在每个样本中最少总 count，低于该阈值不计算 EMD
min_gene_count_per_sample <- 10

# 输出的 topN
top_n_values <- c(20, 30)

# 样本集合
sample_set_list <- list(
  cfRNA_only = list(component_keep = "cfRNA"),
  all_components = list(component_keep = c("cfRNA", "debris", "cell"))
)

# 热图颜色
heat_colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)
mean_heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# 是否对单基因热图按行做 z-score
scale_gene_heatmap_by_row <- TRUE


# ---------------------------- 工具函数 ----------------------------

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

is_qc_pass <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr %in% c("1", "TRUE", "True", "true", "PASS", "Pass", "pass")
}

read_tsv_header_fields <- function(file_path) {
  header_line <- readLines(file_path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (length(header_line) == 0) {
    stop(sprintf("文件为空，无法读取表头: %s", file_path))
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
    if (is.finite(anchor_idx) && anchor_idx < length(clean_names)) {
      length_idx <- seq.int(anchor_idx + 1L, length(clean_names))
    }
  }

  list(
    sample_idx = if (length(sample_idx) == 0 || is.na(sample_idx)) NA_integer_ else sample_idx,
    gene_idx = if (length(gene_idx) == 0 || is.na(gene_idx)) NA_integer_ else gene_idx,
    length_idx = length_idx,
    length_names = if (length(length_idx) > 0) clean_names[length_idx] else character(0)
  )
}

parse_length_positions <- function(length_cols) {
  pos <- suppressWarnings(as.numeric(sub("^len_", "", sub("^X", "", length_cols))))
  if (anyNA(pos)) {
    pos <- seq_along(length_cols)
  }
  pos
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
    stop("p, q, positions 长度不一致。")
  }
  if (sum(p) <= 0 || sum(q) <= 0) {
    return(NA_real_)
  }
  ord <- order(positions)
  p <- as.numeric(p[ord]) / sum(p)
  q <- as.numeric(q[ord]) / sum(q)
  pos <- positions[ord]
  if (length(pos) == 1) return(0)
  cdf_diff <- cumsum(p - q)
  sum(abs(cdf_diff[-length(cdf_diff)]) * diff(pos))
}

make_annotation_colors <- function(annotation_df) {
  out <- list()

  if ("Cell" %in% colnames(annotation_df)) {
    lev <- unique(as.character(annotation_df$Cell))
    pal <- colorRampPalette(brewer.pal(min(8, max(3, length(lev))), "Set2"))(length(lev))
    names(pal) <- lev
    out$Cell <- pal
  }

  if ("Component" %in% colnames(annotation_df)) {
    lev <- unique(as.character(annotation_df$Component))
    base_cols <- c(cfRNA = "#D55E00", debris = "#0072B2", cell = "#009E73")
    pal <- base_cols[lev]
    if (any(is.na(pal))) {
      extra_lev <- lev[is.na(pal)]
      extra_pal <- colorRampPalette(brewer.pal(8, "Dark2"))(length(extra_lev))
      names(extra_pal) <- extra_lev
      pal[is.na(pal)] <- extra_pal[extra_lev]
    }
    names(pal) <- lev
    out$Component <- pal
  }

  if ("feature_set" %in% colnames(annotation_df)) {
    lev <- unique(as.character(annotation_df$feature_set))
    base_cols <- c(K562 = "#E64B35", Placenta = "#4DBBD5", Liver = "#00A087")
    pal <- base_cols[lev]
    if (any(is.na(pal))) {
      extra_lev <- lev[is.na(pal)]
      extra_pal <- colorRampPalette(brewer.pal(8, "Set1"))(length(extra_lev))
      names(extra_pal) <- extra_lev
      pal[is.na(pal)] <- extra_pal[extra_lev]
    }
    names(pal) <- lev
    out$feature_set <- pal
  }

  out
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}


# ---------------------------- Step 1. 读取样本信息 ----------------------------

message(">> 读取样本信息 ...")
sample_info <- fread(cell_info_file)
setnames(sample_info, names(sample_info), trimws(names(sample_info)))

required_info_cols <- c("Cell", "Sample", "Component", "QC")
miss_info_cols <- setdiff(required_info_cols, names(sample_info))
if (length(miss_info_cols) > 0) {
  stop(sprintf("样本信息文件缺少列: %s", paste(miss_info_cols, collapse = ", ")))
}

sample_info[, QC_pass := is_qc_pass(QC)]
sample_info <- sample_info[QC_pass == TRUE]

if (nrow(sample_info) == 0) {
  stop("样本信息文件中没有通过 QC 的样本。")
}

sample_info_out <- copy(sample_info[, .(Sample, Cell, Component, QC)])
fwrite(sample_info_out, file.path(output_dir, "samples_used_qc_pass.tsv"), sep = "\t")


# ---------------------------- Step 2. 读取基因评分并选择 topN ----------------------------

message(">> 读取 one-vs-rest 基因评分 ...")
gene_score_dt <- fread(gene_score_file)

required_score_cols <- c("comparison", "Gene", "wasserstein_specificity_score")
miss_score_cols <- setdiff(required_score_cols, names(gene_score_dt))
if (length(miss_score_cols) > 0) {
  stop(sprintf("gene score 文件缺少列: %s", paste(miss_score_cols, collapse = ", ")))
}

gene_score_dt <- merge(gene_score_dt, target_feature_map, by = "comparison", all.x = FALSE, all.y = FALSE)

# 自动识别显著性列
fdr_col <- NULL
if ("perm_fdr_wasserstein" %in% names(gene_score_dt)) {
  fdr_col <- "perm_fdr_wasserstein"
} else if ("fdr" %in% names(gene_score_dt)) {
  fdr_col <- "fdr"
} else if ("perm_p_wasserstein" %in% names(gene_score_dt)) {
  fdr_col <- "perm_p_wasserstein"
} else {
  stop("gene score 文件中未找到可用于筛选的显著性列（如 perm_fdr_wasserstein / fdr / perm_p_wasserstein）。")
}

gene_score_dt <- gene_score_dt[!is.na(get(fdr_col)) & get(fdr_col) < fdr_cutoff]

if (nrow(gene_score_dt) == 0) {
  stop("没有基因满足显著性筛选条件，请检查 fdr_cutoff 或输入文件。")
}

setorderv(
  gene_score_dt,
  cols = c("feature_set", fdr_col, "wasserstein_specificity_score", "between_wasserstein", "Gene"),
  order = c(1, 1, -1, -1, 1),
  na.last = TRUE
)

selected_gene_list <- list()

for (n_top in top_n_values) {
  top_dt <- gene_score_dt[
    ,
    head(.SD, n_top),
    by = feature_set
  ]
  top_dt[, rank_within_feature := seq_len(.N), by = feature_set]
  selected_gene_list[[paste0("top", n_top)]] <- top_dt

  fwrite(
    top_dt[, .(feature_set, comparison, rank_within_feature, Gene,
               wasserstein_specificity_score,
               between_wasserstein,
               selected_p = get(fdr_col))],
    file.path(output_dir, sprintf("selected_genes_top%d.tsv", n_top)),
    sep = "\t"
  )
}

all_selected_genes <- unique(unlist(lapply(selected_gene_list, function(x) x$Gene)))

message(sprintf("   共纳入目标基因: %d", length(all_selected_genes)))


# ---------------------------- Step 3. 读取 reference 分布 ----------------------------

message(">> 读取参考长度分布 ...")
reference_dt <- fread(reference_file)

required_ref_cols <- c("comparison", "reference_group", "Gene")
miss_ref_cols <- setdiff(required_ref_cols, names(reference_dt))
if (length(miss_ref_cols) > 0) {
  stop(sprintf("reference 文件缺少列: %s", paste(miss_ref_cols, collapse = ", ")))
}

reference_dt <- merge(reference_dt, target_feature_map, by = "comparison", all.x = FALSE, all.y = FALSE)
reference_dt <- reference_dt[reference_group == "target" & Gene %in% all_selected_genes]

ref_len_cols <- grep("^len_", names(reference_dt), value = TRUE)
if (length(ref_len_cols) == 0) {
  stop("reference 文件中未找到 len_ 开头的长度分布列。")
}

reference_dt <- reference_dt[, c("feature_set", "comparison", "Gene", ref_len_cols), with = FALSE]
reference_dt <- unique(reference_dt)

if (nrow(reference_dt) == 0) {
  stop("参考分布表中没有匹配到目标基因。")
}


# ---------------------------- Step 4. 读取长度矩阵并构建样本-基因长度分布 ----------------------------

message(">> 读取长度矩阵并构建样本-基因长度分布 ...")
header_fields <- read_tsv_header_fields(cell_length_file)
layout <- detect_length_file_layout(header_fields)

if (is.na(layout$sample_idx) || is.na(layout$gene_idx) || length(layout$length_idx) == 0) {
  stop("长度矩阵文件列结构无法识别。")
}

select_idx <- c(layout$sample_idx, layout$gene_idx, layout$length_idx)
len_raw <- fread(cell_length_file, sep = "\t", header = FALSE, skip = 1L, select = select_idx)
setnames(len_raw, c("sample", "Gene", layout$length_names))

len_raw <- len_raw[sample %in% sample_info$Sample & Gene %in% all_selected_genes]
if (nrow(len_raw) == 0) {
  stop("长度矩阵中没有匹配到目标样本或目标基因。")
}

for (col in layout$length_names) {
  set(len_raw, j = col, value = as.numeric(len_raw[[col]]))
}

gene_len_dt <- len_raw[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = .(sample, Gene),
  .SDcols = layout$length_names
]

setnames(gene_len_dt, layout$length_names, paste0("len_", parse_length_positions(layout$length_names)))
len_cols <- grep("^len_", names(gene_len_dt), value = TRUE)

gene_len_dt[, total_count := rowSums(.SD, na.rm = TRUE), .SDcols = len_cols]
gene_len_dt <- gene_len_dt[total_count >= min_gene_count_per_sample]

if (nrow(gene_len_dt) == 0) {
  stop("所有样本-基因记录均未通过 min_gene_count_per_sample 过滤。")
}

prop_mat <- normalize_rows(as.matrix(gene_len_dt[, ..len_cols]))
gene_prop_dt <- cbind(gene_len_dt[, .(sample, Gene, total_count)], as.data.table(prop_mat))
setnames(gene_prop_dt, names(gene_prop_dt)[4:ncol(gene_prop_dt)], len_cols)

message(sprintf("   保留样本-基因记录数: %d", nrow(gene_prop_dt)))


# ---------------------------- Step 5. 计算每个样本、每个基因到 reference 的 EMD ----------------------------

message(">> 计算逐基因 EMD 距离 ...")

ref_positions <- parse_length_positions(ref_len_cols)

merged_dt <- merge(
  gene_prop_dt,
  reference_dt,
  by = "Gene",
  suffixes = c("_sample", "_ref"),
  allow.cartesian = TRUE
)

sample_len_cols <- paste0(ref_len_cols, "_sample")
ref_len_cols_merged <- paste0(ref_len_cols, "_ref")

gene_emd_dt <- merged_dt[
  ,
  {
    sample_vec <- as.numeric(.SD[1, ..sample_len_cols])
    ref_vec    <- as.numeric(.SD[1, ..ref_len_cols_merged])
    list(
      D_target = wasserstein_1d(sample_vec, ref_vec, ref_positions),
      total_count = first(total_count)
    )
  },
  by = .(feature_set, comparison, sample, Gene),
  .SDcols = c("total_count", sample_len_cols, ref_len_cols_merged)
]

gene_emd_dt <- merge(
  gene_emd_dt,
  sample_info[, .(Sample, Cell, Component)],
  by.x = "sample",
  by.y = "Sample",
  all.x = TRUE
)

setcolorder(gene_emd_dt, c("sample", "Cell", "Component", "feature_set", "comparison", "Gene", "total_count", "D_target"))

fwrite(gene_emd_dt, file.path(output_dir, "all_selected_gene_level_EMD.tsv.gz"), sep = "\t")

if (nrow(gene_emd_dt) == 0) {
  stop("未计算出任何逐基因 EMD 结果。")
}


# ---------------------------- Step 6. 作图函数 ----------------------------

plot_mean_heatmap <- function(mean_mat,
                              annotation_col,
                              out_prefix,
                              title_text) {
  if (nrow(mean_mat) == 0 || ncol(mean_mat) == 0) return(NULL)

  ann_col <- as.data.frame(annotation_col)
  rownames(ann_col) <- ann_col$sample
  ann_col$sample <- NULL
  ann_col <- ann_col[colnames(mean_mat), , drop = FALSE]

  ann_colors <- make_annotation_colors(ann_col)

  pdf(sprintf("%s.pdf", out_prefix), width = max(10, ncol(mean_mat) * 0.22), height = 4.8)
  pheatmap(
    mean_mat,
    color = mean_heat_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    border_color = NA,
    main = title_text,
    angle_col = 45,
    fontsize_row = 11,
    fontsize_col = 8
  )
  dev.off()

  png(sprintf("%s.png", out_prefix), width = max(1800, ncol(mean_mat) * 38), height = 900, res = 150)
  pheatmap(
    mean_mat,
    color = mean_heat_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    border_color = NA,
    main = title_text,
    angle_col = 45,
    fontsize_row = 11,
    fontsize_col = 8
  )
  dev.off()
}

plot_gene_heatmap <- function(gene_mat,
                              annotation_col,
                              annotation_row,
                              out_prefix,
                              title_text,
                              scale_by_row = TRUE) {
  if (nrow(gene_mat) == 0 || ncol(gene_mat) == 0) return(NULL)

  plot_mat <- gene_mat
  if (scale_by_row) {
    plot_mat <- row_zscore(plot_mat)
  }

  ann_col <- as.data.frame(annotation_col)
  rownames(ann_col) <- ann_col$sample
  ann_col$sample <- NULL
  ann_col <- ann_col[colnames(plot_mat), , drop = FALSE]

  ann_row <- as.data.frame(annotation_row)
  rownames(ann_row) <- ann_row$row_id
  ann_row$row_id <- NULL
  ann_row <- ann_row[rownames(plot_mat), , drop = FALSE]

  ann_colors_col <- make_annotation_colors(ann_col)
  ann_colors_row <- make_annotation_colors(ann_row)
  ann_colors <- c(ann_colors_col, ann_colors_row)

  pdf(sprintf("%s.pdf", out_prefix), width = max(12, ncol(plot_mat) * 0.22), height = max(8, nrow(plot_mat) * 0.16))
  pheatmap(
    plot_mat,
    color = heat_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    border_color = NA,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = title_text,
    angle_col = 45,
    fontsize_row = 7,
    fontsize_col = 8
  )
  dev.off()

  png(sprintf("%s.png", out_prefix), width = max(2200, ncol(plot_mat) * 38), height = max(1600, nrow(plot_mat) * 22), res = 160)
  pheatmap(
    plot_mat,
    color = heat_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    border_color = NA,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = title_text,
    angle_col = 45,
    fontsize_row = 7,
    fontsize_col = 8
  )
  dev.off()
}


# ---------------------------- Step 7. 循环输出 8 套结果 ----------------------------

message(">> 开始输出热图和矩阵 ...")

for (sample_set_name in names(sample_set_list)) {
  sample_cfg <- sample_set_list[[sample_set_name]]
  component_keep <- sample_cfg$component_keep

  sample_sub_info <- sample_info[Component %in% component_keep, .(sample = Sample, Cell, Component)]
  if (nrow(sample_sub_info) == 0) {
    warning(sprintf("[%s] 没有符合条件的样本，跳过。", sample_set_name))
    next
  }

  sample_sub_info <- unique(sample_sub_info)
  sample_sub_info <- sample_sub_info[order(Cell, Component, sample)]

  for (n_top in top_n_values) {
    key_top <- paste0("top", n_top)
    gene_sel_dt <- copy(selected_gene_list[[key_top]])

    out_subdir <- file.path(output_dir, sample_set_name, key_top)
    dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)

    # 保存当前使用的基因列表
    fwrite(
      gene_sel_dt[, .(feature_set, comparison, rank_within_feature, Gene,
                      wasserstein_specificity_score,
                      between_wasserstein,
                      selected_p = get(fdr_col))],
      file.path(out_subdir, sprintf("%s_%s_selected_genes.tsv", sample_set_name, key_top)),
      sep = "\t"
    )

    # 当前样本集合 + 当前 topN 的逐基因 EMD
    sub_gene_emd <- merge(
      gene_emd_dt,
      gene_sel_dt[, .(feature_set, Gene, rank_within_feature)],
      by = c("feature_set", "Gene"),
      all = FALSE
    )
    sub_gene_emd <- sub_gene_emd[sample %in% sample_sub_info$sample]

    if (nrow(sub_gene_emd) == 0) {
      warning(sprintf("[%s - %s] 无逐基因 EMD 数据，跳过。", sample_set_name, key_top))
      next
    }

    # ---------- A. mean EMD matrix ----------
    mean_dt <- sub_gene_emd[
      ,
      .(mean_EMD = mean(D_target, na.rm = TRUE),
        n_gene_used = sum(!is.na(D_target))),
      by = .(feature_set, sample)
    ]

    mean_wide <- dcast(
      mean_dt,
      feature_set ~ sample,
      value.var = "mean_EMD",
      fill = NA_real_
    )

    feature_order <- c("K562", "Placenta", "Liver")
    mean_wide[, feature_set := factor(feature_set, levels = feature_order)]
    setorder(mean_wide, feature_set)

    mean_mat <- as.matrix(mean_wide[, -1, with = FALSE])
    rownames(mean_mat) <- as.character(mean_wide$feature_set)

    # 保证列顺序与注释一致
    common_samples_mean <- intersect(sample_sub_info$sample, colnames(mean_mat))
    mean_mat <- mean_mat[, common_samples_mean, drop = FALSE]
    ann_col_mean <- sample_sub_info[sample %in% common_samples_mean]

    fwrite(
      as.data.table(cbind(feature_set = rownames(mean_mat), mean_mat)),
      file.path(out_subdir, sprintf("%s_%s_mean_EMD_matrix.tsv", sample_set_name, key_top)),
      sep = "\t"
    )

    fwrite(
      mean_dt[order(feature_set, sample)],
      file.path(out_subdir, sprintf("%s_%s_mean_EMD_long.tsv", sample_set_name, key_top)),
      sep = "\t"
    )

    plot_mean_heatmap(
      mean_mat = mean_mat,
      annotation_col = ann_col_mean,
      out_prefix = file.path(out_subdir, sprintf("%s_%s_mean_EMD_heatmap", sample_set_name, key_top)),
      title_text = sprintf("%s | %s | Mean EMD to Reference", sample_set_name, key_top)
    )

    # ---------- B. gene-level EMD matrix ----------
    gene_plot_dt <- copy(sub_gene_emd)
    gene_plot_dt[, row_id := sprintf("%s_%02d_%s", feature_set, rank_within_feature, Gene)]

    gene_wide <- dcast(
      gene_plot_dt,
      row_id + feature_set + rank_within_feature + Gene ~ sample,
      value.var = "D_target",
      fun.aggregate = mean,
      fill = NA_real_
    )

    gene_wide[, feature_set := factor(feature_set, levels = feature_order)]
    setorder(gene_wide, feature_set, rank_within_feature, Gene)

    meta_cols <- c("row_id", "feature_set", "rank_within_feature", "Gene")
    gene_mat <- as.matrix(gene_wide[, !..meta_cols, with = FALSE])
    rownames(gene_mat) <- gene_wide$row_id

    common_samples_gene <- intersect(sample_sub_info$sample, colnames(gene_mat))
    gene_mat <- gene_mat[, common_samples_gene, drop = FALSE]
    ann_col_gene <- sample_sub_info[sample %in% common_samples_gene]

    row_ann <- gene_wide[, .(row_id, feature_set, rank_within_feature, Gene)]

    fwrite(
      gene_wide,
      file.path(out_subdir, sprintf("%s_%s_gene_EMD_matrix.tsv", sample_set_name, key_top)),
      sep = "\t"
    )

    plot_gene_heatmap(
      gene_mat = gene_mat,
      annotation_col = ann_col_gene,
      annotation_row = row_ann,
      out_prefix = file.path(out_subdir, sprintf("%s_%s_gene_EMD_heatmap", sample_set_name, key_top)),
      title_text = sprintf("%s | %s | Gene-level EMD to Reference", sample_set_name, key_top),
      scale_by_row = scale_gene_heatmap_by_row
    )

    # 同时保存当前矩阵用到的样本注释
    fwrite(
      sample_sub_info[order(Cell, Component, sample)],
      file.path(out_subdir, sprintf("%s_%s_sample_annotation.tsv", sample_set_name, key_top)),
      sep = "\t"
    )

    fwrite(
      row_ann,
      file.path(out_subdir, sprintf("%s_%s_gene_row_annotation.tsv", sample_set_name, key_top)),
      sep = "\t"
    )

    message(sprintf("   完成: %s / %s", sample_set_name, key_top))
  }
}

message(">> 全部完成")
message(sprintf("输出目录: %s", output_dir))