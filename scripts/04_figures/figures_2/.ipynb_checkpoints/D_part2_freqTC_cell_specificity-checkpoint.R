###############################################################################
# Part2: freqTC cell-line specificity analysis
#
# 1. Sample-level freqTC boxplot by cell line
#    - raw P-value annotation
#    - BH-adjusted P-value annotation
#
# 2. Gene-level freqTC PCA after count > 30 filter
#
# 3. Single-gene freqTC boxplots for RACK1, GAPDH, RPL23A
#    - raw P-value annotation
#    - BH-adjusted P-value annotation
###############################################################################

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

# ======================== 0. 路径 ========================

info_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"

sample_freq_file <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/mRNA_sample_motif_ratios_summary.txt"

gene_freq_file <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/freqTC.mRNA1.txt"

count_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mRNA_Count.txt"

outdir <- "/data/work/01_2603cell_culture/07_figure/02_figure2add/01_part2_freqTC_cell_specificity"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================== 1. 参数 ========================

cells <- c(
  "Hep3B2.1-7",
  "HepG2",
  "K562",
  "HTR-8/SVneo",
  "HEK293T"
)

cell_colors <- c(
  "Hep3B2.1-7" = "#F8766D",
  "HepG2"      = "#A3A500",
  "K562"       = "#00BF7D",
  "HTR-8/SVneo"= "#00B0F6",
  "HEK293T"    = "#E76BF3"
)

target_genes <- c(
  "RACK1",
  "GAPDH",
  "RPL23A"
)

count_cutoff <- 30


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
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.border = element_rect(
      linewidth = 1,
      fill = NA
    )
  )


# ======================== 2. 通用函数 ========================

sig_label <- function(p) {
  
  fcase(
    is.na(p), "",
    p < 0.0001, "****",
    p < 0.001,  "***",
    p < 0.01,   "**",
    p < 0.05,   "*",
    default = "ns"
  )
}


format_p <- function(p) {
  
  ifelse(
    is.na(p),
    NA_character_,
    formatC(
      p,
      format = "e",
      digits = 3
    )
  )
}


pairwise_wilcox <- function(
    dt,
    value_col,
    group_col = "Cell"
) {
  
  groups <- levels(dt[[group_col]])
  groups <- groups[groups %in% unique(dt[[group_col]])]
  
  combs <- combn(
    groups,
    2,
    simplify = FALSE
  )
  
  res <- rbindlist(
    lapply(
      combs,
      function(cp) {
        
        g1 <- cp[1]
        g2 <- cp[2]
        
        v1 <- dt[get(group_col) == g1, get(value_col)]
        v2 <- dt[get(group_col) == g2, get(value_col)]
        
        if (
          length(na.omit(v1)) >= 2 &&
          length(na.omit(v2)) >= 2
        ) {
          
          p <- wilcox.test(
            v1,
            v2,
            exact = FALSE
          )$p.value
          
        } else {
          
          p <- NA_real_
        }
        
        data.table(
          group1 = g1,
          group2 = g2,
          n1 = length(na.omit(v1)),
          n2 = length(na.omit(v2)),
          p_value = p
        )
      }
    )
  )
  
  res[, p_adjusted := p.adjust(
    p_value,
    method = "BH"
  )]
  
  res[, sig_raw := sig_label(p_value)]
  res[, sig_adj := sig_label(p_adjusted)]
  
  res[, p_value_format := format_p(p_value)]
  res[, p_adjusted_format := format_p(p_adjusted)]
  
  return(res)
}


add_sig_annotation <- function(
    p,
    plot_dt,
    stat_dt,
    value_col,
    group_col = "Cell",
    p_col = "p_adjusted",
    sig_col = "sig_adj"
) {
  
  sig_dt <- stat_dt[
    !is.na(get(p_col)) &
      get(p_col) < 0.05
  ]
  
  if (nrow(sig_dt) == 0) {
    return(p)
  }
  
  groups <- levels(plot_dt[[group_col]])
  x_map <- setNames(
    seq_along(groups),
    groups
  )
  
  y_max <- max(
    plot_dt[[value_col]],
    na.rm = TRUE
  )
  
  y_min <- min(
    plot_dt[[value_col]],
    na.rm = TRUE
  )
  
  y_range <- y_max - y_min
  
  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 0.1
  }
  
  sig_dt[, x1 := x_map[group1]]
  sig_dt[, x2 := x_map[group2]]
  sig_dt[, label := get(sig_col)]
  
  sig_dt[, y := y_max + seq_len(.N) * y_range * 0.10]
  sig_dt[, y_text := y + y_range * 0.025]
  
  p +
    geom_segment(
      data = sig_dt,
      aes(
        x = x1,
        xend = x2,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_segment(
      data = sig_dt,
      aes(
        x = x1,
        xend = x1,
        y = y - y_range * 0.015,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_segment(
      data = sig_dt,
      aes(
        x = x2,
        xend = x2,
        y = y - y_range * 0.015,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_text(
      data = sig_dt,
      aes(
        x = (x1 + x2) / 2,
        y = y_text,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold"
    ) +
    coord_cartesian(
      ylim = c(
        y_min,
        max(sig_dt$y_text, na.rm = TRUE) + y_range * 0.05
      ),
      clip = "off"
    )
}


draw_freqTC_boxplot <- function(
    plot_dt,
    title_text,
    outfile_prefix
) {
  
  plot_dt[, Cell := factor(
    Cell,
    levels = cells
  )]
  
  stat_dt <- pairwise_wilcox(
    plot_dt,
    value_col = "freqTC",
    group_col = "Cell"
  )
  
  fwrite(
    stat_dt,
    file.path(
      outdir,
      paste0(
        outfile_prefix,
        "_pairwise_wilcoxon.tsv"
      )
    ),
    sep = "\t"
  )
  
  base_p <- ggplot(
    plot_dt,
    aes(
      x = Cell,
      y = freqTC,
      color = Cell,
      fill = Cell
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0,
      width = 0.6,
      linewidth = 0.6
    ) +
    geom_jitter(
      width = 0.18,
      size = 1.8,
      alpha = 0.75
    ) +
    scale_color_manual(
      values = cell_colors,
      guide = "none"
    ) +
    scale_fill_manual(
      values = cell_colors,
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = "freqTC",
      title = title_text
    ) +
    theme_pub
  
  # 原始 p 值标记
  p_raw <- add_sig_annotation(
    p = base_p,
    plot_dt = plot_dt,
    stat_dt = stat_dt,
    value_col = "freqTC",
    group_col = "Cell",
    p_col = "p_value",
    sig_col = "sig_raw"
  ) +
    labs(
      subtitle = "Wilcoxon rank-sum test, nominal P value"
    )
  
  ggsave(
    file.path(
      outdir,
      paste0(
        outfile_prefix,
        "_rawP.pdf"
      )
    ),
    p_raw,
    width = 6,
    height = 5.5
  )
  
  ggsave(
    file.path(
      outdir,
      paste0(
        outfile_prefix,
        "_rawP.png"
      )
    ),
    p_raw,
    width = 6,
    height = 5.5,
    dpi = 300
  )
  
  # BH 校正后 p 值标记
  p_adj <- add_sig_annotation(
    p = base_p,
    plot_dt = plot_dt,
    stat_dt = stat_dt,
    value_col = "freqTC",
    group_col = "Cell",
    p_col = "p_adjusted",
    sig_col = "sig_adj"
  ) +
    labs(
      subtitle = "Wilcoxon rank-sum test, BH-adjusted P value"
    )
  
  ggsave(
    file.path(
      outdir,
      paste0(
        outfile_prefix,
        "_adjP.pdf"
      )
    ),
    p_adj,
    width =6,
    height = 6
  )
  
  ggsave(
    file.path(
      outdir,
      paste0(
        outfile_prefix,
        "_adjP.png"
      )
    ),
    p_adj,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  return(
    list(
      plot_raw = p_raw,
      plot_adj = p_adj,
      stat = stat_dt
    )
  )
}


# ======================== 3. 读取样本信息，只保留 cfRNA + QC == 1 ========================

message("[1/7] Reading sample information ...")

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

qc_cfRNA[, Cell := factor(
  Cell,
  levels = cells
)]

fwrite(
  qc_cfRNA,
  file.path(
    outdir,
    "01_qc_pass_cfRNA_samples.tsv"
  ),
  sep = "\t"
)

message("QC-passed cfRNA samples:")
print(
  qc_cfRNA[
    ,
    .N,
    by = Cell
  ]
)


# ======================== 4. 样本层面 freqTC 箱型图 ========================

message("[2/7] Sample-level freqTC boxplot ...")

sample_freq <- fread(sample_freq_file)

if (!"sample_name" %in% names(sample_freq)) {
  stop("sample motif file should contain column: sample_name")
}

sample_freq <- sample_freq[
  ,
  .(
    Sample = sample_name,
    freqA = as.numeric(freqA),
    freqTC = as.numeric(freqTC),
    freqTCA = as.numeric(freqTCA)
  )
]

sample_plot_dt <- merge(
  qc_cfRNA,
  sample_freq,
  by = "Sample",
  all.x = FALSE,
  all.y = FALSE
)

sample_plot_dt[, Cell := factor(
  Cell,
  levels = cells
)]

fwrite(
  sample_plot_dt,
  file.path(
    outdir,
    "02_sample_level_freqTC_data.tsv"
  ),
  sep = "\t"
)

message("Matched sample-level freqTC samples:")
print(
  sample_plot_dt[
    ,
    .N,
    by = Cell
  ]
)

draw_freqTC_boxplot(
  sample_plot_dt,
  title_text = "Sample-level cfRNA freqTC across cell lines",
  outfile_prefix = "03_sample_level_freqTC_boxplot"
)


# ======================== 5. 读取 count 矩阵并筛选基因 ========================

message("[3/7] Reading count matrix and filtering genes ...")

counts <- fread(count_file)

gene_col <- intersect(
  c(
    "sym_id",
    "Gene",
    "gene",
    "gene_id"
  ),
  names(counts)
)[1]

if (is.na(gene_col)) {
  stop("Cannot find gene column in count matrix.")
}

setnames(
  counts,
  gene_col,
  "Gene"
)

sample_cols_count <- intersect(
  qc_cfRNA$Sample,
  names(counts)
)

if (length(sample_cols_count) == 0) {
  stop("No QC cfRNA samples matched count matrix columns.")
}

missing_count_samples <- setdiff(
  qc_cfRNA$Sample,
  sample_cols_count
)

if (length(missing_count_samples) > 0) {
  warning(
    "Some QC cfRNA samples are missing in count matrix: ",
    paste(
      missing_count_samples,
      collapse = ", "
    )
  )
}

counts_sub <- counts[
  ,
  c(
    "Gene",
    sample_cols_count
  ),
  with = FALSE
]

counts_sub[
  ,
  keep_count :=
    rowSums(
      .SD > count_cutoff,
      na.rm = TRUE
    ) == length(sample_cols_count),
  .SDcols = sample_cols_count
]

genes_pass_count <- counts_sub[
  keep_count == TRUE,
  unique(Gene)
]

count_summary <- data.table(
  total_genes_in_count = nrow(counts_sub),
  matched_samples_in_count = length(sample_cols_count),
  genes_count_gt_30_all_samples = length(genes_pass_count)
)

fwrite(
  count_summary,
  file.path(
    outdir,
    "04_count_filter_summary.tsv"
  ),
  sep = "\t"
)

writeLines(
  genes_pass_count,
  file.path(
    outdir,
    "04_genes_count_gt30_all_cfRNA_samples.txt"
  )
)

message(
  sprintf(
    "Genes passing count > %s in all matched cfRNA samples: %d",
    count_cutoff,
    length(genes_pass_count)
  )
)


# ======================== 6. 基因层面 freqTC PCA ========================

message("[4/7] Gene-level freqTC PCA ...")

gene_freq <- fread(gene_freq_file)

gene_col2 <- intersect(
  c(
    "Gene",
    "gene",
    "sym_id",
    "gene_id"
  ),
  names(gene_freq)
)[1]

if (is.na(gene_col2)) {
  stop("Cannot find gene column in freqTC matrix.")
}

setnames(
  gene_freq,
  gene_col2,
  "Gene"
)

sample_cols_freq <- intersect(
  sample_cols_count,
  names(gene_freq)
)

if (length(sample_cols_freq) == 0) {
  stop("No QC cfRNA samples matched gene-level freqTC matrix columns.")
}

missing_freq_samples <- setdiff(
  sample_cols_count,
  sample_cols_freq
)

if (length(missing_freq_samples) > 0) {
  warning(
    "Some count-matched samples are missing in freqTC matrix: ",
    paste(
      missing_freq_samples,
      collapse = ", "
    )
  )
}

gene_freq_sub <- gene_freq[
  Gene %chin% genes_pass_count,
  c(
    "Gene",
    sample_cols_freq
  ),
  with = FALSE
]

for (s in sample_cols_freq) {
  gene_freq_sub[[s]] <- as.numeric(gene_freq_sub[[s]])
}

gene_freq_sub <- gene_freq_sub[
  complete.cases(
    gene_freq_sub[
      ,
      ..sample_cols_freq
    ]
  )
]

gene_sd <- apply(
  as.matrix(
    gene_freq_sub[
      ,
      ..sample_cols_freq
    ]
  ),
  1,
  sd,
  na.rm = TRUE
)

gene_freq_sub <- gene_freq_sub[
  gene_sd > 0
]

message(
  sprintf(
    "Genes used for PCA after count/freqTC/variance filter: %d",
    nrow(gene_freq_sub)
  )
)

writeLines(
  gene_freq_sub$Gene,
  file.path(
    outdir,
    "05_genes_used_for_freqTC_PCA.txt"
  )
)

pca_mat <- t(
  as.matrix(
    gene_freq_sub[
      ,
      ..sample_cols_freq
    ]
  )
)

colnames(pca_mat) <- gene_freq_sub$Gene

pca_res <- prcomp(
  pca_mat,
  center = TRUE,
  scale. = TRUE
)

pca_df <- as.data.table(
  pca_res$x[
    ,
    1:2,
    drop = FALSE
  ]
)

pca_df[, Sample := rownames(pca_res$x)]

pca_df <- merge(
  pca_df,
  qc_cfRNA[
    ,
    .(
      Sample,
      Cell,
      day
    )
  ],
  by = "Sample",
  all.x = TRUE
)

pca_df[, Cell := factor(
  Cell,
  levels = cells
)]

var_explain <- round(
  100 * pca_res$sdev^2 / sum(pca_res$sdev^2),
  2
)

fwrite(
  pca_df,
  file.path(
    outdir,
    "06_gene_level_freqTC_PCA_coordinates.tsv"
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
    x = paste0(
      "PC1 (",
      var_explain[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      var_explain[2],
      "%)"
    ),
    title = "PCA of gene-level cfRNA freqTC",
    color = "Cell line"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

ggsave(
  file.path(
    outdir,
    "06_gene_level_freqTC_PCA.pdf"
  ),
  p_pca,
  width = 7,
  height = 5.5
)

ggsave(
  file.path(
    outdir,
    "06_gene_level_freqTC_PCA.png"
  ),
  p_pca,
  width = 7,
  height = 5.5,
  dpi = 300
)


# ======================== 7. 指定基因 freqTC 箱型图 ========================

message("[5/7] Target-gene freqTC boxplots ...")

target_status <- data.table(
  Gene = target_genes,
  in_freqTC_matrix = target_genes %in% gene_freq$Gene,
  in_count_matrix = target_genes %in% counts$Gene,
  pass_count_filter = target_genes %in% genes_pass_count,
  used_for_plot = target_genes %in% gene_freq_sub$Gene
)

fwrite(
  target_status,
  file.path(
    outdir,
    "07_target_gene_status.tsv"
  ),
  sep = "\t"
)

print(target_status)

target_genes_use <- target_status[
  used_for_plot == TRUE,
  Gene
]

if (length(target_genes_use) == 0) {
  
  warning(
    "None of target genes pass filters. Skip target-gene boxplots."
  )
  
} else {
  
  target_freq <- gene_freq_sub[
    Gene %chin% target_genes_use,
    c(
      "Gene",
      sample_cols_freq
    ),
    with = FALSE
  ]
  
  target_long <- melt(
    target_freq,
    id.vars = "Gene",
    measure.vars = sample_cols_freq,
    variable.name = "Sample",
    value.name = "freqTC"
  )
  
  target_long[
    ,
    freqTC := as.numeric(freqTC)
  ]
  
  target_long <- merge(
    target_long,
    qc_cfRNA[
      ,
      .(
        Sample,
        Cell,
        day
      )
    ],
    by = "Sample",
    all.x = TRUE
  )
  
  target_long[
    ,
    Cell := factor(
      Cell,
      levels = cells
    )
  ]
  
  fwrite(
    target_long,
    file.path(
      outdir,
      "08_target_gene_freqTC_long.tsv"
    ),
    sep = "\t"
  )
  
  all_gene_stats <- list()
  
  for (g in target_genes_use) {
    
    gene_dt <- target_long[
      Gene == g
    ]
    
    res <- draw_freqTC_boxplot(
      gene_dt,
      title_text = paste0(
        g,
        " cfRNA freqTC across cell lines"
      ),
      outfile_prefix = paste0(
        "08_target_gene_",
        g,
        "_freqTC_boxplot"
      )
    )
    
    stat_g <- copy(res$stat)
    stat_g[, Gene := g]
    
    all_gene_stats[[g]] <- stat_g
  }
  
  gene_stat_all <- rbindlist(
    all_gene_stats,
    fill = TRUE
  )
  
  setcolorder(
    gene_stat_all,
    c(
      "Gene",
      "group1",
      "group2",
      "n1",
      "n2",
      "p_value",
      "p_adjusted",
      "p_value_format",
      "p_adjusted_format",
      "sig_raw",
      "sig_adj"
    )
  )
  
  fwrite(
    gene_stat_all,
    file.path(
      outdir,
      "08_target_gene_pairwise_wilcoxon_all.tsv"
    ),
    sep = "\t"
  )
  
  p_target_all <- ggplot(
    target_long,
    aes(
      x = Cell,
      y = freqTC,
      color = Cell,
      fill = Cell
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0,
      width = 0.6,
      linewidth = 0.6
    ) +
    geom_jitter(
      width = 0.18,
      size = 1.6,
      alpha = 0.75
    ) +
    facet_wrap(
      ~ Gene,
      scales = "free_y",
      ncol = 3
    ) +
    scale_color_manual(
      values = cell_colors,
      guide = "none"
    ) +
    scale_fill_manual(
      values = cell_colors,
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = "freqTC",
      title = "Gene-level cfRNA freqTC of selected genes"
    ) +
    theme_pub +
    theme(
      strip.background = element_rect(
        fill = "grey95",
        color = "grey80"
      ),
      strip.text = element_text(
        face = "bold",
        size = 11
      )
    )
  
  ggsave(
    file.path(
      outdir,
      "08_target_genes_freqTC_boxplot_facet.pdf"
    ),
    p_target_all,
    width = 12,
    height = 5
  )
  
  ggsave(
    file.path(
      outdir,
      "08_target_genes_freqTC_boxplot_facet.png"
    ),
    p_target_all,
    width = 12,
    height = 5,
    dpi = 300
  )
}


# ======================== 8. 输出样本匹配汇总 ========================

message("[6/7] Writing summary ...")

summary_dt <- data.table(
  item = c(
    "QC_pass_cfRNA_samples",
    "sample_level_freqTC_matched_samples",
    "count_matrix_matched_samples",
    "gene_freqTC_matrix_matched_samples",
    "genes_pass_count_gt30_all_samples",
    "genes_used_for_PCA",
    "target_genes_used_for_plot"
  ),
  value = c(
    nrow(qc_cfRNA),
    nrow(sample_plot_dt),
    length(sample_cols_count),
    length(sample_cols_freq),
    length(genes_pass_count),
    nrow(gene_freq_sub),
    length(target_genes_use)
  )
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    "00_part2_freqTC_analysis_summary.tsv"
  ),
  sep = "\t"
)

message("[7/7] Done.")
message("Output directory: ", outdir)