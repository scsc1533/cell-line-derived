###############################################################################
# Plasma validation for gene-list-based cfRNA fragmentomic features
#
# Groups:
#   Healthy: group_21 == 1 and health_gravida_end == 1
#   HBV:     group_21 == 1 and HBV == 1
#   ICP:     group_21 == 1 and ICP_preICP == 1
#
# Features:
#   1. Length distribution
#   2. Long-fragment proportion >36 bp
#   3. freqTC
###############################################################################

library(data.table)
library(ggplot2)

# ======================== 1. Paths ========================

gene_list_file <- "/data/work/01_2603cell_culture/07_figure/04_figure4/01_plasma_validation_gene_list/gene_list.txt"

clin_file <- "/data/work/01_2603cell_culture/01_raw_date/01_Plasma/SH/SH_cell_line_clin.txt"

length_file <- "/data/work/01_2603cell_culture/01_raw_date/01_Plasma/SH/SH_samples_length_merged.mlncRNA.txt"

freqtc_file <- "/data/work/01_2603cell_culture/01_raw_date/01_Plasma/SH/SH_gene_freqTC_ratio.txt"

outdir <- "/data/work/01_2603cell_culture/07_figure/04_figure4/01_plasma_validation_gene_list"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

single_gene_plot_dir <- file.path(outdir, "single_gene_plots")
dir.create(single_gene_plot_dir, recursive = TRUE, showWarnings = FALSE)


# ======================== 2. Parameters ========================

group_levels <- c("Healthy", "HBV", "ICP")

group_colors <- c(
  "Healthy" = "#4DBBD5FF",
  "HBV"     = "#E64B35FF",
  "ICP"     = "#00A087FF"
)

long_cutoff <- 36
min_gene_sample_total_count <- 10

# "p_adjusted" or "p_value"
p_label_use <- "p_adjusted"

save_each_single_gene_plot <- TRUE


# ======================== 3. Helper functions ========================

clean_names <- function(x) {
  x <- gsub("\ufeff", "", x, fixed = TRUE)
  trimws(x)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

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
  ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 3))
}

read_gene_list <- function(file) {
  x <- fread(file, header = FALSE)
  genes <- unique(trimws(gsub("\"", "", as.character(x[[1]]))))
  genes <- genes[
    !is.na(genes) &
      genes != "" &
      !tolower(genes) %in% c("gene", "gene_name", "symbol", "sym_id")
  ]
  if (length(genes) == 0) stop("No valid genes found in gene_list.txt.")
  genes
}

read_length_matrix_like_HPA <- function(file) {
  
  x <- fread(file)
  
  if (ncol(x) <= 4) {
    stop("Length matrix should contain at least 5 columns.")
  }
  
  # 参考 05_length_liver_HPA.R：
  # 不依赖原始列名，直接按位置强制命名
  n_len_cols <- ncol(x) - 4L
  
  setnames(
    x,
    seq_along(x),
    c(
      "sample",
      "Transcript",
      "Gene",
      "Type",
      as.character(seq_len(n_len_cols))
    )
  )
  
  length_cols <- as.character(seq_len(n_len_cols))
  
  message(
    "Detected length range: ",
    min(as.integer(length_cols)),
    "-",
    max(as.integer(length_cols)),
    " bp"
  )
  
  list(
    data = x,
    length_cols = length_cols
  )
}

pairwise_wilcox <- function(dt, value_col, group_col = "Group") {
  
  dt <- copy(dt)
  dt <- dt[!is.na(get(value_col)) & !is.na(get(group_col))]
  dt[, (group_col) := factor(as.character(get(group_col)), levels = group_levels)]
  
  combs <- combn(group_levels, 2, simplify = FALSE)
  
  res <- rbindlist(
    lapply(combs, function(cp) {
      
      g1 <- cp[1]
      g2 <- cp[2]
      
      v1 <- dt[get(group_col) == g1, get(value_col)]
      v2 <- dt[get(group_col) == g2, get(value_col)]
      
      p <- if (length(na.omit(v1)) >= 2 && length(na.omit(v2)) >= 2) {
        wilcox.test(v1, v2, exact = FALSE)$p.value
      } else {
        NA_real_
      }
      
      data.table(
        group1 = g1,
        group2 = g2,
        n1 = length(na.omit(v1)),
        n2 = length(na.omit(v2)),
        p_value = p
      )
    }),
    fill = TRUE
  )
  
  res[, p_adjusted := p.adjust(p_value, method = "BH")]
  res[, sig_raw := sig_label(p_value)]
  res[, sig_adj := sig_label(p_adjusted)]
  res[, p_value_format := format_p(p_value)]
  res[, p_adjusted_format := format_p(p_adjusted)]
  
  res
}

add_sig_annotation <- function(
    p,
    plot_dt,
    stat_dt,
    value_col,
    p_col = "p_adjusted"
) {
  
  sig_col <- ifelse(
    p_col == "p_value",
    "sig_raw",
    "sig_adj"
  )
  
  # 不再筛选 p < 0.05
  # 所有两两比较均标注，包括 ns
  sig_dt <- copy(stat_dt)
  
  if (nrow(sig_dt) == 0) {
    return(p)
  }
  
  x_map <- setNames(
    seq_along(group_levels),
    group_levels
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
  
  # 有 p 值时显示星号或 ns；无 p 值时显示 NA
  sig_dt[
    ,
    label := ifelse(
      is.na(get(p_col)),
      "NA",
      get(sig_col)
    )
  ]
  
  sig_dt <- sig_dt[
    !is.na(x1) &
      !is.na(x2)
  ]
  
  if (nrow(sig_dt) == 0) {
    return(p)
  }
  
  sig_dt[
    ,
    y := y_max + seq_len(.N) * y_range * 0.12
  ]
  
  sig_dt[
    ,
    y_text := y + y_range * 0.025
  ]
  
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
        y = y - y_range * 0.02,
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
        y = y - y_range * 0.02,
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
      size = 4,
      fontface = "bold"
    ) +
    coord_cartesian(
      ylim = c(
        y_min,
        max(sig_dt$y_text, na.rm = TRUE) + y_range * 0.08
      ),
      clip = "off"
    )
}

theme_box <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.border = element_rect(linewidth = 1, fill = NA)
  )

theme_length <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.border = element_rect(linewidth = 1, fill = NA)
  )

make_boxplot <- function(dt, value_col, ylab, title_text,
                         outfile_prefix = NULL,
                         outdir_use = outdir) {
  
  plot_dt <- copy(dt)
  plot_dt[, plot_value := as.numeric(get(value_col))]
  plot_dt <- plot_dt[!is.na(plot_value) & !is.na(Group)]
  
  if (nrow(plot_dt) == 0) {
    warning("No valid data for: ", title_text)
    return(NULL)
  }
  
  plot_dt[, Group := factor(as.character(Group), levels = group_levels)]
  
  stat_dt <- pairwise_wilcox(plot_dt, value_col = "plot_value")
  
  p <- ggplot(
    plot_dt,
    aes(x = Group, y = plot_value, color = Group, fill = Group)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0,
      width = 0.6,
      linewidth = 0.6
    ) +
    geom_jitter(width = 0.18, size = 1.8, alpha = 0.75) +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_fill_manual(values = group_colors, guide = "none") +
    labs(x = NULL, y = ylab, title = title_text) +
    theme_box
  
  p <- add_sig_annotation(
    p = p,
    plot_dt = plot_dt,
    stat_dt = stat_dt,
    value_col = "plot_value",
    p_col = p_label_use
  )
  
  if (!is.null(outfile_prefix)) {
    
    safe_prefix <- safe_filename(outfile_prefix)
    
    fwrite(
      stat_dt,
      file.path(outdir_use, paste0(safe_prefix, "_pairwise_wilcoxon.tsv")),
      sep = "\t"
    )
    
    ggsave(
      file.path(outdir_use, paste0(safe_prefix, ".pdf")),
      p,
      width = 5.6,
      height = 5.2
    )
    
    ggsave(
      file.path(outdir_use, paste0(safe_prefix, ".png")),
      p,
      width = 5.6,
      height = 5.2,
      dpi = 300
    )
  }
  
  list(plot = p, stat = stat_dt)
}

make_length_plot <- function(dist_dt, title_text) {
  
  plot_dt <- copy(dist_dt)
  plot_dt <- plot_dt[
    Length >= global_length_x_min &
      Length <= global_length_x_max
  ]
  
  if (nrow(plot_dt) == 0) {
    warning("No valid length data for: ", title_text)
    return(NULL)
  }
  
  plot_dt[
    ,
    Group := factor(as.character(Group), levels = group_levels)
  ]
  
  plot_summary <- plot_dt[
    ,
    .(
      MeanFreq = mean(Proportion, na.rm = TRUE)
    ),
    by = .(Group, Length)
  ]
  
  p <- ggplot(
    plot_summary,
    aes(x = Length, y = MeanFreq, color = Group)
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = group_colors, breaks = names(group_colors)) +
    scale_x_continuous(
      limits = c(global_length_x_min, global_length_x_max),
      breaks = function(x) sort(unique(c(pretty(x), long_cutoff)))
    ) +
    labs(
      x = "Read length (bp)",
      y = "Mean frequency",
      title = title_text,
      color = "Group"
    ) +
    theme_length
  
  list(plot = p, summary = plot_summary)
}

save_length_plot <- function(plot_res, outfile_prefix, outdir_use = outdir) {
  
  if (is.null(plot_res)) return(NULL)
  
  safe_prefix <- safe_filename(outfile_prefix)
  
  fwrite(
    plot_res$summary,
    file.path(outdir_use, paste0(safe_prefix, "_summary.tsv")),
    sep = "\t"
  )
  
  ggsave(
    file.path(outdir_use, paste0(safe_prefix, ".pdf")),
    plot_res$plot,
    width = 8,
    height = 5
  )
  
  ggsave(
    file.path(outdir_use, paste0(safe_prefix, ".png")),
    plot_res$plot,
    width = 8,
    height = 5,
    dpi = 300
  )
}


# ======================== 4. Read gene list ========================

message("[1/8] Reading gene list ...")

gene_list <- read_gene_list(gene_list_file)

fwrite(
  data.table(Gene = gene_list),
  file.path(outdir, "01_gene_list_used_for_plasma_validation.tsv"),
  sep = "\t"
)

message("Genes in list: ", length(gene_list))


# ======================== 5. Read clinical file and define groups ========================

message("[2/8] Reading clinical file ...")

clin <- fread(clin_file)
setnames(clin, clean_names(names(clin)))

required_clin_cols <- c(
  "sample",
  "group_21",
  "health_gravida_end",
  "HBV",
  "ICP_preICP"
)

missing_cols <- setdiff(required_clin_cols, names(clin))

if (length(missing_cols) > 0) {
  stop("Missing columns in clinical file: ", paste(missing_cols, collapse = ", "))
}

clin[, group_21_num := to_num(group_21)]
clin[, health_num := to_num(health_gravida_end)]
clin[, HBV_num := to_num(HBV)]
clin[, ICP_num := to_num(ICP_preICP)]

clin_use <- clin[group_21_num == 1]

clin_use[
  ,
  group_hit_n :=
    as.integer(!is.na(health_num) & health_num == 1) +
    as.integer(!is.na(HBV_num) & HBV_num == 1) +
    as.integer(!is.na(ICP_num) & ICP_num == 1)
]

clin_use[
  ,
  Group := fcase(
    !is.na(health_num) & health_num == 1, "Healthy",
    !is.na(HBV_num) & HBV_num == 1, "HBV",
    !is.na(ICP_num) & ICP_num == 1, "ICP",
    default = NA_character_
  )
]

excluded_samples <- clin_use[group_hit_n != 1]

if (nrow(excluded_samples) > 0) {
  fwrite(
    excluded_samples,
    file.path(outdir, "02_excluded_ambiguous_or_unclassified_samples.tsv"),
    sep = "\t"
  )
}

optional_cols <- intersect(
  c("Sampling_GW", "ALT", "AST", "Alb", "TBA"),
  names(clin_use)
)

sample_meta <- clin_use[
  group_hit_n == 1 & Group %in% group_levels,
  c(
    "sample",
    "Group",
    "group_21",
    "health_gravida_end",
    "HBV",
    "ICP_preICP",
    optional_cols
  ),
  with = FALSE
]

sample_meta[
  ,
  Group := factor(as.character(Group), levels = group_levels)
]

fwrite(
  sample_meta,
  file.path(outdir, "02_selected_plasma_samples.tsv"),
  sep = "\t"
)

message("Selected plasma samples:")
print(sample_meta[, .N, by = Group])


# ======================== 6. Read length matrix ========================

message("[3/8] Reading length matrix ...")

length_obj <- read_length_matrix_like_HPA(length_file)

length_raw <- length_obj$data
length_cols <- length_obj$length_cols

length_sub <- length_raw[
  sample %chin% sample_meta$sample &
    Gene %chin% gene_list
]

if (nrow(length_sub) == 0) {
  stop("No matched rows in length matrix for selected samples and genes.")
}

matched_length_genes <- sort(unique(length_sub$Gene))
missing_length_genes <- setdiff(gene_list, matched_length_genes)

fwrite(
  data.table(Gene = matched_length_genes),
  file.path(outdir, "03_genes_matched_in_length_matrix.tsv"),
  sep = "\t"
)

fwrite(
  data.table(Gene = missing_length_genes),
  file.path(outdir, "03_genes_missing_in_length_matrix.tsv"),
  sep = "\t"
)

message("Genes matched in length matrix: ", length(matched_length_genes))
message("Genes missing in length matrix: ", length(missing_length_genes))

# 同一样本、同一基因的多个 transcript 在每个长度上求和
gene_length_wide <- length_sub[
  ,
  lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
  by = .(sample, Gene),
  .SDcols = length_cols
]

gene_length_wide <- merge(
  gene_length_wide,
  sample_meta[, .(sample, Group)],
  by = "sample",
  all.x = TRUE
)

fwrite(
  gene_length_wide,
  file.path(outdir, "03_gene_sample_length_count_wide.tsv"),
  sep = "\t"
)

length_long <- melt(
  gene_length_wide,
  id.vars = c("sample", "Gene", "Group"),
  measure.vars = length_cols,
  variable.name = "Length",
  value.name = "Count"
)

length_long[, Length := as.integer(as.character(Length))]
length_long[, Count := as.numeric(Count)]

fwrite(
  length_long,
  file.path(outdir, "03_gene_sample_length_count_long.tsv"),
  sep = "\t"
)

# 全局 x 轴范围：基于所有入选基因和所有入选血浆样本
global_length_count <- length_long[
  ,
  .(
    total_count = sum(Count, na.rm = TRUE)
  ),
  by = Length
][order(Length)]

nonzero_length <- global_length_count[
  total_count > 0,
  Length
]

if (length(nonzero_length) == 0) {
  stop("No non-zero length found across all selected genes and samples.")
}

global_length_x_min <- min(nonzero_length)
global_length_x_max <- max(nonzero_length)

fwrite(
  data.table(
    x_min = global_length_x_min,
    x_max = global_length_x_max,
    nonzero_length_n = length(nonzero_length),
    note = "Calculated using all selected genes and all selected plasma samples"
  ),
  file.path(outdir, "03_global_length_distribution_x_axis_range.tsv"),
  sep = "\t"
)

message(
  "Global length x-axis range: ",
  global_length_x_min,
  "-",
  global_length_x_max,
  " bp"
)


# ======================== 7. Gene-set length distribution and Prop >36 ========================

message("[4/8] Gene-set length features ...")

gene_set_length_by_sample <- length_long[
  ,
  .(
    Count = sum(Count, na.rm = TRUE)
  ),
  by = .(sample, Group, Length)
]

gene_set_length_by_sample[
  ,
  TotalCount := sum(Count, na.rm = TRUE),
  by = sample
]

gene_set_length_by_sample[
  ,
  Proportion := ifelse(TotalCount > 0, Count / TotalCount, NA_real_)
]

fwrite(
  gene_set_length_by_sample,
  file.path(outdir, "04_gene_set_length_distribution_by_sample.tsv"),
  sep = "\t"
)

p_gene_set_len <- make_length_plot(
  gene_set_length_by_sample,
  "Gene-set length distribution in plasma"
)

save_length_plot(
  p_gene_set_len,
  "04_gene_set_length_distribution"
)

gene_set_prop36 <- gene_set_length_by_sample[
  ,
  .(
    TotalCount = sum(Count, na.rm = TRUE),
    LongCount = sum(Count[Length > long_cutoff], na.rm = TRUE)
  ),
  by = .(sample, Group)
]

gene_set_prop36[
  ,
  Prop_gt36 := ifelse(TotalCount > 0, LongCount / TotalCount, NA_real_)
]

fwrite(
  gene_set_prop36,
  file.path(outdir, "05_gene_set_prop_gt36_by_sample.tsv"),
  sep = "\t"
)

make_boxplot(
  gene_set_prop36,
  value_col = "Prop_gt36",
  ylab = "Proportion of fragments >36 bp",
  title_text = "Gene-set long-fragment proportion in plasma",
  outfile_prefix = "05_gene_set_prop_gt36_boxplot"
)


# ======================== 8. Single-gene length distribution and Prop >36 ========================

message("[5/8] Single-gene length features ...")

single_gene_length_by_sample <- copy(length_long)

single_gene_length_by_sample[
  ,
  TotalCount := sum(Count, na.rm = TRUE),
  by = .(Gene, sample)
]

single_gene_length_by_sample[
  ,
  Proportion := ifelse(TotalCount > 0, Count / TotalCount, NA_real_)
]

fwrite(
  single_gene_length_by_sample,
  file.path(outdir, "06_single_gene_length_distribution_by_sample.tsv"),
  sep = "\t"
)

single_gene_prop36 <- single_gene_length_by_sample[
  ,
  .(
    TotalCount = sum(Count, na.rm = TRUE),
    LongCount = sum(Count[Length > long_cutoff], na.rm = TRUE)
  ),
  by = .(Gene, sample, Group)
]

single_gene_prop36[
  ,
  Prop_gt36 := ifelse(
    TotalCount >= min_gene_sample_total_count,
    LongCount / TotalCount,
    NA_real_
  )
]

fwrite(
  single_gene_prop36,
  file.path(outdir, "07_single_gene_prop_gt36_by_sample.tsv"),
  sep = "\t"
)

prop_stat_list <- list()

pdf(
  file.path(outdir, "06_single_gene_length_distribution_multipage.pdf"),
  width = 8,
  height = 5
)

for (g in sort(unique(single_gene_length_by_sample$Gene))) {
  
  tmp_len <- single_gene_length_by_sample[
    Gene == g &
      TotalCount >= min_gene_sample_total_count
  ]
  
  if (nrow(tmp_len) == 0) next
  
  p_len <- make_length_plot(
    tmp_len,
    paste0(g, " length distribution in plasma")
  )
  
  if (is.null(p_len)) next
  
  print(p_len$plot)
  
  if (save_each_single_gene_plot) {
    
    g_safe <- safe_filename(g)
    
    ggsave(
      file.path(single_gene_plot_dir, paste0(g_safe, "_length_distribution.pdf")),
      p_len$plot,
      width = 8,
      height = 5
    )
    
    ggsave(
      file.path(single_gene_plot_dir, paste0(g_safe, "_length_distribution.png")),
      p_len$plot,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

dev.off()

pdf(
  file.path(outdir, "07_single_gene_prop_gt36_boxplot_multipage.pdf"),
  width = 5.6,
  height = 5.2
)

for (g in sort(unique(single_gene_prop36$Gene))) {
  
  tmp <- single_gene_prop36[
    Gene == g &
      !is.na(Prop_gt36)
  ]
  
  if (nrow(tmp) == 0) next
  
  res <- make_boxplot(
    tmp,
    value_col = "Prop_gt36",
    ylab = "Proportion of fragments >36 bp",
    title_text = paste0(g, " long-fragment proportion"),
    outfile_prefix = NULL
  )
  
  if (is.null(res)) next
  
  print(res$plot)
  
  stat_g <- copy(res$stat)
  stat_g[, Gene := g]
  prop_stat_list[[g]] <- stat_g
  
  if (save_each_single_gene_plot) {
    
    g_safe <- safe_filename(g)
    
    ggsave(
      file.path(single_gene_plot_dir, paste0(g_safe, "_prop_gt36_boxplot.pdf")),
      res$plot,
      width = 5.6,
      height = 5.2
    )
    
    ggsave(
      file.path(single_gene_plot_dir, paste0(g_safe, "_prop_gt36_boxplot.png")),
      res$plot,
      width = 5.6,
      height = 5.2,
      dpi = 300
    )
  }
}

dev.off()

single_gene_prop_stats <- if (length(prop_stat_list) > 0) {
  rbindlist(prop_stat_list, fill = TRUE)
} else {
  data.table()
}

if (nrow(single_gene_prop_stats) > 0) {
  setcolorder(
    single_gene_prop_stats,
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
}

fwrite(
  single_gene_prop_stats,
  file.path(outdir, "07_single_gene_prop_gt36_pairwise_wilcoxon.tsv"),
  sep = "\t"
)


# ======================== 9. freqTC matrix ========================

message("[6/8] Reading freqTC matrix ...")

freqtc_raw <- fread(freqtc_file)
setnames(freqtc_raw, clean_names(names(freqtc_raw)))
setnames(freqtc_raw, names(freqtc_raw)[1], "Gene")

freqtc_sample_cols <- intersect(sample_meta$sample, names(freqtc_raw))

if (length(freqtc_sample_cols) == 0) {
  stop("No selected plasma samples matched freqTC matrix columns.")
}

freqtc_sub <- freqtc_raw[
  Gene %chin% gene_list,
  c("Gene", freqtc_sample_cols),
  with = FALSE
]

matched_freqtc_genes <- sort(unique(freqtc_sub$Gene))
missing_freqtc_genes <- setdiff(gene_list, matched_freqtc_genes)

fwrite(
  data.table(Gene = matched_freqtc_genes),
  file.path(outdir, "08_genes_matched_in_freqTC_matrix.tsv"),
  sep = "\t"
)

fwrite(
  data.table(Gene = missing_freqtc_genes),
  file.path(outdir, "08_genes_missing_in_freqTC_matrix.tsv"),
  sep = "\t"
)

message("Genes matched in freqTC matrix: ", length(matched_freqtc_genes))
message("Genes missing in freqTC matrix: ", length(missing_freqtc_genes))

freqtc_long <- melt(
  freqtc_sub,
  id.vars = "Gene",
  measure.vars = freqtc_sample_cols,
  variable.name = "sample",
  value.name = "freqTC"
)

freqtc_long[, freqTC := as.numeric(freqTC)]

freqtc_long <- merge(
  freqtc_long,
  sample_meta[, .(sample, Group)],
  by = "sample",
  all.x = TRUE
)

gene_total_count <- single_gene_prop36[
  ,
  .(
    TotalCount = unique(TotalCount)[1]
  ),
  by = .(Gene, sample)
]

freqtc_long <- merge(
  freqtc_long,
  gene_total_count,
  by = c("Gene", "sample"),
  all.x = TRUE
)

freqtc_long[
  is.na(TotalCount) |
    TotalCount < min_gene_sample_total_count,
  freqTC := NA_real_
]

fwrite(
  freqtc_long,
  file.path(outdir, "08_single_gene_freqTC_by_sample.tsv"),
  sep = "\t"
)


# ======================== 10. Gene-set freqTC ========================

message("[7/8] Gene-set freqTC ...")

gene_set_freqtc <- freqtc_long[
  !is.na(freqTC),
  .(
    freqTC = median(freqTC, na.rm = TRUE),
    mean_freqTC = mean(freqTC, na.rm = TRUE),
    n_gene = uniqueN(Gene)
  ),
  by = .(sample, Group)
]

fwrite(
  gene_set_freqtc,
  file.path(outdir, "09_gene_set_freqTC_by_sample.tsv"),
  sep = "\t"
)

make_boxplot(
  gene_set_freqtc,
  value_col = "freqTC",
  ylab = "Gene-set median freqTC",
  title_text = "Gene-set freqTC in plasma",
  outfile_prefix = "09_gene_set_freqTC_boxplot"
)


# ======================== 11. Single-gene freqTC ========================

message("[8/8] Single-gene freqTC ...")

freqtc_stat_list <- list()

pdf(
  file.path(outdir, "10_single_gene_freqTC_boxplot_multipage.pdf"),
  width = 5.6,
  height = 5.2
)

for (g in sort(unique(freqtc_long$Gene))) {
  
  tmp <- freqtc_long[
    Gene == g &
      !is.na(freqTC)
  ]
  
  if (nrow(tmp) == 0) next
  
  res <- make_boxplot(
    tmp,
    value_col = "freqTC",
    ylab = "freqTC",
    title_text = paste0(g, " freqTC in plasma"),
    outfile_prefix = NULL
  )
  
  if (is.null(res)) next
  
  print(res$plot)
  
  stat_g <- copy(res$stat)
  stat_g[, Gene := g]
  freqtc_stat_list[[g]] <- stat_g
  
  if (save_each_single_gene_plot) {
    
    g_safe <- safe_filename(g)
    
    ggsave(
      file.path(single_gene_plot_dir, paste0(g_safe, "_freqTC_boxplot.pdf")),
      res$plot,
      width = 5.6,
      height = 5.2
    )
    
    ggsave(
      file.path(single_gene_plot_dir, paste0(g_safe, "_freqTC_boxplot.png")),
      res$plot,
      width = 5.6,
      height = 5.2,
      dpi = 300
    )
  }
}

dev.off()

single_gene_freqtc_stats <- if (length(freqtc_stat_list) > 0) {
  rbindlist(freqtc_stat_list, fill = TRUE)
} else {
  data.table()
}

if (nrow(single_gene_freqtc_stats) > 0) {
  setcolorder(
    single_gene_freqtc_stats,
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
}

fwrite(
  single_gene_freqtc_stats,
  file.path(outdir, "10_single_gene_freqTC_pairwise_wilcoxon.tsv"),
  sep = "\t"
)


# ======================== 12. Summary ========================

summary_dt <- data.table(
  item = c(
    "Gene list file",
    "Number of genes in list",
    "Selected plasma samples",
    "Healthy samples",
    "HBV samples",
    "ICP samples",
    "Genes matched in length matrix",
    "Genes missing in length matrix",
    "Genes matched in freqTC matrix",
    "Genes missing in freqTC matrix",
    "Global length x-axis range",
    "Long-fragment cutoff",
    "Minimum gene-sample total count",
    "P value used for plot annotation"
  ),
  value = c(
    gene_list_file,
    length(gene_list),
    nrow(sample_meta),
    sample_meta[Group == "Healthy", .N],
    sample_meta[Group == "HBV", .N],
    sample_meta[Group == "ICP", .N],
    length(matched_length_genes),
    length(missing_length_genes),
    length(matched_freqtc_genes),
    length(missing_freqtc_genes),
    paste0(global_length_x_min, "-", global_length_x_max, " bp"),
    paste0(">", long_cutoff, " bp"),
    min_gene_sample_total_count,
    p_label_use
  )
)

fwrite(
  summary_dt,
  file.path(outdir, "00_plasma_validation_summary.tsv"),
  sep = "\t"
)

cat("\nFinished.\n")
cat("Output directory:", outdir, "\n")