###############################################################################
# 一 vs 其余 特异性基因 与 cfRNA 高表达基因的 Venn 分析
# 目标：对 scores 文件中 perm_fdr_wasserstein < 0.05 的基因，
#       分别与对应细胞系 cfRNA 高表达基因做 overlap，输出 Venn 图和表格
###############################################################################

library(data.table)
library(VennDiagram)

# ======================== 文件路径 ========================
scores_file <- "/data/work/01_2603cell_culture/07_figure/03_figure3/01_length_distribution_one_vs_rest/all_comparisons_gene_scores.tsv"
overlap_dir <- "/data/work/01_2603cell_culture/07_figure/02_figure2/01_venn_output"

out_dir <- "02_venn_one_vs_rest"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======================== 配置：comparison → overlap CSV → 标签 ========================
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

# ======================== 1. 读取 scores 文件 ========================
message("[1/3] Reading one-vs-rest scores ...")
scores <- fread(scores_file)
message(sprintf("  Total rows: %d", nrow(scores)))

# 过滤：target comparisons + fdr < 0.05
scores_sig <- scores[comparison %chin% comp_map$comparison &
                     perm_fdr_wasserstein < 0.05]
message(sprintf("  Significant genes (fdr < 0.05): %d", nrow(scores_sig)))

# ======================== 2. 逐组做 Venn ========================
message("[2/3] Processing Venn overlaps ...")

# Venn 绘图函数（固定圆大小，无黑框，左 One vs Rest / 右 cfRNA high）
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

# 遍历每组
for (i in seq_len(nrow(comp_map))) {

  comp    <- comp_map$comparison[i]
  csv_fn  <- comp_map$overlap_csv[i]
  label   <- comp_map$label[i]

  message(sprintf("\n--- %s ---", label))

  # 获取 one-vs-rest 显著基因
  ovr_genes <- scores_sig[comparison == comp, unique(Gene)]
  message(sprintf("  One-vs-rest significant genes: %d", length(ovr_genes)))

  # 读取 overlap CSV，提取 cfRNA 高表达基因
  csv_path <- file.path(overlap_dir, csv_fn)
  if (!file.exists(csv_path)) {
    message(sprintf("  Overlap CSV not found: %s, skipping.", csv_path))
    next
  }

  overlap_dt <- fread(csv_path)
  cfRNA_genes <- overlap_dt[cfRNA == TRUE, unique(Gene)]
  message(sprintf("  cfRNA-high genes: %d", length(cfRNA_genes)))

  # 跳过任一集合为空的组
  if (length(ovr_genes) == 0 && length(cfRNA_genes) == 0) {
    message("  Both sets empty, skipping.")
    next
  }

  # 构建 overlap 表格
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

  # 输出 overlap 表格
  safe_label <- gsub("[ /]", "_", label)
  table_out <- file.path(out_dir, sprintf("overlap_one_vs_rest_%s.csv", safe_label))
  fwrite(dt, table_out)
  message(sprintf("  Overlap table saved: %s (%d genes)", table_out, nrow(dt)))

  # 打印交集统计
  intersect_counts <- dt[, .N, by = Intersection]
  for (j in seq_len(nrow(intersect_counts))) {
    message(sprintf("    %-18s : %d", intersect_counts$Intersection[j], intersect_counts$N[j]))
  }

  # 绘制 Venn 图
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
