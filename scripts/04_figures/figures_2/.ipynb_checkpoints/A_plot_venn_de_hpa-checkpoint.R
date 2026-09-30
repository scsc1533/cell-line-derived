###############################################################################
# DE 基因与 HPA 组织特异性基因 Venn 图
# 目标：对每个细胞系，绘制 cfRNA DEGs、cell_debris DEGs、HPA 组织特异性基因
#       三组 Venn 图，并输出 overlap 详细表格
###############################################################################

# ======================== 0. 加载依赖 ========================
library(data.table)

# 尝试加载 Venn 图包，优先 ggVennDiagram，回退到 VennDiagram
if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
  library(ggVennDiagram)
  use_ggvenn <- TRUE
} else if (requireNamespace("VennDiagram", quietly = TRUE)) {
  library(VennDiagram)
  use_ggvenn <- FALSE
} else {
  stop("Please install either 'ggVennDiagram' or 'VennDiagram' package.")
}

# ======================== 1. 配置 ========================
base_dir <- "/data/work/01_2603cell_culture/06_xinyue/02_DE_cell_specific/mlRNA"
hpa_file <- "/data/work/01_2603cell_culture/01_raw_date/02_database/HPA/tissue_enriched_genes.tsv"
out_dir  <- "01_venn_output"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 细胞系配置：R 中名称 → 文件夹名称 → HPA 组织
cell_config <- list(
  "Hep3B2.1-7"       = list(folder = "Hep3B2.1-7",       tissue = "liver"),
  "HepG2"             = list(folder = "HepG2",             tissue = "liver"),
  "Hep3B2.1-7_HepG2" = list(folder = "Hep3B2.1-7_HepG2", tissue = "liver"),
  "K562"              = list(folder = "K562",              tissue = "bone marrow"),
  "HTR-8/SVneo"       = list(folder = "HTR-8_SVneo",       tissue = "placenta"),
  "HEK293T"           = list(folder = "HEK293T",           tissue = "adrenal gland")
)

# ======================== 2. 辅助函数 ========================

#' 读取 DE 文件，返回 Combination == "DEL" 的基因
read_deg <- function(file_path) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(character(0))
  }
  dt <- fread(file_path)
      dt[Combination %in% c("DEL", "EL", "DL", "DE"), unique(Gene)]
}

#' 构建 overlap 表格
#' @param set_left  左侧集合基因（cell+debris）
#' @param set_right 右侧集合基因（cfRNA）
#' @param set_hpa   HPA 组织特异性基因
#' @param set_labels 三组标签，顺序为 left, right, hpa
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

  # 添加所属交集类别标签
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

#' 绘制并保存 Venn 图
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

# ======================== 3. 读取 HPA 数据 ========================
message("[1/3] Reading HPA tissue-enriched genes ...")
hpa <- fread(hpa_file)
message(sprintf("  HPA entries: %d", nrow(hpa)))

# ======================== 4. 逐细胞系处理 ========================
message("[2/3] Processing each cell line ...")

for (cell_name in names(cell_config)) {

  cfg     <- cell_config[[cell_name]]
  folder  <- cfg$folder
  tissue  <- cfg$tissue

  message(sprintf("\n--- %s (tissue: %s) ---", cell_name, tissue))

  # 读取 DE 基因
  cf_file <- file.path(base_dir, "cfRNA", "single_cell_line", folder,
                       "gene_overlap_detailed_logFC2_adj0.01.csv")
  cd_file <- file.path(base_dir, "cell_debris", "single_cell_line", folder,
                       "gene_overlap_detailed_logFC2_adj0.01.csv")

  cf_genes <- read_deg(cf_file)
  cd_genes <- read_deg(cd_file)

  # 读取 HPA 组织特异性基因
  hpa_genes <- hpa[Tissues == tissue, unique(Gene)]

  message(sprintf("  cfRNA DEGs: %d, cell_debris DEGs: %d, HPA (%s): %d",
                  length(cf_genes), length(cd_genes), tissue, length(hpa_genes)))

  # 跳过三组全空的细胞系
  if (length(cf_genes) == 0 && length(cd_genes) == 0 && length(hpa_genes) == 0) {
    message("  All three sets empty, skipping.")
    next
  }

  # 三组标签（左 cell+debris，右 cfRNA，下 HPA）
  set_labels <- c("cell+debris", "cfRNA", paste0("HPA\n", tissue))

  # 输出 overlap 表格
  overlap_dt <- build_overlap_table(cd_genes, cf_genes, hpa_genes, set_labels)

  table_out <- file.path(out_dir, sprintf("overlap_%s.csv", gsub("/", "_", cell_name)))
  fwrite(overlap_dt, table_out)
  message(sprintf("  Overlap table saved: %s (%d genes)", table_out, nrow(overlap_dt)))

  # 打印交集统计
  intersect_counts <- overlap_dt[, .N, by = Intersection]
  for (i in seq_len(nrow(intersect_counts))) {
    message(sprintf("    %-30s : %d", intersect_counts$Intersection[i], intersect_counts$N[i]))
  }

  # 绘制 Venn 图（有边缘）
  venn_out <- file.path(out_dir,
                        sprintf("venn_%s.pdf", gsub("/", "_", cell_name)))
  draw_venn(cd_genes, cf_genes, hpa_genes,
            labels    = set_labels,
            title     = cell_name,
            out_pdf   = venn_out,
            edge_size = 1.2)
  message(sprintf("  Venn plot saved: %s", venn_out))

  # 绘制 Venn 图（无边缘）
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
