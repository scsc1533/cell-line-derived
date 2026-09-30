###############################################################################
# 片段长度分布 PCA
# 用 RACK1 和 GAPDH 在各样本中的长度比例向量作为特征进行主成分分析
###############################################################################

library(ggplot2)
library(dplyr)
library(data.table)

# ======================== 文件路径 ========================
length_file     <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/all_samples_length_merged.mlncRNA.txt"
sample_info_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
output_dir      <- "/data/work/01_2603cell_culture/07_figure/02_figure2/04_plot_length_PCA"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ======================== 读取样本信息 ========================
sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_info_filtered <- sample_info %>%
  filter(QC == 1) %>%
  mutate(
    Cell      = factor(Cell, levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")),
    Component = factor(Component, levels = c("cfRNA", "debris", "cell"))
  )

# ======================== 读取长度数据 ========================
message("[1/4] Reading length distribution data ...")
length_raw <- fread(length_file)

# 列名：前 4 列为元数据，其余为长度位点 1-100
n_len_cols <- ncol(length_raw) - 4L
setnames(length_raw, seq_along(length_raw),
         c("sample", "Transcript", "Gene", "Type", as.character(seq_len(n_len_cols))))

message(sprintf("  Raw rows: %d, cols: %d", nrow(length_raw), ncol(length_raw)))

# ======================== 过滤与聚合 ========================
message("[2/4] Pre-filtering by samples and target genes ...")

# 在宽表阶段同时过滤 QC 样本和 RACK1/GAPDH（melt 前，大幅减少内存）
target_genes <- c("RACK1", "GAPDH", "RPL23A")
length_filt <- length_raw[sample %in% sample_info_filtered$Sample & Gene %in% target_genes]
message(sprintf("  Filtered rows: %d (was %d)", nrow(length_filt), nrow(length_raw)))

# 合并元信息（使用 data.table join，避免 dplyr 导致类型丢失）
length_filt <- merge(length_filt,
                     sample_info_filtered[, c("Sample", "Cell", "Component")],
                     by.x = "sample", by.y = "Sample",
                     all.x = FALSE)

length_cols <- grep("^[0-9]+$", names(length_filt), value = TRUE)

# 转长表（此时数据量已大幅缩减）
length_long <- melt(length_filt,
                    id.vars       = c("sample", "Cell", "Component", "Gene"),
                    measure.vars  = length_cols,
                    variable.name = "Length",
                    value.name    = "Count")
length_long[, Length := as.integer(as.character(Length))]

# 按 sample × Gene 聚合转录本，计算比例
length_long <- length_long[, .(Count = sum(Count, na.rm = TRUE)),
                           by = .(sample, Cell, Component, Gene, Length)]
length_long[, Proportion := Count / sum(Count), by = .(sample, Gene)]

message(sprintf("  Aggregated rows: %d", nrow(length_long)))
rm(length_raw, length_filt)
gc()

# ======================== PCA 绘图 ========================
message("[3/4] Running PCA ...")

# PCA + 绘图子函数
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

  # 仅 cfRNA 图添加置信椭圆
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

  # 全样本
  run_pca(gene_data, g, "all components", "")

  # 仅 cfRNA
  gene_cf <- gene_data[Component == "cfRNA"]
  run_pca(gene_cf, g, "cfRNA only", "_cfRNA", add_ellipse = TRUE)
}

message("[4/4] Done. Outputs: PCA_length_RACK1.pdf, PCA_length_GAPDH.pdf")
