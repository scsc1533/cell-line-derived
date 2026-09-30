###############################################################################
# All_three 交集基因长度分布图
# 目标：针对 overlap_Hep3B2.1-7_HepG2.csv 中 Intersection == "All_three" 的基因，
#       绘制其肝相关样本中 cfRNA / cell / debris / cell+debris 的长度分布
#       每个组分一张图，Hep3B2.1-7 和 HepG2 不同颜色
#       含整体合并图和单基因图，x 轴含 33 bp 刻度，自动裁剪前后零值区间
# 风格：照搬 length_distribution_by_component.R
###############################################################################

library(data.table)
library(ggplot2)

# ======================== 文件路径 ========================
overlap_csv  <- "01_venn_output/overlap_Hep3B2.1-7_HepG2.csv"
info_file    <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
length_file  <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/all_samples_length_merged.mlncRNA.txt"

out_dir <- "05_length_liver_HPA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 肝相关细胞系及其配色
liver_cells <- c("Hep3B2.1-7", "HepG2")
cell_colors <- c("Hep3B2.1-7" = "#E64B35FF", "HepG2" = "#4DBBD5FF")

# 四个组分
all_comps <- c("cfRNA", "cell", "debris", "cell+debris")

# ======================== 1. 读取 overlap，筛选 All_three 基因 ========================
message("[1/6] Reading overlap CSV and filtering All_three genes ...")
overlap_dt <- fread(overlap_csv)
all_three_genes <- overlap_dt[Intersection == "All_three", unique(Gene)]
message(sprintf("  All_three genes (before count filter): %d", length(all_three_genes)))
if (length(all_three_genes) == 0) stop("No All_three genes found.")

# ======================== 2. 读取样本信息 + 长度数据 ========================
message("[2/6] Reading sample info and length data ...")

qc_raw <- fread(info_file)
setnames(qc_raw, trimws(sub("\ufeff", "", names(qc_raw), fixed = TRUE)))
qc_pass <- qc_raw[QC == 1 & Cell %chin% liver_cells,
                  .(Sample, Cell, Component)]
message(sprintf("  Liver QC-passed samples: %d", nrow(qc_pass)))

# ======================== 2.5. 表达量 count 过滤 ========================
message("[2.5/6] Filtering genes by expression count > 30 in liver cfRNA samples ...")
counts_file <- "/data/work/01_2603cell_culture/01_raw_date/expression_matrix/all_mlRNA_counts.txt"

# 所有 QC==1 的肝细胞系 cfRNA 样本
cfRNA_liver <- qc_pass[Component == "cfRNA", Sample]
message(sprintf("  Liver cfRNA QC-passed samples: %d", length(cfRNA_liver)))

counts_raw <- fread(counts_file)

# 定位基因名列（可能是 gene_id 或 sym_id）
gene_col <- intersect(c("gene_id", "sym_id", "Gene"), names(counts_raw))[1]
if (is.na(gene_col)) stop("Cannot find gene ID column in counts file.")
setnames(counts_raw, gene_col, "Gene")

# 找出 counts 文件中存在的肝 cfRNA 样本列
count_cols <- intersect(names(counts_raw), cfRNA_liver)
message(sprintf("  Matched sample columns in counts file: %d / %d",
                length(count_cols), length(cfRNA_liver)))
if (length(count_cols) == 0) {
  stop("No cfRNA sample columns matched. Counts columns: ",
       paste(head(setdiff(names(counts_raw), "Gene"), 10), collapse = ", "))
}

counts_filt <- counts_raw[Gene %chin% all_three_genes,
                          .SD, .SDcols = c("Gene", count_cols)]
message(sprintf("  Genes found in counts file: %d / %d",
                nrow(counts_filt), length(all_three_genes)))

# 保留在所有肝 cfRNA 样本中 count > 30 的基因
n_samp <- length(count_cols)
counts_filt[, keep := rowSums(.SD > 30, na.rm = TRUE) == n_samp, .SDcols = -"Gene"]
genes_pass_count <- counts_filt[keep == TRUE, unique(Gene)]
all_three_genes <- intersect(all_three_genes, genes_pass_count)
message(sprintf("  Genes passing count > 30 in all %d cfRNA samples: %d",
                n_samp, length(genes_pass_count)))
message(sprintf("  All_three genes (after count filter): %d", length(all_three_genes)))
if (length(all_three_genes) == 0) stop("No genes pass the count filter.")

# ======================== 新增：输出最终通过质控的基因列表 ========================
gene_list_file <- file.path(out_dir, "All_three_genes_QCpassed.txt")
writeLines(all_three_genes, gene_list_file)
message(sprintf("  Gene list written to: %s", gene_list_file))
# <--- NEW 以上三行

length_raw <- fread(length_file)
n_len_cols <- ncol(length_raw) - 4L
setnames(length_raw, seq_along(length_raw),
         c("sample", "Transcript", "Gene", "Type", as.character(seq_len(n_len_cols))))

message(sprintf("  Raw length rows: %d, cols: %d", nrow(length_raw), ncol(length_raw)))

# 预过滤：肝样本 + All_three 基因
length_filt <- length_raw[sample %in% qc_pass$Sample & Gene %chin% all_three_genes]
message(sprintf("  Filtered rows: %d (was %d)", nrow(length_filt), nrow(length_raw)))

# 合并元信息
length_filt <- merge(length_filt,
                     qc_pass[, .(Sample, Cell, Component)],
                     by.x = "sample", by.y = "Sample", all.x = FALSE)
length_filt[, c("Transcript", "Type") := NULL]

# ======================== 3. melt + 聚合 + 计算 per-gene 比例 ========================
message("[3/6] Melting, aggregating, and computing proportions ...")

length_cols <- grep("^[0-9]+$", names(length_filt), value = TRUE)

length_long <- melt(length_filt,
                    id.vars       = c("sample", "Gene", "Cell", "Component"),
                    measure.vars  = length_cols,
                    variable.name = "Length",
                    value.name    = "Count")
length_long[, Length := as.integer(as.character(Length))]

# 按 sample × Gene 聚合所有转录本
gene_agg <- length_long[, .(Count = sum(Count, na.rm = TRUE)),
                        by = .(sample, Gene, Cell, Component, Length)]

# 计算 per-gene 比例（每个 sample × Gene 内部）
gene_agg[, Total := sum(Count), by = .(sample, Gene)]
gene_agg[, Proportion := ifelse(Total > 0, Count / Total, 0)]

message(sprintf("  Aggregated rows: %d", nrow(gene_agg)))
rm(length_raw, length_filt, length_long)
gc()

# ======================== 4. 创建 cell+debris 合并分组 ========================
message("[4/6] Creating cell+debris pooled group ...")

cd_pooled <- gene_agg[Component %in% c("cell", "debris")]
cd_pooled[, Component := "cell+debris"]
gene_agg <- rbind(gene_agg, cd_pooled)

gene_agg[, Cell      := factor(Cell,      levels = liver_cells)]
gene_agg[, Component := factor(Component, levels = all_comps)]

message(sprintf("  Aggregated rows (with cell+debris): %d", nrow(gene_agg)))

# ======================== 5. 预计算均值 ========================

# (A) 按 Gene × Cell × Component × Length：用于单基因图
plot_dt_gene <- gene_agg[, .(
  MeanFreq = mean(Proportion, na.rm = TRUE)
), by = .(Gene, Cell, Component, Length)]

# (B) 按 Cell × Component × Length（所有 All_three 基因聚合）：用于整体图
plot_dt_all <- gene_agg[, .(
  MeanFreq = mean(Proportion, na.rm = TRUE)
), by = .(Cell, Component, Length)]

# ======================== 6. 绘图 ========================
message("[5/5] Generating plots ...")

# 动态计算有数据的 x 轴范围（裁剪前后零值），取所有 Cell 的并集
get_x_range <- function(dt) {
  dt[MeanFreq > 0, .(x_min = min(Length), x_max = max(Length)),
     by = .(Cell, Component)]
}

range_all  <- get_x_range(plot_dt_all)
range_gene <- plot_dt_gene[, get_x_range(.SD), by = Gene]

# 绘图子函数（照搬 length_distribution_by_component.R）
draw_cell_plot <- function(plot_dt, range_dt, title_text) {

  x_min <- range_dt[, min(x_min)]
  x_max <- range_dt[, max(x_max)]

  # 只保留有效区间内的数据
  plot_dt <- plot_dt[Length >= x_min & Length <= x_max]

  ggplot(plot_dt, aes(x = Length, y = MeanFreq, color = Cell)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = cell_colors,
                       breaks = names(cell_colors)) +
    scale_x_continuous(breaks = function(x) sort(unique(c(pretty(x), 33, 35, 37)))) +
    labs(x = "Read length (bp)", y = "Mean frequency", title = title_text) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# ---------- 6a. 整体合并图（每个组分一张）----------
message("  Generating overall combined plots ...")

for (comp_i in all_comps) {
  comp_data <- plot_dt_all[Component == comp_i]
  comp_range <- range_all[Component == comp_i]

  if (nrow(comp_data) == 0 || nrow(comp_range) == 0) next

  p <- draw_cell_plot(comp_data, comp_range,
                      sprintf("All_three genes (n=%d) — %s",
                              length(all_three_genes), comp_i))

  safe_comp <- gsub("[+ ]", "_", comp_i)
  ggsave(file.path(out_dir, paste0("length_all_three_overall_", safe_comp, ".pdf")),
         p, width = 8, height = 5, device = "pdf")
  message(sprintf("  Saved: length_all_three_overall_%s.pdf", safe_comp))
}

# ---------- 6b. 单基因图（每个基因每个组分一张）----------
message("  Generating per-gene plots ...")

for (g in all_three_genes) {

  safe_gene <- gsub("[/\\:*?\"<>| ]", "_", g)
  gene_data <- plot_dt_gene[Gene == g]
  gene_range <- range_gene[Gene == g]

  if (nrow(gene_data) == 0) next

  for (comp_i in all_comps) {
    comp_gene_data <- gene_data[Component == comp_i]
    comp_gene_range <- gene_range[Component == comp_i]
    if (nrow(comp_gene_data) == 0 || nrow(comp_gene_range) == 0) next

    p <- draw_cell_plot(comp_gene_data, comp_gene_range,
                        sprintf("%s — %s", g, comp_i))

    safe_comp <- gsub("[+ ]", "_", comp_i)
    ggsave(file.path(out_dir, paste0("length_all_three_", safe_gene, "_", safe_comp, ".pdf")),
           p, width = 8, height = 5, device = "pdf")
    message(sprintf("  Saved: length_all_three_%s_%s.pdf", safe_gene, safe_comp))
  }
}

message("Done. All outputs in: ", out_dir)