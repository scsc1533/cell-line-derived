###############################################################################
# 持家基因片段长度分布绘图脚本
# 目标：绘制 RACK1, EEF1A1, HNRNPA1, RPL23A, RPL31, GAPDH 在不同
#       Component（cfRNA, debris, cell, debris+cell）中的长度分布
# 输出：每组 Component 一张图，不同 Cell 用不同颜色区分，附带误差阴影
# 优化：在 melt 前先按样本 + 基因过滤，大幅降低内存占用
###############################################################################

# ======================== 0. 加载依赖 ========================
library(data.table)
library(ggplot2)
library(scales)

# ======================== 1. 读取数据 ========================
message("[1/6] Reading input files ...")

qc_raw     <- fread("/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt")
length_raw <- fread("/data/work/01_2603cell_culture/01_raw_date/length_motif/all_samples_length_merged.mlncRNA.txt")

# 清理列名：去除可能携带的 UTF-8 BOM 和首尾空白
setnames(qc_raw, trimws(sub("\ufeff", "", names(qc_raw), fixed = TRUE)))

# 长度表列名：前4列为元数据，其余为长度位点 1-100
# 如果 fread 未识别表头（列名为 V1, V2, ...），则手动指定
n_len_cols <- ncol(length_raw) - 4L
length_col_names <- c("sample", "Transcript", "Gene", "Type", as.character(seq_len(n_len_cols)))
setnames(length_raw, seq_along(length_raw), length_col_names)

message(sprintf("  Raw length table rows: %d, cols: %d", nrow(length_raw), ncol(length_raw)))
message(sprintf("  Length table columns: %s", paste(head(names(length_raw), 6), collapse = ", ")))

# ======================== 2. 样本 + 基因双维度预过滤 ========================
message("[2/6] Pre-filtering by QC samples and housekeeping genes ...")

# 只保留 QC == 1 的样本
qc_pass <- qc_raw[QC == 1, .(Sample, Cell, Component)]
message(sprintf("  QC-passed samples: %d", nrow(qc_pass)))

# 持家基因列表
hk_genes <- c("SPN", "DRAXIN", "NUPR1", "WDR74", "MIR99AHG", "TPT1-AS1")
hk_dt <- data.table(Gene = hk_genes)

# 先用 data.table join 过滤样本（避免 %in% 在大表上的类型歧义）
length_filt <- length_raw[qc_pass[, .(Sample)], on = .(sample = Sample), nomatch = 0]
# 再用 %chin% 过滤基因（%chin% 是 data.table 专门为字符列优化的 fast %in%）
length_filt <- length_filt[Gene %chin% hk_genes]

message(sprintf("  Length table rows after pre-filtering: %d (was %d)",
                nrow(length_filt), nrow(length_raw)))

# 检查基因覆盖情况
found_genes <- unique(length_filt$Gene)
missing_genes <- setdiff(hk_genes, found_genes)
if (length(missing_genes) > 0) {
  warning("Genes NOT found in data: ", paste(missing_genes, collapse = ", "))
}
message(sprintf("  Genes retained: %s", paste(found_genes, collapse = ", ")))

# ======================== 3. 合并样本元信息 ========================
message("[3/6] Merging sample metadata ...")

# 用 data.table join 追加 Cell 和 Component 信息
length_filt <- qc_pass[length_filt, on = .(Sample = sample), nomatch = 0]
setnames(length_filt, "Sample", "sample")
length_filt[, c("Transcript", "Type") := NULL]  # 不再需要转录本信息

# Cell 因子顺序（控制图例）
length_filt[, Cell := factor(Cell, levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T"))]

# ======================== 4. 数据整形与基因聚合 ========================
message("[4/6] Melting to long format and aggregating ...")

length_cols <- grep("^[0-9]+$", names(length_filt), value = TRUE)
message(sprintf("  Length positions: %d (range %s - %s)",
                length(length_cols), length_cols[1], length_cols[length(length_cols)]))

# 宽表 -> 长表（此时数据量已大幅缩减）
length_long <- melt(length_filt,
                    id.vars       = c("sample", "Gene", "Cell", "Component"),
                    measure.vars  = length_cols,
                    variable.name = "Length",
                    value.name    = "Count")
length_long[, Length := as.integer(as.character(Length))]

# 按 sample × Gene 聚合所有转录本，求和每个长度位点的 count
# （由于宽表阶段已过滤基因且每行就是一个转录本，这里直接按 sample × Gene 求和即可）
gene_agg <- length_long[, .(Count = sum(Count, na.rm = TRUE)),
                        by = .(sample, Gene, Cell, Component, Length)]

# 计算每个 sample × Gene 内部的比例
gene_agg[, Total := sum(Count), by = .(sample, Gene)]
gene_agg[, Proportion := ifelse(Total > 0, Count / Total, 0)]

message(sprintf("  Aggregated rows: %d", nrow(gene_agg)))
rm(length_raw, length_filt, length_long)  # 释放中间变量
gc()

# ======================== 4.5. 创建 cell+debris 合并分组 ========================
message("[4.5/6] Creating cell+debris pooled group ...")

# 将 cell 和 debris 样本的 Proportion 合并，Component 统一标记为 cell+debris
cd_pooled <- gene_agg[Component %in% c("cell", "debris")]
cd_pooled[, Component := "cell+debris"]
# 不重新聚合 Count，直接复用已算好的 Proportion（下游按 Cell × Gene × Component 计算 Mean/SE）
gene_agg <- rbind(gene_agg, cd_pooled)

message(sprintf("  Aggregated rows (with cell+debris): %d", nrow(gene_agg)))

# ======================== 5. 绘图 ========================
message("[5/6] Generating plots ...")

# 主题设置（学术期刊风格，无背景网格线）
theme_academic <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = "grey80"),
    strip.text         = element_text(size = 11, face = "bold"),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y        = element_text(size = 9),
    axis.title         = element_text(size = 12),
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.margin        = margin(10, 15, 10, 10)
  )

# 预计算全量汇总统计量（Cell × Gene × Component × Length）
message("  Pre-computing summary statistics ...")
plot_dt_all <- gene_agg[, .(
  Mean = mean(Proportion, na.rm = TRUE),
  SE   = sd(Proportion, na.rm = TRUE) / sqrt(.N),
  N    = .N
), by = .(Cell, Gene, Component, Length)]
plot_dt_all[, `:=`(
  ymin = Mean - SE,
  ymax = Mean + SE
)]
plot_dt_all[N <= 1L, `:=`(ymin = NA_real_, ymax = NA_real_)]

# 按基因 × 组分依次出图（一个基因一个组分一张图）
for (gene in hk_genes) {
  for (comp in sort(unique(plot_dt_all$Component))) {

    sub_data <- plot_dt_all[Gene == gene & Component == comp]

    if (nrow(sub_data) == 0) {
      message(sprintf("  Skipping '%s' × '%s': no data", gene, comp))
      next
    }

    # 动态确定 x 轴范围
    x_min <- sub_data[Mean > 0, min(Length, na.rm = TRUE)]
    x_max <- sub_data[Mean > 0, max(Length, na.rm = TRUE)]
    x_range <- seq(x_min, x_max, by = 1)

    has_se <- sub_data[, any(!is.na(ymin))]

    message(sprintf("  %s | %s: x-axis %d - %d bp, SE: %s",
                    gene, comp, x_min, x_max, has_se))

    p <- ggplot(sub_data,
                aes(x     = Length,
                    y     = Mean,
                    color = Cell,
                    fill  = Cell,
                    group = Cell))

    if (has_se) {
      p <- p + geom_ribbon(data = sub_data[N > 1L],
                           aes(ymin = ymin, ymax = ymax),
                           alpha    = 0.15,
                           colour   = NA)
    }

    p <- p +
      geom_line(linewidth = 0.8, na.rm = TRUE) +

      scale_x_continuous(breaks = x_range,
                         limits = c(min(x_range) - 0.5, max(x_range) + 0.5),
                         expand = c(0, 0)) +

      scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                         labels = label_number(accuracy = 0.001)) +

      labs(x     = "Fragment Length (bp)",
           y     = "Proportion",
           title = sprintf("%s — %s", gene, comp)) +

      theme_academic

    out_name <- sprintf("length_distribution_%s_%s.pdf",
                        gene, gsub("[+ ]", "_", comp))
    ggsave(out_name, plot = p, width = 8, height = 6, device = "pdf")
    message(sprintf("  Saved: %s", out_name))
  }
}

message("[6/6] Done.")
