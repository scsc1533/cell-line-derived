library(tidyverse)
library(corrplot)
library(RColorBrewer)

# ================== 读取数据 ==================
exp <- read.table("/data/work/1/01_quality_control/all_mRNA_TPM.txt",
                  row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)

info <- read_tsv("/data/work/1/cellcult_Sample_Information.txt")

# ================== QC筛选 ==================
info_filtered <- info %>%
  filter(QC %in% c(1, 2) & !Component %in% c("blank", "water"))

# ================== 按 Cell × Component 合并取均值 ==================
exp_t <- as.data.frame(t(exp))
exp_t$Sample <- rownames(exp_t)
rownames(exp_t) <- NULL

exp_long <- exp_t %>%
  pivot_longer(-Sample, names_to = "Gene", values_to = "TPM") %>%
  left_join(info_filtered %>% select(Sample, Cell, Component), by = "Sample") %>%
  group_by(Cell, Component, Gene) %>%
  summarise(TPM_mean = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  unite("Group", Cell, Component, sep = "_") %>%
  pivot_wider(names_from = "Gene", values_from = "TPM_mean")

exp_mat <- as.data.frame(t(exp_long[, -1]))
colnames(exp_mat) <- exp_long$Group
exp_mat <- as.matrix(sapply(exp_mat, as.numeric))
rownames(exp_mat) <- exp_long$Gene

exp_log <- log2(exp_mat + 1)

# ================== 排序 ==================
cell_order <- c("Hep3B2.1-7", "HepG2", "HTR-8/SVneo", "K562", "HEK293T")
comp_order <- c("cfRNA", "debris", "cell")

ordered_groups <- c()
for (cl in cell_order) {
  for (co in comp_order) {
    ordered_groups <- c(ordered_groups, paste0(cl, "_", co))
  }
}
ordered_groups <- ordered_groups[ordered_groups %in% colnames(exp_log)]

exp_log <- exp_log[, ordered_groups, drop = FALSE]

# ================== Spearman 相关 ==================
cor_mat <- cor(exp_log, method = "spearman")

# ================== 画右上三角热图 ==================
pdf("All_Cell_Component_Correlation_mRNA_grouped.pdf", width = 14, height = 14)

# 红→白→蓝配色（-1红，0白，1蓝）- 200个渐变
col <- colorRampPalette(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
                          "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200)

# 计算实际范围用于调试图
cat("相关性范围:", round(min(cor_mat, na.rm = TRUE), 2), "到", 
    round(max(cor_mat, na.rm = TRUE), 2), "\n")

corrplot(cor_mat, 
         method = "color",
         type = "upper",
         col = col,
         addCoef.col = "white",      # 数值改成白色
         number.cex = 0.7,
         number.font = 2,            # 加粗数值（可选）
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,
         cl.pos = "r",
         cl.lim = c(-1, 1),
         is.corr = TRUE)

dev.off()

cat("\nDone! 输出: All_Cell_Component_Correlation_mRNA_grouped.pdf\n")
