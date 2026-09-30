# ============================================================
# Script: Cell x Component correlation heatmap for mRNA
# Purpose:
#   Read mRNA TPM matrix and sample information, keep QC 1/2
#   samples excluding blank/water, average TPM by Cell x
#   Component group, compute Spearman correlations between
#   groups on log2(TPM+1) values, and draw an upper-triangular
#   correlation heatmap with coefficients, ordered by cell line
#   and component.
# ============================================================

library(tidyverse)
library(corrplot)
library(RColorBrewer)

# ================== Read data ==================
exp <- read.table("./01_quality_control/all_mRNA_TPM.txt",
                  row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)

info <- read_tsv("./cellcult_Sample_Information.txt")

# ================== QC filtering ==================
info_filtered <- info %>%
  filter(QC %in% c(1, 2) & !Component %in% c("blank", "water"))

# ================== Average by Cell x Component ==================
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

# ================== Ordering ==================
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

# ================== Spearman correlation ==================
cor_mat <- cor(exp_log, method = "spearman")

# ================== Draw upper-triangular heatmap ==================
pdf("All_Cell_Component_Correlation_mRNA_grouped.pdf", width = 14, height = 14)

# Red -> White -> Blue palette (-1 red, 0 white, 1 blue) - 200 gradients
col <- colorRampPalette(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
                          "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200)

# Compute actual range for debugging
cat("Correlation range:", round(min(cor_mat, na.rm = TRUE), 2), "to", 
    round(max(cor_mat, na.rm = TRUE), 2), "\n")

corrplot(cor_mat, 
         method = "color",
         type = "upper",
         col = col,
         addCoef.col = "white",      # Change coefficient text to white
         number.cex = 0.7,
         number.font = 2,            # Bold coefficients (optional)
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,
         cl.pos = "r",
         cl.lim = c(-1, 1),
         is.corr = TRUE)

dev.off()

cat("\nDone! Output: All_Cell_Component_Correlation_mRNA_grouped.pdf\n")