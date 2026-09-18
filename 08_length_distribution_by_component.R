# =============================================================================
# 按组分绘制长度分布图（cfRNA / debris / cell 各一条聚合线）
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------ 文件路径 -------------------------------------
info_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
freq_file <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/sample_mRNA_length_fre.txt"

out_dir <- "./08_length_distribution_by_component"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------ 1. 读取数据 ----------------------------------
info <- read.delim(info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
info_filtered <- info %>% filter(QC == 1)

freq <- read.delim(freq_file, header = TRUE, sep = "\t", check.names = FALSE)
freq_long <- freq %>%
  pivot_longer(cols = -readLength, names_to = "Sample", values_to = "Frequency") %>%
  mutate(readLength = as.numeric(readLength))

data_merged <- freq_long %>%
  inner_join(info_filtered, by = "Sample") %>%
  select(Sample, Cell, Component, readLength, Frequency) %>%
  filter(!is.na(Cell), Cell != "")

# ------------------------------ 2. 按组分聚合均值 -----------------------------
comp_order <- c("cell", "debris", "cfRNA")
comp_colors <- c("cfRNA" = "#E64B35", "debris" = "#00A087", "cell" = "#4DBBD5")

line_data <- data_merged %>%
  filter(readLength >= 17 & readLength <= 90) %>%
  group_by(Component, readLength) %>%
  summarise(MeanFreq = mean(Frequency, na.rm = TRUE), .groups = "drop") %>%
  mutate(Component = factor(Component, levels = comp_order))

# ------------------------------ 3. 绘图 --------------------------------------
p <- ggplot(line_data, aes(x = readLength, y = MeanFreq, color = Component)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = comp_colors, breaks = comp_order) +
  labs(x = "Read length (bp)", y = "Mean frequency") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "length_distribution_by_component.pdf"),
       p, width = 8, height = 5, device = "pdf")
message("Saved: length_distribution_by_component.pdf")
