# =============================================================================
# 按组分绘制长度分布图（不同细胞系用不同颜色）
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------ 文件路径 -------------------------------------
info_file <- "/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt"
freq_file <- "/data/work/01_2603cell_culture/01_raw_date/length_motif/sample_mRNA_length_fre.txt"

out_dir <- "./08_sample_length_distribution"
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

# ------------------------------ 2. 均值计算 ----------------------------------
# Cell 因子顺序（控制图例）
cell_order <- c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T")

line_data <- data_merged %>%
  filter(readLength >= 17 & readLength <= 90) %>%
  group_by(Cell, Component, readLength) %>%
  summarise(MeanFreq = mean(Frequency, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Cell      = factor(Cell,      levels = cell_order),
    Component = factor(Component, levels = c("cfRNA", "debris", "cell"))
  )

# ------------------------------ 3. 按组分分别绘图 -----------------------------
components <- levels(line_data$Component)

for (comp_i in components) {
  p_comp <- line_data %>%
    filter(Component == comp_i) %>%
    ggplot(aes(x = readLength, y = MeanFreq, color = Cell)) +
    geom_line(linewidth = 0.8) +
    labs(title = paste("Length distribution -", comp_i),
         x = "Read length (bp)", y = "Mean frequency") +
    theme_bw() +
    theme(legend.position = "bottom")

  safe_name <- gsub("[/\\:*?\"<>| ]", "_", comp_i)
  ggsave(filename = file.path(out_dir, paste0("length_distribution_", safe_name, ".pdf")),
         plot = p_comp, width = 6, height = 5, device = "pdf")
  message("Saved: length_distribution_", safe_name, ".pdf")
}

# ------------------------------ 4. 合并图（分面）------------------------------
p_all <- ggplot(line_data, aes(x = readLength, y = MeanFreq, color = Cell)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Component, scales = "free_y", ncol = 3) +
  labs(x = "Read length (bp)", y = "Mean frequency") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "length_distribution_all_components.pdf"),
       p_all, width = 14, height = 5, device = "pdf")
message("Saved: length_distribution_all_components.pdf")
