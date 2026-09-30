# 加载必要的包
library(tidyverse)
library(ggpubr)
library(gridExtra)

# 定义分组顺序
group_order <- c("cfRNA", "debris", "cell", "blank", "water")

# ------------------------------
# 1. 读取数据并整理
# ------------------------------

genenum <- read_tsv("/data/work/01_2603cell_culture/02_qc_batch/01_qc/genenum_all.txt",
                    col_select = c(Sample, 
                                   `mRNA(TPM>0)`, 
                                   `lncRNA(TPM>0)`,
                                   `miRNA(RPM>0)`,
                                   `tRNA(RPM>0)`,
                                   `piRNA(RPM>0)`)) %>%
  rename(mRNA   = `mRNA(TPM>0)`,
         lncRNA = `lncRNA(TPM>0)`,
         miRNA  = `miRNA(RPM>0)`,
         tRNA   = `tRNA(RPM>0)`,
         piRNA  = `piRNA(RPM>0)`)

info <- read_tsv("/data/work/01_2603cell_culture/01_raw_date/information_table/cellcult_Sample_Information.txt") %>%
  select(Sample, Component)

data <- left_join(genenum, info, by = "Sample") %>%
  mutate(mlncRNA = mRNA + lncRNA,
         Component = factor(Component, levels = group_order))

print(unique(data$Component))

# ------------------------------
# 2. 定义函数：分组两两显著性检验并输出表格
# ------------------------------
pairwise_test <- function(df, y_var, y_name, output_prefix) {
  groups <- levels(df$Component)
  groups <- groups[groups %in% unique(df$Component)]
  combs <- combn(groups, 2, simplify = FALSE)
  
  results <- map_dfr(combs, function(g) {
    x <- df %>% filter(Component == g[1]) %>% pull({{ y_var }})
    y <- df %>% filter(Component == g[2]) %>% pull({{ y_var }})
    test <- wilcox.test(x, y, exact = FALSE)
    tibble(
      group1 = g[1],
      group2 = g[2],
      statistic = test$statistic,
      p_value = test$p.value,
      method = test$method
    )
  }) %>%
    mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  filename <- paste0(output_prefix, "_pairwise_test.txt")
  write_tsv(results, filename)
  message("Saved pairwise test results to ", filename)
  return(results)
}

# ------------------------------
# 3. 定义绘图函数（添加blank背景上限和实验组IQR界限，带数值标注）
# ------------------------------
plot_box_with_iqr <- function(df, y_var, y_label, comparisons = NULL, filename = NULL) {
  
  # 获取分组水平
  group_levels <- levels(df$Component)
  n_groups <- length(group_levels)
  
  # 计算每组的统计量
  iqr_stats <- df %>%
    group_by(Component) %>%
    summarise(
      q1 = quantile({{ y_var }}, 0.25, na.rm = TRUE),
      q3 = quantile({{ y_var }}, 0.75, na.rm = TRUE),
      median = median({{ y_var }}, na.rm = TRUE),
      sd = sd({{ y_var }}, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      IQR = q3 - q1,
      lower_bound = q1 - 1.5 * IQR,
      upper_bound = q3 + 1.5 * IQR,
      background_ceiling = median + 3 * sd,
      # 为每个组分配x位置（1, 2, 3...）
      x_pos = match(Component, group_levels)
    )
  
  # 基础箱线图
  p <- ggplot(df, aes(x = Component, y = {{ y_var }})) +
    geom_boxplot(fill = "skyblue", outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1)
  
  # 标注blank组的背景上限（红色虚线，仅覆盖blank组宽度）
  blank_stats <- iqr_stats %>% filter(Component == "blank")
  if (nrow(blank_stats) > 0) {
    x_blank <- blank_stats$x_pos[1]
    ceiling_val <- blank_stats$background_ceiling[1]
    
    # 短线：仅覆盖blank组宽度（左右各0.3）
    p <- p +
      annotate("segment", 
               x = x_blank - 0.3, xend = x_blank + 0.3,
               y = ceiling_val, yend = ceiling_val,
               linetype = "dashed", color = "red", linewidth = 1) +
      # 数值标签
      annotate("text",
               x = x_blank + 0.35,  # 稍微偏右
               y = ceiling_val,
               label = paste0("Ceiling=", round(ceiling_val, 0)),
               hjust = 0, vjust = 0.5, size = 3, color = "red", fontface = "bold")
  }
  
  # 标注cfRNA、debris、cell组的上下界（仅覆盖各自组宽度）
  exp_groups <- c("cfRNA", "debris", "cell", "blank")
  exp_stats <- iqr_stats %>% filter(Component %in% exp_groups)
  
  for (i in seq_len(nrow(exp_stats))) {
    group_name <- exp_stats$Component[i]
    x_pos <- exp_stats$x_pos[i]
    lb <- exp_stats$lower_bound[i]
    ub <- exp_stats$upper_bound[i]
    q1_val <- exp_stats$q1[i]
    q3_val <- exp_stats$q3[i]
    
    # 下界线（绿色虚线，仅该组宽度）
    p <- p + 
      annotate("segment",
               x = x_pos - 0.3, xend = x_pos + 0.3,
               y = lb, yend = lb,
               linetype = "dashed", color = "darkgreen", linewidth = 0.8, alpha = 0.8) +
      # 下界数值标签（左侧）
      annotate("text",
               x = x_pos - 0.35,
               y = lb,
               label = paste0("L=", round(lb, 0)),
               hjust = 1, vjust = 0.5, size = 2.8, color = "darkgreen")
    
    # 上界线（紫色虚线，仅该组宽度）
    p <- p + 
      annotate("segment",
               x = x_pos - 0.3, xend = x_pos + 0.3,
               y = ub, yend = ub,
               linetype = "dashed", color = "purple", linewidth = 0.8, alpha = 0.8) +
      # 上界数值标签（右侧）
      annotate("text",
               x = x_pos + 0.35,
               y = ub,
               label = paste0("U=", round(ub, 0)),
               hjust = 0, vjust = 0.5, size = 2.8, color = "purple")
  }
  
  # 添加Q1和Q3数值标签（箱线图上）
  p <- p +
    geom_text(data = iqr_stats, aes(x = Component, y = q1, label = round(q1, 0)),
              vjust = 1.5, size = 3, color = "darkred", fontface = "bold") +
    geom_text(data = iqr_stats, aes(x = Component, y = q3, label = round(q3, 0)),
              vjust = -0.8, size = 3, color = "darkred", fontface = "bold")
  
  # 添加图例说明
  p <- p +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = "Red: Blank Ceiling (M+3SD)\nGreen: Lower Fence (Q1-1.5IQR)\nPurple: Upper Fence (Q3+1.5IQR)",
      hjust = 1, vjust = 1, size = 2.5, color = "gray30",
      xjust = 1, yjust = 1
    ) +
    labs(title = paste(y_label, "detection number by group"),
         x = "Group", y = y_label) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 100, 10, 60)  # 为左右标签留空间
    )
  
  # 如果提供了比较列表，添加显著性标记
  if (!is.null(comparisons) && length(comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons,
                                 method = "wilcox.test",
                                 label = "p.signif",
                                 tip.length = 0.02,
                                 bracket.size = 0.3,
                                 size = 4,
                                 step.increase = 0.05)
  }
  
  # 如果指定了文件名，保存图片（PNG + PDF）
  if (!is.null(filename)) {
    # 保存 PNG（原行为）
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
    # ---- 修改开始：同时保存 PDF ----
    pdf_filename <- sub("\\.png$", ".pdf", filename)
    ggsave(pdf_filename, p, width = 12, height = 8)
    # ---- 修改结束 ----
  }
  
  return(p)
}

# ------------------------------
# 4. 构建所需的比较列表
# ------------------------------
actual_groups <- levels(data$Component)[levels(data$Component) %in% unique(data$Component)]

# water 与其他组的比较
if ("water" %in% actual_groups) {
  water_comparisons <- map(setdiff(actual_groups, "water"), ~ c("water", .x))
} else {
  water_comparisons <- list()
}

# cell、debris、cfRNA 之间的两两比较
cell_debris_cfrna <- intersect(c("cell", "debris", "cfRNA"), actual_groups)
if (length(cell_debris_cfrna) >= 2) {
  cell_debris_comparisons <- combn(cell_debris_cfrna, 2, simplify = FALSE)
} else {
  cell_debris_comparisons <- list()
}

# 合并所有需要标记的比较（用于 mRNA、lncRNA、mlncRNA）
all_comparisons <- c(water_comparisons, cell_debris_comparisons)

# 对于其他 RNA 类型，只使用 water 与其他组的比较
other_comparisons <- water_comparisons

# ------------------------------
# 5. 对每种 RNA 类型执行分析和绘图
# ------------------------------

rna_list <- list(
  list(var = "mRNA",    label = "mRNA count"),
  list(var = "lncRNA",  label = "lncRNA count"),
  list(var = "mlncRNA", label = "mlncRNA count"),
  list(var = "miRNA",   label = "miRNA count"),
  list(var = "tRNA",    label = "tRNA count"),
  list(var = "piRNA",   label = "piRNA count")
)

plot_list <- list()

for (rna in rna_list) {
  var_name <- rna$var
  y_label  <- rna$label
  file_prefix <- var_name
  
  # 根据 RNA 类型选择比较列表
  if (var_name %in% c("mRNA", "lncRNA", "mlncRNA")) {
    comps <- all_comparisons
  } else {
    comps <- other_comparisons
  }
  
  # 生成图形并保存单个文件（函数内已同时保存 PNG 和 PDF）
  p <- plot_box_with_iqr(data, !!sym(var_name), y_label, 
                         comparisons = comps, 
                         filename = paste0(file_prefix, "_boxplot.png"))
  
  # 如果是前三种，保存到列表用于合并
  if (var_name %in% c("mRNA", "lncRNA", "mlncRNA")) {
    plot_list[[var_name]] <- p
  }
  
  # 进行两两显著性检验并输出表格
  pairwise_test(data, !!sym(var_name), var_name, file_prefix)
}

# ------------------------------
# 6. 合并 mRNA、lncRNA、mlncRNA 的三张图
# ------------------------------
g <- arrangeGrob(plot_list[["mRNA"]], 
                 plot_list[["lncRNA"]], 
                 plot_list[["mlncRNA"]], 
                 ncol = 3)

# 保存合并图（PNG）
ggsave("combined_mlRNA_boxplot.png", g, width = 28, height = 9, dpi = 300)
# ---- 修改开始：同时保存 PDF ----
ggsave("combined_mlRNA_boxplot.pdf", g, width = 28, height = 9)
# ---- 修改结束 ----