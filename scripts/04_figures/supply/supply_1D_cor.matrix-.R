setwd("/data/work/01_2603cell_culture/02_qc_batch/02_batch_PCA")
library(corrplot)
library(psych)
library(pheatmap)
library(dplyr)

rna <- "mlRNA"

# 添加错误处理机制，确保文件存在
if(!file.exists(paste0("all_", rna, "_TPM.txt"))) {
  stop(paste("错误：文件 all_", rna, "_TPM.txt 不存在，请确认文件路径和文件名！"))
}

# 读取表达数据
exp <- read.table(paste0("all_", rna, "_TPM.txt"),
                  row.names = 1, header = TRUE, sep = "\t")

# 数据处理
exp <- data.frame(t(exp))
exp1 <- arrange(exp, row.names(exp))
dat <- data.frame(t(exp1))

# 计算相关性矩阵
cor1 <- cor(dat, method = "spearman")

# 修复文件名问题：原来的文件名有错误的连字符
base_filename <- paste0("cor-", rna)

# 动态计算高度：根据样本数量自适应调整高度
n_samples <- ncol(dat)
height_pdf <- max(8, min(10, n_samples * 0.15)) # 每增加10个样本增加1英寸高度
height_png <- height_pdf * 100 # PNG的高度（像素）

# 导出PDF - 修复文件名并添加图形设备安全措施
pdf(file = paste0(base_filename, ".pdf"), 
    width = 8, height = height_pdf)

pheatmap(cor1, 
         scale = "none", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = colorRampPalette(colors = c("blue", "white", "red"))(100),
         legend = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         treeheight_row = 40,
         treeheight_col = 40,
         border_color = NA,
         fontsize = 15,
         fontsize_row = 8)

# 确保正确关闭PDF设备
dev.off()

# 导出PNG格式 - 提供另一种输出格式
png(file = paste0(base_filename, ".png"), 
    width = 2000, 
    height = height_png, 
    res = 150) # 150 DPI分辨率

pheatmap(cor1, 
         scale = "none", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = colorRampPalette(colors = c("blue", "white", "red"))(100),
         legend = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         treeheight_row = 40,
         treeheight_col = 40,
         border_color = NA,
         fontsize = 15, # PNG中增加字体大小
         fontsize_row = 8) # 行文字稍大

# 确保关闭PNG设备
dev.off()

# 完成提示
cat("相关性热图已成功导出为两种格式：\n")
cat("1. PDF: ", paste0(base_filename, ".pdf"), "\n")
cat("2. PNG: ", paste0(base_filename, ".png"), "\n")
cat("文件高度已根据样本数量自动调整为：", height_pdf, "英寸 (PDF) / ", height_png, "像素 (PNG)\n")