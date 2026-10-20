############################################################
# 40 个免疫轨迹基因热图（横向布局版本）- 按指定顺序排列
# 适配edgeR处理的bulk测序结果
############################################################

library(edgeR)
library(limma)
library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(tidyr)

## 1. 读取和处理表达矩阵 --------------------------------------------------------

# 读取counts文件
counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/8bc_counts_gene", 
                     header = TRUE, sep = "\t", 
                     comment.char = "#", row.names = 1)

cat("原始数据维度:", dim(counts), "\n")
cat("列名:\n")
print(colnames(counts))

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:18]




# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix))

# 创建样本信息表
sample_info <- data.frame(
  sample = simplified_names,
  group = sapply(strsplit(simplified_names, "_"), function(x) paste(x[2:3], collapse = "_")),
  replicate = sapply(strsplit(simplified_names, "_"), function(x) x[4]),
  stringsAsFactors = FALSE
)

cat("\n=== 样本信息 ===\n")
print(sample_info)

# 检查count数据的统计信息
cat("\n=== Count数据统计 ===\n")
cat("总基因数:", nrow(count_matrix), "\n")
cat("总测序reads数:", sum(count_matrix), "\n")
cat("每个样本的reads数:\n")
print(colSums(count_matrix))

# 检查是否有零表达基因
zero_genes <- rowSums(count_matrix) == 0
cat("零表达基因数:", sum(zero_genes), "\n")
cat("占比:", round(mean(zero_genes) * 100, 2), "%\n")











