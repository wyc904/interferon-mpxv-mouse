


# 加载包
library(edgeR)
library(limma)
library(ggplot2)


# # 检查文件的前几行
# Y_A_B_counts<-readLines("D:/Desktop/Project/生信/3-差异基因分析/out/A-B-counts.txt", n = 10)
# # 读取 featureCounts 生成的文件
# counts <- Y_A_B_counts
# 
# head(Y_A_B_counts)
# 
# 
# 
# count_data <- counts[, 7:ncol(counts)] 
# 
# 
# D:\Desktop\10-干扰素小鼠\1-源
# 

counts<- read.delim("D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# print(dim(counts))
counts2 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/AGB6_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# print(dim(counts))
counts3 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# print(dim(counts))






colnames(counts)






# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:20]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix))





# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix2 <- counts2[, 6:19]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix2)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix2) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix2))





# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix3 <- counts3[, 6:20]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix3)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix3) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix3))







