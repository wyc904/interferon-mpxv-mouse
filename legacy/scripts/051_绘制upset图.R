# 加载必要的包
library(tidyverse)
library(UpSetR)
library(ComplexUpset)  # 更现代的UpSet图包
library(openxlsx)
library(viridis)
library(RColorBrewer)

data <- read.xlsx("绘制.xlsx", sheet = 1)

# 检查数据结构
cat("数据维度:", dim(data), "\n")
cat("列名:", colnames(data), "\n")
cat("time唯一值:", unique(data$time), "\n")
cat("mice唯一值:", unique(data$mice), "\n")

# 查看前几行
print(head(data))

# 1. 数据预处理 - 提取显著差异基因
# 定义显著性阈值
padj_threshold <- 0.05
log2fc_threshold <- 1

# 筛选显著差异基因
sig_data <- data %>%
  filter(!is.na(padj)) %>%
  mutate(
    direction = case_when(
      padj < padj_threshold & log2FoldChange > log2fc_threshold ~ "Up",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down",
      TRUE ~ "Not_sig"
    )
  ) %>%
  filter(direction != "Not_sig")  # 只保留显著差异基因

cat("\n显著差异基因总数:", nrow(sig_data), "\n")
cat("各条件显著基因数:\n")
print(table(sig_data$mice, sig_data$time, sig_data$direction))

# 2. 准备UpSet图数据
# 为每个条件创建基因列表
conditions <- c("c57_3d_up", "c57_5d_up", "8bc_3d_up", "8bc_5d_up", 
                "c57_3d_down", "c57_5d_down", "8bc_3d_down", "8bc_5d_down")

# 创建空的基因列表
gene_lists <- list()

# 填充各个条件的基因列表
for(cond in conditions) {
  parts <- str_split(cond, "_", simplify = TRUE)
  mice_type <- parts[1]
  time_point <- parts[2]
  direction_type <- parts[3]
  
  # 筛选符合条件的基因
  cond_genes <- sig_data %>%
    filter(mice == mice_type, 
           time == time_point, 
           direction == str_to_title(direction_type)) %>%
    pull(gene) %>%
    unique()
  
  gene_lists[[cond]] <- cond_genes
  cat(paste(cond, "基因数:", length(cond_genes), "\n"))
}

# 3. 创建二元矩阵用于UpSet图
# 获取所有显著差异基因
all_sig_genes <- unique(sig_data$gene)
cat("\n所有显著差异基因总数:", length(all_sig_genes), "\n")

# 创建二元矩阵
binary_matrix <- data.frame(Gene = all_sig_genes)

for(cond in conditions) {
  binary_matrix[[cond]] <- as.integer(binary_matrix$Gene %in% gene_lists[[cond]])
}

# 查看二元矩阵结构
cat("\n二元矩阵维度:", dim(binary_matrix), "\n")
print(head(binary_matrix[, 1:5]))

# 4. 使用UpSetR绘制基础UpSet图
# 4.1 上调基因的UpSet图
up_conditions <- conditions[str_detect(conditions, "_up$")]

cat("\n上调基因条件:", paste(up_conditions, collapse = ", "), "\n")

# 准备上调基因的二元矩阵
up_binary <- binary_matrix[, c("Gene", up_conditions)]

# 计算每个基因所在的条件数
up_binary$n_sets <- rowSums(up_binary[, up_conditions])

# 只保留至少在1个条件中上调的基因
up_binary <- up_binary[up_binary$n_sets > 0, ]

cat("上调基因矩阵维度:", dim(up_binary), "\n")

# 绘制上调基因UpSet图
pdf("UpSet_Upregulated_Genes.pdf", width = 12, height = 8)
upset(
  up_binary[, up_conditions],  # 只包含二元列
  nsets = length(up_conditions),
  nintersects = 20,  # 显示前20个交集
  mb.ratio = c(0.6, 0.4),  # 矩阵和条形图的比例
  order.by = "freq",  # 按频率排序
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.5,
  line.size = 1,
  mainbar.y.label = "基因数量",
  sets.x.label = "各条件基因数",
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1),  # 调整文本大小
  set_size.show = TRUE,  # 显示集合大小
  set_size.scale_max = max(colSums(up_binary[, up_conditions])) * 1.1
)
dev.off()

# 4.2 下调基因的UpSet图
down_conditions <- conditions[str_detect(conditions, "_down$")]

cat("\n下调基因条件:", paste(down_conditions, collapse = ", "), "\n")

# 准备下调基因的二元矩阵
down_binary <- binary_matrix[, c("Gene", down_conditions)]

# 计算每个基因所在的条件数
down_binary$n_sets <- rowSums(down_binary[, down_conditions])

# 只保留至少在1个条件中下调的基因
down_binary <- down_binary[down_binary$n_sets > 0, ]

cat("下调基因矩阵维度:", dim(down_binary), "\n")

# 绘制下调基因UpSet图
pdf("UpSet_Downregulated_Genes.pdf", width = 12, height = 8)
upset(
  down_binary[, down_conditions],
  nsets = length(down_conditions),
  nintersects = 20,
  mb.ratio = c(0.6, 0.4),
  order.by = "freq",
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.5,
  line.size = 1,
  mainbar.y.label = "基因数量",
  sets.x.label = "各条件基因数",
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1),
  set_size.show = TRUE,
  set_size.scale_max = max(colSums(down_binary[, down_conditions])) * 1.1
)
dev.off()

# 5. 使用ComplexUpset绘制更高级的UpSet图
# 5.1 准备所有条件的UpSet图（上下调分开显示）
# 创建绘图数据
upset_data <- binary_matrix %>%
  mutate(Gene = as.character(Gene)) %>%
  column_to_rownames("Gene")

# 修复ComplexUpset绘图部分
pdf("UpSet_All_Conditions.pdf", width = 14, height = 9)
ComplexUpset::upset(
  upset_data,
  intersect = conditions,
  name = "",
  width_ratio = 0.2,
  height_ratio = 0.8,
  min_size = 5,  # 最小交集大小
  sort_sets = FALSE
)
dev.off()

# 5.2 分别按时间和品系绘制UpSet图
# 5.2.1 按时间点绘制
time_conditions <- list(
  "3d_up" = intersect(gene_lists[["c57_3d_up"]], gene_lists[["8bc_3d_up"]]),
  "3d_down" = intersect(gene_lists[["c57_3d_down"]], gene_lists[["8bc_3d_down"]]),
  "5d_up" = intersect(gene_lists[["c57_5d_up"]], gene_lists[["8bc_5d_up"]]),
  "5d_down" = intersect(gene_lists[["c57_5d_down"]], gene_lists[["8bc_5d_down"]])
)

# 创建时间点的二元矩阵
time_genes <- unique(unlist(time_conditions))
time_binary <- data.frame(Gene = time_genes)
for(cond in names(time_conditions)) {
  time_binary[[cond]] <- as.integer(time_binary$Gene %in% time_conditions[[cond]])
}

cat("\n时间点交集基因数:\n")
for(cond in names(time_conditions)) {
  cat(paste(cond, ":", length(time_conditions[[cond]]), "\n"))
}

# 绘制时间点UpSet图
pdf("UpSet_By_Time.pdf", width = 10, height = 7)
upset(
  time_binary[, -1],  # 移除Gene列
  nsets = length(time_conditions),
  nintersects = 10,
  mb.ratio = c(0.6, 0.4),
  order.by = "freq",
  mainbar.y.label = "基因数量",
  sets.x.label = "各时间点基因数",
  text.scale = c(1.2, 1.2, 1, 1, 1.2, 1)
)
dev.off()

# 5.2.2 按品系绘制
strain_conditions <- list(
  "C57_up" = union(gene_lists[["c57_3d_up"]], gene_lists[["c57_5d_up"]]),
  "C57_down" = union(gene_lists[["c57_3d_down"]], gene_lists[["c57_5d_down"]]),
  "8BC_up" = union(gene_lists[["8bc_3d_up"]], gene_lists[["8bc_5d_up"]]),
  "8BC_down" = union(gene_lists[["8bc_3d_down"]], gene_lists[["8bc_5d_down"]])
)

# 创建品系的二元矩阵
strain_genes <- unique(unlist(strain_conditions))
strain_binary <- data.frame(Gene = strain_genes)
for(cond in names(strain_conditions)) {
  strain_binary[[cond]] <- as.integer(strain_binary$Gene %in% strain_conditions[[cond]])
}

cat("\n品系基因数:\n")
for(cond in names(strain_conditions)) {
  cat(paste(cond, ":", length(strain_conditions[[cond]]), "\n"))
}

# 绘制品系UpSet图
pdf("UpSet_By_Strain.pdf", width = 10, height = 7)
upset(
  strain_binary[, -1],
  nsets = length(strain_conditions),
  nintersects = 10,
  mb.ratio = c(0.6, 0.4),
  order.by = "freq",
  mainbar.y.label = "基因数量",
  sets.x.label = "各品系基因数",
  text.scale = c(1.2, 1.2, 1, 1, 1.2, 1)
)
dev.off()

# 6. 统计分析函数
# 计算各种交集的统计信息
calculate_intersection_stats <- function(gene_lists) {
  # 获取所有条件
  all_conditions <- names(gene_lists)
  
  # 初始化结果数据框
  stats_df <- data.frame(
    Intersection = character(),
    Genes = integer(),
    Gene_List = character()
  )
  
  # 计算所有可能的交集
  for(i in 1:length(all_conditions)) {
    # 单个条件的组合
    combs <- combn(all_conditions, i, simplify = FALSE)
    
    for(comb in combs) {
      # 计算交集
      if(length(comb) == 1) {
        intersect_genes <- gene_lists[[comb]]
      } else {
        intersect_genes <- Reduce(intersect, gene_lists[comb])
      }
      
      # 保存结果
      intersection_name <- paste(comb, collapse = " ∩ ")
      stats_df <- rbind(stats_df, data.frame(
        Intersection = intersection_name,
        Genes = length(intersect_genes),
        Gene_List = paste(intersect_genes, collapse = ";")
      ))
    }
  }
  
  return(stats_df)
}

# 计算上调基因的交集统计
up_stats <- calculate_intersection_stats(gene_lists[str_detect(names(gene_lists), "_up$")])
cat("\n上调基因交集统计（前10个）:\n")
print(head(up_stats[order(-up_stats$Genes), ], 10))

# 计算下调基因的交集统计
down_stats <- calculate_intersection_stats(gene_lists[str_detect(names(gene_lists), "_down$")])
cat("\n下调基因交集统计（前10个）:\n")
print(head(down_stats[order(-down_stats$Genes), ], 10))

# 保存统计结果
write.csv(up_stats, "Upregulated_Intersection_Stats.csv", row.names = FALSE)
write.csv(down_stats, "Downregulated_Intersection_Stats.csv", row.names = FALSE)
write.csv(binary_matrix, "Binary_Matrix_All_Genes.csv", row.names = FALSE)

# 7. 创建总结报告
cat("\n", strrep("=", 60), "\n", sep = "")
cat("UpSet分析总结报告\n")
cat(strrep("=", 60), "\n", sep = "")
cat("分析条件: C57和8BC小鼠，3d和5d时间点\n")
cat("显著性阈值: padj <", padj_threshold, ", |log2FC| >", log2fc_threshold, "\n")
cat("总显著差异基因数:", nrow(sig_data), "\n")
cat("上调基因总数:", sum(sig_data$direction == "Up"), "\n")
cat("下调基因总数:", sum(sig_data$direction == "Down"), "\n")
cat("\n各条件基因数:\n")
for(cond in conditions) {
  cat(paste(cond, ":", length(gene_lists[[cond]]), "\n"))
}
cat("\n关键交集发现:\n")
cat("1. 所有条件共同上调基因:", 
    length(Reduce(intersect, gene_lists[str_detect(names(gene_lists), "_up$")])), "\n")
cat("2. 所有条件共同下调基因:", 
    length(Reduce(intersect, gene_lists[str_detect(names(gene_lists), "_down$")])), "\n")
cat("3. C57特有上调基因:", 
    length(setdiff(
      union(gene_lists[["c57_3d_up"]], gene_lists[["c57_5d_up"]]),
      union(gene_lists[["8bc_3d_up"]], gene_lists[["8bc_5d_up"]])
    )), "\n")
cat("4. 8BC特有上调基因:", 
    length(setdiff(
      union(gene_lists[["8bc_3d_up"]], gene_lists[["8bc_5d_up"]]),
      union(gene_lists[["c57_3d_up"]], gene_lists[["c57_5d_up"]])
    )), "\n")

cat("\n生成的文件:\n")
cat("1. UpSet_Upregulated_Genes.pdf - 上调基因UpSet图\n")
cat("2. UpSet_Downregulated_Genes.pdf - 下调基因UpSet图\n")
cat("3. UpSet_All_Conditions.pdf - 所有条件UpSet图\n")
cat("4. UpSet_By_Time.pdf - 按时间点UpSet图\n")
cat("5. UpSet_By_Strain.pdf - 按品系UpSet图\n")
cat("6. Upregulated_Intersection_Stats.csv - 上调基因交集统计\n")
cat("7. Downregulated_Intersection_Stats.csv - 下调基因交集统计\n")
cat("8. Binary_Matrix_All_Genes.csv - 二元矩阵数据\n")

cat(strrep("=", 60), "\n", sep = "")