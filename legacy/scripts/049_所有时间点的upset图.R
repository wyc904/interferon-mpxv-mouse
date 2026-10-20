# 加载必要的包
library(tidyverse)
library(UpSetR)
library(openxlsx)

# 读取数据
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

# 2. 准备总的差异基因数据（不区分上调下调）
# 定义所有时间点和小鼠品系
time_points <- c("3d", "5d", "7d", "10d")
mouse_strains <- c("c57", "8bc")

# 创建总的差异基因条件（不区分上调下调）
total_conditions <- c()
for(mice_type in mouse_strains) {
  for(time_point in time_points) {
    total_conditions <- c(total_conditions, paste0(mice_type, "_", time_point))
  }
}

cat("\n总条件数:", length(total_conditions), "\n")
cat("条件列表:", paste(total_conditions, collapse = ", "), "\n")

# 创建总的基因列表（不区分上调下调）
total_gene_lists <- list()

# 填充各个条件的基因列表（合并上下调）
for(cond in total_conditions) {
  parts <- str_split(cond, "_", simplify = TRUE)
  mice_type <- parts[1]
  time_point <- parts[2]
  
  # 筛选该条件下所有显著差异基因（上调和下调的并集）
  up_genes <- sig_data %>%
    filter(mice == mice_type, 
           time == time_point, 
           direction == "Up") %>%
    pull(gene) %>%
    unique()
  
  down_genes <- sig_data %>%
    filter(mice == mice_type, 
           time == time_point, 
           direction == "Down") %>%
    pull(gene) %>%
    unique()
  
  # 合并上下调基因
  total_genes <- union(up_genes, down_genes)
  total_gene_lists[[cond]] <- total_genes
  cat(paste(cond, "总差异基因数:", length(total_genes), "\n"))
}

# 3. 创建总的二元矩阵用于UpSet图
# 获取所有显著差异基因
all_sig_genes <- unique(sig_data$gene)
cat("\n所有显著差异基因总数:", length(all_sig_genes), "\n")

# 创建总的二元矩阵
total_binary <- data.frame(Gene = all_sig_genes)

for(cond in total_conditions) {
  total_binary[[cond]] <- as.integer(total_binary$Gene %in% total_gene_lists[[cond]])
}

# 查看二元矩阵结构
cat("\n总差异基因二元矩阵维度:", dim(total_binary), "\n")

# 4. 绘制总的差异基因UpSet图（不区分上调下调）
pdf("UpSet_Total_DEGs_All_Timepoints.pdf", width = 14, height = 10)
upset(
  total_binary[, -1],  # 移除Gene列，只保留二元数据
  nsets = length(total_conditions),
  nintersects = 30,  # 显示前30个交集
  mb.ratio = c(0.6, 0.4),  # 矩阵和条形图的比例
  order.by = "freq",  # 按频率排序
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.5,
  line.size = 1,
  mainbar.y.label = "基因数量",
  sets.x.label = "各条件总差异基因数",
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1),
  set_size.show = TRUE,
  set_size.scale_max = max(colSums(total_binary[, -1])) * 1.1
)
dev.off()

# 5. 创建简要总结报告
cat("\n", strrep("=", 70), "\n", sep = "")
cat("总差异基因UpSet图分析报告\n")
cat(strrep("=", 70), "\n", sep = "")
cat("分析条件: C57和8BC小鼠，3d、5d、7d、10d四个时间点\n")
cat("显著性阈值: padj <", padj_threshold, ", |log2FC| >", log2fc_threshold, "\n")
cat("总显著差异基因数:", nrow(sig_data), "\n")

cat("\n各条件总差异基因数（不区分上下调）:\n")
for(cond in total_conditions) {
  cat(paste(cond, ":", length(total_gene_lists[[cond]]), "\n"))
}

cat("\n关键交集发现:\n")
# 所有条件共同差异基因
all_common <- Reduce(intersect, total_gene_lists)
cat("所有8个条件共同差异基因数:", length(all_common), "\n")

# 各品系内所有时间点共同差异基因
for(mice_type in mouse_strains) {
  strain_conditions <- total_conditions[str_detect(total_conditions, paste0("^", mice_type))]
  strain_common <- Reduce(intersect, total_gene_lists[strain_conditions])
  cat(paste0(str_to_upper(mice_type), "所有时间点共同差异基因数: ", length(strain_common), "\n"))
}

# 各时间点品系间共同差异基因
for(time_point in time_points) {
  time_conditions <- total_conditions[str_detect(total_conditions, paste0("_", time_point, "$"))]
  time_common <- Reduce(intersect, total_gene_lists[time_conditions])
  cat(paste0(time_point, "时间点品系间共同差异基因数: ", length(time_common), "\n"))
}

cat("\n生成的文件:\n")
cat("UpSet_Total_DEGs_All_Timepoints.pdf - 总差异基因UpSet图（不区分上下调）\n")

cat(strrep("=", 70), "\n", sep = "")