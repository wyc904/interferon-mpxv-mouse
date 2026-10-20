# 载入必要的库
library(openxlsx)  # 用于读取Excel文件
library(ggplot2)  # 用于绘图
library(dplyr)  # 数据处理
library(tidyr)  # 数据整理
library(stringr)  # 字符串处理
library(viridis)  # 颜色方案
library(forcats)  # 因子处理

# 读取GO分析数据
file_path <- "D:/Desktop/8bc_MPXV_GO_Analysis.xlsx"
go_data <- read.xlsx(file_path, sheet = 1)

# 检查数据结构
cat("数据维度:", dim(go_data), "\n")
cat("列名:", colnames(go_data), "\n")

# 检查Cluster列的内容
cat("Cluster列的唯一值:\n")
if("Cluster" %in% colnames(go_data)) {
  print(table(go_data$Cluster))
} else {
  cat("错误: 数据中不存在Cluster列!\n")
  # 检查实际的列名
  cat("实际列名:", colnames(go_data), "\n")
  stop("请检查数据文件，确保包含Cluster列")
}

# 选择与目标GO通路相关的描述
target_go_terms <- c(
  "antiviral innate immune response", "defense response to virus", "type I interferon production", 
  "cellular response to interferon-beta", "chemokine production", "chemokine (C-X-C motif) ligand 2 production", 
  "positive regulation of chemokine production", "immune response-activating signaling pathway", 
  "negative regulation of viral genome replication", "mitotic nuclear division", "chromosome segregation", 
  "positive regulation of cell cycle", "negative regulation of cytokine-mediated signaling pathway", 
  "negative regulation of response to cytokine stimulus", "negative regulation of viral process", 
  "wound healing", "positive regulation of DNA repair", "T cell differentiation", 
  "positive regulation of lymphocyte activation", "regulation of actin cytoskeleton organization"
)

# 从数据中筛选出目标GO通路
filtered_go_data <- go_data %>%
  filter(Description %in% target_go_terms)

# 确保数值列正确
numeric_columns <- c("RichFactor", "p.adjust", "Count", "FoldEnrichment")
for(col in numeric_columns) {
  if(col %in% colnames(filtered_go_data)) {
    filtered_go_data[[col]] <- as.numeric(as.character(filtered_go_data[[col]]))
  }
}

# 排除包含NA或非数值的行
filtered_go_data <- filtered_go_data %>%
  filter(if_all(all_of(numeric_columns), ~!is.na(.)))

# 设置显著性阈值
filtered_go_data <- filtered_go_data %>%
  filter(p.adjust < 0.05)

# 计算-log10(p.adjust)
filtered_go_data$LogP <- -log10(filtered_go_data$p.adjust)

# 确保Cluster列正确
available_clusters <- as.character(unique(filtered_go_data$Cluster))
available_clusters <- available_clusters[!is.na(available_clusters)]

# 按期望的顺序设置Cluster水平
expected_clusters <- c("3d", "5d", "7d")
actual_clusters <- intersect(expected_clusters, available_clusters)

if(length(actual_clusters) > 0) {
  filtered_go_data$Cluster <- factor(filtered_go_data$Cluster, levels = actual_clusters)
} else {
  filtered_go_data$Cluster <- factor(filtered_go_data$Cluster)
}

# 检查数据
cat("最终数据维度:", dim(filtered_go_data), "\n")
cat("可用的Cluster:", levels(filtered_go_data$Cluster), "\n")
cat("显著通路数量:", nrow(filtered_go_data), "\n")
# 
# # 简化版1: 基本点图，横坐标为时间点
# p_simple <- ggplot(filtered_go_data, aes(x = Cluster, y = fct_reorder(Description, RichFactor))) +
#   geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
#   scale_color_viridis_c(name = "-Log10(Adj. P)") +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     title = "GO Enrichment Analysis Across Time Points",
#     x = "Time Point", 
#     y = "GO Terms"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#     axis.text.y = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, size = 14),
#     legend.position = "right"
#   )
# 
# print(p_simple)

# 简化版2: 去除所有复杂标识的基本点图
p_clean <- ggplot(filtered_go_data, aes(x = Cluster, y = fct_reorder(Description, RichFactor))) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  scale_color_viridis_c(name = "-Log10(P.adj)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(
    x = "Time Point", 
    y = "GO Terms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p_clean)



# 
# # 简化版3: 使用自定义颜色方案
# p_custom_colors <- ggplot(filtered_go_data, aes(x = Cluster, y = fct_reorder(Description, RichFactor))) +
#   geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
#   scale_color_gradientn(
#     colors = c("#E1F5FE", "#81D4FA", "#29B6F6", "#0288D1", "#01579B"),
#     name = "-Log10(P.adj)"
#   ) +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     x = "Time Point", 
#     y = "GO Terms"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#     axis.text.y = element_text(size = 10),
#     panel.grid.major = element_line(color = "grey90")
#   )
# 
# print(p_custom_colors)




# 简化版4: 水平布局，更适合长通路名称
p_horizontal <- ggplot(filtered_go_data, aes(y = Cluster, x = fct_reorder(Description, RichFactor))) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  scale_color_viridis_c(name = "-Log10(P.adj)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(
    x = "GO Terms", 
    y = "Time Point"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  ) +
  coord_flip()

print(p_horizontal)


# 
# 
# # 简化版5: 最简风格，适合发表
# p_minimal <- ggplot(filtered_go_data, aes(x = Cluster, y = fct_reorder(Description, RichFactor))) +
#   geom_point(aes(size = Count, fill = LogP), alpha = 0.8, shape = 21, color = "black") +
#   scale_fill_viridis_c(name = "-Log10(P.adj)") +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     x = "Days Post Infection", 
#     y = "GO Biological Processes"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 10),
#     axis.title = element_text(size = 12, face = "bold"),
#     panel.grid.major = element_line(color = "grey80"),
#     legend.position = "right"
#   )
# 
# print(p_minimal)

# 保存详细结果到Excel（简化版）
wb <- createWorkbook()

# 添加所有筛选结果
addWorksheet(wb, "Filtered_GO_Terms")
writeData(wb, sheet = "Filtered_GO_Terms", filtered_go_data, rowNames = FALSE)

# 按时间点添加工作表
for(cluster in levels(filtered_go_data$Cluster)) {
  cluster_data <- filtered_go_data %>% filter(Cluster == cluster)
  if(nrow(cluster_data) > 0) {
    addWorksheet(wb, paste("GO_Terms", cluster, sep = "_"))
    writeData(wb, sheet = paste("GO_Terms", cluster, sep = "_"), cluster_data, rowNames = FALSE)
  }
}

# 添加汇总统计
addWorksheet(wb, "Summary")
writeData(wb, sheet = "Summary", timepoint_summary, rowNames = FALSE)

# 保存Excel文件
saveWorkbook(wb, "D:/Desktop/GO_Enrichment_Results_Simplified.xlsx", overwrite = TRUE)

# 打印分析摘要
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("GO富集分析完成！\n")
cat("分析的时间点:", paste(levels(filtered_go_data$Cluster), collapse = ", "), "\n")
cat("显著富集通路总数:", nrow(filtered_go_data), "\n")
cat("各时间点通路数:\n")
for(cluster in levels(filtered_go_data$Cluster)) {
  count <- sum(filtered_go_data$Cluster == cluster)
  cat("  ", cluster, ":", count, "个通路\n")
}














