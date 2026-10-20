# 载入必要的库
library(openxlsx)  # 用于读取Excel文件
library(ggplot2)  # 用于绘图
library(dplyr)  # 数据处理
library(tidyr)  # 数据整理
library(stringr)  # 字符串处理
library(viridis)  # 颜色方案
library(forcats)  # 因子处理
library(patchwork)  # 用于多图组合

# 读取GO分析数据
file_path <- "D:/Desktop/go分析总.xlsx"
go_data <- read.xlsx(file_path, sheet = 1)

# 检查数据结构
cat("数据维度:", dim(go_data), "\n")
cat("列名:", colnames(go_data), "\n")

# 检查Cluster列和mice列的内容
cat("Cluster列的唯一值:\n")
if("Cluster" %in% colnames(go_data)) {
  print(table(go_data$Cluster))
} else {
  cat("错误: 数据中不存在Cluster列!\n")
  stop("请检查数据文件，确保包含Cluster列")
}

cat("Mice列的唯一值:\n")
if("mice" %in% colnames(go_data)) {
  print(table(go_data$mice))
} else {
  cat("错误: 数据中不存在mice列!\n")
  stop("请检查数据文件，确保包含mice列")
}

# 指定GO通路的顺序
target_go_terms <- c(
  "antiviral innate immune response", "defense response to virus", "type I interferon production", 
  "cellular response to interferon-beta", "chemokine production", "chemokine (C-X-C motif) ligand 2 production", 
  "positive regulation of chemokine production", "immune response-activating signaling pathway", 
  "negative regulation of viral genome replication","type I interferon-mediated signaling pathway", "toll-like receptor 3 signaling pathway","regulation of viral life cycle",
  
  
  
  "mitotic nuclear division", "chromosome segregation", 
  "positive regulation of cell cycle",
  
  "negative regulation of cytokine-mediated signaling pathway", 
  "negative regulation of response to cytokine stimulus", "negative regulation of viral process", 
  "wound healing", "positive regulation of DNA repair"
  
  
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

# 确保Cluster列和mice列正确
available_clusters <- as.character(unique(filtered_go_data$Cluster))
available_clusters <- available_clusters[!is.na(available_clusters)]

# 按期望的顺序设置Cluster水平
expected_clusters <- c("3d", "5d", "7d", "10d")
actual_clusters <- intersect(expected_clusters, available_clusters)

if(length(actual_clusters) > 0) {
  filtered_go_data$Cluster <- factor(filtered_go_data$Cluster, levels = actual_clusters)
} else {
  filtered_go_data$Cluster <- factor(filtered_go_data$Cluster)
}

# 确保mice列正确
filtered_go_data$mice <- factor(filtered_go_data$mice, levels = c("8bc", "c57"))

# 按指定顺序设置GO通路顺序
# 首先确保所有在数据中出现的GO通路都在指定的列表中
existing_go_terms <- intersect(target_go_terms, unique(filtered_go_data$Description))
# 按指定的顺序设置水平
filtered_go_data$Description <- factor(filtered_go_data$Description, 
                                       levels = rev(existing_go_terms))  # 使用rev反转顺序，使第一个在顶部

# 检查数据
cat("最终数据维度:", dim(filtered_go_data), "\n")
cat("可用的Cluster:", levels(filtered_go_data$Cluster), "\n")
cat("可用的Mice:", levels(filtered_go_data$mice), "\n")
cat("显著通路数量:", nrow(filtered_go_data), "\n")
cat("GO通路顺序（从上到下）:\n")
cat(paste(rev(levels(filtered_go_data$Description)), collapse = "\n"), "\n")











# 定义红色渐变色 (浅红 -> 深红)
red_colors <- c("#FFE5E5", "#FF9999", "#FF6666", "#FF3333", "#CC0000", "#990000", "#660000")

# 方法1: 分面显示两种品系
p_faceted <- ggplot(filtered_go_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  facet_grid(. ~ mice, scales = "free", space = "free") +
  scale_color_gradientn(
    colors = red_colors,
    name = "-Log10(P.adj)",
    guide = guide_colorbar(
      title.position = "top",  # 图例标题在顶部
      barwidth = unit(5, "cm"),  # 设置颜色条宽度
      barheight = unit(0.5, "cm"),  # 设置颜色条高度
      direction = "horizontal"  # 颜色条水平方向
    )
  ) +
  scale_size_continuous(
    range = c(3, 8), 
    name = "Gene Count",
    guide = guide_legend(
      title.position = "top",  # 图例标题在顶部
      direction = "horizontal"  # 图例项水平排列
    )
  ) +
  labs(
    title = "",
    x = "Time Point", 
    y = "GO Terms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",  # 图例在底部
    legend.box = "horizontal",  # 水平排列图例
    legend.box.just = "center",  # 图例居中
    legend.direction = "horizontal",  # 水平方向
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1, "cm"),
    legend.margin = margin(t = 5, r = 0, b = 0, l = 0)  # 调整图例间距
  )

print(p_faceted)




















# 
# 
# # 方法3: 水平布局，更适合长通路名称
# p_horizontal <- ggplot(filtered_go_data, aes(y = Cluster, x = Description)) +
#   geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
#   facet_grid(mice ~ ., scales = "free", space = "free") +
#   scale_color_viridis_c(name = "-Log10(P.adj)") +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     title = "",
#     x = "GO Terms", 
#     y = "Time Point"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#     axis.text.y = element_text(size = 10),
#     strip.text = element_text(face = "bold", size = 12),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
#   ) +
#   coord_flip()
# 
# print(p_horizontal)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Comparison_Horizontal.pdf", p_horizontal, width = 16, height = 10)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Comparison_Horizontal.png", p_horizontal, width = 16, height = 10, dpi = 300)
# 














# 
# 
# # 方法4: 分别绘制两种品系的图，便于详细比较
# # 筛选8bc品系数据
# p_8bc <- ggplot(filtered_go_data %>% filter(mice == "8bc"), 
#                 aes(x = Cluster, y = Description)) +
#   geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
#   scale_color_viridis_c(name = "-Log10(P.adj)") +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     title = "8bc Mice: GO Enrichment Analysis",
#     x = "Time Point", 
#     y = "GO Terms"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#     axis.text.y = element_text(size = 9),
#     plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
#   )
# 
# # 筛选C57品系数据
# p_c57 <- ggplot(filtered_go_data %>% filter(mice == "c57"), 
#                 aes(x = Cluster, y = Description)) +
#   geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
#   scale_color_viridis_c(name = "-Log10(P.adj)") +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     title = "C57 Mice: GO Enrichment Analysis",
#     x = "Time Point", 
#     y = "GO Terms"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#     axis.text.y = element_text(size = 9),
#     plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
#   )
# 
# # 使用patchwork包组合两个图
# if(require(patchwork)) {
#   p_combined_side <- p_8bc + p_c57 + 
#     plot_layout(ncol = 2, guides = "collect") +
#     plot_annotation(title = "Comparison of GO Enrichment in 8bc vs C57 Mice",
#                     theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
#   
#   print(p_combined_side)
#   ggsave("D:/Desktop/GO_Enrichment_Mice_Side_by_Side.pdf", p_combined_side, width = 20, height = 10)
#   ggsave("D:/Desktop/GO_Enrichment_Mice_Side_by_Side.png", p_combined_side, width = 20, height = 10, dpi = 300)
# } else {
#   print(p_8bc)
#   ggsave("D:/Desktop/GO_Enrichment_8bc_Mice.pdf", p_8bc, width = 10, height = 8)
#   ggsave("D:/Desktop/GO_Enrichment_8bc_Mice.png", p_8bc, width = 10, height = 8, dpi = 300)
#   
#   print(p_c57)
#   ggsave("D:/Desktop/GO_Enrichment_C57_Mice.pdf", p_c57, width = 10, height = 8)
#   ggsave("D:/Desktop/GO_Enrichment_C57_Mice.png", p_c57, width = 10, height = 8, dpi = 300)
# }
# 
# 
# 
# 



















# 
# 
# # 方法5: 按生物学功能分组的分面图
# # 添加生物学功能分组
# filtered_go_data <- filtered_go_data %>%
#   mutate(
#     Functional_Group = case_when(
#       Description %in% c("antiviral innate immune response", "defense response to virus", 
#                          "type I interferon production", "cellular response to interferon-beta") ~ "1. Antiviral Response",
#       Description %in% c("chemokine production", "chemokine (C-X-C motif) ligand 2 production",
#                          "positive regulation of chemokine production", "immune response-activating signaling pathway",
#                          "negative regulation of viral genome replication") ~ "2. Inflammatory Response",
#       Description %in% c("mitotic nuclear division", "chromosome segregation", 
#                          "positive regulation of cell cycle") ~ "3. Cell Cycle",
#       Description %in% c("negative regulation of cytokine-mediated signaling pathway",
#                          "negative regulation of response to cytokine stimulus", 
#                          "negative regulation of viral process") ~ "4. Negative Regulation",
#       Description %in% c("wound healing", "positive regulation of DNA repair", "T cell differentiation",
#                          "positive regulation of lymphocyte activation", 
#                          "regulation of actin cytoskeleton organization") ~ "5. Tissue Repair & Immune Reconstitution",
#       TRUE ~ "6. Other"
#     ),
#     # 确保功能分组按指定顺序
#     Functional_Group = factor(Functional_Group, 
#                               levels = c("1. Antiviral Response", "2. Inflammatory Response", 
#                                          "3. Cell Cycle", "4. Negative Regulation", 
#                                          "5. Tissue Repair & Immune Reconstitution", "6. Other"))
#   )
# 
# 





# # 按功能分组的分面图
# p_grouped <- ggplot(filtered_go_data, aes(x = Cluster, y = Description)) +
#   geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
#   facet_grid(Functional_Group ~ mice, scales = "free_y", space = "free_y") +
#   scale_color_viridis_c(name = "-Log10(P.adj)") +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     title = "GO Enrichment Analysis by Functional Group: 8bc vs C57 Mice",
#     x = "Time Point", 
#     y = "GO Terms"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#     axis.text.y = element_text(size = 8),
#     strip.text.y = element_text(angle = 0, hjust = 0, size = 10, face = "bold"),
#     strip.text.x = element_text(size = 10, face = "bold"),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
#   )
# 
# print(p_grouped)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Grouped.pdf", p_grouped, width = 16, height = 12)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Grouped.png", p_grouped, width = 16, height = 12, dpi = 300)
# 











# 
# 
# 
# 
# # 方法2: 用颜色区分品系，合并显示
# p_combined <- ggplot(filtered_go_data, aes(x = Cluster, y = Description)) +
#   geom_point(aes(size = Count, color = mice, alpha = LogP), position = position_dodge(width = 0.7)) +
#   scale_color_manual(values = c("8bc" = "#E69F00", "c57" = "#56B4E9"), name = "Mouse Strain") +
#   scale_alpha_continuous(name = "-Log10(P.adj)", range = c(0.6, 1)) +
#   scale_size_continuous(range = c(3, 8), name = "Gene Count") +
#   labs(
#     title = "GO Enrichment Analysis: Comparison of 8bc and C57 Mice",
#     x = "Time Point", 
#     y = "GO Terms"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#     axis.text.y = element_text(size = 9),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#     legend.position = "right"
#   )
# 
# print(p_combined)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Comparison_Combined.pdf", p_combined, width = 12, height = 10)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Comparison_Combined.png", p_combined, width = 12, height = 10, dpi = 300)
# 
# 
# 
# 



















# 创建汇总表格
timepoint_summary <- filtered_go_data %>%
  group_by(mice, Cluster) %>%
  summarise(
    n_pathways = n(),
    total_genes = sum(Count),
    avg_rich_factor = mean(RichFactor),
    .groups = "drop"
  )

# 详细比较表格
detailed_comparison <- filtered_go_data %>%
  group_by(Description, mice) %>%
  summarise(
    n_timepoints = n_distinct(Cluster),
    timepoints = paste(sort(unique(Cluster)), collapse = ", "),
    avg_count = mean(Count),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = mice,
    values_from = c(n_timepoints, timepoints, avg_count),
    names_sep = "_"
  )

print("按品系和时间点汇总:")
print(timepoint_summary)

# 保存详细结果到Excel
wb <- createWorkbook()

# 添加所有筛选结果
addWorksheet(wb, "All_Filtered_GO_Terms")
writeData(wb, sheet = "All_Filtered_GO_Terms", filtered_go_data, rowNames = FALSE)

# 按品系和时间点添加工作表
for(mouse_strain in levels(filtered_go_data$mice)) {
  for(cluster in levels(filtered_go_data$Cluster)) {
    strain_cluster_data <- filtered_go_data %>% 
      filter(mice == mouse_strain, Cluster == cluster)
    if(nrow(strain_cluster_data) > 0) {
      sheet_name <- paste("GO_Terms", mouse_strain, cluster, sep = "_")
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet = sheet_name, strain_cluster_data, rowNames = FALSE)
    }
  }
}

# 添加汇总统计
addWorksheet(wb, "Timepoint_Summary")
writeData(wb, sheet = "Timepoint_Summary", timepoint_summary, rowNames = FALSE)

addWorksheet(wb, "Detailed_Comparison")
writeData(wb, sheet = "Detailed_Comparison", detailed_comparison, rowNames = FALSE)

# 保存Excel文件
saveWorkbook(wb, "D:/Desktop/GO_Enrichment_Mice_Comparison_Results.xlsx", overwrite = TRUE)

# 打印分析摘要
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("两种品系小鼠GO富集分析完成！\n")
cat("分析的时间点:", paste(levels(filtered_go_data$Cluster), collapse = ", "), "\n")
cat("分析的品系:", paste(levels(filtered_go_data$mice), collapse = " vs "), "\n")
cat("显著富集通路总数:", nrow(filtered_go_data), "\n")
cat("GO通路排列顺序（从顶部到底部）:\n")
for(i in seq_along(existing_go_terms)) {
  cat(i, ":", existing_go_terms[i], "\n")
}

cat("\n各品系通路数:\n")
for(mouse_strain in levels(filtered_go_data$mice)) {
  count <- sum(filtered_go_data$mice == mouse_strain)
  cat("  ", mouse_strain, ":", count, "个通路\n")
}

cat("\n各时间点通路数:\n")
for(cluster in levels(filtered_go_data$Cluster)) {
  count <- sum(filtered_go_data$Cluster == cluster)
  cat("  ", cluster, ":", count, "个通路\n")
}

cat("\n输出文件已保存到 D:/Desktop/\n")
cat("生成的文件:\n")
cat("1. GO_Enrichment_Mice_Comparison_Faceted.pdf/png - 分面比较图\n")
cat("2. GO_Enrichment_Mice_Comparison_Combined.pdf/png - 合并比较图\n")
cat("3. GO_Enrichment_Mice_Comparison_Horizontal.pdf/png - 水平布局图\n")
cat("4. GO_Enrichment_Mice_Side_by_Side.pdf/png - 并排比较图\n")
cat("5. GO_Enrichment_Mice_Grouped.pdf/png - 按功能分组图\n")
cat("6. GO_Enrichment_Mice_Comparison_Results.xlsx - 详细结果表格\n")
cat(paste(rep("=", 60), collapse = ""), "\n")