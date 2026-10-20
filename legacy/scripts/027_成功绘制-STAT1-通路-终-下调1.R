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
file_path <- "D:/Desktop/GO分析下调.xlsx"
go_data <- read.xlsx(file_path, sheet = 1)

# 检查数据结构
cat("数据维度:", dim(go_data), "\n")
cat("列名:", colnames(go_data), "\n")

# 检查Cluster列和Comparison列的内容
cat("Cluster列的唯一值:\n")
if("Cluster" %in% colnames(go_data)) {
  cluster_values <- unique(go_data$Cluster)
  cluster_values <- cluster_values[!is.na(cluster_values)]
  print(cluster_values)
} else {
  cat("错误: 数据中不存在Cluster列!\n")
  stop("请检查数据文件，确保包含Cluster列")
}

# 指定GO通路的顺序
target_go_terms <- c(
  "immune response-activating signaling pathway",
  "immune response-regulating signaling pathway", 
  "innate immune response-activating signaling pathway",
  "activation of innate immune response", 
  "regulation of innate immune response", 
  "cellular response to interferon-beta", 
  "response to interferon-beta",
  "pattern recognition receptor signaling pathway", 
  "cytoplasmic pattern recognition receptor signaling pathway",
  "cytokine-mediated signaling pathway",
  "positive regulation of macrophage migration", 
  "regulation of macrophage migration", 
  "antigen processing and presentation", "antigen processing and presentation of peptide antigen via MHC class II",
  "lymphocyte differentiation", 
  "lymphocyte proliferation",
  "regulation of lymphocyte proliferation", 
  "B cell activation", 
  "T cell differentiation", 
  "response to peptide", 
  "cellular response to peptide", 
  "regulation of GTPase activity",
  "peptidyl-serine phosphorylation", 
  "intracellular calcium ion homeostasis",
  "hydrogen peroxide-mediated programmed cell death")

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

# 重新定义Cluster列为比较组


filtered_go_data$Cluster


# 假设Cluster列中已经有"MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock", "MPXV_Ilb_vs_lb"等比较组
expected_comparisons <- c("MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock")

# 检查实际存在的比较组
available_comparisons <- as.character(unique(filtered_go_data$Cluster))
available_comparisons <- available_comparisons[!is.na(available_comparisons)]
cat("可用的比较组:\n", paste(available_comparisons, collapse = ", "), "\n")

# 按期望的顺序设置比较组水平
actual_comparisons <- intersect(expected_comparisons, available_comparisons)
if(length(actual_comparisons) > 0) {
  filtered_go_data$Cluster <- factor(filtered_go_data$Cluster, levels = actual_comparisons)
} else {
  filtered_go_data$Cluster <- factor(filtered_go_data$Cluster)
}

# 按指定顺序设置GO通路顺序
# 首先确保所有在数据中出现的GO通路都在指定的列表中
existing_go_terms <- intersect(target_go_terms, unique(filtered_go_data$Description))
# 按指定的顺序设置水平
filtered_go_data$Description <- factor(filtered_go_data$Description, 
                                       levels = rev(existing_go_terms))  # 使用rev反转顺序，使第一个在顶部

# 检查数据
cat("最终数据维度:", dim(filtered_go_data), "\n")
cat("可用的比较组:", levels(filtered_go_data$Cluster), "\n")
cat("显著通路数量:", nrow(filtered_go_data), "\n")
cat("GO通路顺序（从上到下）:\n")
cat(paste(rev(levels(filtered_go_data$Description)), collapse = "\n"), "\n")

# 定义红色渐变色 (浅红 -> 深红)
red_colors <- c("#A1D99B", "#74C476", "#41AB5D",
                "#238B45", "#006D2C", "#00441B", "#003314")




# 如果需要重命名x轴标签
p_basic <- ggplot(filtered_go_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  scale_color_gradientn(colors = red_colors, name = "-Log10(P.adj)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  scale_x_discrete(
    breaks = actual_comparisons,
    labels = c(
      "MPXV_lb_vs_Mock" = "MPXV_lb vs Mock",
      "MPXV_Ilb_vs_Mock" = "MPXV_Ilb vs Mock"
    ),
    drop = FALSE
  ) +
  labs(
    title = "GO Enrichment Analysis for Different Comparisons",
    x = "Comparison Groups", 
    y = "GO Terms"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(p_basic)
# ggsave("D:/Desktop/GO_Enrichment_Comparisons_Basic.pdf", p_basic, width = 10, height = 8)
# ggsave("D:/Desktop/GO_Enrichment_Comparisons_Basic.png", p_basic, width = 10, height = 8, dpi = 300)

