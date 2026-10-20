# 载入必要的库
library(openxlsx)  # 用于读取Excel文件
library(clusterProfiler)  # 用于GO分析
library(ggplot2)  # 用于绘图
library(dplyr)  # 数据处理

# 读取GO分析数据
file_path <- "D:/Desktop/GO_analysis_3d_vs_mock.xlsx"
go_data <- read.xlsx(file_path, sheet = 1)

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

# 确保 RichFactor 为数值类型
filtered_go_data$RichFactor <- as.numeric(filtered_go_data$RichFactor)

# 排除包含 NA 或非数值的行
filtered_go_data <- filtered_go_data %>%
  filter(!is.na(RichFactor))

# 设置临界值，选择显著性较强的GO通路
filtered_go_data <- filtered_go_data %>%
  filter(p.adjust < 0.05)

# 绘制点图，展示所有20条GO通路
ggplot(filtered_go_data, aes(x = reorder(Description, -RichFactor), y = RichFactor, color = ONTOLOGY)) +
  geom_point(size = 4) +  # 点图
  coord_flip() +  # 横向排列
  theme_minimal() +
  labs(
    title = "GO Enrichment Analysis for MPXV Infected Mice",
    x = "GO Terms",
    y = "Rich Factor"
  ) +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转X轴标签
