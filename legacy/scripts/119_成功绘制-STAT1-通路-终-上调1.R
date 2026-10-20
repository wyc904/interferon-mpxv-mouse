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
file_path <- "D:/Desktop/GO分析上调.xlsx"
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
  "chromosome segregation",  # 染色体分离
  "nuclear chromosome segregation",  # 核染色体分离
  "double-strand break repair",  # 双链断裂修复
  "DNA recombination",  # DNA重组
  "T cell receptor signaling pathway",  # T细胞受体信号通路
  "immune response-activating signaling pathway",  # 免疫反应激活信号通路
  "immune response-regulating signaling pathway",  # 免疫反应调节信号通路
  "response to radiation",  # 对辐射反应
  "double-strand break repair via homologous recombination",  # 同源重组修复
  "positive regulation of type II interferon production",  # II型干扰素产生正调控
  "regulation of actin cytoskeleton organization",  # 肌动蛋白细胞骨架调节
  "lymphocyte differentiation",  # 淋巴细胞分化
  "lymphocyte proliferation",  # 淋巴细胞增殖
  "response to UV",  # UV反应
  "regulation of DNA repair",  # DNA修复调节
  "protein kinase C signaling",  # PKC信号通路
  "response to endoplasmic reticulum stress",  # ER应激反应
  "nucleocytoplasmic transport",  # 核质运输
  "mitotic nuclear division",  # 有丝分裂核分裂
  "miRNA processing"  # miRNA加工
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

# 重新定义Cluster列为比较组


filtered_go_data$Cluster


# 假设Cluster列中已经有"MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock", "MPXV_Ilb_vs_lb"等比较组
expected_comparisons <- c("MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock", "MPXV_Ilb_vs_lb")

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
red_colors <- c("#FFE5E5", "#FF9999", "#FF6666", "#FF3333", "#CC0000", "#990000", "#660000")




# 如果需要重命名x轴标签
p_basic <- ggplot(filtered_go_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  scale_color_gradientn(colors = red_colors, name = "-Log10(P.adj)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  scale_x_discrete(
    breaks = actual_comparisons,
    labels = c(
      "MPXV_lb_vs_Mock" = "MPXV_lb vs Mock",
      "MPXV_Ilb_vs_Mock" = "MPXV_Ilb vs Mock",
      "MPXV_lb_vs_MPXV_Ilb" = "MPXV_lb vs MPXV_Ilb"
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

