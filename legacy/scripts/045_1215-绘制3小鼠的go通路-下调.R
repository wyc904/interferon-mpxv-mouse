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
file_path <- "D:/Desktop/3种小鼠go-下调.xlsx"
go_data <- read.xlsx(file_path, sheet = 1)

# 检查数据结构
cat("数据维度:", dim(go_data), "\n")
cat("列名:", colnames(go_data), "\n")
cat("\n前6行数据:\n")
print(head(go_data))

# 检查MICE列的内容
cat("\nMICE列的唯一值:\n")
if("MICE" %in% colnames(go_data)) {
  mice_values <- unique(go_data$MICE)
  mice_values <- mice_values[!is.na(mice_values)]
  print(mice_values)
} else {
  cat("错误: 数据中不存在MICE列!\n")
  stop("请检查数据文件，确保包含MICE列")
}

# 检查统计值
cat("\n统计值检查:\n")
cat("p.adjust的范围:", range(go_data$p.adjust, na.rm = TRUE), "\n")
cat("Count的范围:", range(go_data$Count, na.rm = TRUE), "\n")
cat("FoldEnrichment的范围:", range(go_data$FoldEnrichment, na.rm = TRUE), "\n")
cat("EnrichmentRatio的范围:", range(go_data$EnrichmentRatio, na.rm = TRUE), "\n")

# 选择有意义的GO通路进行可视化
# 选择的15条下调GO通路
target_go_terms <- c(
  # 1. 细胞周期相关通路（C57特异性下调）
  "mitotic cell cycle",                    # 有丝分裂细胞周期
  "spindle organization",                  # 纺锤体组织
  "mitotic nuclear division",              # 有丝分裂核分裂
  
  # 2. MAPK信号通路（C57和AGB6共同下调）
  "positive regulation of MAPK cascade",   # MAPK级联正调控
  "regulation of MAPK cascade",            # MAPK级联调控
  "MAPK cascade",                          # MAPK级联
  
  # 3. 抗原呈递通路（AGB6和STAT1共同下调）
  "positive regulation of antigen processing and presentation",  # 抗原加工呈递正调控
  "antigen processing and presentation",                         # 抗原加工呈递
  
  # 4. 免疫激活通路
  "positive regulation of immune response",  # 免疫反应正调控
  "regulation of immune response",           # 免疫反应调控
  
  # 5. STAT1特异性下调通路
  "immune response-activating signaling pathway",   # 免疫反应激活信号通路
  "innate immune response-activating signaling pathway",  # 天然免疫激活信号通路
  "cellular response to interferon-beta",   # 对干扰素-β的细胞反应
  
  # 6. 细胞迁移和吞噬相关
  "phagocytosis",                          # 吞噬作用
  "leukocyte migration",                   # 白细胞迁移
  "macrophage migration"                   # 巨噬细胞迁移
  ,"regulation of apoptotic signaling pathway",  # 凋亡信号通路调控
  "intrinsic apoptotic signaling pathway",     # 内源性凋亡信号通路
  "extrinsic apoptotic signaling pathway"  ,"cellular response to starvation","cellular response to nutrient levels", # 细胞对营养水平的反应
  "hexose metabolic process"    
)

# 从数据中筛选出目标GO通路
filtered_go_data <- go_data %>%
  filter(Description %in% target_go_terms)

# 检查筛选结果
cat("\n筛选后的数据维度:", dim(filtered_go_data), "\n")
cat("每个品系的通路数量:\n")
print(table(filtered_go_data$MICE))

# 计算-log10(p.adjust)，处理p.adjust为0的情况
# 添加一个小值避免取log(0)的无限大
min_padj <- min(filtered_go_data$p.adjust[filtered_go_data$p.adjust > 0], na.rm = TRUE)
filtered_go_data <- filtered_go_data %>%
  mutate(
    p.adjust_adj = ifelse(p.adjust == 0, min_padj/10, p.adjust),
    LogP = -log10(p.adjust_adj)
  )

# 设置MICE的顺序
filtered_go_data$MICE <- factor(filtered_go_data$MICE, 
                                levels = c("C57", "AGB6", "STAT1"))

# 为通路添加功能类别标签
filtered_go_data <- filtered_go_data %>%
  mutate(
    Pathway_Category = case_when(
      str_detect(Description, "chromosome|cell cycle|segregation|heterochromatin") ~ 
        "Genome Stability/Cell Cycle",
      str_detect(Description, "interferon|IFN") ~ 
        "Interferon Response",
      str_detect(Description, "innate|toll-like|receptor|NK") ~ 
        "Innate Immunity",
      str_detect(Description, "T cell|B cell|immune system|differentiation") ~ 
        "Adaptive Immunity",
      TRUE ~ "Other"
    )
  )

# 设置类别顺序
category_order <- c("Genome Stability/Cell Cycle", "Interferon Response", 
                    "Innate Immunity", "Adaptive Immunity", "Other")
filtered_go_data$Pathway_Category <- factor(filtered_go_data$Pathway_Category, 
                                            levels = category_order)

# 按指定顺序设置GO通路顺序
# 首先确保所有在数据中出现的GO通路都在指定的列表中
existing_go_terms <- intersect(target_go_terms, unique(filtered_go_data$Description))

# 按类别和富集程度排序
term_order <- filtered_go_data %>%
  group_by(Description, Pathway_Category) %>%
  summarise(
    Mean_LogP = mean(LogP, na.rm = TRUE),
    Max_FoldEnrich = max(FoldEnrichment, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(Pathway_Category, desc(Mean_LogP), desc(Max_FoldEnrich)) %>%
  pull(Description)

# 设置因子水平
filtered_go_data$Description <- factor(filtered_go_data$Description, 
                                       levels = rev(term_order))

# 检查数据
cat("\n最终数据维度:", dim(filtered_go_data), "\n")
cat("MICE顺序:", levels(filtered_go_data$MICE), "\n")
cat("显著通路数量:", nrow(filtered_go_data), "\n")
cat("\n通路类别分布:\n")
print(table(filtered_go_data$Pathway_Category))

# 定义红色渐变色 (浅红 -> 深红)
red_colors <- c("#A1D99B", "#74C476", "#41AB5D",
                "#238B45", "#006D2C", "#00441B", "#003314")


# 创建基础点图 - 与示例图保持一致
p_basic <- ggplot(filtered_go_data, aes(x = MICE, y = Description)) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  scale_color_gradientn(
    colors = red_colors, 
    name = "-Log10(P.adj)",
    guide = guide_colorbar(
      barwidth = 10,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  scale_size_continuous(
    range = c(3, 12), 
    breaks = c(5, 10, 15, 20),
    name = "Gene Count",
    guide = guide_legend(
      keywidth = unit(0.5, "cm"),
      keyheight = unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    title = "GO Enrichment Analysis of downregulated Pathways in MPXV-IB Infected Mice",
    subtitle = "Day 7 Post-Infection (Weight Nadir)",
    x = "Mouse Strain", 
    y = "GO Biological Process Terms"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0, 
      hjust = 0.5, 
      size = 12, 
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 12, 
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 12, 
      face = "bold",
      margin = margin(r = 10)
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.spacing.x = unit(0.5, "cm"),
    plot.title = element_text(
      hjust = 0.5, 
      face = "bold", 
      size = 16, 
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, 
      size = 12, 
      margin = margin(b = 15)
    ),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # 添加水平参考线
  geom_hline(
    yintercept = seq(0.5, length(unique(filtered_go_data$Description)) + 0.5, 1), 
    color = "grey90", 
    size = 0.2
  )

print(p_basic)


# print(p_with_annotation)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Comparison_WithAnnotation.pdf", p_with_annotation, 
#        width = 10, height = 8.5, dpi = 300)
# ggsave("D:/Desktop/GO_Enrichment_Mice_Comparison_WithAnnotation.png", p_with_annotation, 
#        width = 10, height = 8.5, dpi = 300)

