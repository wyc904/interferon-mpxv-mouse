# 载入必要的库
library(openxlsx)  # 用于读取Excel文件
library(ggplot2)   # 用于绘图
library(dplyr)     # 数据处理
library(tidyr)     # 数据整理
library(stringr)   # 字符串处理
library(viridis)   # 颜色方案
library(forcats)   # 因子处理 (用于排序)

# ================= 1. 数据读取与清洗 =================

# 读取GO分析数据

go_data <-read.csv("D:/Desktop/compareCluster_Down_GO_BP.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 确保数值列正确
numeric_columns <- c("RichFactor", "p.adjust", "Count", "FoldEnrichment")
for(col in numeric_columns) {
  if(col %in% colnames(go_data)) {
    go_data[[col]] <- as.numeric(as.character(go_data[[col]]))
  }
}

# 排除包含NA或非数值的行，并初筛显著性
filtered_go_data <- go_data %>%
  filter(if_all(all_of(numeric_columns), ~!is.na(.))) %>%
  filter(p.adjust < 0.05) # 初步过滤不显著的

# 计算-log10(p.adjust) 用于绘图颜色
filtered_go_data$LogP <- -log10(filtered_go_data$p.adjust)

# ================= 2. 比较组（X轴）处理 =================

# 期望的比较组顺序
expected_comparisons <- c("DMSO_vs_mock","Tecovirimat_vs_mock","Tecovirimat_vs_DMSO")


# 检查实际存在的比较组
available_comparisons <- unique(filtered_go_data$Cluster)
available_comparisons <- available_comparisons[!is.na(available_comparisons)]

# 确定最终用于绘图的比较组顺序（取交集）
actual_comparisons <- intersect(expected_comparisons, available_comparisons)

# 如果交集为空（说明列名不匹配），则使用数据中原有的所有组
if(length(actual_comparisons) == 0) {
  actual_comparisons <- available_comparisons
  cat("警告：未找到预期的比较组，将使用数据中存在的所有组。\n")
}

# 设置Cluster因子的水平（控制X轴顺序）
filtered_go_data$Cluster <- factor(filtered_go_data$Cluster, levels = actual_comparisons)
# 过滤掉不在选定比较组中的数据
filtered_go_data <- filtered_go_data %>% filter(!is.na(Cluster))


# ================= 3. 筛选 Top 15 GO Terms (核心修改) =================

# 设置要筛选的数量
top_n_per_cluster <- 15

# 逻辑：
# 1. 按Cluster分组
# 2. 按p.adjust升序排列（越小越显著）
# 3. 取每组的前15个
# 4. 提取Description名称并去重
top_terms_list <- filtered_go_data %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>% 
  slice_head(n = top_n_per_cluster) %>%
  ungroup() %>%
  pull(Description) %>%
  unique()

# 从原始过滤数据中，只保留这部分Top Terms
plot_data <- filtered_go_data %>%
  filter(Description %in% top_terms_list)

cat("筛选完成：\n")
cat("共筛选出唯一的GO通路数量:", length(top_terms_list), "\n")

# ================= 4. 设置绘图顺序 (Y轴) =================

# 为了让图好看，我们需要对Y轴进行排序。
# 方法：计算每个Term在所有组中的平均LogP或最小p值，或者使用聚类顺序。
# 这里使用简单的处理：按总体的显著性排序，或者按它们在第一个Cluster中的出现顺序。
# 也可以使用 forcats::fct_reorder 根据 Cluster 和 LogP 进行排序

# 这里我们将 Description 转换为因子，并根据 Cluster 和 LogP 重新排序
# 这样可以让在同一个Cluster显著的通路聚在一起
# 修改前（报错）：
# mutate(Description = fct_reorder(Description, Cluster, .fun = min, .desc = TRUE))

# 修改后（正确）：加入 as.numeric()
plot_data <- plot_data %>%
  mutate(Description = fct_reorder(Description, as.numeric(Cluster), .fun = min, .desc = TRUE))
# 如果想单纯按显著性排序（最显著的在最上面），可以用下面这行代替上面那行：
# plot_data$Description <- fct_reorder(plot_data$Description, plot_data$LogP, .fun = max)


# ================= 5. 绘图 =================

# 定义红色渐变色
# 定义红色渐变色 (浅红 -> 深红)
red_colors <- c("#EAF3FA", "#BCD7EA", "#86B4D4",  "#315F8A", "#1D3F5F")


# 标签映射
label_mapping <- c(
  "DMSO_vs_mock" = "DMSO_vs_mock",
  "Tecovirimat_vs_mock" = "Tecovirimat_vs_mock", 
  "Tecovirimat_vs_DMSO" = "Tecovirimat_vs_DMSO"
)



actual_labels <- label_mapping[levels(plot_data$Cluster)]
# 如果映射中有NA（也就是有未定义的组），用原始名称填补
actual_labels[is.na(actual_labels)] <- names(actual_labels[is.na(actual_labels)])
p_basic <- ggplot(plot_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = LogP), alpha = 0.8) +
  scale_color_gradientn(colors = red_colors, name = "-Log10(P.adj)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  scale_x_discrete(labels = actual_labels) +
  labs(
    title = paste("Top", top_n_per_cluster, "Enriched GO Terms per Group"),
    x = NULL, 
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",  # 将图例放在底部
    legend.direction = "horizontal",  # 图例水平排列
    legend.box = "horizontal",  # 图例框水平排列
    legend.box.just = "center",  # 图例框居中对齐
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.8, "cm"),  # 调整图例键高度
    legend.key.width = unit(1.5, "cm"),  # 调整颜色图例宽度
    legend.spacing.x = unit(1, "cm"),  # 图例之间的水平间距
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),  # 减小图例边距
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major.x = element_blank(),  # 去掉X轴网格线
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt")  # 增加底部边距，为图例留空间
  ) +
  # 单独调整每个图例的指南
  guides(
    color = guide_colorbar(
      title.position = "top",  # 标题放在顶部
      title.hjust = 0.5,  # 标题水平居中
      barwidth = unit(3, "cm"),  # 颜色条宽度
      barheight = unit(0.4, "cm"),  # 颜色条高度
      direction = "horizontal"  # 颜色条水平
    ),
    size = guide_legend(
      title.position = "top",  # 标题放在顶部
      title.hjust = 0.5,  # 标题水平居中
      direction = "horizontal",  # 图例水平排列
      nrow = 1,  # 一行显示
      label.position = "bottom",  # 标签放在底部
      override.aes = list(alpha = 1)  # 覆盖透明度
    )
  )

# 打印图形
print(p_basic)
# 打印图形
print(p_basic)


# ggsave("D:/Desktop/GO_Top15_Enrichment.pdf", p_basic, width = 10, height = dynamic_height)
# ggsave("D:/Desktop/GO_Top15_Enrichment.png", p_basic, width = 10, height = dynamic_height, dpi = 300)