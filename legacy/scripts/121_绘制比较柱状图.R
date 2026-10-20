# 加载必要的包
library(ggplot2)
library(reshape2)

# 创建数据框
data <- data.frame(
  Group = rep(c("C57", "AGB6", "STAT1"), each = 2),
  Condition = rep(c("UP", "DOWN"), times = 3),
  Mock = c(462, 589, 472, 584, 242, 71),
  Vpox = c(NA, NA, NA, NA, NA, NA)  # 占位符，实际数据需补充
)

# 转换数据为长格式
melted_data <- melt(data, id.vars = c("Group", "Condition"),
                    variable.name = "Treatment",
                    value.name = "Count")

# 删除NA值
melted_data <- melted_data[!is.na(melted_data$Count), ]

# 设置JAMA风格主题
jama_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "top",  # 将图例放在顶部
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
  
    legend.box.margin = margin(2, 2, 2, 2),
    text = element_text(family = "sans", size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# 绘制柱状图
ggplot(melted_data, aes(x = Group, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name = "Regulation Direction",  # 图例标题
    values = c("UP" = "#4F94CD", "DOWN" = "#FF6B6B"),  # JAMA常用蓝红配色
    labels = c("UP" = "Up-regulated", "DOWN" = "Down-regulated")  # 修改图例标签
  ) +
  labs(
    title = "",
    x = "Mouse Strain",
    y = "Number of Differentially Expressed Genes"
  ) +
  jama_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    # 图例布局调整
    legend.direction = "horizontal",  # 水平排列
    legend.key.size = unit(0.5, "cm"),  # 图例键大小
    legend.margin = margin(0, 0, 0, 0)  # 图例外边距
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, max(melted_data$Count) * 1.1)  # 为顶部图例留出空间
  ) +
  geom_text(
    aes(label = Count), 
    position = position_dodge(width = 0.8), 
    vjust = -0.5, 
    size = 3.5,
    fontface = "bold"  # 数值标签加粗
  ) 

# 如果想要将图例放在右侧，可以使用以下设置：
# 将上面的 legend.position = "top" 改为 "right"
# 将 legend.direction = "horizontal" 改为 "vertical"