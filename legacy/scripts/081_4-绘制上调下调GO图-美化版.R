library(ggplot2)
library(dplyr)

# 定义一个辅助函数：将 GeneRatio (如 "5/100") 转换为数值 (0.05)
parse_ratio <- function(ratio) {
  sapply(strsplit(ratio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
}

# ==============================================================================
# 1. 上调基因 (UP): 红-蓝主色调
#    设定：蓝色 = 不显著 (LogP小) -> 红色 = 显著 (LogP大)
# ==============================================================================
if(nrow(ego_up) > 0) {
  
  # 1. 提取数据并预处理 - 按GeneRatio大小排序（从大到小改为从小到大）
  up_data <- as.data.frame(ego_up) %>%
    mutate(
      LogP = -log10(p.adjust),              # 计算 -Log10
      GeneRatioVal = parse_ratio(GeneRatio)  # 转换比例
    ) %>%
    arrange(p.adjust) %>%        # 先按显著性排序
    slice(1:10) %>%              # 取前10个
    arrange(GeneRatioVal) %>%     # 改为按GeneRatio从小到大排序
    mutate(
      # 将Description转为因子，按GeneRatio大小排序（从小到大）
      Description = factor(Description, levels = Description) 
    )
  
  # 2. 定义红蓝渐变色 (冷 -> 热)
  red_blue_colors <- c("#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C",
                       "#D7301F", "#67000D")
  
  # 3. 绘图
  p_up_custom <- ggplot(up_data, aes(x = GeneRatioVal, y = Description)) +
    geom_point(aes(size = Count, color = LogP)) + # 核心：点的大小和颜色
    scale_color_gradientn(
      colors = red_blue_colors,
      name = "-Log10(Adj. P)", 
      guide = guide_colorbar(
        title.position = "top", 
        title.hjust = 0.5, 
        barwidth = unit(5, "cm"),
        direction = "horizontal",  # 水平方向
        order = 1  # 颜色图例排在第一位
      )
    ) +
    scale_size_continuous(
      range = c(3, 8), # 调整点的大小范围
      name = "Gene Count",
      breaks = function(x) {
        # 创建4个均匀分布的断点
        if (length(unique(x)) > 1) {
          seq(min(x), max(x), length.out = 4) %>% round()
        } else {
          x
        }
      },
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",  # 水平方向
        nrow = 1,                 # 单行排列
        ncol = 4,                 # 4列（即4个点）
        byrow = TRUE,             # 按行填充
        order = 2                 # 大小图例排在第二位
      )
    ) +
    labs(
      title = "Top Enriched Terms: Up-Regulated Genes",
      x = "Gene Ratio", 
      y = NULL
    ) +
    theme_minimal() +
    theme(
      # 增加图形周围留白
      plot.margin = margin(20, 20, 20, 20),  # 上、右、下、左的边距
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
      # 图例设置 - 专业杂志风格
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center",  # 居中
      legend.box.just = "center",      # 图例框居中
      legend.margin = margin(t = 10, r = 0, b = 5, l = 0),
      legend.spacing.x = unit(0.3, "cm"),  # 图例项之间的水平间距
      legend.spacing.y = unit(0.2, "cm"),  # 图例项之间的垂直间距
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 9),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.major.y = element_line(color = "grey90")
    )
  
  print(p_up_custom)
  ggsave("GO_UP_RedBlue_SortedByRatio_Ascending.pdf", p_up_custom, width = 10, height = 7)
}


# ==============================================================================
# 2. 下调基因 (DOWN): 黄-绿主色调
#    设定：深绿 = 不显著 (LogP小) -> 亮黄 = 显著 (LogP大)
# ==============================================================================
if(nrow(ego_down) > 0) {
  
  # 1. 提取数据并预处理 - 按GeneRatio大小排序（从大到小改为从小到大）
  down_data <- as.data.frame(ego_down) %>%
    mutate(
      LogP = -log10(p.adjust),
      GeneRatioVal = parse_ratio(GeneRatio)
    ) %>%
    arrange(p.adjust) %>%        # 先按显著性排序
    slice(1:10) %>%              # 取前10个
    arrange(GeneRatioVal) %>%     # 改为按GeneRatio从小到大排序
    mutate(
      Description = factor(Description, levels = Description)
    )
  
  # 2. 定义黄绿渐变色 (深绿 -> 亮黄)
  yellow_green_colors <- c("#A1D99B", "#74C476", "#41AB5D",
                           "#238B45", "#006D2C", "#00441B", "#003314")
  
  # 3. 绘图
  p_down_custom <- ggplot(down_data, aes(x = GeneRatioVal, y = Description)) +
    geom_point(aes(size = Count, color = LogP)) +
    scale_color_gradientn(
      colors = yellow_green_colors,
      name = "-Log10(Adj. P)", 
      guide = guide_colorbar(
        title.position = "top", 
        title.hjust = 0.5, 
        barwidth = unit(5, "cm"),
        direction = "horizontal",  # 水平方向
        order = 1  # 颜色图例排在第一位
      )
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name = "Gene Count",
      breaks = function(x) {
        # 创建4个均匀分布的断点
        if (length(unique(x)) > 1) {
          seq(min(x), max(x), length.out = 4) %>% round()
        } else {
          x
        }
      },
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",  # 水平方向
        nrow = 1,                 # 单行排列
        ncol = 4,                 # 4列（即4个点）
        byrow = TRUE,             # 按行填充
        order = 2                 # 大小图例排在第二位
      )
    ) +
    labs(
      title = "Top Enriched Terms: Down-Regulated Genes",
      x = "Gene Ratio", 
      y = NULL
    ) +
    theme_minimal() +
    theme(
      # 增加图形周围留白
      plot.margin = margin(20, 20, 20, 20),  # 上、右、下、左的边距
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
      # 图例设置 - 专业杂志风格
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center",  # 居中
      legend.box.just = "center",      # 图例框居中
      legend.margin = margin(t = 10, r = 0, b = 5, l = 0),
      legend.spacing.x = unit(0.3, "cm"),  # 图例项之间的水平间距
      legend.spacing.y = unit(0.2, "cm"),  # 图例项之间的垂直间距
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 9),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.major.y = element_line(color = "grey90")
    )
  
  print(p_down_custom)
  ggsave("GO_DOWN_YellowGreen_SortedByRatio_Ascending.pdf", p_down_custom, width = 10, height = 7)
}