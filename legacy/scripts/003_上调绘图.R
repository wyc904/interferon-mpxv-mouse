# ================= 绘图 =================

actual_labels <- label_mapping[levels(plot_data$Cluster)]

# 如果存在未定义label，则保留原始名称
actual_labels[is.na(actual_labels)] <- 
  names(actual_labels[is.na(actual_labels)])

p_basic <- ggplot(
  plot_data,
  aes(
    x = Cluster,
    y = Description
  )
) +
  
  # 气泡图
  geom_point(
    aes(
      size = Count,
      color = LogP
    ),
    alpha = 0.85
  ) +
  
  # 颜色渐变
  scale_color_gradientn(
    
    colors = red_colors,
    
    name = "-Log10(P.adj)",
    
    guide = guide_colorbar(
      
      title.position = "top",
      
      title.hjust = 0.5,
      
      barwidth = unit(5, "cm"),
      
      barheight = unit(0.5, "cm"),
      
      direction = "horizontal",
      
      order = 1
    )
  ) +
  
  # 点大小
  scale_size_continuous(
    
    range = c(3, 8),
    
    name = "Gene Count",
    
    guide = guide_legend(
      
      title.position = "top",
      
      title.hjust = 0.5,
      
      direction = "horizontal",
      
      nrow = 1,
      
      label.position = "bottom",
      
      order = 2,
      
      override.aes = list(alpha = 1)
    )
  ) +
  
  # X轴标签
  scale_x_discrete(
    labels = actual_labels
  ) +
  
  labs(
    
    title = paste(
      "Top",
      top_n_per_cluster,
      "Enriched GO Terms per Group"
    ),
    
    x = NULL,
    
    y = NULL
  ) +
  
  # 使用更干净的主题
  theme_minimal() +
  
  theme(
    
    # ================= X轴 =================
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "black"
    ),
    
    # ================= Y轴 =================
    
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    # ================= 标题 =================
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    # ================= 网格线 =================
    
    # 竖直参考线
    panel.grid.major.x = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    
    panel.grid.minor.x = element_blank(),
    
    # 水平参考线
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    
    panel.grid.minor.y = element_blank(),
    
    # ================= 图例 =================
    
    legend.position = "bottom",
    
    legend.box = "horizontal",
    
    legend.box.just = "center",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.width = unit(1, "cm"),
    
    legend.margin = margin(
      t = 5,
      r = 0,
      b = 0,
      l = 0
    ),
    
    legend.spacing.x = unit(0.5, "cm"),
    
    legend.box.margin = margin(
      t = 10,
      r = 0,
      b = 5,
      l = 0
    ),
    
    # ================= 外边距 =================
    
    plot.margin = unit(
      c(1, 1, 1.5, 1),
      "cm"
    )
  )

# ================= 打印 =================

print(p_basic)



ggsave("D:/Desktop/GO_Top15_上调.pdf", p_basic, width =16, height = 10)




