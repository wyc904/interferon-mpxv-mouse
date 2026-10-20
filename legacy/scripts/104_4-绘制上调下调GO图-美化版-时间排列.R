# --------------------------------------------------------------------------
# 模块 E: 跨时间点 GO 富集结果比较可视化
# --------------------------------------------------------------------------

if (run_enrichment) {
  
  # 检查是否收集到结果
  if (length(all_ego_up_list) == 0 && length(all_ego_down_list) == 0) {
    cat("\n警告: 全局 GO 结果列表 (all_ego_up_list / all_ego_down_list) 为空，无法进行跨时间点比较分析。\n")
  } else {
    
    cat("\n======================================================\n")
    cat(">>> 正在生成跨时间点 GO 富集结果比较图 (Top 20 BP 条目)...\n")
    cat("======================================================\n")
    
    # 定义时间点顺序 (按照实验时间顺序)
    timepoint_order <- c("mock", "3d", "5d", "7d", "10d")
    
    # 定义颜色方案
    red_blue_colors <- c("#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C",
                         "#D7301F", "#67000D")
    
    yellow_green_colors <- c("#A1D99B", "#74C476", "#41AB5D",
                             "#238B45", "#006D2C", "#00441B", "#003314")
    
    # --- 1. 上调基因 (Upregulated Genes) GO-BP 比较图 ---
    if (length(all_ego_up_list) > 0) {
      
      # 整合所有 enrichResult 对象为一个 data.frame 
      up_df_list <- lapply(names(all_ego_up_list), function(tp) {
        df <- as.data.frame(all_ego_up_list[[tp]])
        df$Cluster <- tp # 添加时间点作为分组列
        return(df)
      })
      
      up_df_combined <- do.call(rbind, up_df_list)
      
      if (nrow(up_df_combined) > 0) {
        
        # 筛选跨时间点 p.adjust 最小的 Top 20 Description
        top_descriptions_up <- up_df_combined %>%
          group_by(Description) %>%
          summarise(min_p.adjust = min(p.adjust, na.rm = TRUE)) %>%
          arrange(min_p.adjust) %>%
          slice(1:20) %>% # 取前20个
          pull(Description)
        
        # 筛选出要绘图的数据，并按照时间点顺序排序
        plot_data_up <- up_df_combined %>%
          filter(Description %in% top_descriptions_up) %>%
          mutate(
            LogP = -log10(p.adjust),  # 计算 -Log10
            GeneRatioVal = sapply(strsplit(GeneRatio, "/"), 
                                  function(x) as.numeric(x[1])/as.numeric(x[2])),
            # 将 Cluster 转换为因子，按照时间点顺序排序
            Cluster = factor(Cluster, levels = timepoint_order)
          ) %>%
          # 为 Description 排序，基于最小 p.adjust
          mutate(Description = factor(Description, levels = rev(top_descriptions_up))) %>%
          # 按照时间点顺序排序数据
          arrange(Cluster)
        
        # 绘图 - 使用与示例相同的样式
        p_up_compare <- ggplot(plot_data_up, 
                               aes(x = Cluster, y = Description, size = Count, color = LogP)) +
          geom_point() +
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
            title = "GO-BP Enrichment Comparison (Upregulated Genes - Top 20)",
            x = "Timepoint", 
            y = "GO Biological Process"
          ) +
          theme_minimal() +
          theme(
            # 增加图形周围留白
            plot.margin = margin(20, 20, 20, 20),  # 上、右、下、左的边距
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
            axis.text.y = element_text(color = "black", size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
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
        
        # 保存图表
        ggsave(file.path(output_dir, "GO_Upregulated_BP_AllTimepoints_Comparison.pdf"), 
               p_up_compare, width = 12, height = 10)
        cat("上调基因 GO-BP 比较图已保存。\n")
        
      } else {
        cat("警告: 上调基因 GO-BP 比较数据为空。\n")
      }
    }
    
    # --- 2. 下调基因 (Downregulated Genes) GO-BP 比较图 ---
    if (length(all_ego_down_list) > 0) {
      
      down_df_list <- lapply(names(all_ego_down_list), function(tp) {
        df <- as.data.frame(all_ego_down_list[[tp]])
        df$Cluster <- tp # 添加时间点作为分组列
        return(df)
      })
      
      down_df_combined <- do.call(rbind, down_df_list)
      
      if (nrow(down_df_combined) > 0) {
        
        # 筛选跨时间点 p.adjust 最小的 Top 20 Description
        top_descriptions_down <- down_df_combined %>%
          group_by(Description) %>%
          summarise(min_p.adjust = min(p.adjust, na.rm = TRUE)) %>%
          arrange(min_p.adjust) %>%
          slice(1:20) %>% # 取前20个
          pull(Description)
        
        # 筛选出要绘图的数据，并按照时间点顺序排序
        plot_data_down <- down_df_combined %>%
          filter(Description %in% top_descriptions_down) %>%
          mutate(
            LogP = -log10(p.adjust),  # 计算 -Log10
            GeneRatioVal = sapply(strsplit(GeneRatio, "/"), 
                                  function(x) as.numeric(x[1])/as.numeric(x[2])),
            # 将 Cluster 转换为因子，按照时间点顺序排序
            Cluster = factor(Cluster, levels = timepoint_order)
          ) %>%
          # 为 Description 排序，基于最小 p.adjust
          mutate(Description = factor(Description, levels = rev(top_descriptions_down))) %>%
          # 按照时间点顺序排序数据http://127.0.0.1:44321/graphics/plot_zoom_png?width=731&height=721
          arrange(Cluster)
        
        # 绘图 - 使用与示例相同的样式
        p_down_compare <- ggplot(plot_data_down, 
                                 aes(x = Cluster, y = Description, size = Count, color = LogP)) +
          geom_point() +
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
            title = "GO-BP Enrichment Comparison (Downregulated Genes - Top 20)",
            x = "Timepoint", 
            y = "GO Biological Process"
          ) +
          theme_minimal() +
          theme(
            # 增加图形周围留白
            plot.margin = margin(20, 20, 20, 20),  # 上、右、下、左的边距
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
            axis.text.y = element_text(color = "black", size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
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
        
        # 保存图表
        ggsave(file.path(output_dir, "GO_Downregulated_BP_AllTimepoints_Comparison.pdf"), 
               p_down_compare, width = 12, height = 10)
        cat("下调基因 GO-BP 比较图已保存。\n")
        
      } else {
        cat("警告: 下调基因 GO-BP 比较数据为空。\n")
      }
    }
  }
}