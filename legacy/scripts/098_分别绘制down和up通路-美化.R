# --------------------------------------------------------------------------
# Module E: Visualization of GO enrichment results for merging all comparison groups (optional)
# --------------------------------------------------------------------------

# 确保必要的包已加载
# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
# if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
# library(dplyr)
# library(stringr)

# 定义自定义颜色（已在您的原代码块末尾提供）
red_blue_colors <- c("#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C",
                     "#D7301F", "#67000D")
yellow_green_colors <- c("#A1D99B", "#74C476", "#41AB5D",
                         "#238B45", "#006D2C", "#00441B", "#003314")

# 定义一个辅助函数来解析 GeneRatio 字符串
parse_ratio <- function(ratio_string) {
  sapply(str_split(ratio_string, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}


# --- 1. Combined Upregulated GO BP Dotplot (ALL Groups, Top 20, Custom Styles) ---
if(length(all_ego_up_list) > 1) {
  cat("\n>>> Generate a comparison dotplot for Upregulated GO BP results (All Groups, Top 20, Custom Styles)...\n")
  
  # 1. 合并 GO 结果
  ego_up_combined <- merge_result(all_ego_up_list)
  
  if(nrow(ego_up_combined) > 0) {
    
    # 2. 提取数据并预处理
    up_data <- as.data.frame(ego_up_combined) %>%
      mutate(
        LogP = -log10(p.adjust),
        GeneRatioVal = parse_ratio(GeneRatio)
      ) %>%
      # 根据 p.adjust 排序，取前 20 个独特的 GO Description
      arrange(p.adjust)
    
    # 获取前 20 个 Description 的因子顺序
    top_descriptions <- unique(up_data$Description)[1:20]
    
    # 筛选数据并设置 Description 的因子水平
    up_data_top20 <- up_data %>%
      filter(Description %in% top_descriptions) %>%
      mutate(
        Description = factor(Description, levels = rev(top_descriptions)) # Y轴排序：将最显著的放在最上面
      )
    
    # 3. 绘图 (基于 ggplot2)
    p_compare_up_custom <- ggplot(up_data_top20, 
                                  aes(x = Cluster, y = Description)) +
      
      # 使用 GeneRatioVal 映射点的大小 (size)，LogP 映射颜色 (color)
      geom_point(aes(size = GeneRatioVal, color = LogP)) +
      
      # 颜色标尺（使用自定义的红蓝色）
      scale_color_gradientn(
        colors = red_blue_colors,
        name = "-Log10(Adj. P)",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(5, "cm"),
          direction = "horizontal",
          order = 1
        )
      ) +
      
      # 大小标尺（映射 GeneRatioVal）
      scale_size_continuous(
        range = c(3, 8), # 调整点的大小范围
        name = "Gene Ratio",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          direction = "horizontal",
          nrow = 1,
          byrow = TRUE,
          order = 2
        )
      ) +
      
      # 标签和主题
      labs(
        title = "Comparison of GO BP Enrichment - Upregulated Genes (All Groups)",
        x = NULL, # X轴为 Cluster，不需要标题
        y = NULL
      ) +
      theme_minimal() +
      theme(
        plot.margin = margin(20, 20, 20, 20),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
        axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11), # X轴标签倾斜
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.margin = margin(t = 10, r = 0, b = 5, l = 0),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 9),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "grey90")
      )
    
    # 4. 保存图
    ggsave(file.path(output_dir, "GO_BP_Upregulated_Comparison_ALL_Groups_CustomStyle.pdf"),
           p_compare_up_custom, width = 14, height = 10)
    
    cat("Combined dotplot for upregulated genes (All Groups) saved as GO_BP_Upregulated_Comparison_ALL_Groups_CustomStyle.pdf\n")
  } else {
    cat("Warning: Merged upregulated GO BP results (All Groups) are empty. Skipping combined plot.\n")
  }
} else {
  cat("\nWarning: Not enough upregulated GO BP results (All Groups) to create the combined plot.\n")
}


# --- 2. Combined Downregulated GO BP Dotplot (ALL Groups, Top 20, Custom Styles) ---
# 同样对 Downregulated 组进行相同的定制（使用绿/黄色调）

if(length(all_ego_down_list) > 1) {
  cat("\n>>> Generate a comparison dotplot for Downregulated GO BP results (All Groups, Top 20, Custom Styles)...\n")
  
  # 1. 合并 GO 结果
  ego_down_combined <- merge_result(all_ego_down_list)
  
  if(nrow(ego_down_combined) > 0) {
    
    # 2. 提取数据并预处理
    down_data <- as.data.frame(ego_down_combined) %>%
      mutate(
        LogP = -log10(p.adjust),
        GeneRatioVal = parse_ratio(GeneRatio)
      ) %>%
      # 根据 p.adjust 排序
      arrange(p.adjust)
    
    # 获取前 20 个 Description 的因子顺序
    top_descriptions_down <- unique(down_data$Description)[1:20]
    
    # 筛选数据并设置 Description 的因子水平
    down_data_top20 <- down_data %>%
      filter(Description %in% top_descriptions_down) %>%
      mutate(
        Description = factor(Description, levels = rev(top_descriptions_down))
      )
    
    # 3. 绘图 (基于 ggplot2)
    p_compare_down_custom <- ggplot(down_data_top20, 
                                    aes(x = Cluster, y = Description)) +
      
      geom_point(aes(size = GeneRatioVal, color = LogP)) +
      
      # 颜色标尺（使用自定义的绿黄色）
      scale_color_gradientn(
        colors = yellow_green_colors,
        name = "-Log10(Adj. P)",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(5, "cm"),
          direction = "horizontal",
          order = 1
        )
      ) +
      
      # 大小标尺（映射 GeneRatioVal）
      scale_size_continuous(
        range = c(3, 8),
        name = "Gene Ratio",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          direction = "horizontal",
          nrow = 1,
          byrow = TRUE,
          order = 2
        )
      ) +
      
      # 标签和主题
      labs(
        title = "Comparison of GO BP Enrichment - Downregulated Genes (All Groups)",
        x = NULL, 
        y = NULL
      ) +
      theme_minimal() +
      theme(
        plot.margin = margin(20, 20, 20, 20),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
        axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.margin = margin(t = 10, r = 0, b = 5, l = 0),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 9),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "grey90")
      )
    
    # 4. 保存图
    ggsave(file.path(output_dir, "GO_BP_Downregulated_Comparison_ALL_Groups_CustomStyle.pdf"),
           p_compare_down_custom, width = 14, height = 10)
    
    cat("Combined dotplot for downregulated genes (All Groups) saved as GO_BP_Downregulated_Comparison_ALL_Groups_CustomStyle.pdf\n")
  } else {
    cat("Warning: Merged downregulated GO BP results (All Groups) are empty. Skipping combined plot.\n")
  }
} else {
  cat("\nWarning: Not enough downregulated GO BP results (All Groups) to create the combined plot.\n")
}