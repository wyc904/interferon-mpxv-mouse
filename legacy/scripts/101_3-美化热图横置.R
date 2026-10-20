# 8. 差异表达基因热图分析 -------------------------------------------------
# 选择最有差异的前30个基因（跨所有时间点）

# 8.1 合并所有时间点的差异分析结果
all_degs <- data.frame()
for(tp in timepoints) {
  if(tp %in% names(deg_list)) {
    deg_df <- deg_list[[tp]] %>%
      mutate(timepoint = tp) %>%
      filter(change != "Not Significant")  # 只保留显著差异基因
    all_degs <- bind_rows(all_degs, deg_df)
  }
}

# 8.2 选择前30个最显著的基因（按最小padj排序）
if(nrow(all_degs) > 0) {
  # 对每个基因，取其最小padj（跨时间点）
  top_genes <- all_degs %>%
    group_by(gene) %>%
    summarize(min_padj = min(padj, na.rm = TRUE),
              max_abs_log2fc = max(abs(log2FoldChange), na.rm = TRUE)) %>%
    arrange(min_padj) %>%
    head(30) %>%
    pull(gene)
  
  cat("选择前30个最显著差异基因进行热图分析\n")
  cat("基因数量:", length(top_genes), "\n")
  
  if(length(top_genes) > 2) {
    
    # --------------------------------------------------------------------------
    # --- 核心修改部分: 样本重命名和重排序 ---
    # --------------------------------------------------------------------------
    
    # 1. 定义期望的排序顺序
    ordered_timepoints <- c("mock", "3d", "5d", "7d", "10d")
    
    # 2. 生成新的样本标签 (例如: "mock-1", "mock-2", "3d-1", ...)
    sample_info$new_label <- ave(
      seq_along(sample_info$timepoint), 
      sample_info$timepoint, 
      FUN = function(x) paste0(sample_info$timepoint[x], "-", seq_along(x))
    )
    
    # 3. 对样本信息表进行排序
    sample_info_ordered <- sample_info %>%
      # 将 timepoint 转换为因子，确保排序顺序正确
      mutate(timepoint_fct = factor(timepoint, levels = ordered_timepoints)) %>%
      # 按时间点因子和新标签排序，确保时间组内顺序稳定
      arrange(timepoint_fct, new_label)
    
    # 提取新的列名和原始的样本名（用于匹配矩阵）
    ordered_cols_new_labels <- sample_info_ordered$new_label
    ordered_cols_original_names <- rownames(sample_info_ordered) 
    
    # 8.3 提取表达矩阵并标准化
    vsd <- vst(dds, blind = FALSE)
    mat_original <- assay(vsd)[top_genes, ]
    
    # 4. 根据排序后的样本名对矩阵进行重排序和重命名
    mat_ordered <- mat_original[, ordered_cols_original_names] 
    colnames(mat_ordered) <- ordered_cols_new_labels # 设置新列名
    
    mat <- t(scale(t(mat_ordered)))  # Z-score标准化
    
    # --------------------------------------------------------------------------
    # --- 修改: 转置矩阵，使热图横向放置 ---
    # --------------------------------------------------------------------------
    mat_transposed <- t(mat)  # 转置矩阵，使样本在行，基因在列
    
    # 8.4 创建样本注释
    # 注释表也需要使用新的样本名和排序
    annotation_row <- data.frame(
      Timepoint = sample_info_ordered$timepoint,
      Condition = sample_info_ordered$group_label,
      row.names = sample_info_ordered$new_label # 确保行名是新标签
    )
    
    # --------------------------------------------------------------------------
    # --- 核心修改部分结束 ---
    # --------------------------------------------------------------------------
    
    # 8.5 创建颜色映射
    # 确保颜色映射的名称与 Timepoint/Condition 列中的值相匹配
    timepoint_colors <- c("mock" = "#708090", "3d" = "#1B9E77", 
                          "5d" = "#D95F02", "7d" = "#7570B3", "10d" = "#E7298A")
    condition_colors <- c("MOCK" = "#708090", "MPXV_lb" = "#A52A2A")
    
    annotation_colors <- list(
      Timepoint = timepoint_colors,
      Condition = condition_colors
    )
    
    # 8.6 绘制横向热图
    p_heat <- pheatmap(mat_transposed,
                       color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                       annotation_row = annotation_row,  # 使用行注释而非列注释
                       annotation_colors = annotation_colors,
                       cluster_rows = FALSE,  # 不聚类行（样本），保持时间顺序
                       cluster_cols = TRUE,   # 聚类列（基因）
                       show_colnames = TRUE,  # 显示基因名
                       show_rownames = TRUE,  # 显示样本名
                       border_color = NA,
                       fontsize_row = 9,       # 样本名字体大小
                       fontsize_col = 8,       # 基因名字体大小
                       fontsize = 10,
                       main = "Top 30 Most Significant DEGs\nAcross All Timepoints",
                       clustering_distance_cols = "euclidean",
                       clustering_method = "complete",
                       cutree_cols = 2,       # 将基因分为2个簇
                       angle_col = 45,         # 基因名倾斜45度
                       treeheight_col = 20,     # 基因聚类树高度
                       treeheight_row = 0)      # 样本不聚类，不显示树
    
    # 8.7 保存横向热图
    # 调整PDF尺寸以适应横向布局
    pdf("Top30_DEGs_Heatmap_Horizontal.pdf", width = 14, height = 8)  # 宽度大于高度
    grid::grid.newpage()
    grid::grid.draw(p_heat$gtable)
    dev.off()
    
    print(paste("成功绘制横向热图"))
    
    # 8.8 保存基因列表
    top_genes_info <- all_degs %>%
      filter(gene %in% top_genes) %>%
      arrange(min_padj = padj) %>%
      select(gene, log2FoldChange, padj, change, timepoint) %>%
      distinct(gene, .keep_all = TRUE)
    
    write.csv(top_genes_info, "Top30_Significant_Genes_List.csv", row.names = FALSE)
    
  } else {
    cat("显著差异基因数量不足，无法绘制热图\n")
  }
} else {
  cat("未发现显著差异基因\n")
}