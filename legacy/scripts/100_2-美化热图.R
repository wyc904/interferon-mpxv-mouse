# 8. 差异表达基因热图分析 -------------------------------------------------
# 确保已加载必要的库
library(dplyr)
library(pheatmap)
library(RColorBrewer)
# 假设 vsd, dds, deg_list, timepoints, sample_info 已在环境中定义

# 8.1 合并所有时间点的差异分析结果
all_degs <- data.frame()
for(tp in timepoints) {
  if(tp %in% names(deg_list)) {
    deg_df <- deg_list[[tp]] %>%
      mutate(timepoint = tp) %>%
      filter(change != "Not Significant") # 只保留显著差异基因
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
    
    # 8.4 创建样本注释
    # 注释表也需要使用新的样本名和排序
    annotation_col <- data.frame(
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
    
    # 8.6 绘制热图
    p_heat <- pheatmap(mat,
                       color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE, # 关闭列聚类，以保持时间顺序
                       show_colnames = TRUE,
                       show_rownames = TRUE,
                       border_color = NA,
                       fontsize_row = 8,
                       fontsize_col = 9,
                       fontsize = 10,
                    
                       clustering_distance_rows = "euclidean",
                       clustering_method = "complete",
                       cutree_rows = 2, # 将基因分为2个簇
                       cutree_cols = 5, # 可尝试根据时间点数量分簇
                       angle_col = 45,
                       treeheight_row = 20,
                       treeheight_col = 0) # 移除列树高，因为我们不聚类
    
    # 8.7 保存热图
    # 注意：pheatmap对象不能直接用pdf()保存，需要使用print或grid.draw
    pdf("Top30_DEGs_Heatmap_TimeOrdered.pdf", width = 12, height = 10)
    print(p_heat)
    dev.off()
    
    print(paste("成功绘制前30个最显著差异基因热图"))
    
    # 8.8 保存基因列表
    # (这部分代码保持不变，因为它只处理基因列表)
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