# 7.3 合并所有比较组的差异基因热图
# 提取所有比较组的显著差异基因
all_sig_genes <- unique(unlist(lapply(deg_list, function(x) {
  x %>% filter(change != "Not Significant") %>% pull(gene)
})))

if(length(all_sig_genes) > 10) {
  # 取前30个最显著的基因
  gene_min_padj <- sapply(all_sig_genes, function(gene) {
    min(sapply(deg_list, function(df) {
      idx <- which(df$gene == gene)
      if(length(idx) > 0) df$padj[idx] else 1
    }), na.rm = TRUE)
  })
  
  top_genes <- names(sort(gene_min_padj))[1:min(30, length(gene_min_padj))]
  
  # 检查基因是否存在于vsd对象中
  genes_in_vsd <- top_genes[top_genes %in% rownames(assay(vsd))]
  
  if(length(genes_in_vsd) > 5) {
    # 提取表达矩阵并标准化
    mat_original <- assay(vsd)[genes_in_vsd, ]
    
    # 按照指定顺序重新排列样本：MOCK, MPXV_Ib, MPXV_IIb
    # 获取每个分组的样本名
    mock_samples <- rownames(sample_info)[sample_info$group == "MOCK"]
    mpvx_Ib_samples <- rownames(sample_info)[sample_info$group == "MPXV_Ib"]
    mpvx_IIb_samples <- rownames(sample_info)[sample_info$group == "MPXV_IIb"]
    
    # 按指定顺序组合样本: MOCK, MPXV_Ib, MPXV_IIb
    sample_order <- c(mock_samples, mpvx_Ib_samples, mpvx_IIb_samples)
    
    # 重新排列矩阵列
    mat_ordered <- mat_original[, sample_order]
    
    # 生成新的样本标签
    sample_info$new_label <- paste0(sample_info$group, "-", 
                                    gsub(".*_(\\d)$", "\\1", rownames(sample_info)))
    
    # 对样本信息表进行排序
    sample_info_ordered <- sample_info[sample_order, ]
    ordered_cols_new_labels <- sample_info_ordered$new_label
    
    # 设置新列名
    colnames(mat_ordered) <- ordered_cols_new_labels
    
    # 标准化矩阵
    mat <- t(scale(t(mat_ordered)))  # Z-score标准化
    
    # 转置矩阵，使样本在行，基因在列
    mat_transposed <- t(mat)
    
    # 创建行注释（样本注释）
    annotation_row <- data.frame(
      Group = sample_info_ordered$group,
      row.names = sample_info_ordered$new_label
    )
    
    # 创建颜色映射
    annotation_colors <- list(
      Group = c("MOCK" = color_mock, "MPXV_Ib" = color_Ib, "MPXV_IIb" = color_IIb)
    )
    
    # 绘制横向热图
    p_heat <- pheatmap(mat_transposed,
                       color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                       annotation_row = annotation_row,  # 使用行注释
                       annotation_colors = annotation_colors,
                       cluster_rows = FALSE,  # 不聚类行（样本），保持指定顺序
                       cluster_cols = TRUE,   # 聚类列（基因）
                       show_colnames = TRUE,  # 显示基因名
                       show_rownames = TRUE,  # 显示样本名
                       border_color = NA,
                       fontsize_row = 9,       # 样本名字体大小
                       fontsize_col = 8,       # 基因名字体大小
                       fontsize = 10,
                       main = "Top 30 Most Significant DEGs\nMPXV Groups Comparison",
                       clustering_distance_cols = "euclidean",
                       clustering_method = "complete",
                       cutree_cols = 2,       # 将基因分为2个簇
                       angle_col = 45,         # 基因名倾斜45度
                       treeheight_col = 20,    # 基因聚类树高度
                       treeheight_row = 0)     # 样本不聚类，不显示树
    
    # 保存横向热图
    pdf("Heatmap_MPXV_Groups_DEGs_Horizontal.pdf", width = 14, height = 8)
    grid::grid.newpage()
    grid::grid.draw(p_heat$gtable)
    dev.off()
    
    cat("已生成按指定顺序排列的横向热图：MOCK -> MPXV_Ib -> MPXV_IIb\n")
    
    # 保存基因列表
    top_genes_info <- data.frame(
      gene = genes_in_vsd,
      min_padj = sapply(genes_in_vsd, function(g) {
        min(sapply(deg_list, function(df) {
          idx <- which(df$gene == g)
          if(length(idx) > 0) df$padj[idx] else NA
        }), na.rm = TRUE)
      })
    ) %>%
      arrange(min_padj)
    
    write.csv(top_genes_info, "Top30_Significant_Genes_MPXV_Groups.csv", row.names = FALSE)
    
  } else {
    cat("警告：只有", length(genes_in_vsd), "个基因在vsd对象中找到，不足以绘制热图\n")
  }
} else {
  cat("显著差异基因数量不足，无法绘制热图\n")
}