# 8. 干扰素通路基因分析（按时间顺序排列）-----------------------------------------
# 定义小鼠干扰素相关基因
ifn_genes <- c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Mx2", 
               "Oas1a", "Oas1b", "Oas2", "Oas3",
               "Stat1", "Stat2", "Irf7", "Irf9", "Ifnb1", "Ifna1", "Ifna2",
               "Irf1", "Irf3", "Irf5", "Jak1", "Jak2", "Tyk2",
               "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2")

# 检查哪些基因在我们的数据中
ifn_genes_present <- ifn_genes[ifn_genes %in% rownames(vsd)]

if(length(ifn_genes_present) > 0) {
  # 提取干扰素基因表达矩阵
  ifn_expr <- assay(vsd)[ifn_genes_present, ]
  ifn_expr_scaled <- t(scale(t(ifn_expr)))  # 行标准化
  
  # 按照时间顺序重新排列样本
  ifn_expr_ordered <- ifn_expr_scaled[, sample_order]
  
  # 创建注释数据框（按时间顺序）
  annotation_col <- data.frame(
    Timepoint = sample_info[sample_order, "timepoint"],
    Condition = sample_info[sample_order, "condition"],
    row.names = sample_order
  )
  
  # 定义注释颜色
  ann_colors <- list(
    Timepoint = c(
      "mock" = "#808080",  # 灰色
      "3d" = "#800080",    # 紫色
      "5d" = "#FFC0CB",    # 粉色
      "7d" = "#008000",    # 绿色
      "10d" = "#0000FF"    # 蓝色
    ),
    Condition = c(
      "mock" = "#FFFF00",  # 黄色
      "MPXV1b" = "#00FFFF" # 青色
    )
  )
  
  # 绘制热图（按时间顺序排列，不进行聚类）
  pheatmap::pheatmap(ifn_expr_ordered,
                     cluster_rows = TRUE,  # 对基因进行聚类
                     cluster_cols = FALSE, # 不对样本进行聚类，保持时间顺序
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     main = "干扰素通路基因表达热图（按时间顺序）",
                     fontsize_row = 8,
                     fontsize_col = 8,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     gaps_col = c(3, 6, 9, 12))  # 在时间点之间添加间隔
  
  # 可选：绘制按时间点分面的箱线图，显示关键干扰素基因表达趋势
  ifn_expr_long <- as.data.frame(ifn_expr_ordered) %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
    left_join(sample_info, by = "sample")
  
  # 选择几个关键的干扰素基因进行可视化
  key_ifn_genes <- c("Isg15", "Mx1", "Stat1", "Irf7", "Ifnb1")
  key_ifn_genes <- key_ifn_genes[key_ifn_genes %in% ifn_genes_present]
  
  if(length(key_ifn_genes) > 0) {
    ifn_expr_long_key <- ifn_expr_long %>% filter(gene %in% key_ifn_genes)
    
    ggplot(ifn_expr_long_key, aes(x = timepoint, y = expression, fill = condition)) +
      geom_boxplot() +
      facet_wrap(~ gene, scales = "free_y") +
      labs(title = "关键干扰素基因表达趋势",
           x = "时间点", y = "标准化表达量") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  print(paste("找到", length(ifn_genes_present), "个干扰素相关基因"))
} else {
  print("未找到干扰素相关基因，请检查基因名")
}

# 9. 干扰素刺激基因(ISG)的专门分析 -----------------------------------------
# 提取已知的干扰素刺激基因列表（更全面的列表）
isg_genes <- c("Ifit1", "Ifit2", "Ifit3", "Ifit3b", "Isg15", "Isg20", 
               "Mx1", "Mx2", "Oas1a", "Oas1b", "Oas2", "Oas3", "Oasl1", "Oasl2",
               "Rsad2", "Iigp1", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Gbp6", "Gbp7",
               "Irf1", "Irf7", "Irf9", "Stat1", "Stat2", "Usp18", "Samhd1", "Trim25")

isg_genes_present <- isg_genes[isg_genes %in% rownames(vsd)]

if(length(isg_genes_present) > 0) {
  # 提取ISG表达矩阵并按时间顺序排列
  isg_expr <- assay(vsd)[isg_genes_present, ]
  isg_expr_scaled <- t(scale(t(isg_expr)))
  isg_expr_ordered <- isg_expr_scaled[, sample_order]
  
  # 绘制ISG热图
  pheatmap::pheatmap(isg_expr_ordered,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     main = "干扰素刺激基因(ISG)表达热图",
                     fontsize_row = 7,
                     gaps_col = c(3, 6, 9, 12),
                     color = colorRampPalette(c("blue", "white", "red"))(100))
  
  # 计算每个样本的ISG得分（平均表达）
  isg_scores <- colMeans(isg_expr_scaled, na.rm = TRUE)
  isg_score_df <- data.frame(
    sample = names(isg_scores),
    ISG_score = isg_scores
  ) %>%
    left_join(sample_info, by = "sample")
  
  # 绘制ISG得分箱线图
  ggplot(isg_score_df, aes(x = timepoint, y = ISG_score, fill = condition)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
    labs(title = "干扰素刺激基因(ISG)得分",
         x = "时间点", y = "ISG得分（平均标准化表达）") +
    theme_minimal()
  
  print(paste("找到", length(isg_genes_present), "个ISG基因"))
}