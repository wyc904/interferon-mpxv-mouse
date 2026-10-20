# 优化后的热图代码
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# 按照时间顺序排列样本
time_order <- c("mock", "3d", "5d", "7d", "10d")
sample_order <- unlist(lapply(time_order, function(tp) {
  rownames(sample_info)[sample_info$timepoint == tp]
}))

# 按照时间顺序重新排列矩阵
sampleDistMatrix_ordered <- sampleDistMatrix[sample_order, sample_order]

# 创建注释信息
annotation_df <- data.frame(
  Timepoint = sample_info[sample_order, "timepoint"],
  Condition = sample_info[sample_order, "condition"],
  row.names = sample_order
)

# 定义与图中一致的配色方案
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

# 绘制优化后的热图
pheatmap::pheatmap(sampleDistMatrix_ordered,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = "样本间距离热图",
                   annotation_col = annotation_df,
                   annotation_row = annotation_df,
                   annotation_colors = ann_colors,  # 应用自定义颜色
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   fontsize_row = 8,
                   fontsize_col = 8,
                   color = colorRampPalette(c("white", "#008000", "red"))(100),
                   border_color = "grey60",  # 添加边框使单元格更清晰
                   display_numbers = FALSE,  # 不显示数字，避免过于拥挤
                   gaps_row = c(3, 6, 9, 12),  # 在mock/3d/5d/7d/10d之间添加间隔
                   gaps_col = c(3, 6, 9, 12),  # 同上
                   cutree_rows = 5,  # 按5个时间点切割聚类树
                   cutree_cols = 5)  # 同上