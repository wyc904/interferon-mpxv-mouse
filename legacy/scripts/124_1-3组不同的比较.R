# 加载必要的包
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(tidyr)
library(grid)
library(ggVennDiagram)
library(VennDiagram)
library(gridExtra)

# 1. 全局配色设置 -------------------------------------------------
color_mock <- "#708090"    # 灰色
color_mpxv <- "#A52A2A"    # 深棕色
color_up <- "#C75D69"       # 浅红
color_down <- "#4682B4"     # 蓝色
color_ns <- "#D3D3D3"       # 浅灰

# 2. 整合数据并创建表达矩阵 -------------------------------------------------
cat("=== 开始整合数据 ===\n")

# 假设您有三个不同的count矩阵
# count_matrix: 包含C57_7d_1_P, C57_7d_2_P, C57_7d_3_P
# count_matrix2: 包含AGB6的样本
# 假设还有其他包含"JZ25055844.4.1.P.4.1.P_combined"等样本的矩阵

# 首先，从不同的矩阵中提取所需的样本
# 假设count_matrix已经加载
cat("从count_matrix中提取C57_7d样本...\n")
c57_7d_samples <- c("C57_7d_1_P", "C57_7d_2_P", "C57_7d_3_P")

# 检查这些样本是否存在于count_matrix中
available_c57_7d <- intersect(c57_7d_samples, colnames(count_matrix))
cat("找到的C57_7d样本:", paste(available_c57_7d, collapse = ", "), "\n")

# 提取C57_7d样本的表达数据
if(length(available_c57_7d) > 0) {
  c57_7d_expr <- count_matrix[, available_c57_7d, drop = FALSE]
} else {
  stop("在count_matrix中未找到C57_7d样本")
}

# 从count_matrix2中提取AGB6.MPXV_lb样本
cat("从count_matrix2中提取AGB6.MPXV_lb样本...\n")
agb6_mpxv_samples <- c(
  "JZ25070169.AGB6.MPXV_lb.1.P.AGB6.MPXV_lb.1.P_combined",
  "JZ25070170.AGB6.MPXV_lb.2.P.AGB6.MPXV_lb.2.P_combined", 
  "JZ25070171.AGB6.MPXV_lb.3.P.AGB6.MPXV_lb.3.P_combined"
)

available_agb6_mpxv <- intersect(agb6_mpxv_samples, colnames(count_matrix2))
cat("找到的AGB6.MPXV_lb样本:", paste(available_agb6_mpxv, collapse = ", "), "\n")

if(length(available_agb6_mpxv) > 0) {
  agb6_mpxv_expr <- count_matrix2[, available_agb6_mpxv, drop = FALSE]
} else {
  stop("在count_matrix2中未找到AGB6.MPXV_lb样本")
}

# 从count_matrix3中提取4系列样本
cat("从count_matrix3中提取4系列样本...\n")
four_series_samples <- c(
  "JZ25055844.4.1.P.4.1.P_combined",
  "JZ25055845.4.2.P.4.2.P_combined",
  "JZ25055846.4.3.P.4.3.P_combined"
)

available_four_series <- intersect(four_series_samples, colnames(count_matrix3))
cat("找到的4系列样本:", paste(available_four_series, collapse = ", "), "\n")

if(length(available_four_series) > 0) {
  four_series_expr <- count_matrix3[, available_four_series, drop = FALSE]
} else {
  stop("在count_matrix3中未找到4系列样本")
}

# 3. 合并所有样本的表达矩阵 -------------------------------------------------
cat("合并所有样本的表达矩阵...\n")

# 获取共同的基因（行名）
common_genes <- Reduce(intersect, list(
  rownames(c57_7d_expr),
  rownames(agb6_mpxv_expr),
  rownames(four_series_expr)
))

cat("共同基因数量:", length(common_genes), "\n")

# 合并矩阵
combined_expr <- cbind(
  c57_7d_expr[common_genes, ],
  agb6_mpxv_expr[common_genes, ],
  four_series_expr[common_genes, ]
)

# 简化列名以便处理
colnames(combined_expr) <- c(
  "C57_7d_1", "C57_7d_2", "C57_7d_3",
  "AGB6_MPXV_1", "AGB6_MPXV_2", "AGB6_MPXV_3",
  "Four_Series_1", "Four_Series_2", "Four_Series_3"
)

cat("合并后的表达矩阵维度:", dim(combined_expr), "\n")

# 4. 创建样本信息 -------------------------------------------------
cat("创建样本信息...\n")

sample_info <- data.frame(
  sample = colnames(combined_expr),
  group = factor(c(
    rep("C57_7d", 3),
    rep("AGB6_MPXV", 3),
    rep("Four_Series", 3)
  ), levels = c("C57_7d", "AGB6_MPXV", "Four_Series")),
  strain = factor(c(
    rep("C57", 3),
    rep("AGB6", 3),
    rep("Four_Series", 3)
  )),
  condition = factor(c(
    rep("MPXV", 3),
    rep("MPXV", 3),
    rep("Unknown", 3)
  )),
  stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$sample
print(sample_info)

# 5. 创建DESeq2对象 -------------------------------------------------
cat("创建DESeq2对象...\n")

# 确保表达矩阵是整数
combined_expr_int <- round(as.matrix(combined_expr))

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = combined_expr_int,
  colData = sample_info,
  design = ~ group
)

cat("原始基因数量:", nrow(dds), "\n")

# 6. 过滤低表达基因 -------------------------------------------------
cat("过滤低表达基因...\n")

# 过滤在至少3个样本中counts >= 10的基因
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("过滤后基因数量:", nrow(dds), "\n")

# 7. 差异表达分析 -------------------------------------------------
cat("进行差异表达分析...\n")

# 设置参考水平
dds$group <- relevel(dds$group, ref = "C57_7d")

# 运行DESeq2
dds <- DESeq(dds)

# 获取所有比较的结果
comparisons <- list(
  "AGB6_MPXV_vs_C57_7d" = c("group", "AGB6_MPXV", "C57_7d"),
  "Four_Series_vs_C57_7d" = c("group", "Four_Series", "C57_7d"),
  "AGB6_MPXV_vs_Four_Series" = c("group", "AGB6_MPXV", "Four_Series")
)

# 存储所有结果
results_list <- list()
deg_list <- list()

# 定义阈值
FC_cutoff <- 1.5
P_cutoff <- 0.05

# 处理每个比较
for(comp_name in names(comparisons)) {
  cat("处理比较:", comp_name, "\n")
  
  comp <- comparisons[[comp_name]]
  res <- results(dds, contrast = comp, alpha = P_cutoff)
  
  # 存储完整结果
  results_list[[comp_name]] <- res
  
  # 转换为数据框并添加基因名
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # 标记差异表达状态
  res_df$change <- "Not Significant"
  res_df$change[res_df$padj < P_cutoff & res_df$log2FoldChange > log2(FC_cutoff)] <- "Up-regulated"
  res_df$change[res_df$padj < P_cutoff & res_df$log2FoldChange < -log2(FC_cutoff)] <- "Down-regulated"
  res_df$change <- factor(res_df$change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
  
  # 存储DEGs
  deg_list[[comp_name]] <- res_df
  
  # 输出统计信息
  up_count <- sum(res_df$change == "Up-regulated", na.rm = TRUE)
  down_count <- sum(res_df$change == "Down-regulated", na.rm = TRUE)
  cat(comp_name, ": 上调 =", up_count, ", 下调 =", down_count, "\n")
}

# 8. 可视化分析 -------------------------------------------------
cat("开始可视化分析...\n")

# 创建输出目录
output_dir <- "Multi_Group_Comparison_Analysis"
dir.create(output_dir, showWarnings = FALSE)

# 8.1 PCA分析
cat("进行PCA分析...\n")
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA图
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(aes(label = name), size = 3.5, fontface = "bold", max.overlaps = 20) +
  scale_color_manual(values = c("C57_7d" = "#1B9E77", 
                                "AGB6_MPXV" = "#D95F02", 
                                "Four_Series" = "#7570B3")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA Analysis of Three Groups",
       subtitle = "C57_7d vs AGB6_MPXV vs Four_Series") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

print(p_pca)
ggsave(file.path(output_dir, "PCA_Three_Groups.pdf"), p_pca, width = 9, height = 7)

# 8.2 火山图
cat("绘制火山图...\n")
for(comp_name in names(comparisons)) {
  res_df <- deg_list[[comp_name]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(20) %>% 
      pull(gene)
    
    p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down-regulated" = color_down, 
                                    "Not Significant" = color_ns, 
                                    "Up-regulated" = color_up)) +
      geom_vline(xintercept = c(-log2(FC_cutoff), log2(FC_cutoff)), 
                 linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(P_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_text_repel(data = subset(res_df, gene %in% top_genes),
                      aes(label = gene, color = change),
                      size = 2.5, show.legend = FALSE, max.overlaps = 15) +
      labs(title = paste("Volcano Plot:", comp_name), 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      theme_bw() + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "bottom"
      )
    
    print(p_vol)
    ggsave(file.path(output_dir, paste0(comp_name, "_Volcano.pdf")), p_vol, width = 8, height = 7)
  }
}

# 8.3 热图分析
cat("绘制热图...\n")
# 提取所有比较的显著差异基因
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
  
  # 提取表达矩阵并标准化
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  
  # 创建样本注释
  annotation_col <- data.frame(
    Group = sample_info$group,
    Strain = sample_info$strain
  )
  rownames(annotation_col) <- colnames(mat)
  
  # 创建颜色映射
  annotation_colors <- list(
    Group = c("C57_7d" = "#1B9E77", "AGB6_MPXV" = "#D95F02", "Four_Series" = "#7570B3"),
    Strain = c("C57" = "#66C2A5", "AGB6" = "#FC8D62", "Four_Series" = "#8DA0CB")
  )
  
  # 绘制热图
  p_heat <- pheatmap(mat,
                     color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                     annotation_col = annotation_col,
                     annotation_colors = annotation_colors,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     fontsize_row = 8,
                     fontsize_col = 9,
                     main = "Top 30 Differentially Expressed Genes Across Groups")
  
  # 保存热图
  pdf(file.path(output_dir, "Heatmap_Top30_DEGs.pdf"), width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  dev.off()
}

# 8.4 DEGs数量统计
cat("统计DEGs数量...\n")
deg_counts <- data.frame(
  Comparison = names(comparisons),
  Up = sapply(names(comparisons), function(x) {
    sum(deg_list[[x]]$change == "Up-regulated", na.rm = TRUE)
  }),
  Down = sapply(names(comparisons), function(x) {
    sum(deg_list[[x]]$change == "Down-regulated", na.rm = TRUE)
  })
)

deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count")

p_bar <- ggplot(deg_counts_long, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Differentially Expressed Genes Count", 
       x = "Comparison", y = "Number of DEGs") +
  scale_fill_manual(values = c("Down" = color_down, "Up" = color_up)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar)
ggsave(file.path(output_dir, "DEGs_Count_By_Comparison.pdf"), p_bar, width = 8, height = 6)

# 8.5 Venn图分析
cat("绘制Venn图...\n")
# 提取每个比较的DEGs
deg_sets <- lapply(deg_list, function(df) {
  df %>% filter(change != "Not Significant") %>% pull(gene)
})

# 检查是否有足够的基因
if(all(sapply(deg_sets, length) > 0)) {
  # 使用ggVennDiagram
  p_venn <- ggVennDiagram(deg_sets, label_alpha = 0) +
    scale_fill_gradient(low = "#F4F4F4", high = "#2E86AB") +
    labs(title = "Venn Diagram of DEGs Across Comparisons") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  print(p_venn)
  ggsave(file.path(output_dir, "Venn_Diagram_DEGs.pdf"), p_venn, width = 8, height = 7)
}

# 8.6 干扰素通路基因分析
cat("进行干扰素通路基因分析...\n")

# 定义小鼠干扰素相关基因
ifn_genes <- c(
  "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2",
  "Stat1", "Stat2", "Irp9", "Irf7", "Irf3",
  "Ifi44", "Ifi44l", "Ifit1", "Ifit2", "Ifit3", "Ifitm1", "Ifitm3",
  "Isg15", "Isg20", "Mx1", "Mx2", "Oas1a", "Oas1b", "Oas2", "Oas3"
)

# 检查哪些基因在我们的数据中
ifn_genes_present <- ifn_genes[ifn_genes %in% rownames(vsd)]

if(length(ifn_genes_present) > 0) {
  cat("找到的干扰素相关基因:", length(ifn_genes_present), "\n")
  
  # 提取干扰素基因表达矩阵
  ifn_expr <- assay(vsd)[ifn_genes_present, ]
  ifn_expr_scaled <- t(scale(t(ifn_expr)))  # 行标准化
  
  # 创建注释数据框
  annotation_col <- data.frame(
    Group = sample_info$group,
    row.names = colnames(ifn_expr_scaled)
  )
  
  # 绘制热图
  p_ifn <- pheatmap(ifn_expr_scaled,
                    annotation_col = annotation_col,
                    show_rownames = TRUE,
                    cluster_cols = TRUE,
                    main = "Interferon Pathway Genes Expression",
                    fontsize_row = 8)
  
  # 保存干扰素热图
  pdf(file.path(output_dir, "Interferon_Pathway_Heatmap.pdf"), width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_ifn$gtable)
  dev.off()
} else {
  cat("未找到干扰素相关基因\n")
}

# 9. 保存结果 -------------------------------------------------
cat("保存分析结果...\n")

# 保存每个比较的完整结果和DEGs
for(comp_name in names(comparisons)) {
  if(comp_name %in% names(results_list)) {
    write.csv(as.data.frame(results_list[[comp_name]]), 
              file.path(output_dir, paste0(comp_name, "_full_results.csv")),
              row.names = TRUE)
    
    if(nrow(deg_list[[comp_name]]) > 0) {
      write.csv(deg_list[[comp_name]], 
                file.path(output_dir, paste0(comp_name, "_DEGs.csv")),
                row.names = FALSE)
    }
  }
}

# 保存样本信息
write.csv(sample_info, file.path(output_dir, "sample_information.csv"))

# 保存基因表达矩阵
write.csv(counts(dds, normalized = TRUE), 
          file.path(output_dir, "normalized_counts.csv"))

# 保存工作空间
save.image(file.path(output_dir, "Multi_Group_Comparison_Analysis.RData"))

# 10. 结果总结 -------------------------------------------------
cat("=== 多组比较分析完成 ===\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总样本数:", ncol(dds), "\n")
cat("分析基因数:", nrow(dds), "\n\n")

cat("各比较组DEGs统计:\n")
for(comp_name in names(comparisons)) {
  if(comp_name %in% names(deg_list)) {
    deg_df <- deg_list[[comp_name]]
    total_degs <- sum(deg_df$change != "Not Significant", na.rm = TRUE)
    up_degs <- sum(deg_df$change == "Up-regulated", na.rm = TRUE)
    down_degs <- sum(deg_df$change == "Down-regulated", na.rm = TRUE)
    cat(comp_name, ": 总DEGs =", total_degs, 
        "(上调 =", up_degs, ", 下调 =", down_degs, ")\n")
  }
}

cat("\n分析完成！结果保存在:", output_dir, "目录中\n")
cat("包含以下文件:\n")
cat("1. PCA图\n")
cat("2. 火山图\n")
cat("3. 热图（前30个DEGs）\n")
cat("4. DEGs数量统计图\n")
cat("5. Venn图\n")
cat("6. 干扰素通路基因热图\n")
cat("7. 完整的差异表达分析结果（CSV格式）\n")
cat("8. 样本信息文件\n")
cat("9. 标准化表达矩阵\n")
cat("10. 完整的R工作空间文件\n")