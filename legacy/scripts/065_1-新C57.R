# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop/Project/10-干扰素小鼠/1-源/C57-counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:20]

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("C57_10d_1_P", "C57_10d_2_P", "C57_10d_3_P",
                            "C57_3d_1_P", "C57_3d_2_P", "C57_3d_3_P",
                            "C57_5d_1_P", "C57_5d_2_P", "C57_5d_3_P",
                            "C57_7d_1_P", "C57_7d_2_P", "C57_7d_3_P",
                            "c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P")

print("重命名后的列名:")
print(colnames(count_matrix))

# 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  timepoint = factor(c(rep("10d", 3), rep("3d", 3), rep("5d", 3), 
                       rep("7d", 3), rep("mock", 3)),
                     levels = c("mock", "3d", "5d", "7d", "10d")),
  condition = factor(c(rep("MPXV_lb", 12), rep("mock", 3)),
                     levels = c("mock", "MPXV_lb")),
  group_label = c(rep("MPXV_lb", 12), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

print("样本信息:")
print(sample_info)

# 3. 创建DESeq2对象 ----------------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ condition
)

print("DESeq2对象创建成功！")
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行DESeq2分析 ----------------------------------------------------
dds <- DESeq(dds)

# 6. 按时间点进行差异分析 ----------------------------------------------
timepoints <- c("3d", "5d", "7d", "10d")
results_list <- list()
deg_list <- list()

# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05

for(tp in timepoints) {
  print(paste("分析时间点:", tp))
  
  # 选择当前时间点和mock的样本
  tp_samples <- sample_info$sample[sample_info$timepoint %in% c(tp, "mock")]
  dds_sub <- dds[, tp_samples]
  
  # 确保condition因子水平正确（mock作为参考）
  dds_sub$condition <- relevel(dds_sub$condition, ref = "mock")
  
  # 运行DESeq2
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "MPXV_lb", "mock"), alpha = 0.05)
  
  # 保存完整结果
  results_list[[tp]] <- res
  
  # 提取显著差异基因并标记上下调
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      change = case_when(
        padj < P_cutoff & log2FoldChange > FC_cutoff ~ "Up-regulated",
        padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      change = factor(change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
    ) %>%
    arrange(padj)
  
  deg_list[[tp]] <- res_df
  print(paste(tp, "vs mock: 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# 7. 可视化 --------------------------------------------------------

# 全局配色设置
color_mock <- "#708090"    # 灰色
color_mpxv <- "#A52A2A"    # 深棕色
color_up <- "#C75D69"       # 浅红
color_down <- "#4682B4"     # 蓝色
color_ns <- "#D3D3D3"       # 浅灰

# 7.1 PCA图（所有样本）
vsd <- vst(dds, blind = FALSE)
# 版本2优化：改进10d时间点的显示效果
pcaData$timepoint_label <- paste0(pcaData$timepoint, "_", 
                                  gsub(".*_(\\d)_P.*", "\\1", pcaData$name))

# 改进的形状和颜色设置
timepoint_shapes <- c("mock" = 16,   # 实心圆
                      "3d" = 16,     # 实心三角
                      "5d" = 16,     # 实心方块
                      "7d" = 16,     # 实心菱形
                      "10d" = 16)     # 星号

timepoint_colors <- c("mock" = "#708090",    # 石板灰
                      "3d" = "#1B9E77",     # 蓝绿色
                      "5d" = "#D95F02",      # 橙色
                      "7d" = "#7570B3",      # 紫色
                      "10d" = "#E7298A")     # 粉红色

# 为10d时间点使用更大的点尺寸和不同的透明度
point_sizes <- c("mock" = 4, "3d" = 4, "5d" = 4, "7d" = 4, "10d" = 5)
point_alphas <- c("mock" = 0.8, "3d" = 0.8, "5d" = 0.8, "7d" = 0.8, "10d" = 0.9)

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  # 使用映射的大小和透明度
  geom_point(aes(size = timepoint, alpha = timepoint)) +
  scale_size_manual(values = point_sizes) +
  scale_alpha_manual(values = point_alphas) +
  scale_shape_manual(values = timepoint_shapes) +
  scale_color_manual(values = timepoint_colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "MPXV lb vs MOCK") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 14, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  geom_text_repel(aes(label = timepoint_label), 
                  size = 3.5, show.legend = FALSE, max.overlaps = 20,
                  fontface = "bold", box.padding = 0.3)

print(p_pca)
# ggsave("PCA_TimeSeries_All_Samples.pdf", p_pca, width = 8, height = 6)













# 7.2 各时间点火山图
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Down-regulated" = color_down, 
                                    "Not Significant" = color_ns, 
                                    "Up-regulated" = color_up)) +
      geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(P_cutoff), 
                 linetype = "dashed", color = "grey50") +
      geom_text_repel(data = subset(res_df, gene %in% top_genes),
                      aes(label = gene, color = change),
                      size = 2.5, show.legend = FALSE, max.overlaps = 15) +
      labs(title = paste("Volcano:", tp, "vs Mock"), 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      theme_bw() + 
      theme(legend.position = "bottom")
    
    print(p_vol)
    # 修改ggsave调用：添加dpi=600和指定路径
    ggsave(paste0("D:/Desktop/1/", tp, "_vs_Mock.pdf"), p_vol, width = 4, height = 4.5, dpi = 800)
  }
}






# 7.3 合并所有时间点的热图
# 提取所有时间点的显著差异基因
all_sig_genes <- unique(unlist(lapply(deg_list, function(x) {
  x %>% filter(change != "Not Significant") %>% pull(gene)
})))

if(length(all_sig_genes) > 10) {
  # 取前30个最显著的基因（跨所有时间点）
  # 计算每个基因在不同时间点出现的最小p值
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
    Timepoint = sample_info$timepoint,
    Condition = sample_info$group_label
  )
  rownames(annotation_col) <- colnames(mat)
  
  # 创建颜色映射
  annotation_colors <- list(
    Timepoint = c("mock" = color_mock, "3d" = "#1B9E77", "5d" = "#D95F02", 
                  "7d" = "#7570B3", "10d" = "#E7298A"),
    Condition = c("MOCK" = color_mock, "MPXV_lb" = color_mpxv)
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
                     main = "MPXV_lb 时间序列差异表达基因热图")
  
  # 保存热图
  pdf("Heatmap_TimeSeries_All_DEGs.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  dev.off()
}

# 7.4 各时间点DEGs数量统计
deg_counts <- data.frame(
  timepoint = timepoints,
  up = sapply(timepoints, function(x) {
    if(x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$change == "Up-regulated")
    } else {0}
  }),
  down = sapply(timepoints, function(x) {
    if(x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$change == "Down-regulated")
    } else {0}
  })
)

deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count")

p_bar <- ggplot(deg_counts_long, aes(x = timepoint, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "各时间点差异基因数量统计", 
       x = "时间点", y = "DEGs数量") +
  scale_fill_manual(values = c("down" = color_down, "up" = color_up)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_bar)
ggsave("DEGs_Count_By_Timepoint.pdf", p_bar, width = 8, height = 6)

# 8. 干扰素通路基因分析 -------------------------------------------------
# 定义小鼠干扰素相关基因
选择最有差异的前30基因

# 检查哪些基因在我们的数据中
ifn_genes_present <- ifn_genes[ifn_genes %in% rownames(vsd)]

if(length(ifn_genes_present) > 0) {
  # 提取干扰素基因表达矩阵
  ifn_expr <- assay(vsd)[ifn_genes_present, ]
  ifn_expr_scaled <- t(scale(t(ifn_expr)))  # 行标准化
  
  # 创建注释数据框
  annotation_col <- data.frame(
    Timepoint = sample_info$timepoint,
    Condition = sample_info$group_label,
    row.names = sample_info$sample
  )
  
  # 绘制热图
  p_ifn <- pheatmap(ifn_expr_scaled,
                    annotation_col = annotation_col,
                    show_rownames = TRUE,
                    cluster_cols = TRUE,
                    main = "干扰素通路基因表达热图",
                    fontsize_row = 8)
  
  # 保存干扰素热图
  pdf("Interferon_Pathway_Heatmap.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_ifn$gtable)
  dev.off()
  
  print(paste("找到", length(ifn_genes_present), "个干扰素相关基因"))
} else {
  print("未找到干扰素相关基因，请检查基因名")
}

# 9. 保存结果 ---------------------------------------------------------
output_dir <- "DESeq2_TimeSeries_Analysis_Results"
dir.create(output_dir, showWarnings = FALSE)

# 保存每个时间点的完整结果和DEGs
for(tp in timepoints) {
  if(tp %in% names(results_list)) {
    write.csv(as.data.frame(results_list[[tp]]), 
              file.path(output_dir, paste0("C57_", tp, "_vs_mock_full.csv")),
              row.names = TRUE)
    
    if(nrow(deg_list[[tp]]) > 0) {
      write.csv(deg_list[[tp]], 
                file.path(output_dir, paste0("C57_", tp, "_vs_mock_DEGs.csv")),
                row.names = FALSE)
    }
  }
}

# 保存样本信息
write.csv(sample_info, file.path(output_dir, "sample_information.csv"))

# 10. 结果总结 ---------------------------------------------------------
cat("=== MPXV_lb 时间序列分析完成 ===\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总样本数:", ncol(dds), "\n")
cat("分析基因数:", nrow(dds), "\n\n")

cat("各时间点DEGs统计:\n")
for(tp in timepoints) {
  if(tp %in% names(deg_list)) {
    deg_df <- deg_list[[tp]]
    total_degs <- sum(deg_df$change != "Not Significant")
    up_degs <- sum(deg_df$change == "Up-regulated")
    down_degs <- sum(deg_df$change == "Down-regulated")
    cat(tp, "vs mock: 总DEGs =", total_degs, 
        "(上调 =", up_degs, ", 下调 =", down_degs, ")\n")
  }
}

# 保存工作空间
save.image(file.path(output_dir, "MPXV_lb_timecourse_analysis_complete.RData"))
cat("\n分析完成！结果保存在:", output_dir, "目录中\n")
















