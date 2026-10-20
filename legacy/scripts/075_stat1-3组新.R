# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop//Project/10-干扰素小鼠/1-源/STA1_counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 查看原始列名，确认样本顺序
print("原始列名（第7-15列）:")
print(colnames(counts)[1:20])

# 提取count矩阵（第7-15列是我们要分析的样本）
count_matrix <- counts[, 12:20]

# 根据您提供的信息，正确重命名列名
# 第7-9列: MPXV_IIb样本
# 第10-12列: MPXV_Ib样本  
# 第13-15列: MOCK样本

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("MPXV_IIb_1", "MPXV_IIb_2", "MPXV_IIb_3",
                            "MPXV_Ib_1", "MPXV_Ib_2", "MPXV_Ib_3",
                            "MOCK_1", "MOCK_2", "MOCK_3")

print("重命名后的列名:")
print(colnames(count_matrix))

# 验证样本对应关系
sample_mapping <- data.frame(
  Original_Name = colnames(counts)[12:20],
  New_Name = colnames(count_matrix),
  Group = c(rep("MPXV_IIb", 3), rep("MPXV_Ib", 3), rep("MOCK", 3))
)

print("样本对应关系:")
print(sample_mapping)

# 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  group = factor(c(rep("MPXV_IIb", 3), rep("MPXV_Ib", 3), rep("MOCK", 3)),
                 levels = c("MOCK", "MPXV_Ib", "MPXV_IIb")),
  replicate = rep(1:3, 3)
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
  design = ~ group
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
print("DESeq2分析完成！")

# 后续可以继续添加差异分析等代码...








# 6. 进行差异分析比较 ----------------------------------------------
comparisons <- list(
  c("MPXV_IIb", "MOCK"),
  c("MPXV_Ib", "MOCK"),
  c("MPXV_IIb", "MPXV_Ib")
)

results_list <- list()
deg_list <- list()

# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05

for(comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  comp_name <- paste0(group1, "_vs_", group2)
  
  print(paste("分析比较:", comp_name))
  
  # 确保参考水平正确
  dds$group <- relevel(dds$group, ref = group2)
  
  # 重新运行DESeq2
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", group1, group2), alpha = 0.05)
  
  # 保存完整结果
  results_list[[comp_name]] <- res
  
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
  
  deg_list[[comp_name]] <- res_df
  print(paste(comp_name, ": 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# 7. 可视化 --------------------------------------------------------

# 全局配色设置
color_mock <- "#708090"    # 灰色
color_IIb <- "#A52A2A"     # 深棕色 - MPXV IIb
color_Ib <- "#8B4513"      # 马鞍棕色 - MPXV Ib
color_up <- "#C75D69"      # 浅红
color_down <- "#4682B4"    # 蓝色
color_ns <- "#D3D3D3"      # 浅灰

# 7.1 PCA图（所有样本）
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("MOCK" = color_mock, "MPXV_Ib" = color_Ib, "MPXV_IIb" = color_IIb)) +
  scale_shape_manual(values = c("MOCK" = 16, "MPXV_Ib" = 17, "MPXV_IIb" = 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "MPXV IIb vs MPXV Ib vs MOCK") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  geom_text_repel(aes(label = name), 
                  size = 3.5, show.legend = FALSE, max.overlaps = 20,
                  fontface = "bold", box.padding = 0.3)

print(p_pca)
# ggsave("PCA_MPXV_Groups.pdf", p_pca, width = 8, height = 6)




# 7.2 各比较组火山图
for(comp_name in names(deg_list)) {
  res_df <- deg_list[[comp_name]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(30) %>% 
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
      labs(title = paste("Volcano:", comp_name), 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      theme_bw() + 
      theme(legend.position = "bottom")
    
    print(p_vol)
    
    ggsave(paste0("Volcano_", comp_name, ".pdf"), p_vol, width = 7, height = 6)
  }
}










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
  
  # 提取表达矩阵并标准化
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  
  # 创建样本注释
  annotation_col <- data.frame(
    Group = sample_info$group
  )
  rownames(annotation_col) <- colnames(mat)
  
  # 创建颜色映射
  annotation_colors <- list(
    Group = c("MOCK" = color_mock, "MPXV_Ib" = color_Ib, "MPXV_IIb" = color_IIb)
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
                     main = "MPXV各组差异表达基因热图")
  
  
  
  # 保存热图
  pdf("Heatmap_MPXV_Groups_DEGs.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  dev.off()
}






# 7.4 各比较组DEGs数量统计
deg_counts <- data.frame(
  comparison = names(deg_list),
  up = sapply(names(deg_list), function(x) {
    if(x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$change == "Up-regulated")
    } else {0}
  }),
  down = sapply(names(deg_list), function(x) {
    if(x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$change == "Down-regulated")
    } else {0}
  })
)

deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count")

p_bar <- ggplot(deg_counts_long, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "各比较组差异基因数量统计", 
       x = "比较组", y = "DEGs数量") +
  scale_fill_manual(values = c("down" = color_down, "up" = color_up)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_bar)
ggsave("DEGs_Count_By_Comparison.pdf", p_bar, width = 8, height = 6)

# 8. 干扰素通路基因分析 -------------------------------------------------
# 定义小鼠干扰素相关基因
ifn_genes <- c("Ifna1", "Ifna2", "Ifna4", "Ifna5", "Ifna6", "Ifna7", "Ifnb1", 
               "Ifng", "Ifnl2", "Ifnl3", "Irf1", "Irf3", "Irf7", "Irf9",
               "Stat1", "Stat2", "Isg15", "Isg20", "Mx1", "Mx2", "Oas1a", 
               "Oas1b", "Oas1g", "Oas2", "Oas3", "Oasl1", "Oasl2")

# 检查哪些基因在我们的数据中
ifn_genes_present <- ifn_genes[ifn_genes %in% rownames(vsd)]

if(length(ifn_genes_present) > 0) {
  # 提取干扰素基因表达矩阵
  ifn_expr <- assay(vsd)[ifn_genes_present, ]
  ifn_expr_scaled <- t(scale(t(ifn_expr)))  # 行标准化
  
  # 创建注释数据框
  annotation_col <- data.frame(
    Group = sample_info$group,
    row.names = sample_info$sample
  )
  
  # 创建颜色映射
  annotation_colors <- list(
    Group = c("MOCK" = color_mock, "MPXV_Ib" = color_Ib, "MPXV_IIb" = color_IIb)
  )
  
  # 绘制热图
  p_ifn <- pheatmap(ifn_expr_scaled,
                    annotation_col = annotation_col,
                    annotation_colors = annotation_colors,
                    show_rownames = TRUE,
                    cluster_cols = TRUE,
                    main = "干扰素通路基因表达热图",
                    fontsize_row = 8)
  
  # 保存干扰素热图
  pdf("Interferon_Pathway_Heatmap_MPXV.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_ifn$gtable)
  dev.off()
  
  print(paste("找到", length(ifn_genes_present), "个干扰素相关基因"))
} else {
  print("未找到干扰素相关基因，请检查基因名")
}

# 9. 保存结果 ---------------------------------------------------------
output_dir <- "DESeq2_MPXV









