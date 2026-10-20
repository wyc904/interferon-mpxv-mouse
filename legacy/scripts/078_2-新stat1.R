
library(dplyr)

library(tidyverse)

library(DESeq2)
library(ggrepel)






# 1.1 构建样本信息表 (根据你提供的新样本列表)
sample_info_raw <- data.frame(
  sample = c(
    # MPXV_lb 组
    "JZ25055844.4.1.P.4.1.P_combined",
    "JZ25055845.4.2.P.4.2.P_combined",
    "JZ25055846.4.3.P.4.3.P_combined",
    # mock 组
    "JZ25055847.5.1.P.5.1.P_combined",
    "JZ25055848.5.2.P.5.2.P_combined",
    "JZ25055849.5.3.P.5.3.P_combined"
  ),
  group_clear = c(
    "MPXV_lb", "MPXV_lb", "MPXV_lb",
    "mock", "mock", "mock"
  ),
  replicate = c(1, 2, 3, 1, 2, 3)
)

# 1.2 使用原始的group_clear名称，不进行映射
coldata <- sample_info_raw %>%
  mutate(
    condition = factor(group_clear, levels = c("mock", "MPXV_lb")) # 设定mock为参考水平
  ) %>%
  column_to_rownames("sample")

# 1.3 匹配 Count Matrix
if(!exists("count_matrix")) stop("错误：请先加载 count_matrix！")

# 检查样本匹配情况
target_samples <- rownames(coldata)
if(!all(target_samples %in% colnames(count_matrix))) {
  missing <- target_samples[!target_samples %in% colnames(count_matrix)]
  stop(paste("错误：以下样本在 count_matrix 中未找到:", paste(missing, collapse=", ")))
}

counts_subset <- count_matrix[, target_samples]
counts_subset <- counts_subset[, rownames(coldata)] # 确保顺序一致

cat(">>> 数据准备完成。分析组别: MPXV_lb vs mock\n")

# ==============================================================================
# 2. DESeq2 差异分析
# ==============================================================================

# 2.1 构建对象与过滤
dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3 # 预过滤
dds <- dds[keep, ]

# 2.2 运行 DESeq2
dds <- DESeq(dds)

# 2.3 提取结果 (MPXV_lb / mock)
# log2FC > 0 表示在 MPXV_lb 中上调
res <- results(dds, contrast = c("condition", "MPXV_lb", "mock"))

res_df <- as.data.frame(res) %>%
  mutate(gene_id = rownames(.)) %>%
  arrange(padj)

# 2.4 定义阈值与标记
FC_cutoff <- 1.0  # log2FoldChange > 1
P_cutoff <- 0.05  # padj < 0.05

res_df$change <- "Not Significant"
res_df$change[res_df$padj < P_cutoff & res_df$log2FoldChange > FC_cutoff] <- "Up-regulated"
res_df$change[res_df$padj < P_cutoff & res_df$log2FoldChange < -FC_cutoff] <- "Down-regulated"
res_df$change <- factor(res_df$change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))

# ==============================================================================
# 3. 可视化 (PCA, Volcano, Heatmap)
# ==============================================================================

# 全局配色
color_mock <- "#708090"  # 灰色
color_mpxv <- "#A52A2A"  # 褐色
color_up   <- "#C75D69"  # 浅红
color_down <- "#4682B4"  # 蓝色
color_ns   <- "#D3D3D3"  # 浅灰

# --- 图 A: PCA Plot ---
vsd <- vst(dds, blind = FALSE)
# 从样本名称中提取重复编号
pcaData$replicate <- gsub(".*\\.(\\d)\\..*", "\\1", pcaData$name)
# 创建新的标签：分组名称+重复编号
pcaData$new_label <- paste0(pcaData$condition, "_", pcaData$replicate)

# 然后用new_label作为标签
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_shape_manual(values = c("mock" = 16, "MPXV_lb" = 17)) + 
  scale_color_manual(values = c("mock" = color_mock, "MPXV_lb" = color_mpxv)) +
  geom_text_repel(aes(label = new_label), size = 3, show.legend = FALSE, nudge_y = 0.5) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "MPXV lb vs mock") +
  theme_bw()

print(p_pca)
# ggsave("Figure_A_PCA.pdf", p_pca, width = 6, height = 5)

# --- 图 B: Volcano Plot ---
# 选取最显著的 Top 15 基因做标签
top_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% head(15) %>% pull(gene_id)

p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = color_down, 
                                "Not Significant" = color_ns, 
                                "Up-regulated" = color_up)) +
  geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = subset(res_df, gene_id %in% top_genes),
                  aes(label = gene_id, color = change),
                  size = 2.5, show.legend = FALSE, max.overlaps = 15) +
  labs(title = "MPXV lb vs mock", 
       x = expression(Log[2]~Fold~Change), 
       y = expression(-Log[10]~Adjusted~P-value)) +
  theme_bw() + theme(legend.position = "bottom")

print(p_vol)
# ggsave("Figure_B_Volcano.pdf", p_vol, width = 7, height = 6)











# --- 图 C: Heatmap ---
# 选取 Top 30 显著基因
# --- 图 C: Heatmap ---
# 选取 Top 30 显著基因
top30_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% head(30) %>% pull(gene_id)

if(length(top30_genes) > 2) {
  # 提取并 Z-score 标准化
  mat <- assay(vsd)[top30_genes, ]
  mat <- t(scale(t(mat)))
  
  # 创建与PCA图相同的简洁样本标签
  # 从原始样本名中提取重复编号
  sample_replicates <- gsub(".*\\.(\\d)\\..*", "\\1", colnames(mat))
  
  # 获取对应的分组信息
  sample_conditions <- coldata[colnames(mat), "condition"]
  
  # 创建简洁的样本标签：分组名称+重复编号
  new_colnames <- paste0(sample_conditions, "_", sample_replicates)
  colnames(mat) <- new_colnames
  
  # 注释条设置 - 使用原始group_clear名称
  anno_col <- data.frame(Group = sample_conditions)
  rownames(anno_col) <- new_colnames
  
  anno_colors <- list(Group = c("mock" = color_mock, "MPXV_lb" = color_mpxv))
  
  # 绘制热图
  p_heat <- pheatmap(mat,
                     color = colorRampPalette(c("#36648B", "white", "#E9967A"))(100),
                     annotation_col = anno_col,
                     annotation_colors = anno_colors,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, # 保持样本顺序
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     border_color = NA,
                     fontsize_row = 9,
                     fontsize_col = 10,  # 增加列名字体大小
                     main = "Top 30 DEGs: MPXV lb vs mock")
  
  # 保存热图
  save_pheatmap_pdf <- function(x, filename, width=7, height=9) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  save_pheatmap_pdf(p_heat, "Figure_C_Heatmap.pdf")
  
  # 同时保存一个更宽的热图，确保标签完全显示
  p_heat_wide <- pheatmap(mat,
                          color = colorRampPalette(c("#36648B", "white", "#E9967A"))(100),
                          annotation_col = anno_col,
                          annotation_colors = anno_colors,
                          cluster_rows = TRUE,
                          cluster_cols = FALSE,
                          show_colnames = TRUE,
                          show_rownames = TRUE,
                          border_color = NA,
                          fontsize_row = 9,
                          fontsize_col = 10,
                          main = "Top 30 DEGs: MPXV lb vs mock")
  
  save_pheatmap_pdf(p_heat_wide, "Figure_C_Heatmap_wide.pdf", width=9, height=9)
  
} else {
  cat("显著差异基因过少，无法绘制热图。\n")
}

# 保存结果
write.csv(res_df, "DESeq2_results_MPXV_lb_vs_mock.csv", row.names = FALSE)

cat(">>> 分析完成！结果文件已保存。\n")












