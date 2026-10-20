# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
packages <- c("DESeq2", "dplyr", "ggplot2", "pheatmap", "ggrepel", "RColorBrewer", "tibble")
lapply(packages, library, character.only = TRUE)

# ==============================================================================
# 1. 数据准备与筛选 (MPXV-Ilb vs MPXV-lb)
# ==============================================================================

# 1.1 构建样本信息表 (3组 vs 4组)
sample_info_raw <- data.frame(
  sample = c(
    # MPXV-Ilb 组 (Group 3)
    "JZ25055841.3.1.P.3.1.P_combined",
    "JZ25055842.3.2.P.3.2.P_combined", 
    "JZ25055843.3.3.P.3.3.P_combined",
    # MPXV-lb 组 (Group 4)
    "JZ25055844.4.1.P.4.1.P_combined",
    "JZ25055845.4.2.P.4.2.P_combined",
    "JZ25055846.4.3.P.4.3.P_combined"
  ),
  group_clear = c(
    "MPXV_Ilb", "MPXV_Ilb", "MPXV_Ilb",
    "MPXV_lb", "MPXV_lb", "MPXV_lb"
  ),
  replicate = c(1, 2, 3, 1, 2, 3)
)

# 1.2 映射分组标签 (MPXV_lb -> Control, MPXV_Ilb -> Treatment)
coldata <- sample_info_raw %>%
  mutate(
    condition = case_when(
      group_clear == "MPXV_lb" ~ "MPXV_lb",
      group_clear == "MPXV_Ilb" ~ "MPXV_Ilb"
    ),
    condition = factor(condition, levels = c("MPXV_lb", "MPXV_Ilb")) # 设定 Control 为参考水平
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

cat(">>> 数据准备完成。分析组别: MPXV_Ilb (Treatment) vs MPXV_lb (Control)\n")
cat(">>> 样本数量:", ncol(counts_subset), "\n")

# ==============================================================================
# 2. DESeq2 差异分析
# ==============================================================================

# 2.1 构建对象与过滤
dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 2 # 调整为至少2个样本有10个counts
dds <- dds[keep, ]

cat(">>> 过滤后保留基因数:", nrow(dds), "\n")

# 2.2 运行 DESeq2
dds <- DESeq(dds)

# 2.3 提取结果 (Treatment / Control)
# log2FC > 0 表示在 MPXV_Ilb 中上调
res <- results(dds, contrast = c("condition", "Treatment", "Control"))

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

# 2.5 输出统计信息
sig_summary <- table(res_df$change)
cat(">>> 差异基因统计:\n")
cat("上调基因:", sig_summary["Up-regulated"], "\n")
cat("下调基因:", sig_summary["Down-regulated"], "\n")
cat("不显著基因:", sig_summary["Not Significant"], "\n")

# ==============================================================================
# 3. 可视化 (PCA, Volcano, Heatmap)
# ==============================================================================

# 全局配色
color_control <- "#708090"  # 灰色
color_treat   <- "#8B4513"  # 深棕色
color_up      <- "#C75D69"  # 浅红
color_down    <- "#4682B4"  # 蓝色
color_ns      <- "#D3D3D3"  # 浅灰

# --- 图 A: PCA Plot ---
vsd <- vst(dds, blind = FALSE)
pcaData$replicate <- gsub(".*\\.(\\d)\\..*", "\\1", pcaData$name)
pcaData$new_label <- paste0(pcaData$condition, "_", pcaData$replicate)

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_shape_manual(values = c("Control" = 16, "Treatment" = 17)) + 
  scale_color_manual(values = c("Control" = color_control, "Treatment" = color_treat)) +
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE, nudge_y = 0.5) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "PCA: MPXV Ilb vs MPXV lb") +
  theme_bw()

print(p_pca)
# ggsave("Figure_A_PCA_MPXV_Ilb_vs_lb.pdf", p_pca, width = 6, height = 5)

# --- 图 B: Volcano Plot ---
# 选取最显著的 Top 15 基因做标签
top_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% 
  head(15) %>% 
  pull(gene_id)

p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = color_down, 
                                "Not Significant" = color_ns, 
                                "Up-regulated" = color_up)) +
  geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = subset(res_df, gene_id %in% top_genes),
                  aes(label = gene_id, color = change),
                  size = 3, show.legend = FALSE, max.overlaps = 20) +
  labs(title = "Volcano: MPXV_Ilb vs MPXV_lb", 
       x = expression(Log[2]~Fold~Change), 
       y = expression(-Log[10]~Adjusted~P-value)) +
  theme_bw() + 
  theme(legend.position = "bottom")

print(p_vol)
# ggsave("Figure_B_Volcano_MPXV_Ilb_vs_lb.pdf", p_vol, width = 7, height = 6)

# --- 图 C: Heatmap ---
# 选取 Top 30 显著基因
top30_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% 
  head(30) %>% 
  pull(gene_id)

if(length(top30_genes) > 2) {
  # 提取并 Z-score 标准化
  mat <- assay(vsd)[top30_genes, ]
  mat <- t(scale(t(mat)))
  
  # 创建唯一的列名用于热图显示
  unique_colnames <- paste0(coldata$condition, "_", coldata$group_clear, ".", 1:ncol(mat))
  colnames(mat) <- unique_colnames
  
  # 注释条设置
  anno_col <- data.frame(
    Condition = coldata$condition,
    Group = coldata$group_clear
  )
  rownames(anno_col) <- unique_colnames
  
  anno_colors <- list(
    Condition = c("Control" = color_control, "Treatment" = color_treat),
    Group = c("MPXV_lb" = "#1b9e77", "MPXV_Ilb" = "#d95f02")
  )
  
  p_heat <- pheatmap(mat,
                     color = colorRampPalette(c("#36648B", "white", "#E9967A"))(100),
                     annotation_col = anno_col,
                     annotation_colors = anno_colors,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     border_color = NA,
                     fontsize_row = 8,
                     fontsize_col = 8,
                     main = "Top 30 DEGs: MPXV_Ilb vs MPXV_lb")
  
  # 保存热图
  save_pheatmap_pdf <- function(x, filename, width=7, height=10) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  save_pheatmap_pdf(p_heat, "Figure_C_Heatmap_MPXV_Ilb_vs_lb.pdf")
  
} else {
  cat("显著差异基因过少，无法绘制热图。\n")
}

# ==============================================================================
# 4. 结果保存
# ==============================================================================

# 4.1 保存完整结果
write.csv(res_df, "DESeq2_results_MPXV_Ilb_vs_lb.csv", row.names = FALSE)

# 4.2 保存显著差异基因
sig_genes <- res_df %>% filter(change != "Not Significant")
write.csv(sig_genes, "Significant_DEGs_MPXV_Ilb_vs_lb.csv", row.names = FALSE)

# 4.3 保存上调和下调基因
up_genes <- sig_genes %>% filter(change == "Up-regulated")
down_genes <- sig_genes %>% filter(change == "Down-regulated")

write.csv(up_genes, "Upregulated_genes_MPXV_Ilb_vs_lb.csv", row.names = FALSE)
write.csv(down_genes, "Downregulated_genes_MPXV_Ilb_vs_lb.csv", row.names = FALSE)

# 4.4 保存样本信息
write.csv(coldata, "Sample_info_MPXV_Ilb_vs_lb.csv")

cat(">>> 分析完成！所有结果文件已保存。\n")
cat(">>> 比较组: MPXV_Ilb (Treatment) vs MPXV_lb (Control)\n")












