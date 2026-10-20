# ==============================================================================
# 1. 数据准备与筛选 (只保留 AGB6 组)
# ==============================================================================

# 1.1 构建完整的样本信息表
sample_info_raw <- data.frame(
  sample = c(
    "JZ25070166.AGB6.MOCK.1.P.AGB6.MOCK.1.P_combined",
    "JZ25070167.AGB6.MOCK.2.P.AGB6.MOCK.2.P_combined", 
    "JZ25070168.AGB6.MOCK.3.P.AGB6.MOCK.3.P_combined",
    "JZ25070169.AGB6.MPXV_lb.1.P.AGB6.MPXV_lb.1.P_combined",
    "JZ25070170.AGB6.MPXV_lb.2.P.AGB6.MPXV_lb.2.P_combined",
    "JZ25070171.AGB6.MPXV_lb.3.P.AGB6.MPXV_lb.3.P_combined"
  ),
  group_clear = c(
    "AGB6.MOCK", "AGB6.MOCK", "AGB6.MOCK",
    "AGB6.MPXV_lb", "AGB6.MPXV_lb", "AGB6.MPXV_lb"
  ),
  replicate = c(1, 2, 3, 1, 2, 3)
)

# 1.2 添加绘图用的分组标签 (改为 MOCK vs MPXV_lb)
coldata <- sample_info_raw %>%
  mutate(
    condition = ifelse(grepl("MOCK", group_clear), "MOCK", "MPXV_lb"),
    condition = factor(condition, levels = c("MOCK", "MPXV_lb")) # MOCK为基准
  ) %>%
  column_to_rownames("sample")

# 1.3 匹配 Count Matrix
# 注意：请确保你的 count_matrix 已经加载到 R 环境中
if(!exists("count_matrix")) stop("错误：请先加载 count_matrix！")

# 只保留选中的这6个样本
target_samples <- rownames(coldata)
# 检查样本是否存在于矩阵中
if(!all(target_samples %in% colnames(count_matrix))) stop("错误：部分样本在 count_matrix 中找不到！")

counts_subset <- count_matrix[, target_samples]

# 再次确认顺序一致
counts_subset <- counts_subset[, rownames(coldata)]

cat(">>> 数据准备完成。样本数:", ncol(counts_subset), " (3 MOCK, 3 MPXV_lb)\n")

# ==============================================================================
# 2. DESeq2 差异分析
# ==============================================================================

# 2.1 构建 DESeq对象
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = coldata,
                              design = ~ condition)

# 2.2 过滤低表达基因 (至少3个样本count>10)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# 2.3 运行分析
dds <- DESeq(dds)

# 2.4 提取结果 (MPXV_lb vs MOCK)
res <- results(dds, contrast = c("condition", "MPXV_lb", "MOCK"))
res_df <- as.data.frame(res) %>%
  mutate(gene_id = rownames(.)) %>%
  arrange(padj)

# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05

# 标记上下调
res_df$change <- "Not Significant"
res_df$change[res_df$padj < P_cutoff & res_df$log2FoldChange > FC_cutoff] <- "Up-regulated"
res_df$change[res_df$padj < P_cutoff & res_df$log2FoldChange < -FC_cutoff] <- "Down-regulated"
res_df$change <- factor(res_df$change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))

# ==============================================================================
# 3. 可视化复刻 (图 A, B, C)
# ==============================================================================

# --- 全局配色设置 ---
# 更新配色以匹配新的分组名称
color_mock <- "#708090"    # 灰色 (SlateGray)
color_mpxv_lb <- "#A52A2A" # 红色/褐色 (Brown)
color_up      <- "#C75D69"  # 浅红
color_down    <- "#4682B4"  # 蓝色
color_ns      <- "#D3D3D3"  # 浅灰

# ------------------------------------------------------------------------------
# 图 A: PCA Plot (散点图)
# ------------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# 自定义绘图
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  # 使用新的分组名称
  scale_shape_manual(values = c("MOCK" = 16, "MPXV_lb" = 17)) + 
  scale_color_manual(values = c("MOCK" = color_mock, "MPXV_lb" = color_mpxv_lb)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "MPXV lb vs MOCK") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face="bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0, face = "bold", margin = margin(b = 15))  # 标题左对齐
  ) +
  # 添加样本标签 - 使用自定义格式
  geom_text_repel(
    aes(label = ifelse(condition == "MOCK", 
                       paste0("MOCK_", gsub(".*\\.(\\d)\\..*", "\\1", name)),
                       paste0("MPXV_lb_", gsub(".*\\.(\\d)\\..*", "\\1", name)))),
    size = 3, 
    show.legend = FALSE, 
    nudge_y = 0.5
  )

print(p_pca)
# ggsave("Figure_A_PCA_MPXV_lb_vs_MOCK.pdf", p_pca, width = 6, height = 5)

# ------------------------------------------------------------------------------
# 图 B: Volcano Plot (火山图)
# ------------------------------------------------------------------------------
# 选取最显著的基因进行标记
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
  labs(title = "MPXV lb vs MOCK", 
       x = expression(Log[2]~Fold~Change), 
       y = expression(-Log[10]~Adjusted~P-value)) +
  theme_bw() + 
  theme(legend.position = "bottom")

print(p_vol)
# ggsave("Figure_B_Volcano_MPXV_lb_vs_MOCK.pdf", p_vol, width = 7, height = 6)

# ------------------------------------------------------------------------------
# 图 C: Heatmap (热图)
# ------------------------------------------------------------------------------
# 选取Top 30显著差异基因
top30_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% 
  head(30) %>% 
  pull(gene_id)

if(length(top30_genes) > 2) {
  # 提取并Z-score标准化
  mat <- assay(vsd)[top30_genes, ]
  mat <- t(scale(t(mat)))
  
  # 创建样本标签
  sample_replicates <- gsub(".*\\.(\\d)\\..*", "\\1", colnames(mat))
  sample_conditions <- coldata[colnames(mat), "condition"]
  new_colnames <- paste0(sample_conditions, "_", sample_replicates)
  colnames(mat) <- new_colnames
  
  # 注释条设置
  anno_col <- data.frame(Condition = sample_conditions)
  rownames(anno_col) <- new_colnames
  
  anno_colors <- list(
    Condition = c("MOCK" = color_mock, "MPXV_lb" = color_mpxv_lb)
  )
  
  # 绘制热图
  p_heat <- pheatmap(mat,
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
                     main = "MPXV lb vs MOCK")
  
  # 保存热图
  save_pheatmap_pdf <- function(x, filename, width=7, height=9) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  save_pheatmap_pdf(p_heat, "Figure_C_Heatmap_MPXV_lb_vs_MOCK.pdf")
}

# ==============================================================================
# 4. 结果保存
# ==============================================================================

# 保存完整结果
write.csv(res_df, "DESeq2_results_MPXV_lb_vs_MOCK.csv", row.names = FALSE)

# 保存显著差异基因
sig_genes <- res_df %>% filter(change != "Not Significant")
write.csv(sig_genes, "Significant_DEGs_MPXV_lb_vs_MOCK.csv", row.names = FALSE)

# 保存样本信息
write.csv(coldata, "Sample_info_MPXV_lb_vs_MOCK.csv")

cat(">>> 分析完成！所有结果文件已保存。\n")
cat(">>> 比较组: MPXV_lb vs MOCK\n")