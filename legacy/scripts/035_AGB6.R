# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
# 检查并加载必要的包
packages <- c("DESeq2", "dplyr", "ggplot2", "pheatmap", "ggrepel", "RColorBrewer", "tibble")
# 如果需要做富集分析(图D)，需要以下包:
# packages <- c(packages, "clusterProfiler", "org.Mm.eg.db", "enrichplot") 

lapply(packages, library, character.only = TRUE)

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

# 1.2 添加绘图用的分组标签 (模仿图A图例: Control vs Infection)
coldata <- sample_info_raw %>%
  mutate(
    condition = ifelse(grepl("MOCK", group_clear), "Control", "Infection"),
    condition = factor(condition, levels = c("Control", "Infection")) # Control为基准
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

cat(">>> 数据准备完成。样本数:", ncol(counts_subset), " (3 Control, 3 Infection)\n")

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

# 2.4 提取结果 (Infection vs Control)
res <- results(dds, contrast = c("condition", "Infection", "Control"))
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
# 模仿图中配色
color_control <- "#708090"  # 灰色 (SlateGray)
color_infect  <- "#A52A2A"  # 红色/褐色 (Brown)
color_up      <- "#C75D69"  # 浅红 (参考图B)
color_down    <- "#4682B4"  # 蓝色 (SteelBlue, 参考图B)
color_ns      <- "#D3D3D3"  # 浅灰

# ------------------------------------------------------------------------------
# 图 A: PCA Plot (散点图)
# ------------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition", "replicate"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 自定义绘图
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  # 模仿图A：Control用圆点，Infection用三角形
  scale_shape_manual(values = c("Control" = 16, "Infection" = 17)) + 
  scale_color_manual(values = c("Control" = color_control, "Infection" = color_infect)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face="bold"),
    axis.title = element_text(size = 12)
  ) +
  # 添加样本标签 (可选，类似图中文字)
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE, nudge_y = 0.5)

# ggsave("Figure_A_PCA.pdf", p_pca, width = 6, height = 5)
print(p_pca)



dev.off()


# ------------------------------------------------------------------------------
# 图 B: Volcano Plot (火山图)
# ------------------------------------------------------------------------------
# 筛选 top 基因用于标注 (例如 padj 最小的前10个 + log2FC 最大的前10个)
top_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% 
  head(15) %>% 
  pull(gene_id)

p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  # 添加散点
  geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
  # 颜色映射
  scale_color_manual(values = c("Down-regulated" = color_down, 
                                "Not Significant" = color_ns, 
                                "Up-regulated" = color_up)) +
  # 添加阈值线
  geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "grey50") +
  # 添加基因标签
  geom_text_repel(data = subset(res_df, gene_id %in% top_genes),
                  aes(label = gene_id, color = change),
                  size = 3, show.legend = FALSE, max.overlaps = 20) +
  labs(title = "Infection vs. Control", x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~Adjusted~P-value)) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("Figure_B_Volcano.pdf", p_vol, width = 7, height = 6)
print(p_vol)

# ------------------------------------------------------------------------------
# 图 C: Heatmap (热图)
# ------------------------------------------------------------------------------
# 选取 top 30 显著差异基因 (按 padj 排序)
top30_genes <- res_df %>% 
  filter(change != "Not Significant") %>%
  arrange(padj) %>% 
  head(30) %>% 
  pull(gene_id)

# 提取标准化后的表达量
mat <- assay(vsd)[top30_genes, ]
# Z-score 标准化 (按行)
mat <- t(scale(t(mat)))

# 列注释
anno_col <- data.frame(Group = coldata$condition)
rownames(anno_col) <- rownames(coldata$group_clear)
anno_colors <- list(Group = c("Control" = color_control, "Infection" = color_infect))

# 绘图 (模仿图C配色: 橘色高表达，深青色低表达)
# 这里的 break 需要手动调整以匹配 Z-score 范围 (-2 到 2)
p <-pheatmap(mat,
         color = colorRampPalette(c("#36648B", "white", "#E9967A"))(100), # 深青-白-深橘
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE, # 保持样本顺序
         show_colnames = TRUE,
         show_rownames = TRUE,
         border_color = NA,
         fontsize_row = 9,
         main = "Top 30 DEGs",
         filename = "Figure_C_Heatmap.pdf", width = 5, height = 7)
print(p)
p

# dev.off()
# ------------------------------------------------------------------------------
# 图 D: Enrichment Plot (富集分析) - 代码框架
# ------------------------------------------------------------------------------
# 注意：这部分需要你的基因名为 Symbol 或 EntrezID，并且安装了 org.Mm.eg.db
# 下面代码仅为示例，如果不需要运行富集分析可忽略





# GO分析转到GO分析代码













