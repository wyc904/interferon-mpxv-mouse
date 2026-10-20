# ==============================================================================
# 0. 环境准备与包加载
# ==============================================================================
# 如果没有安装这些包，请先运行 BiocManager::install()
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)

# ==============================================================================
# 1. 数据清洗与构建 DESeq2 对象
# ==============================================================================

# --- 关键检查 ---
# 假设你的表达矩阵变量名为 count_matrix
if(!exists("count_matrix")) stop("请先加载你的基因Counts矩阵，并命名为 count_matrix")

# 1.1 整理样本信息 (colData)
coldata <- as.data.frame(sample_info_cleaned)
rownames(coldata) <- coldata$sample  # 将样本名设为行名

# 将分组变量转换为因子 (Factor)，这对于后续分析至关重要
coldata$group_clear <- factor(coldata$group_clear)
coldata$strain      <- factor(coldata$strain)


# 1.2 确保样本顺序一致性
# DESeq2要求：colData的行名必须与count_matrix的列名顺序严格一致
common_samples <- intersect(rownames(coldata), colnames(count_matrix))
if(length(common_samples) < nrow(coldata)) warning("注意：部分样本在Counts矩阵中未找到！")

coldata <- coldata[common_samples, ]
count_matrix_ordered <- count_matrix[, rownames(coldata)]

# 再次检查
if(!all(rownames(coldata) == colnames(count_matrix_ordered))) {
  stop("错误：Counts矩阵列名与样本信息行名不匹配，请检查数据！")
} else {
  cat(">>> 数据匹配成功，包含", ncol(count_matrix_ordered), "个样本，", nrow(count_matrix_ordered), "个基因。\n")
}

# 1.3 构建 DESeqDataSet 对象
# 使用 group_clear 作为设计公式，解决非平衡设计问题
dds <- DESeqDataSetFromMatrix(countData = count_matrix_ordered,
                              colData = coldata,
                              design = ~ group_clear)

# 1.4 预过滤
# 策略：至少在3个样本中（最小分组样本数）Count >= 10
# 这一步可以去除低表达噪音，提高统计效能
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat(">>> 过滤低表达基因后，剩余基因数:", nrow(dds), "\n")




# ==============================================================================
# 2. 差异表达分析 (Differential Expression Analysis)
# ==============================================================================

cat(">>> 正在运行 DESeq2 核心算法...\n")
dds <- DESeq(dds)

# 定义感兴趣的比较组合 (Contrasts)
# 格式: c("设计列名", "实验组(分子)", "对照组(分母)")
contrasts_list <- list(
  # 1. AGB6品系：病毒感染 vs Mock (查看感染引起的基因变化)
  "AGB6_Inf_vs_Mock" = c("group_clear", "AGB6.MPXV_lb", "AGB6.MOCK")
)




# 批量提取结果并保存
results_storage <- list()

for(comp_name in names(contrasts_list)){
  contrast_vec <- contrasts_list[[comp_name]]
  
  # 提取结果，默认 alpha(padj) < 0.05
  res <- results(dds, contrast = contrast_vec, alpha = 0.05)
  
  # LFC收缩 (可选，推荐用于可视化，特别是Volcano Plot，需要安装 apeglm 或 ashr)
  # res <- lfcShrink(dds, contrast = contrast_vec, res = res, type = 'ashr') 
  
  # 整理为Dataframe
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # 添加显著性标记
  res_df$regulation <- "NS"
  res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
  res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
  
  results_storage[[comp_name]] <- res_df
  
  # 打印简报
  cat(paste0("\n[", comp_name, "] 统计结果:\n"))
  print(table(res_df$regulation))
}

# ==============================================================================
# 3. 可视化 (Visualization)
# ==============================================================================

# 3.1 数据标准化 (VST转换) - 用于PCA和热图
vsd <- vst(dds, blind = FALSE)

# --- A. PCA 主成分分析图 ---
pcaData <- plotPCA(vsd, intgroup = c("strain", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = treatment, shape = strain)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") + # 使用鲜明的配色
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("PCA Plot: Strain & Treatment Effect") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(pca_plot)
ggsave("PCA_plot.pdf", pca_plot, width = 6, height = 5)



# --- B. 样本间距离热图 (QC) ---
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distance",
         filename = "Sample_Distance_Heatmap.pdf") # 直接保存

# --- C. 差异基因热图 (Top 50 DEGs) ---
# 选取 AGB6_Inf_vs_Mock 组中最显著的50个基因展示
target_comp <- "AGB6_Inf_vs_Mock"
res_target <- results_storage[[target_comp]]

# 排序并取前50
top_genes <- res_target %>% 
  arrange(padj) %>% 
  head(50) %>% 
  pull(gene_id)

# 提取表达矩阵并标准化 (Z-score)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

# 设置注释条
anno_col <- coldata[, c("strain", "treatment")]

pheatmap(mat, 
         annotation_col = anno_col, 
         show_rownames = TRUE, 
         show_colnames = FALSE,
         fontsize_row = 8,
         cluster_cols = TRUE,
         main = paste("Top 50 DEGs:", target_comp),
         filename = "Top50_DEGs_Heatmap.pdf")

# --- D. 火山图 (Volcano Plot) ---
# 为每个比较组画火山图
for(comp_name in names(results_storage)){
  res_df <- results_storage[[comp_name]]
  
  # 去除NA值
  res_df <- na.omit(res_df)
  
  p <- EnhancedVolcano(res_df,
                       lab = res_df$gene_id,
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = paste('Volcano Plot:', comp_name),
                       subtitle = 'Cutoff: pAdj<0.05, |Log2FC|>1',
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       pointSize = 1.5,
                       labSize = 3.0,
                       colAlpha = 0.6,
                       legendPosition = 'right'
  )
  
  print(p)
  ggsave(paste0("Volcano_", comp_name, ".png"), p, width = 8, height = 8)
}

# ==============================================================================
# 4. 结果导出
# ==============================================================================
dir.create("DESeq2_Results", showWarnings = FALSE)

for(comp_name in names(results_storage)){
  write.csv(results_storage[[comp_name]], 
            file = paste0("DESeq2_Results/", comp_name, "_All_Genes.csv"), 
            row.names = FALSE)
  
  # 仅导出显著基因
  sig_genes <- results_storage[[comp_name]] %>% 
    filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  write.csv(sig_genes, 
            file = paste0("DESeq2_Results/", comp_name, "_Sig_DEGs.csv"), 
            row.names = FALSE)
}

cat(">>> 分析完成！所有结果已保存至 DESeq2_Results 文件夹。\n")