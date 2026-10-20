# 加载必要的包
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# ==========================================
# 1. 数据准备与整理
# ==========================================

# 确保 sample_info_cleaned 是数据框
coldata <- as.data.frame(sample_info_cleaned)

# 设置行名为样本名（DESeq2要求）
rownames(coldata) <- coldata$sample

# 将核心分组变量转换为因子 (Factor)
# 建议将 MOCK 设为参考水平 (Reference Level)
coldata$group_clear <- factor(coldata$group_clear)
coldata$strain <- factor(coldata$strain)
coldata$treatment <- factor(coldata$treatment, levels = c("MOCK", "MPXV_lb", "Tecovirimat"))

# 检查 count_matrix 是否存在
if(!exists("count_matrix")) stop("错误：环境中未找到 'count_matrix'，请先加载你的表达矩阵。")

# 确保 count_matrix 列顺序与 coldata 行顺序严格一致
# 这一步至关重要，否则结果会完全错乱
common_samples <- intersect(rownames(coldata), colnames(count_matrix))
coldata <- coldata[common_samples, ]
count_matrix_ordered <- count_matrix[, rownames(coldata)]

# 检查匹配情况
if(all(rownames(coldata) == colnames(count_matrix_ordered))) {
  cat("数据匹配成功：样本名顺序一致。\n")
} else {
  stop("错误：Count矩阵的列名与样本信息的行名不匹配！")
}

# ==========================================
# 2. 构建 DESeq2 对象与运行分析
# ==========================================

# 使用 group_clear 作为设计公式，这样可以灵活对比任意组
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_ordered,
  colData = coldata,
  design = ~ group_clear
)

# 预处理：过滤低表达基因
# 策略：至少在 3 个样本（最小分组样本数）中 count >= 10
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("过滤后保留基因数:", nrow(dds), "\n")

# 运行 DESeq2 标准流程
cat("正在运行 DESeq2 分析...\n")
dds <- DESeq(dds)

# ==========================================
# 3. 定义差异比较组合 (Contrasts)
# ==========================================
# 根据你的数据结构，定义以下具有生物学意义的比较：

comparisons <- list(
  # 1. AGB6 品系：感染 vs 对照
  "AGB6_Infection" = c("group_clear", "AGB6.MPXV_lb", "AGB6.MOCK"),
  
  # 2. C57 品系：感染 vs 对照
  "C57_Infection"  = c("group_clear", "C57.MPXV_lb", "C57.MOCK"),
  
  # 3. AGB6 品系：药物治疗 vs 病毒感染 (查看药物的回补/治疗效果)
  "AGB6_Treatment" = c("group_clear", "AGB6.Tecovirimat", "AGB6.MPXV_lb"),
  
  # 4. (可选) AGB6 品系：药物治疗 vs MOCK (查看药物毒性或是否恢复到正常水平)
  "AGB6_Teco_vs_Mock" = c("group_clear", "AGB6.Tecovirimat", "AGB6.MOCK")
)

# ==========================================
# 4. 提取结果与统计
# ==========================================

results_list <- list()

cat("\n=== 差异分析统计 ===\n")
for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  
  # 提取结果，设置阈值
  res <- results(dds, contrast = comp, alpha = 0.05)
  
  # 收缩 Log2FC (可选，推荐用于可视化，特别是样本少时，使用 ashr 或 apeglm)
  # res <- lfcShrink(dds, contrast = comp, res=res, type = "ashr") 
  
  # 添加基因列并排序
  res$gene <- rownames(res)
  res <- res[order(res$padj), ]
  
  results_list[[comp_name]] <- res
  
  # 统计显著基因
  up_genes <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
  down_genes <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)
  
  cat(sprintf("%-20s | 上调: %-5d | 下调: %-5d | 总计显著: %d\n", 
              comp_name, up_genes, down_genes, up_genes + down_genes))
}

# ==========================================
# 5. 数据可视化
# ==========================================

# 数据标准化 (Variance Stabilizing Transformation) 用于绘图
vsd <- vst(dds, blind = FALSE)

# --- A. PCA 图 (按品系和处理着色) ---
pca_data <- plotPCA(vsd, intgroup = c("strain", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = treatment, shape = strain)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = name), vjust = 2, size = 2.5, check_overlap = TRUE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Analysis: Strain & Treatment") +
  theme_bw() +
  theme(legend.position = "right")

print(pca_plot)

# --- B. 样本间相关性热图 ---
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample Distance Heatmap")

# --- C. 差异基因热图 (Top 30 差异最显著基因 - 示例: AGB6感染组) ---
select_comp <- "AGB6_Infection" # 可以修改为你关注的对比组
res_sig <- results_list[[select_comp]]

# 选取 padj 最小且 FC > 1 的前30个基因
top_genes_df <- as.data.frame(res_sig) %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(30)

if (nrow(top_genes_df) > 0) {
  top_genes <- rownames(top_genes_df)
  mat <- assay(vsd)[top_genes, ]
  mat <- mat - rowMeans(mat) # Z-score标准化
  
  # 设置热图注释条
  annotation_col <- coldata[, c("strain", "treatment")]
  
  pheatmap(mat,
           annotation_col = annotation_col,
           show_rownames = TRUE,
           show_colnames = FALSE,
           cluster_cols = TRUE, 
           main = paste("Top DEGs Heatmap -", select_comp))
} else {
  cat("该对比组没有显著差异基因，跳过热图绘制。\n")
}

# --- D. 火山图 (批量绘制) ---
for (comp_name in names(comparisons)) {
  res <- results_list[[comp_name]]
  
  # 处理可能的 NA 值
  res_clean <- res[!is.na(res$padj), ]
  
  p <- EnhancedVolcano(res_clean,
                       lab = rownames(res_clean),
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = paste('Volcano Plot:', comp_name),
                       subtitle = 'pCutoff=0.05, FC=2',
                       pCutoff = 0.05,
                       FCcutoff = 1, # log2FC=1 等于 FC=2
                       pointSize = 1.5,
                       labSize = 3.0,
                       legendPosition = 'right',
                       drawConnectors = TRUE,
                       widthConnectors = 0.5
  )
  print(p)
}

# ==========================================
# 6. 保存结果
# ==========================================

# 创建输出目录
if(!dir.exists("DESeq2_Results")) dir.create("DESeq2_Results")

for (comp_name in names(comparisons)) {
  res_df <- as.data.frame(results_list[[comp_name]])
  
  # 保存所有结果
  write.csv(res_df, 
            file = paste0("DESeq2_Results/", comp_name, "_all_genes.csv"), 
            row.names = TRUE)
  
  # 保存显著差异基因 (padj < 0.05 & |log2FC| > 1)
  sig_res <- res_df %>% 
    filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  write.csv(sig_res, 
            file = paste0("DESeq2_Results/", comp_name, "_significant_DEGs.csv"), 
            row.names = TRUE)
}

cat("\n分析全部完成！结果文件已保存在 'DESeq2_Results' 文件夹中。\n")