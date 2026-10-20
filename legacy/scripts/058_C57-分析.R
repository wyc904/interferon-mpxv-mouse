# 加载必要的包
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# 修正分组信息 - 从group列中提取正确的时间点信息
sample_info$condition <- sapply(strsplit(sample_info$group, "_"), function(x) x[1])
sample_info$condition <- factor(sample_info$condition, 
                                levels = c("mock", "3d", "5d", "7d", "10d"))

# 检查修正后的分组
cat("修正后的分组信息:\n")
print(table(sample_info$condition))

# 创建coldata
coldata <- data.frame(
  condition = sample_info$condition,
  row.names = sample_info$sample
)

# 确保count矩阵列顺序与coldata一致
count_matrix_ordered <- count_matrix[, rownames(coldata)]

# 再次检查确保没有NA值
cat("condition列中的NA值数量:", sum(is.na(coldata$condition)), "\n")

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_ordered,
  colData = coldata,
  design = ~ condition
)


# 预处理：过滤低表达基因（至少在3个样本中count>10）
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("过滤后保留基因数:", nrow(dds), "\n")

# 运行DESeq2分析
cat("开始DESeq2分析...\n")
dds <- DESeq(dds)

# 获取所有比较组合的结果
comparisons <- list(
  "3d_vs_mock" = c("condition", "3d", "mock"),
  "5d_vs_mock" = c("condition", "5d", "mock"),
  "7d_vs_mock" = c("condition", "7d", "mock"),
  "10d_vs_mock" = c("condition", "10d", "mock")
)

# 存储所有结果
results_list <- list()

# 进行每组的差异分析
for (comp_name in names(comparisons)) {
  cat("处理比较:", comp_name, "\n")
  
  comp <- comparisons[[comp_name]]
  res <- results(dds, contrast = comp, alpha = 0.05)
  
  # 添加基因名
  res$gene <- rownames(res)
  
  # 按p值排序
  res <- res[order(res$padj), ]
  
  # 保存结果
  results_list[[comp_name]] <- res
  
  # 输出显著差异基因统计
  sig_genes <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
  cat(comp_name, "显著差异基因数:", sig_genes, "\n")
}



# 数据标准化用于可视化
vsd <- vst(dds, blind = FALSE)  # 方差稳定变换

# 1. PCA图
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = 1.5, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot - 样本聚类")

print(pca_plot)





# 2. 样本相关性热图
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "样本间距离热图")

# 3. 差异基因热图（选取前20个显著基因）
# 选择一个比较组来展示
select_comparison <- "3d_vs_mock"
res_sig <- results_list[[select_comparison]]
top_genes <- res_sig[res_sig$padj < 0.05 & abs(res_sig$log2FoldChange) > 1, ]
if (nrow(top_genes) > 20) {
  top_genes <- top_genes[1:20, ]
}

if (nrow(top_genes) > 0) {
  mat <- assay(vsd)[rownames(top_genes), ]
  mat <- mat - rowMeans(mat)  # 中心化
  
  pheatmap(mat,
           annotation_col = coldata,
           main = paste("差异表达基因热图 -", select_comparison))
}

# 4. 火山图
for (comp_name in names(comparisons)) {
  res <- results_list[[comp_name]]
  
  volcano_plot <- EnhancedVolcano(res,
                                  lab = rownames(res),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = paste('火山图 -', comp_name),
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 1.0,
                                  labSize = 3.0
  )
  print(volcano_plot)
}

# 5. 保存结果到文件
for (comp_name in names(comparisons)) {
  res_df <- as.data.frame(results_list[[comp_name]])
  
  # 保存完整结果
  write.csv(res_df, 
            file = paste0(comp_name, "_full_results.csv"),
            row.names = FALSE)
  
  # 保存显著差异基因
  sig_res <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
  write.csv(sig_res,
            file = paste0(comp_name, "_significant_genes.csv"),
            row.names = FALSE)
}

# 6. 差异基因统计汇总
summary_table <- data.frame(
  Comparison = names(comparisons),
  Total_Genes = sapply(results_list, nrow),
  Upregulated = sapply(results_list, function(x) {
    sum(x$padj < 0.05 & x$log2FoldChange > 1, na.rm = TRUE)
  }),
  Downregulated = sapply(results_list, function(x) {
    sum(x$padj < 0.05 & x$log2FoldChange < -1, na.rm = TRUE)
  })
)

cat("\n=== 差异表达基因统计 ===\n")
print(summary_table)

# 保存统计汇总
write.csv(summary_table, "DEG_summary_statistics.csv", row.names = FALSE)

cat("\n分析完成！结果文件已保存。\n")