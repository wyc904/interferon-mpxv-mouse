# 加载必要的包
library(DESeq2)
library(tidyverse)

# 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)


print("原始列名（第7-15列）:")
print(colnames(counts)[1:20])



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
  condition = factor(c(rep("MPXV1b", 12), rep("mock", 3)),
                     levels = c("mock", "MPXV1b"))
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

# 5. 按时间点进行差异分析 ----------------------------------------------
timepoints <- c("3d", "5d", "7d", "10d")
results_list <- list()
deg_list <- list()

for(tp in timepoints) {
  print(paste("分析时间点:", tp))
  
  # 选择当前时间点和mock的样本
  tp_samples <- sample_info$sample[sample_info$timepoint %in% c(tp, "mock")]
  dds_sub <- dds[, tp_samples]
  
  # 确保condition因子水平正确（mock作为参考）
  dds_sub$condition <- relevel(dds_sub$condition, ref = "mock")
  
  # 运行DESeq2
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "MPXV1b", "mock"), alpha = 0.05)
  
  # 保存完整结果
  results_list[[tp]] <- res
  
  # 提取显著差异基因
  res_sig <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    arrange(padj)
  
  deg_list[[tp]] <- res_sig
  print(paste(tp, "vs mock: 发现", nrow(res_sig), "个显著DEGs"))
}

# 6. 保存结果 ---------------------------------------------------------
output_dir <- "DESeq2_results_C57"
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

# 7. 基本可视化 --------------------------------------------------------
# 7.1 各时间点DEGs数量统计
deg_counts <- data.frame(
  timepoint = timepoints,
  up = sapply(timepoints, function(x) {
    if(x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$log2FoldChange > 1)
    } else {0}
  }),
  down = sapply(timepoints, function(x) {
    if(x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$log2FoldChange < -1)
    } else {0}
  })
)

deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count")

ggplot(deg_counts_long, aes(x = timepoint, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "C57野生型各时间点差异基因数量", 
       x = "时间点", y = "DEGs数量") +
  scale_fill_manual(values = c("down" = "blue", "up" = "red")) +
  theme_minimal()

# 7.2 样本相关性热图
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap::pheatmap(sampleDistMatrix,
                   clustering_distance_rows = sampleDists,
                   clustering_distance_cols = sampleDists,
                   main = "样本间距离热图",
                   annotation_col = data.frame(
                     Timepoint = sample_info$timepoint,
                     Condition = sample_info$condition,
                     row.names = sample_info$sample
                   ))

# 8. 干扰素通路基因分析 -------------------------------------------------
# 定义小鼠干扰素相关基因
ifn_genes <- c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Mx2", 
               "Oas1a", "Oas1b", "Oas2", "Oas3",
               "Stat1", "Stat2", "Irf7", "Irf9", "Ifnb1", "Ifna1", "Ifna2",
               "Irf1", "Irf3", "Irf5", "Jak1", "Jak2", "Tyk2",
               "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2")

# 检查哪些基因在我们的数据中
ifn_genes_present <- ifn_genes[ifn_genes %in% rownames(vsd)]

if(length(ifn_genes_present) > 0) {
  # 提取干扰素基因表达矩阵
  ifn_expr <- assay(vsd)[ifn_genes_present, ]
  ifn_expr_scaled <- t(scale(t(ifn_expr)))  # 行标准化
  
  # 创建注释数据框
  annotation_col <- data.frame(
    Timepoint = sample_info$timepoint,
    Condition = sample_info$condition,
    row.names = sample_info$sample
  )
  
  # 绘制热图
  pheatmap::pheatmap(ifn_expr_scaled,
                     annotation_col = annotation_col,
                     show_rownames = TRUE,
                     cluster_cols = TRUE,
                     main = "干扰素通路基因表达热图",
                     fontsize_row = 8)
  
  print(paste("找到", length(ifn_genes_present), "个干扰素相关基因"))
} else {
  print("未找到干扰素相关基因，请检查基因名")
}

# 9. 结果总结 ---------------------------------------------------------
cat("=== C57野生型MPXV 1b感染时序分析完成 ===\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总样本数:", ncol(dds), "\n")
cat("分析基因数:", nrow(dds), "\n\n")

cat("各时间点DEGs统计:\n")
for(tp in timepoints) {
  if(tp %in% names(deg_list)) {
    deg_count <- nrow(deg_list[[tp]])
    up_count <- sum(deg_list[[tp]]$log2FoldChange > 1)
    down_count <- sum(deg_list[[tp]]$log2FoldChange < -1)
    cat(tp, "vs mock: 总DEGs =", deg_count, 
        "(上调 =", up_count, ", 下调 =", down_count, ")\n")
  }
}

# 保存工作空间
save.image(file.path(output_dir, "C57_timecourse_analysis_complete.RData"))
cat("\n分析完成！结果保存在:", output_dir, "目录中\n")