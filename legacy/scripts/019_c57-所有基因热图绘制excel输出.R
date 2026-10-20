# 导出所有基因的标准化表达矩阵Excel文件 -------------------------------------------------

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(writexl)
library(openxlsx)  # 如果需要导出到xlsx格式

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
print(paste("总基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行DESeq2分析 ----------------------------------------------------
dds <- DESeq(dds)

# 6. 进行方差稳定化转换（VST）-----------------------------------------
vsd <- vst(dds, blind = FALSE)

# 7. 提取所有基因的标准化表达矩阵 -------------------------------------
# 获取VST转换后的表达矩阵
vsd_matrix <- assay(vsd)
print(paste("VST矩阵维度:", dim(vsd_matrix)[1], "个基因 x", dim(vsd_matrix)[2], "个样本"))

# 8. 创建并导出所有基因的表达矩阵Excel文件 ----------------------------
# 转换为数据框，基因名作为第一列
vsd_df <- as.data.frame(vsd_matrix) %>%
  rownames_to_column("Gene")

# 查看前几行
cat("\n标准化表达矩阵前5行预览:\n")
print(head(vsd_df))

# 8.1 创建包含多个工作表的Excel文件
excel_data <- list()

# 工作表1: 所有基因的标准化表达矩阵
excel_data[["All_Genes_VST_Expression"]] <- vsd_df

# 工作表2: 样本信息
sample_info_export <- sample_info %>%
  dplyr::select(sample, timepoint, condition, group_label) %>%
  dplyr::rename(
    Sample = sample,
    Timepoint = timepoint,
    Condition = condition,
    Group = group_label
  )
excel_data[["Sample_Information"]] <- sample_info_export

# 工作表3: 原始计数矩阵（过滤后）
count_df <- as.data.frame(counts(dds)) %>%
  rownames_to_column("Gene")
excel_data[["All_Genes_Raw_Counts"]] <- count_df

# 8.2 按时间点进行差异分析并添加到Excel
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
  
  # 添加到Excel工作表
  res_df_export <- res_df %>%
    dplyr::rename(
      Gene = gene,
      BaseMean = baseMean,
      Log2FoldChange = log2FoldChange,
      Log2FC_SE = lfcSE,
      Statistic = stat,
      PValue = pvalue,
      AdjustedPValue = padj,
      Regulation = change
    ) %>%
    mutate(
      Log2FoldChange = round(Log2FoldChange, 3),
      Log2FC_SE = round(Log2FC_SE, 3),
      Statistic = round(Statistic, 3),
      PValue = format(PValue, scientific = TRUE, digits = 3),
      AdjustedPValue = format(AdjustedPValue, scientific = TRUE, digits = 3)
    )
  
  excel_data[[paste0("DEGs_", tp, "_vs_Mock")]] <- res_df_export
  
  print(paste(tp, "vs mock: 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# 9. 导出到Excel文件 ----------------------------------------------------
excel_filename <- "All_Genes_Expression_Analysis.xlsx"

# 使用write_xlsx导出
write_xlsx(
  excel_data,
  path = excel_filename
)

print(paste("Excel文件已保存为:", excel_filename))
print(paste("包含", length(excel_data), "个工作表:"))

# 打印工作表信息
for(sheet_name in names(excel_data)) {
  print(paste("-", sheet_name, ":", nrow(excel_data[[sheet_name]]), "行 x", ncol(excel_data[[sheet_name]]), "列"))
}

# 10. 可选：导出单独的热图数据CSV文件 ------------------------------------
# 10.1 导出所有基因的标准化表达矩阵（用于热图绘制）
write.csv(vsd_df, "All_Genes_VST_Expression_Matrix.csv", row.names = FALSE)

# 10.2 导出Z-score标准化后的矩阵（用于热图）
# 对行（基因）进行Z-score标准化
zscore_matrix <- t(scale(t(vsd_matrix)))
zscore_df <- as.data.frame(zscore_matrix) %>%
  rownames_to_column("Gene")
write.csv(zscore_df, "All_Genes_Zscore_Matrix.csv", row.names = FALSE)

# 10.3 导出差异分析汇总表
# 创建一个汇总表，显示每个基因在不同时间点的log2FC和调节方向
all_genes <- rownames(dds)
summary_df <- data.frame(Gene = all_genes)

for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  # 添加log2FoldChange
  log2fc_col <- paste0("log2FC_", tp)
  padj_col <- paste0("padj_", tp)
  regulation_col <- paste0("Regulation_", tp)
  
  # 匹配基因并提取数据
  matched_data <- res_df[match(all_genes, res_df$gene), ]
  
  summary_df[[log2fc_col]] <- round(matched_data$log2FoldChange, 3)
  summary_df[[padj_col]] <- format(matched_data$padj, scientific = TRUE, digits = 3)
  summary_df[[regulation_col]] <- as.character(matched_data$change)
}

write.csv(summary_df, "All_Genes_DEG_Summary.csv", row.names = FALSE)

# 11. 生成统计信息报告 ---------------------------------------------------
cat("\n=== 数据分析统计报告 ===\n")
cat("总样本数:", ncol(dds), "\n")
cat("过滤后基因数:", nrow(dds), "\n")

# 统计每个时间点的差异基因数量
cat("\n各时间点差异基因统计:\n")
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  up_count <- sum(res_df$change == "Up-regulated", na.rm = TRUE)
  down_count <- sum(res_df$change == "Down-regulated", na.rm = TRUE)
  total_count <- up_count + down_count
  
  cat(tp, "vs Mock: 总DEGs =", total_count, 
      "(上调 =", up_count, "下调 =", down_count, ")\n")
}

# 12. 生成示例热图 -------------------------------------------------------
# 12.1 提取前50个最显著差异基因（跨所有时间点）
# 计算每个基因在不同时间点出现的最小padj
all_sig_genes <- unique(unlist(lapply(deg_list, function(x) {
  x %>% filter(change != "Not Significant") %>% pull(gene)
})))

if(length(all_sig_genes) > 0) {
  gene_min_padj <- sapply(all_sig_genes, function(gene) {
    min(sapply(deg_list, function(df) {
      idx <- which(df$gene == gene)
      if(length(idx) > 0) df$padj[idx] else 1
    }), na.rm = TRUE)
  })
  
  top_genes <- names(sort(gene_min_padj))[1:min(50, length(gene_min_padj))]
  
  # 提取表达矩阵并标准化
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  
  # 创建样本注释
  annotation_col <- data.frame(
    Timepoint = sample_info$timepoint,
    Condition = sample_info$group_label
  )
  rownames(annotation_col) <- colnames(mat)
  
  # 定义颜色
  color_mock <- "#708090"    # 灰色
  color_mpxv <- "#A52A2A"    # 深棕色
  
  # 创建颜色映射
  annotation_colors <- list(
    Timepoint = c("mock" = color_mock, "3d" = "#1B9E77", "5d" = "#D95F02", 
                  "7d" = "#7570B3", "10d" = "#E7298A"),
    Condition = c("MOCK" = color_mock, "MPXV_lb" = color_mpxv)
  )
  
  # 绘制热图
  pheatmap(mat,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 9,
           angle_col = 45,
           main = "Top 50 显著差异基因热图",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_rows = TRUE,
           cluster_cols = TRUE)
  
  # 保存热图
  png("Top50_DEGs_Heatmap.png", width = 3200, height = 2400, res = 300)
  pheatmap(mat,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 9,
           angle_col = 45,
           main = "Top 50 显著差异基因热图",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_rows = TRUE,
           cluster_cols = TRUE)
  dev.off()
  
  print("热图已保存为: Top50_DEGs_Heatmap.png")
}

cat("\n=== 文件输出总结 ===\n")
cat("1. Excel文件:", excel_filename, "\n")
cat("2. CSV文件:\n")
cat("   - All_Genes_VST_Expression_Matrix.csv\n")
cat("   - All_Genes_Zscore_Matrix.csv\n")
cat("   - All_Genes_DEG_Summary.csv\n")
cat("3. 热图文件: Top50_DEGs_Heatmap.png\n")
cat("\n分析完成！\n")