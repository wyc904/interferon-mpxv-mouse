# 导出Top30差异基因的Excel文件 -------------------------------------------------

  library(writexl)

# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt", 
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
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行DESeq2分析 ----------------------------------------------------
dds <- DESeq(dds)

# 6. 按时间点进行差异分析 ----------------------------------------------
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
  print(paste(tp, "vs mock: 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# 7. 可视化 --------------------------------------------------------

# 全局配色设置
color_mock <- "#708090"    # 灰色
color_mpxv <- "#A52A2A"    # 深棕色
color_up <- "#C75D69"       # 浅红
color_down <- "#4682B4"     # 蓝色
color_ns <- "#D3D3D3"       # 浅灰




# 7.3 合并所有时间点的热图
# 提取所有时间点的显著差异基因
all_sig_genes <- unique(unlist(lapply(deg_list, function(x) {
  x %>% filter(change != "Not Significant") %>% pull(gene)
})))

if(length(all_sig_genes) > 10) {
  # 取前30个最显著的基因（跨所有时间点）
  # 计算每个基因在不同时间点出现的最小p值
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
    Timepoint = sample_info$timepoint,
    Condition = sample_info$group_label
  )
  rownames(annotation_col) <- colnames(mat)
  
  # 创建颜色映射
  annotation_colors <- list(
    Timepoint = c("mock" = color_mock, "3d" = "#1B9E77", "5d" = "#D95F02", 
                  "7d" = "#7570B3", "10d" = "#E7298A"),
    Condition = c("MOCK" = color_mock, "MPXV_lb" = color_mpxv)
  )
  
  
  








# 1. 提取Top30基因的完整信息
# 计算每个基因在不同时间点出现的最小padj
gene_min_padj <- sapply(all_sig_genes, function(gene) {
  min(sapply(deg_list, function(df) {
    idx <- which(df$gene == gene)
    if (length(idx) > 0) df$padj[idx] else 1
  }), na.rm = TRUE)
})

# 获取Top30基因
top_30_genes <- names(sort(gene_min_padj))[1:min(500, length(gene_min_padj))]

# 2. 创建包含所有时间点数据的Excel文件
excel_data <- list()

# 2.1 创建汇总表 - 包含所有时间点的log2FC和padj
summary_df <- data.frame(
  Gene = top_30_genes,
  stringsAsFactors = FALSE
)

# 为每个时间点添加log2FC和padj列
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  # 添加log2FoldChange
  log2fc_col <- paste0("log2FC_", tp)
  padj_col <- paste0("padj_", tp)
  regulation_col <- paste0("Regulation_", tp)
  
  # 匹配基因并提取数据
  matched_data <- res_df[match(top_30_genes, res_df$gene), ]
  
  summary_df[[log2fc_col]] <- round(matched_data$log2FoldChange, 3)
  summary_df[[padj_col]] <- format(matched_data$padj, scientific = TRUE, digits = 3)
  summary_df[[regulation_col]] <- as.character(matched_data$change)
}

# 添加到Excel数据列表
excel_data[["Summary_All_Timepoints"]] <- summary_df

# 2.2 为每个时间点创建详细表
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  # 筛选Top30基因的详细结果
  detailed_df <- res_df[res_df$gene %in% top_30_genes, ] %>%
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, change) %>%
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
    ) %>%
    arrange(AdjustedPValue)
  
  excel_data[[paste0("Detailed_", tp, "_vs_Mock")]] <- detailed_df
}

# 2.3 创建基因表达量矩阵表
# 提取Top30基因的标准化表达量
expression_matrix <- as.data.frame(assay(vsd)[top_30_genes, ]) %>%
  round(3) %>%
  rownames_to_column("Gene")

excel_data[["Expression_Matrix"]] <- expression_matrix

# 2.4 创建样本信息表
sample_info_export <- sample_info %>%
  dplyr::select(sample, timepoint, condition, group_label) %>%
  dplyr::rename(
    Sample = sample,
    Timepoint = timepoint,
    Condition = condition,
    Group = group_label
  )

excel_data[["Sample_Information"]] <- sample_info_export

# 3. 导出到Excel文件
excel_filename <- "Top30_DEGs_TimeSeries_Analysis.xlsx"

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

# 4. 额外导出CSV文件（可选）
# 4.1 导出Top30基因的汇总表
write.csv(summary_df, "Top30_DEGs_Summary.csv", row.names = FALSE)

# 4.2 导出热图中使用的标准化表达矩阵（Z-score）
heatmap_matrix <- t(scale(t(assay(vsd)[top_30_genes, ])))
write.csv(heatmap_matrix, "Top30_DEGs_Zscore_Matrix.csv")

print("CSV文件也已生成:")
print("- Top30_DEGs_Summary.csv")
print("- Top30_DEGs_Zscore_Matrix.csv")

# 5. 生成统计信息报告
cat("\n=== Top30差异基因统计报告 ===\n")
cat("总基因数量:", length(top_30_genes), "\n")
cat("基因列表:\n")
print(top_30_genes)

# 统计每个时间点的上下调基因数量
cat("\n各时间点Top30基因的调控情况:\n")
for(tp in timepoints) {
  tp_data <- summary_df %>% 
    dplyr::select(Gene, contains(paste0("Regulation_", tp)))
  colnames(tp_data)[2] <- "Regulation"
  
  up_count <- sum(tp_data$Regulation == "Up-regulated", na.rm = TRUE)
  down_count <- sum(tp_data$Regulation == "Down-regulated", na.rm = TRUE)
  ns_count <- sum(tp_data$Regulation == "Not Significant", na.rm = TRUE)
  
  cat(tp, "vs Mock: 上调=", up_count, "下调=", down_count, "不显著=", ns_count, "\n")
}

# 6. 创建基因注释信息（如果可用）
# 如果需要添加基因功能注释，可以在这里扩展
gene_annotation <- data.frame(
  Gene = top_30_genes,
  Description = "待补充基因功能描述",  # 可以手动添加或从数据库获取
  stringsAsFactors = FALSE
)

# 将注释信息添加到Excel（可选）
excel_data_with_annotation <- excel_data
excel_data_with_annotation[["Gene_Annotations"]] <- gene_annotation

# 导出带注释的Excel文件
write_xlsx(
  excel_data_with_annotation,
  path = "Top30_DEGs_TimeSeries_Analysis_With_Annotations.xlsx"
)

print(paste("带注释的Excel文件已保存为: Top30_DEGs_TimeSeries_Analysis_With_Annotations.xlsx"))