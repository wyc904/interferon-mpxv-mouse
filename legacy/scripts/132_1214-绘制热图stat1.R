# =========================
# 0) 依赖包
# =========================
library(DESeq2)
library(tidyverse)
library(writexl)
library(openxlsx)  # 用于写入Excel

# =========================
# 1) 读取counts并整理样本名
# =========================
counts <- read.delim(
  "D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt",
  header = TRUE, sep = "\t", comment.char = "#", row.names = 1
)

# 提取count矩阵（从第6列开始是count数据：按你的文件结构）
count_matrix <- counts[, 6:20]

# 简化列名（按你原先规则）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix))

# =========================
# 2) 样本信息表（按你给的9个样本）
# =========================
sample_info <- data.frame(
  sample = c(
    "JZ25055841.3.1.P.3.1.P_combined",
    "JZ25055842.3.2.P.3.2.P_combined",
    "JZ25055843.3.3.P.3.3.P_combined",
    "JZ25055844.4.1.P.4.1.P_combined",
    "JZ25055845.4.2.P.4.2.P_combined",
    "JZ25055846.4.3.P.4.3.P_combined",
    "JZ25055847.5.1.P.5.1.P_combined",
    "JZ25055848.5.2.P.5.2.P_combined",
    "JZ25055849.5.3.P.5.3.P_combined"
  ),
  condition = c(
    "MPXV_Ilb","MPXV_Ilb","MPXV_Ilb",
    "MPXV_lb","MPXV_lb","MPXV_lb",
    "MOCK","MOCK","MOCK"
  ),
  replicate = c(1,2,3, 1,2,3, 1,2,3),
  stringsAsFactors = FALSE
)

# condition设为因子
sample_info$condition <- factor(sample_info$condition, levels = c("MOCK", "MPXV_lb", "MPXV_Ilb"))
rownames(sample_info) <- sample_info$sample

# 检查样本名是否在count_matrix里
missing_samples <- setdiff(sample_info$sample, colnames(count_matrix))
if (length(missing_samples) > 0) {
  stop("这些样本名在count_matrix列名中找不到，请检查命名是否一致：\n",
       paste(missing_samples, collapse = "\n"))
}

# 只保留这9个样本
count_matrix_filtered <- count_matrix[, sample_info$sample]
count_matrix_numeric <- as.matrix(count_matrix_filtered)
mode(count_matrix_numeric) <- "numeric"

# =========================
# 3) 构建DESeq2对象 + 过滤低表达
# =========================
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ condition
)

# 过滤低表达：至少3个样本count>=10
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

cat("过滤后保留基因数:", nrow(dds), "\n")

# =========================
# 4) 运行DESeq（只跑一次）
# =========================
dds <- DESeq(dds)

# 用于表达矩阵导出（vst更适合下游热图/聚类）
vsd <- vst(dds, blind = FALSE)

# =========================
# 5) 定义"组间比较"列表
# =========================
contrasts <- list(
  MPXV_lb_vs_MOCK     = c("condition", "MPXV_lb",  "MOCK"),
  MPXV_Ilb_vs_MOCK    = c("condition", "MPXV_Ilb", "MOCK"),
  MPXV_Ilb_vs_MPXV_lb = c("condition", "MPXV_Ilb", "MPXV_lb")
)

# 阈值
FC_cutoff <- 1.0   # abs(log2FC) > 1
P_cutoff  <- 0.05  # padj < 0.05
top_n     <- 30    # Top30 (用于单独的工作表)

# =========================
# 6) 对每个contrast做结果表
# =========================
res_tables <- list()        # 每个对比的完整结果df
sig_counts <- list()        # 每个对比的显著基因统计
top_gene_union <- c()       # 所有对比Top30基因的并集
all_genes <- rownames(dds)  # 所有基因

# 获取所有基因
all_genes <- rownames(dds)
cat("总基因数:", length(all_genes), "\n")

for (nm in names(contrasts)) {
  cat("分析对比:", nm, "\n")
  
  res <- results(dds, contrast = contrasts[[nm]], alpha = P_cutoff)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Gene") %>%
    mutate(
      Regulation = case_when(
        !is.na(padj) & padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up-regulated",
        !is.na(padj) & padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      Significance = ifelse(Regulation != "Not Significant", "Yes", "No")
    ) %>%
    arrange(padj)
  
  res_tables[[nm]] <- res_df
  
  # 统计显著基因
  sig_up <- sum(res_df$Regulation == "Up-regulated", na.rm = TRUE)
  sig_down <- sum(res_df$Regulation == "Down-regulated", na.rm = TRUE)
  sig_counts[[nm]] <- c(Up = sig_up, Down = sig_down, Total = sig_up + sig_down)
  
  # Top30：先筛显著，再按padj排序取前top_n
  top_df <- res_df %>%
    filter(Regulation != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = min(top_n, nrow(.)))
  
  top_gene_union <- union(top_gene_union, top_df$Gene)
  
  cat("  ", nm, "显著DEG数:", sig_up + sig_down, 
      "(上调:", sig_up, "下调:", sig_down, ")\n")
  cat("  Top30基因数:", length(top_df$Gene), "\n")
}

# 创建所有基因的汇总表
cat("\n创建所有基因的汇总表...\n")
all_genes_summary <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)

for (nm in names(contrasts)) {
  df <- res_tables[[nm]]
  m <- df[match(all_genes, df$Gene), ]
  
  all_genes_summary[[paste0("log2FC_", nm)]] <- round(m$log2FoldChange, 4)
  all_genes_summary[[paste0("padj_", nm)]] <- format(m$padj, scientific = TRUE, digits = 4)
  all_genes_summary[[paste0("Regulation_", nm)]] <- m$Regulation
  all_genes_summary[[paste0("BaseMean_", nm)]] <- round(m$baseMean, 2)
}

# 统计每个基因在多少个比较中显著
all_genes_summary$Significant_in_N_Comparisons <- rowSums(
  sapply(names(contrasts), function(nm) {
    ifelse(all_genes_summary[[paste0("Regulation_", nm)]] != "Not Significant", 1, 0)
  })
)

# =========================
# 7) 组织Excel输出
# =========================
cat("\n组织Excel输出...\n")
excel_data <- list()

# 7.1 统计摘要
stats_summary <- data.frame(
  Contrast = names(contrasts),
  Upregulated = sapply(names(contrasts), function(x) sig_counts[[x]]["Up"]),
  Downregulated = sapply(names(contrasts), function(x) sig_counts[[x]]["Down"]),
  Total_DEGs = sapply(names(contrasts), function(x) sig_counts[[x]]["Total"]),
  Percent_Total = round(100 * sapply(names(contrasts), function(x) sig_counts[[x]]["Total"]) / nrow(dds), 2)
)

excel_data[["Statistics_Summary"]] <- stats_summary

# 7.2 所有基因的汇总表
cat("  添加所有基因汇总表 (", nrow(all_genes_summary), " 行)...\n")
excel_data[["All_Genes_Summary"]] <- all_genes_summary

# 7.3 每个对比的完整结果表
for (nm in names(contrasts)) {
  cat("  添加对比", nm, "的完整结果表...\n")
  df <- res_tables[[nm]]
  
  detailed_df <- df %>%
    transmute(
      Gene,
      BaseMean = round(baseMean, 3),
      Log2FoldChange = round(log2FoldChange, 4),
      Log2FC_SE = round(lfcSE, 4),
      Statistic = round(stat, 4),
      PValue = pvalue,
      AdjustedPValue = padj,
      Regulation,
      Significance
    ) %>%
    arrange(AdjustedPValue)
  
  # 创建工作表名称（确保不超过31个字符）
  sheet_name <- gsub("MPXV_", "", nm)  # 移除MPXV_前缀
  sheet_name <- gsub("_vs_", "_vs_", sheet_name)
  sheet_name <- paste0("Full_", sheet_name)
  
  # 如果名称太长，截断
  if (nchar(sheet_name) > 31) {
    sheet_name <- substr(sheet_name, 1, 31)
  }
  
  excel_data[[sheet_name]] <- detailed_df
}

# 7.4 每个对比的显著基因表
for (nm in names(contrasts)) {
  cat("  添加对比", nm, "的显著基因表...\n")
  df <- res_tables[[nm]]
  
  sig_df <- df %>%
    filter(Regulation != "Not Significant") %>%
    transmute(
      Gene,
      BaseMean = round(baseMean, 3),
      Log2FoldChange = round(log2FoldChange, 4),
      Log2FC_SE = round(lfcSE, 4),
      Statistic = round(stat, 4),
      PValue = pvalue,
      AdjustedPValue = padj,
      Regulation
    ) %>%
    arrange(AdjustedPValue)
  
  # 创建工作表名称
  sheet_name <- gsub("MPXV_", "", nm)
  sheet_name <- paste0("Sig_", sheet_name)
  
  if (nchar(sheet_name) > 31) {
    sheet_name <- substr(sheet_name, 1, 31)
  }
  
  excel_data[[sheet_name]] <- sig_df
}

# 7.5 Top30基因的汇总表（与原始代码保持兼容）
if (length(top_gene_union) > 0) {
  cat("  添加Top30基因汇总表...\n")
  top_summary_df <- data.frame(Gene = top_gene_union, stringsAsFactors = FALSE)
  
  for (nm in names(contrasts)) {
    df <- res_tables[[nm]]
    m <- df[match(top_gene_union, df$Gene), ]
    
    top_summary_df[[paste0("log2FC_", nm)]] <- round(m$log2FoldChange, 4)
    top_summary_df[[paste0("padj_", nm)]] <- format(m$padj, scientific = TRUE, digits = 4)
    top_summary_df[[paste0("Regulation_", nm)]] <- m$Regulation
  }
  
  excel_data[["Top30_Genes_Summary"]] <- top_summary_df
}

# 7.6 表达矩阵
cat("  添加表达矩阵...\n")
expr_mat <- assay(vsd) %>%
  as.data.frame() %>%
  round(3) %>%
  rownames_to_column("Gene")

excel_data[["Expression_Matrix_VST"]] <- expr_mat

# 7.7 原始计数矩阵
cat("  添加原始计数矩阵...\n")
raw_counts <- counts(dds, normalized = FALSE) %>%
  as.data.frame() %>%
  rownames_to_column("Gene")

excel_data[["Raw_Count_Matrix"]] <- raw_counts

# 7.8 标准化计数矩阵
cat("  添加标准化计数矩阵...\n")
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  round(3) %>%
  rownames_to_column("Gene")

excel_data[["Normalized_Count_Matrix"]] <- norm_counts

# 7.9 样本信息
excel_data[["Sample_Information"]] <- sample_info %>%
  dplyr::rename(
    Sample = sample,
    Condition = condition,
    Replicate = replicate
  )

# 7.10 分析参数
params_df <- data.frame(
  Parameter = c("P-value cutoff", "Log2FC cutoff", "Low-expression filter", 
                "Total genes", "Filtered genes", "Comparisons"),
  Value = c(P_cutoff, FC_cutoff, "Count >= 10 in >= 3 samples", 
            nrow(counts(dds, normalized = FALSE)), nrow(dds), 
            paste(names(contrasts), collapse = "; "))
)

excel_data[["Analysis_Parameters"]] <- params_df

# =========================
# 8) 写入Excel
# =========================
excel_filename <- "DEGs_All_Genes_Analysis.xlsx"

cat("\n正在写入Excel文件:", excel_filename, "...\n")

# 使用openxlsx写入Excel，可以处理更多行
wb <- createWorkbook()

for (sheet_name in names(excel_data)) {
  # 确保工作表名称不超过31个字符
  safe_sheet_name <- ifelse(nchar(sheet_name) > 31, 
                            substr(sheet_name, 1, 31), 
                            sheet_name)
  
  addWorksheet(wb, safe_sheet_name)
  writeData(wb, sheet = safe_sheet_name, excel_data[[sheet_name]], rowNames = FALSE)
  
  # 为大数据表添加筛选器
  if (nrow(excel_data[[sheet_name]]) > 0 && ncol(excel_data[[sheet_name]]) > 0) {
    addFilter(wb, sheet = safe_sheet_name, row = 1, cols = 1:ncol(excel_data[[sheet_name]]))
  }
  
  # 设置列宽
  setColWidths(wb, sheet = safe_sheet_name, cols = 1:ncol(excel_data[[sheet_name]]), widths = "auto")
}

# 保存工作簿
saveWorkbook(wb, excel_filename, overwrite = TRUE)

# =========================
# 9) 打印统计摘要
# =========================
cat("\n✅ Excel文件已生成：", excel_filename, "\n")
cat("工作表数量：", length(excel_data), "\n")
cat("工作表名称：\n")
for (sheet in names(excel_data)) {
  cat("  - ", sheet, ":", nrow(excel_data[[sheet]]), "行 x", ncol(excel_data[[sheet]]), "列\n")
}

cat("\n=== 差异基因统计摘要 ===\n")
print(stats_summary)

cat("\n=== 各比较详细统计 ===\n")
for (nm in names(contrasts)) {
  cat(nm, ":\n")
  cat("  总基因数:", nrow(res_tables[[nm]]), "\n")
  cat("  显著DEGs:", sig_counts[[nm]]["Total"], 
      "(上调:", sig_counts[[nm]]["Up"], 
      "下调:", sig_counts[[nm]]["Down"], ")\n")
  cat("  百分比:", round(100 * sig_counts[[nm]]["Total"] / nrow(res_tables[[nm]]), 2), "%\n\n")
}

# 可选：保存为CSV文件便于其他软件处理
cat("\n✅ 同时生成CSV文件...\n")
write.csv(all_genes_summary, "All_Genes_Summary.csv", row.names = FALSE, na = "")

# 保存每个比较的显著基因
for (nm in names(contrasts)) {
  sig_df <- res_tables[[nm]] %>% filter(Regulation != "Not Significant")
  if (nrow(sig_df) > 0) {
    write.csv(sig_df, paste0("Significant_Genes_", gsub(" ", "_", nm), ".csv"), 
              row.names = FALSE, na = "")
  }
}

cat("✅ 所有文件生成完成！\n")