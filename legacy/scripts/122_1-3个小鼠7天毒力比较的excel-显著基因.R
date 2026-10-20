library(dplyr)
library(DESeq2)
library(tidyr)
library(openxlsx)  # 用于写入多个sheet的Excel文件

# 1. 读取数据 -------------------------------------------------
cat("正在读取数据...\n")
c57_counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt", 
                         header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
cat("C57数据维度:", dim(c57_counts), "\n")

agb6_counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/AGB6_counts_gene.txt", 
                          header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
cat("AGB6数据维度:", dim(agb6_counts), "\n")

stat1_counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt", 
                           header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
cat("STAT1数据维度:", dim(stat1_counts), "\n")

# 2. 处理列名 -------------------------------------------------
# 从第6列开始是真正的count数据
c57_count_matrix <- c57_counts[, 6:20]
agb6_count_matrix <- agb6_counts[, 6:ncol(agb6_counts)]
stat1_count_matrix <- stat1_counts[, 6:ncol(stat1_counts)]

# 简化列名
simplify_colnames <- function(colnames_vec) {
  colnames_vec <- gsub("STAR_aligned\\.trimmomatic\\.", "", colnames_vec)
  colnames_vec <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam", "", colnames_vec)
  colnames_vec <- gsub("_combined$", "", colnames_vec)
  colnames_vec <- gsub("\\.", "_", colnames_vec)
  colnames_vec <- gsub("_+", "_", colnames_vec)  # 移除连续的下划线
  colnames_vec <- gsub("^_|_$", "", colnames_vec)  # 移除开头和结尾的下划线
  return(colnames_vec)
}

colnames(c57_count_matrix) <- simplify_colnames(colnames(c57_count_matrix))
colnames(agb6_count_matrix) <- simplify_colnames(colnames(agb6_count_matrix))
colnames(stat1_count_matrix) <- simplify_colnames(colnames(stat1_count_matrix))

# 3. 手动定义样本信息 -------------------------------------------------
# C57 (Ccs-11) 样本
c57_control_samples <- c("c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P")
c57_infected_samples <- c("C57_7d_1_P", "C57_7d_2_P", "C57_7d_3_P")

# AGB6 (8BC) 样本
agb6_control_samples <- c(
  "JZ25070166_AGB6_MOCK_1_P_AGB6_MOCK_1_P",
  "JZ25070167_AGB6_MOCK_2_P_AGB6_MOCK_2_P", 
  "JZ25070168_AGB6_MOCK_3_P_AGB6_MOCK_3_P"
)
agb6_infected_samples <- c(
  "JZ25070169_AGB6_MPXV_lb_1_P_AGB6_MPXV_lb_1_P",
  "JZ25070170_AGB6_MPXV_lb_2_P_AGB6_MPXV_lb_2_P",
  "JZ25070171_AGB6_MPXV_lb_3_P_AGB6_MPXV_lb_3_P"
)

# STAT1 样本
stat1_control_samples <- c(
  "JZ25055847_5_1_P_5_1_P",
  "JZ25055848_5_2_P_5_2_P",
  "JZ25055849_5_3_P_5_3_P"
)
stat1_infected_samples <- c(
  "JZ25055844_4_1_P_4_1_P",
  "JZ25055845_4_2_P_4_2_P", 
  "JZ25055846_4_3_P_4_3_P"
)

# 4. 创建样本信息表 -------------------------------------------------
create_sample_info <- function(control_samples, infected_samples, strain_name, actual_samples) {
  # 检查样本是否存在
  valid_control <- control_samples[control_samples %in% actual_samples]
  valid_infected <- infected_samples[infected_samples %in% actual_samples]
  
  if(length(valid_control) == 0 && length(valid_infected) == 0) {
    stop(paste("没有找到", strain_name, "的有效样本"))
  }
  
  cat(strain_name, " - Control samples:", length(valid_control), "found:", paste(valid_control, collapse = ", "), "\n")
  cat(strain_name, " - Infected samples:", length(valid_infected), "found:", paste(valid_infected, collapse = ", "), "\n")
  
  samples <- c(valid_control, valid_infected)
  conditions <- c(rep("Control", length(valid_control)), rep("Infected", length(valid_infected)))
  strains <- rep(strain_name, length(samples))
  groups <- paste(strain_name, conditions, sep = "_")
  
  sample_info <- data.frame(
    sample = samples,
    strain = strains,
    condition = conditions,
    group = groups,
    stringsAsFactors = FALSE
  )
  
  return(sample_info)
}

# 创建样本信息
cat("\n=== 创建样本信息 ===\n")
c57_info <- create_sample_info(c57_control_samples, c57_infected_samples, "C57", colnames(c57_count_matrix))
agb6_info <- create_sample_info(agb6_control_samples, agb6_infected_samples, "8BC", colnames(agb6_count_matrix))
stat1_info <- create_sample_info(stat1_control_samples, stat1_infected_samples, "STAT1", colnames(stat1_count_matrix))

# 5. 标准化表达数据 -------------------------------------------------
cat("\n=== 标准化表达数据 ===\n")

# 标准化表达矩阵
normalize_expression <- function(count_matrix, sample_info) {
  # 只使用样本信息中的样本
  samples <- sample_info$sample
  valid_samples <- intersect(samples, colnames(count_matrix))
  
  if(length(valid_samples) < 2) {
    stop("没有足够的样本")
  }
  
  sample_info <- sample_info[sample_info$sample %in% valid_samples, ]
  expr_matrix <- as.matrix(count_matrix[, valid_samples])
  
  # 过滤低表达基因
  keep <- rowSums(expr_matrix >= 10) >= 2
  expr_matrix <- expr_matrix[keep, ]
  
  # 计算CPM (Counts Per Million)
  cpm_matrix <- t(t(expr_matrix) * 1e6 / colSums(expr_matrix))
  # 对数转换
  log_cpm_matrix <- log2(cpm_matrix + 1)
  
  return(log_cpm_matrix)
}

cat("标准化表达数据...\n")
c57_expr <- normalize_expression(c57_count_matrix, c57_info)
agb6_expr <- normalize_expression(agb6_count_matrix, agb6_info)
stat1_expr <- normalize_expression(stat1_count_matrix, stat1_info)

# 6. 计算各组的平均表达量 -------------------------------------------------
cat("\n=== 计算各组的平均表达量 ===\n")

# 计算每个品系Control和Infected组的平均表达量
calculate_group_mean <- function(expr_matrix, sample_info, strain_name, condition) {
  # 获取该品系和条件下的样本
  samples <- sample_info$sample[sample_info$strain == strain_name & sample_info$condition == condition]
  valid_samples <- intersect(samples, colnames(expr_matrix))
  
  if(length(valid_samples) == 0) {
    cat("警告: ", strain_name, "的", condition, "组没有找到样本\n")
    return(NULL)
  }
  
  # 提取这些样本的表达数据
  expr_subset <- expr_matrix[, valid_samples, drop = FALSE]
  
  # 计算平均表达量
  if(length(valid_samples) > 1) {
    mean_expr <- rowMeans(expr_subset, na.rm = TRUE)
  } else {
    mean_expr <- expr_subset[, 1]
  }
  
  cat(strain_name, condition, "组: 使用", length(valid_samples), "个样本计算平均表达量\n")
  return(mean_expr)
}

# 为每个品系和条件计算平均表达量
cat("计算C57 Control组平均表达量...\n")
c57_control_mean <- calculate_group_mean(c57_expr, c57_info, "C57", "Control")

cat("计算C57 Infected组平均表达量...\n")
c57_infected_mean <- calculate_group_mean(c57_expr, c57_info, "C57", "Infected")

cat("计算8BC Control组平均表达量...\n")
agb6_control_mean <- calculate_group_mean(agb6_expr, agb6_info, "8BC", "Control")

cat("计算8BC Infected组平均表达量...\n")
agb6_infected_mean <- calculate_group_mean(agb6_expr, agb6_info, "8BC", "Infected")

cat("计算STAT1 Control组平均表达量...\n")
stat1_control_mean <- calculate_group_mean(stat1_expr, stat1_info, "STAT1", "Control")

cat("计算STAT1 Infected组平均表达量...\n")
stat1_infected_mean <- calculate_group_mean(stat1_expr, stat1_info, "STAT1", "Infected")

# 7. 获取所有数据集中存在的基因 -------------------------------------------------
cat("\n=== 获取所有数据集中存在的基因 ===\n")

# 获取各数据集的基因名
c57_genes <- rownames(c57_expr)
agb6_genes <- rownames(agb6_expr)
stat1_genes <- rownames(stat1_expr)

# 找出在任意数据集中存在的基因
all_genes <- unique(c(c57_genes, agb6_genes, stat1_genes))
cat("在任意数据集中存在的总基因数:", length(all_genes), "\n")

# 8. 创建完整的平均表达量矩阵（包含所有基因） -------------------------------------------------
cat("\n=== 创建完整的平均表达量矩阵 ===\n")

# 创建空矩阵
mean_expr_matrix <- matrix(0, nrow = length(all_genes), ncol = 6)
rownames(mean_expr_matrix) <- all_genes
colnames(mean_expr_matrix) <- c("C57_Control", "8BC_Control", "STAT1_Control", 
                                "C57_Infected", "8BC_Infected", "STAT1_Infected")

# 填充C57 Control组
cat("填充C57 Control组数据...\n")
if(!is.null(c57_control_mean)) {
  for(gene in names(c57_control_mean)) {
    if(gene %in% rownames(mean_expr_matrix)) {
      mean_expr_matrix[gene, "C57_Control"] <- c57_control_mean[gene]
    }
  }
}

# 填充C57 Infected组
cat("填充C57 Infected组数据...\n")
if(!is.null(c57_infected_mean)) {
  for(gene in names(c57_infected_mean)) {
    if(gene %in% rownames(mean_expr_matrix)) {
      mean_expr_matrix[gene, "C57_Infected"] <- c57_infected_mean[gene]
    }
  }
}

# 填充8BC Control组
cat("填充8BC Control组数据...\n")
if(!is.null(agb6_control_mean)) {
  for(gene in names(agb6_control_mean)) {
    if(gene %in% rownames(mean_expr_matrix)) {
      mean_expr_matrix[gene, "8BC_Control"] <- agb6_control_mean[gene]
    }
  }
}

# 填充8BC Infected组
cat("填充8BC Infected组数据...\n")
if(!is.null(agb6_infected_mean)) {
  for(gene in names(agb6_infected_mean)) {
    if(gene %in% rownames(mean_expr_matrix)) {
      mean_expr_matrix[gene, "8BC_Infected"] <- agb6_infected_mean[gene]
    }
  }
}

# 填充STAT1 Control组
cat("填充STAT1 Control组数据...\n")
if(!is.null(stat1_control_mean)) {
  for(gene in names(stat1_control_mean)) {
    if(gene %in% rownames(mean_expr_matrix)) {
      mean_expr_matrix[gene, "STAT1_Control"] <- stat1_control_mean[gene]
    }
  }
}

# 填充STAT1 Infected组
cat("填充STAT1 Infected组数据...\n")
if(!is.null(stat1_infected_mean)) {
  for(gene in names(stat1_infected_mean)) {
    if(gene %in% rownames(mean_expr_matrix)) {
      mean_expr_matrix[gene, "STAT1_Infected"] <- stat1_infected_mean[gene]
    }
  }
}

# 转换为数据框
mean_expr_df <- as.data.frame(mean_expr_matrix)
mean_expr_df$Gene <- rownames(mean_expr_matrix)
mean_expr_df <- mean_expr_df[, c("Gene", "C57_Control", "8BC_Control", "STAT1_Control", 
                                 "C57_Infected", "8BC_Infected", "STAT1_Infected")]

# 9. 差异表达分析 -------------------------------------------------
cat("\n=== 开始差异表达分析 ===\n")

# 定义差异分析函数
run_deseq2_analysis <- function(count_matrix, sample_info, strain_name) {
  cat("\n=== ", strain_name, "品系差异分析 ===\n")
  
  # 确保样本顺序一致
  valid_samples <- intersect(colnames(count_matrix), sample_info$sample)
  if(length(valid_samples) < 2) {
    cat(strain_name, "没有足够的样本进行差异分析\n")
    return(NULL)
  }
  
  # 过滤样本
  count_subset <- count_matrix[, valid_samples]
  sample_info_subset <- sample_info[sample_info$sample %in% valid_samples, ]
  
  # 确保样本顺序一致
  sample_info_subset <- sample_info_subset[match(colnames(count_subset), sample_info_subset$sample), ]
  
  # 过滤低表达基因
  keep <- rowSums(count_subset >= 10) >= 2
  count_subset <- count_subset[keep, ]
  
  cat(strain_name, " - 过滤后基因数:", nrow(count_subset), "\n")
  cat(strain_name, " - 样本数:", ncol(count_subset), "\n")
  
  # 创建DESeq2对象
  dds <- DESeqDataSetFromMatrix(
    countData = count_subset,
    colData = sample_info_subset,
    design = ~ condition
  )
  
  # 设置参考水平
  dds$condition <- relevel(dds$condition, ref = "Control")
  
  # 运行DESeq2
  dds <- DESeq(dds)
  
  # 获取结果
  res <- results(dds, contrast = c("condition", "Infected", "Control"))
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  
  # 添加差异表达标签
  res_df$DE_Status <- "Not Significant"
  res_df$DE_Status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up-regulated"
  res_df$DE_Status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down-regulated"
  
  # 重新排列列
  res_df <- res_df[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "DE_Status")]
  
  # 排序：先按显著性，再按log2FoldChange绝对值
  res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange), na.last = TRUE), ]
  
  # 统计
  n_up <- sum(res_df$DE_Status == "Up-regulated", na.rm = TRUE)
  n_down <- sum(res_df$DE_Status == "Down-regulated", na.rm = TRUE)
  n_total <- n_up + n_down
  
  cat(strain_name, "差异基因统计:\n")
  cat("  上调基因:", n_up, "\n")
  cat("  下调基因:", n_down, "\n")
  cat("  总计:", n_total, "\n")
  
  return(list(
    results = res_df,
    dds = dds,
    stats = c(up = n_up, down = n_down, total = n_total)
  ))
}

# 对每个品系进行差异分析
cat("\n对C57品系进行差异分析...\n")
c57_deseq_result <- run_deseq2_analysis(c57_count_matrix, c57_info, "C57")

cat("\n对8BC品系进行差异分析...\n")
agb6_deseq_result <- run_deseq2_analysis(agb6_count_matrix, agb6_info, "8BC")

cat("\n对STAT1品系进行差异分析...\n")
stat1_deseq_result <- run_deseq2_analysis(stat1_count_matrix, stat1_info, "STAT1")

# 10. 提取显著差异基因 -------------------------------------------------
cat("\n=== 提取显著差异基因 ===\n")

# 提取每个品系的显著差异基因
extract_significant_genes <- function(deseq_result, strain_name) {
  if(is.null(deseq_result)) {
    cat(strain_name, "没有差异分析结果\n")
    return(list(
      up_genes = character(0),
      down_genes = character(0),
      all_sig_genes = character(0)
    ))
  }
  
  res_df <- deseq_result$results
  
  # 提取上调基因
  up_genes <- res_df$Gene[res_df$DE_Status == "Up-regulated"]
  up_genes <- up_genes[!is.na(up_genes)]
  
  # 提取下调基因
  down_genes <- res_df$Gene[res_df$DE_Status == "Down-regulated"]
  down_genes <- down_genes[!is.na(down_genes)]
  
  # 所有显著差异基因
  all_sig_genes <- c(up_genes, down_genes)
  
  cat(strain_name, "显著差异基因:\n")
  cat("  上调基因数:", length(up_genes), "\n")
  cat("  下调基因数:", length(down_genes), "\n")
  cat("  总计:", length(all_sig_genes), "\n")
  
  if(length(all_sig_genes) > 0) {
    cat("  前10个显著差异基因:", paste(head(all_sig_genes, 10), collapse = ", "), "\n")
  }
  
  return(list(
    up_genes = up_genes,
    down_genes = down_genes,
    all_sig_genes = all_sig_genes
  ))
}

# 提取各品系显著差异基因
c57_sig_genes <- extract_significant_genes(c57_deseq_result, "C57")
agb6_sig_genes <- extract_significant_genes(agb6_deseq_result, "8BC")
stat1_sig_genes <- extract_significant_genes(stat1_deseq_result, "STAT1")

# 11. 分析显著差异基因的交集和并集 -------------------------------------------------
cat("\n=== 分析显著差异基因的交集和并集 ===\n")

# 获取各品系显著差异基因
c57_all_sig <- c57_sig_genes$all_sig_genes
agb6_all_sig <- agb6_sig_genes$all_sig_genes
stat1_all_sig <- stat1_sig_genes$all_sig_genes

# 计算交集
intersect_all <- Reduce(intersect, list(c57_all_sig, agb6_all_sig, stat1_all_sig))
intersect_c57_8bc <- intersect(c57_all_sig, agb6_all_sig)
intersect_c57_stat1 <- intersect(c57_all_sig, stat1_all_sig)
intersect_8bc_stat1 <- intersect(agb6_all_sig, stat1_all_sig)

# 计算并集
union_all <- unique(c(c57_all_sig, agb6_all_sig, stat1_all_sig))

cat("显著差异基因统计:\n")
cat("  C57显著差异基因数:", length(c57_all_sig), "\n")
cat("  8BC显著差异基因数:", length(agb6_all_sig), "\n")
cat("  STAT1显著差异基因数:", length(stat1_all_sig), "\n")
cat("  三个品系共有显著差异基因数:", length(intersect_all), "\n")
cat("  C57和8BC共有显著差异基因数:", length(intersect_c57_8bc), "\n")
cat("  C57和STAT1共有显著差异基因数:", length(intersect_c57_stat1), "\n")
cat("  8BC和STAT1共有显著差异基因数:", length(intersect_8bc_stat1), "\n")
cat("  所有品系显著差异基因并集数:", length(union_all), "\n")

# 12. 筛选并保存显著差异基因的平均表达量 -------------------------------------------------
cat("\n=== 筛选并保存显著差异基因的平均表达量 ===\n")

# 筛选平均表达量矩阵中的显著差异基因
# 使用三个品系的并集
if(length(union_all) > 0) {
  # 筛选在平均表达量矩阵中存在的基因
  sig_genes_in_matrix <- union_all[union_all %in% mean_expr_df$Gene]
  cat("在平均表达量矩阵中找到的显著差异基因数:", length(sig_genes_in_matrix), "\n")
  
  if(length(sig_genes_in_matrix) > 0) {
    # 筛选平均表达量
    sig_mean_expr_df <- mean_expr_df[mean_expr_df$Gene %in% sig_genes_in_matrix, ]
    
    # 添加差异表达信息
    sig_mean_expr_df$C57_DE_Status <- "Not Significant"
    sig_mean_expr_df$C57_log2FC <- NA
    sig_mean_expr_df$C57_padj <- NA
    
    sig_mean_expr_df$`8BC_DE_Status` <- "Not Significant"
    sig_mean_expr_df$`8BC_log2FC` <- NA
    sig_mean_expr_df$`8BC_padj` <- NA
    
    sig_mean_expr_df$STAT1_DE_Status <- "Not Significant"
    sig_mean_expr_df$STAT1_log2FC <- NA
    sig_mean_expr_df$STAT1_padj <- NA
    
    # 填充C57差异信息
    if(!is.null(c57_deseq_result)) {
      c57_res <- c57_deseq_result$results
      for(i in 1:nrow(sig_mean_expr_df)) {
        gene <- sig_mean_expr_df$Gene[i]
        if(gene %in% c57_res$Gene) {
          sig_mean_expr_df$C57_DE_Status[i] <- c57_res$DE_Status[c57_res$Gene == gene]
          sig_mean_expr_df$C57_log2FC[i] <- c57_res$log2FoldChange[c57_res$Gene == gene]
          sig_mean_expr_df$C57_padj[i] <- c57_res$padj[c57_res$Gene == gene]
        }
      }
    }
    
    # 填充8BC差异信息
    if(!is.null(agb6_deseq_result)) {
      agb6_res <- agb6_deseq_result$results
      for(i in 1:nrow(sig_mean_expr_df)) {
        gene <- sig_mean_expr_df$Gene[i]
        if(gene %in% agb6_res$Gene) {
          sig_mean_expr_df$`8BC_DE_Status`[i] <- agb6_res$DE_Status[agb6_res$Gene == gene]
          sig_mean_expr_df$`8BC_log2FC`[i] <- agb6_res$log2FoldChange[agb6_res$Gene == gene]
          sig_mean_expr_df$`8BC_padj`[i] <- agb6_res$padj[agb6_res$Gene == gene]
        }
      }
    }
    
    # 填充STAT1差异信息
    if(!is.null(stat1_deseq_result)) {
      stat1_res <- stat1_deseq_result$results
      for(i in 1:nrow(sig_mean_expr_df)) {
        gene <- sig_mean_expr_df$Gene[i]
        if(gene %in% stat1_res$Gene) {
          sig_mean_expr_df$STAT1_DE_Status[i] <- stat1_res$DE_Status[stat1_res$Gene == gene]
          sig_mean_expr_df$STAT1_log2FC[i] <- stat1_res$log2FoldChange[stat1_res$Gene == gene]
          sig_mean_expr_df$STAT1_padj[i] <- stat1_res$padj[stat1_res$Gene == gene]
        }
      }
    }
    
    # 计算差异倍数
    sig_mean_expr_df$C57_FC <- 2^sig_mean_expr_df$C57_log2FC
    sig_mean_expr_df$`8BC_FC` <- 2^sig_mean_expr_df$`8BC_log2FC`
    sig_mean_expr_df$STAT1_FC <- 2^sig_mean_expr_df$STAT1_log2FC
    
    # 重新排列列
    sig_mean_expr_df <- sig_mean_expr_df %>%
      select(Gene, 
             C57_Control, C57_Infected, C57_FC, C57_log2FC, C57_padj, C57_DE_Status,
             `8BC_Control`, `8BC_Infected`, `8BC_FC`, `8BC_log2FC`, `8BC_padj`, `8BC_DE_Status`,
             STAT1_Control, STAT1_Infected, STAT1_FC, STAT1_log2FC, STAT1_padj, STAT1_DE_Status)
    
    # 按C57的log2FC绝对值排序
    sig_mean_expr_df <- sig_mean_expr_df[order(-abs(sig_mean_expr_df$C57_log2FC), na.last = TRUE), ]
    
    cat("显著差异基因平均表达量矩阵维度:", dim(sig_mean_expr_df), "\n")
  } else {
    cat("警告: 没有显著差异基因在平均表达量矩阵中找到\n")
    sig_mean_expr_df <- data.frame()
  }
} else {
  cat("警告: 没有找到任何显著差异基因\n")
  sig_mean_expr_df <- data.frame()
}

# 13. 保存结果 -------------------------------------------------
cat("\n=== 保存结果 ===\n")

# 创建一个Excel工作簿
wb <- createWorkbook()

# 保存所有基因的平均表达量
addWorksheet(wb, "All_Genes_Expression")
writeData(wb, sheet = "All_Genes_Expression", mean_expr_df)

# 保存显著差异基因的平均表达量
if(nrow(sig_mean_expr_df) > 0) {
  addWorksheet(wb, "Significant_Genes_Expression")
  writeData(wb, sheet = "Significant_Genes_Expression", sig_mean_expr_df)
}

# 保存C57差异分析结果
if(!is.null(c57_deseq_result)) {
  addWorksheet(wb, "C57_Differential_Expression")
  writeData(wb, sheet = "C57_Differential_Expression", c57_deseq_result$results)
}

# 保存8BC差异分析结果
if(!is.null(agb6_deseq_result)) {
  addWorksheet(wb, "8BC_Differential_Expression")
  writeData(wb, sheet = "8BC_Differential_Expression", agb6_deseq_result$results)
}

# 保存STAT1差异分析结果
if(!is.null(stat1_deseq_result)) {
  addWorksheet(wb, "STAT1_Differential_Expression")
  writeData(wb, sheet = "STAT1_Differential_Expression", stat1_deseq_result$results)
}

# 保存差异基因统计
de_stats <- data.frame(
  Strain = c("C57", "8BC", "STAT1"),
  Up_Regulated = c(
    if(!is.null(c57_deseq_result)) c57_deseq_result$stats["up"] else 0,
    if(!is.null(agb6_deseq_result)) agb6_deseq_result$stats["up"] else 0,
    if(!is.null(stat1_deseq_result)) stat1_deseq_result$stats["up"] else 0
  ),
  Down_Regulated = c(
    if(!is.null(c57_deseq_result)) c57_deseq_result$stats["down"] else 0,
    if(!is.null(agb6_deseq_result)) agb6_deseq_result$stats["down"] else 0,
    if(!is.null(stat1_deseq_result)) stat1_deseq_result$stats["down"] else 0
  ),
  Total = c(
    if(!is.null(c57_deseq_result)) c57_deseq_result$stats["total"] else 0,
    if(!is.null(agb6_deseq_result)) agb6_deseq_result$stats["total"] else 0,
    if(!is.null(stat1_deseq_result)) stat1_deseq_result$stats["total"] else 0
  )
)

addWorksheet(wb, "DE_Statistics")
writeData(wb, sheet = "DE_Statistics", de_stats)

# 保存交集统计
intersect_stats <- data.frame(
  Comparison = c("All_Three_Strains", "C57_8BC", "C57_STAT1", "8BC_STAT1", "Union_All"),
  Count = c(
    length(intersect_all),
    length(intersect_c57_8bc),
    length(intersect_c57_stat1),
    length(intersect_8bc_stat1),
    length(union_all)
  )
)

addWorksheet(wb, "Intersect_Statistics")
writeData(wb, sheet = "Intersect_Statistics", intersect_stats)

# 保存Excel文件
output_excel_file <- "Differential_Expression_Analysis_Results.xlsx"
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
cat("Excel结果文件已保存为:", output_excel_file, "\n")

# 也保存为单独的CSV文件
write.csv(mean_expr_df, "All_Genes_Average_Expression.csv", row.names = FALSE)
cat("所有基因平均表达量已保存为: All_Genes_Average_Expression.csv\n")

if(nrow(sig_mean_expr_df) > 0) {
  write.csv(sig_mean_expr_df, "Significant_Genes_Average_Expression.csv", row.names = FALSE)
  cat("显著差异基因平均表达量已保存为: Significant_Genes_Average_Expression.csv\n")
}

# 14. 创建汇总报告 -------------------------------------------------
cat("\n=== 分析完成 ===\n")
cat("总分析基因数:", nrow(mean_expr_df), "\n")

if(!is.null(c57_deseq_result)) {
  cat("\nC57差异分析结果:\n")
  cat("  上调基因:", c57_deseq_result$stats["up"], "\n")
  cat("  下调基因:", c57_deseq_result$stats["down"], "\n")
  cat("  总计:", c57_deseq_result$stats["total"], "\n")
}

if(!is.null(agb6_deseq_result)) {
  cat("\n8BC差异分析结果:\n")
  cat("  上调基因:", agb6_deseq_result$stats["up"], "\n")
  cat("  下调基因:", agb6_deseq_result$stats["down"], "\n")
  cat("  总计:", agb6_deseq_result$stats["total"], "\n")
}

if(!is.null(stat1_deseq_result)) {
  cat("\nSTAT1差异分析结果:\n")
  cat("  上调基因:", stat1_deseq_result$stats["up"], "\n")
  cat("  下调基因:", stat1_deseq_result$stats["down"], "\n")
  cat("  总计:", stat1_deseq_result$stats["total"], "\n")
}

cat("\n显著差异基因交集统计:\n")
cat("  三个品系共有显著差异基因:", length(intersect_all), "\n")
cat("  所有品系显著差异基因并集:", length(union_all), "\n")

cat("\n输出文件:\n")
cat("1. Differential_Expression_Analysis_Results.xlsx - 包含所有结果的Excel文件\n")
cat("2. All_Genes_Average_Expression.csv - 所有基因的平均表达量\n")
if(nrow(sig_mean_expr_df) > 0) {
  cat("3. Significant_Genes_Average_Expression.csv - 显著差异基因的平均表达量\n")
}