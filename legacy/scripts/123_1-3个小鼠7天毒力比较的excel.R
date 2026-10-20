library(dplyr)

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

# 辅助函数：填充平均表达量
fill_mean_values <- function(mean_vector, col_name, gene_names) {
  if(!is.null(mean_vector)) {
    for(gene in gene_names) {
      if(gene %in% names(mean_vector)) {
        mean_expr_matrix[gene, col_name] <<- mean_vector[gene]
      } else {
        # 尝试忽略大小写匹配
        matched_names <- names(mean_vector)[tolower(names(mean_vector)) == tolower(gene)]
        if(length(matched_names) > 0) {
          mean_expr_matrix[gene, col_name] <<- mean_vector[matched_names[1]]
        }
      }
    }
  }
}

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

# 检查矩阵
cat("平均表达量矩阵维度:", dim(mean_expr_matrix), "\n")
cat("前3行6列:\n")
print(mean_expr_matrix[1:3, ])

# 9. 保存为Excel文件 -------------------------------------------------
cat("\n=== 保存结果 ===\n")

# 转换为数据框
mean_expr_df <- as.data.frame(mean_expr_matrix)
mean_expr_df$Gene <- rownames(mean_expr_matrix)

# 重新排列列，将Gene列放在第一列
mean_expr_df <- mean_expr_df[, c("Gene", "C57_Control", "8BC_Control", "STAT1_Control", 
                                 "C57_Infected", "8BC_Infected", "STAT1_Infected")]

# 保存为CSV文件（Excel可以打开）
output_file <- "Average_Expression_All_Genes.csv"
write.csv(mean_expr_df, output_file, row.names = FALSE)
cat("平均表达量数据已保存为:", output_file, "\n")

# 保存为TSV文件
output_file_tsv <- "Average_Expression_All_Genes.tsv"
write.table(mean_expr_df, output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat("平均表达量数据已保存为:", output_file_tsv, "\n")

# 10. 创建汇总统计 -------------------------------------------------
cat("\n=== 创建汇总统计 ===\n")

# 计算各组的非零基因数
non_zero_counts <- apply(mean_expr_matrix, 2, function(x) sum(x > 0))

cat("\n各组的非零表达基因数:\n")
for(i in 1:length(non_zero_counts)) {
  cat(sprintf("%-15s: %d个基因\n", names(non_zero_counts)[i], non_zero_counts[i]))
}

# 计算各组的平均表达水平
mean_values <- apply(mean_expr_matrix, 2, function(x) mean(x[x > 0]))

cat("\n各组的平均表达水平（log2CPM，仅非零基因）:\n")
for(i in 1:length(mean_values)) {
  cat(sprintf("%-15s: %.2f\n", names(mean_values)[i], mean_values[i]))
}

# 创建汇总统计表
summary_df <- data.frame(
  Group = colnames(mean_expr_matrix),
  NonZero_Genes = non_zero_counts,
  Mean_Expression = round(mean_values, 2)
)

# 保存汇总统计
write.csv(summary_df, "Expression_Summary_Statistics.csv", row.names = FALSE)
cat("汇总统计已保存为: Expression_Summary_Statistics.csv\n")

cat("\n=== 分析完成 ===\n")
cat("总基因数:", nrow(mean_expr_matrix), "\n")
cat("输出文件:\n")
cat("1. Average_Expression_All_Genes.csv - 所有基因的平均表达量数据\n")
cat("2. Average_Expression_All_Genes.tsv - 所有基因的平均表达量数据（制表符分隔）\n")
cat("3. Expression_Summary_Statistics.csv - 汇总统计信息\n")