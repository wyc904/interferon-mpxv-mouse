library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(grid)
library(ggplot2)
library(tidyr)
library(reshape2)

# 设置全局配色方案
color_palette <- list(
  control = "#708090",    # 灰色
  infected = "#A52A2A",    # 深棕色
  up = "#C75D69",       # 浅红
  down = "#4682B4",     # 蓝色
  ns = "#D3D3D3",       # 浅灰
  strain_colors = c("Ccs-11" = "#1F77B4", "BALB/c" = "#FF7F0E", "STS" = "#2CA02C")
)

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

cat("\n处理后的样本名称:\n")
cat("C57 (前5个):", head(colnames(c57_count_matrix), 5), "\n")
cat("AGB6 (前5个):", head(colnames(agb6_count_matrix), 5), "\n")
cat("STAT1 (前5个):", head(colnames(stat1_count_matrix), 5), "\n")

# 3. 手动定义样本信息（基于图片） -------------------------------------------------
# 图片显示：3个品系，每个品系有Control和Infected两组
# 根据之前的输出，我们知道样本名称，需要手动指定

# Ccs-11 (C57) 样本
c57_control_samples <- c("c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P")
c57_infected_samples <- c("C57_7d_1_P", "C57_7d_2_P", "C57_7d_3_P")

# BALB/c (AGB6) 样本
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

# STS (STAT1) 样本
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
c57_info <- create_sample_info(c57_control_samples, c57_infected_samples, "Ccs-11", colnames(c57_count_matrix))
agb6_info <- create_sample_info(agb6_control_samples, agb6_infected_samples, "BALB/c", colnames(agb6_count_matrix))
stat1_info <- create_sample_info(stat1_control_samples, stat1_infected_samples, "STS", colnames(stat1_count_matrix))





# 5. 从图片中提取基因列表 -------------------------------------------------
genes_from_image <- c(
  "Il6", "Tnf", "Ifnb1", "Ifna1", "Cxcl1", "Cxcl2", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
  "Irf1", "Irf3", "Irf7", "Stat1", "Stat2", "Stat3", "Jak1", "Jak2", "Tyk2",
  "Tlr2", "Tlr3", "Tlr4", "Tlr7", "Tlr8", "Tlr9", "Myd88", "Traf6", "Nfkb1", "Rela",
  "Cd4", "Cd8a", "Cd8b1", "Cd3e", "Cd3g", "Cd3d", "Cd19", "Cd20", "Cd40", "Cd40lg",
  "Il2", "Il4", "Il10", "Il12a", "Il12b", "Il17a", "Il21", "Il23a", "Ifng", "Tgfb1",
  "B2m", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "H2-DMb2"
)

cat("\n从图片中提取的基因列表 (", length(genes_from_image), "个):\n")
cat(paste(genes_from_image, collapse = ", "), "\n")

# 6. 在数据集中查找这些基因 -------------------------------------------------
# 获取各数据集的基因名
c57_genes <- rownames(c57_count_matrix)
agb6_genes <- rownames(agb6_count_matrix)
stat1_genes <- rownames(stat1_count_matrix)

cat("\n=== 数据集基因数量 ===\n")
cat("Ccs-11 (C57) 数据集基因数:", length(c57_genes), "\n")
cat("BALB/c (AGB6) 数据集基因数:", length(agb6_genes), "\n")
cat("STS (STAT1) 数据集基因数:", length(stat1_genes), "\n")

# 查找基因
find_genes_in_dataset <- function(target_genes, dataset_genes, dataset_name) {
  found_genes <- c()
  for(gene in target_genes) {
    # 精确匹配
    if(gene %in% dataset_genes) {
      found_genes <- c(found_genes, gene)
    } else {
      # 尝试忽略大小写匹配
      gene_lower <- tolower(gene)
      dataset_genes_lower <- tolower(dataset_genes)
      if(gene_lower %in% dataset_genes_lower) {
        idx <- which(dataset_genes_lower == gene_lower)[1]
        found_genes <- c(found_genes, dataset_genes[idx])
      } else {
        # 尝试部分匹配
        matches <- grep(paste0("^", gene, "$"), dataset_genes, ignore.case = TRUE, value = TRUE)
        if(length(matches) > 0) {
          found_genes <- c(found_genes, matches[1])
        }
      }
    }
  }
  cat(dataset_name, "中找到图片基因:", length(found_genes), "/", length(target_genes), "个\n")
  return(unique(found_genes))
}

c57_found <- find_genes_in_dataset(genes_from_image, c57_genes, "Ccs-11")
agb6_found <- find_genes_in_dataset(genes_from_image, agb6_genes, "BALB/c")
stat1_found <- find_genes_in_dataset(genes_from_image, stat1_genes, "STS")

# 找出三个数据集共有的基因
common_genes <- Reduce(intersect, list(c57_found, agb6_found, stat1_found))
cat("\n三个数据集共有的图片基因:", length(common_genes), "个\n")
if(length(common_genes) > 0) {
  cat("共有基因:", paste(common_genes, collapse = ", "), "\n")
}

# 7. 如果没有共同基因，使用每个数据集中都存在的基因
if(length(common_genes) < 5) {
  cat("\n共同基因太少，使用在每个数据集中存在的图片基因...\n")
  all_present_genes <- c()
  for(gene in genes_from_image) {
    in_c57 <- gene %in% c57_genes || tolower(gene) %in% tolower(c57_genes)
    in_agb6 <- gene %in% agb6_genes || tolower(gene) %in% tolower(agb6_genes)
    in_stat1 <- gene %in% stat1_genes || tolower(gene) %in% tolower(stat1_genes)
    if(in_c57 && in_agb6 && in_stat1) {
      all_present_genes <- c(all_present_genes, gene)
    }
  }
  common_genes <- all_present_genes
  cat("在三个数据集中都存在的图片基因:", length(common_genes), "个\n")
}

# 8. 准备表达数据 -------------------------------------------------
cat("\n=== 准备表达数据 ===\n")

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

# 9. 创建热图函数 -------------------------------------------------
create_heatmap_for_image <- function(c57_expr, agb6_expr, stat1_expr, 
                                     c57_info, agb6_info, stat1_info, 
                                     gene_list, title) {
  # 检查数据是否可用
  if(is.null(c57_expr) || is.null(agb6_expr) || is.null(stat1_expr)) {
    stop("缺少表达矩阵")
  }
  
  if(length(gene_list) < 5) {
    stop("基因列表太短")
  }
  
  cat("\n准备创建热图...\n")
  cat("使用基因数量:", length(gene_list), "\n")
  cat("Ccs-11 样本数:", nrow(c57_info), "\n")
  cat("BALB/c 样本数:", nrow(agb6_info), "\n")
  cat("STS 样本数:", nrow(stat1_info), "\n")
  
  # 查找基因在数据集中的实际名称
  find_actual_gene_names <- function(gene_list, dataset_genes) {
    actual_genes <- c()
    for(gene in gene_list) {
      if(gene %in% dataset_genes) {
        actual_genes <- c(actual_genes, gene)
      } else {
        # 忽略大小写匹配
        matches <- grep(paste0("^", gene, "$"), dataset_genes, ignore.case = TRUE, value = TRUE)
        if(length(matches) > 0) {
          actual_genes <- c(actual_genes, matches[1])
        }
      }
    }
    return(unique(actual_genes))
  }
  
  
  
  
  # 在每个数据集中查找基因
  c57_genes_found <- find_actual_gene_names(gene_list, rownames(c57_expr))
  agb6_genes_found <- find_actual_gene_names(gene_list, rownames(agb6_expr))
  stat1_genes_found <- find_actual_gene_names(gene_list, rownames(stat1_expr))
  
  # 找出在三个数据集中都存在的基因
  all_genes <- Reduce(intersect, list(c57_genes_found, agb6_genes_found, stat1_genes_found))
  cat("在三个数据集中都找到的基因:", length(all_genes), "个\n")
  
  if(length(all_genes) < 5) {
    # 如果在三个数据集中都存在的基因太少，使用在任意数据集中存在的基因
    all_genes <- unique(c(c57_genes_found, agb6_genes_found, stat1_genes_found))
    cat("使用在任意数据集中存在的基因:", length(all_genes), "个\n")
  }
  
  cat("前20个基因:", paste(head(all_genes, 20), collapse = ", "), "\n")
  
  # 提取表达数据
  c57_expr_sub <- c57_expr[all_genes, c57_info$sample, drop = FALSE]
  agb6_expr_sub <- agb6_expr[all_genes, agb6_info$sample, drop = FALSE]
  stat1_expr_sub <- stat1_expr[all_genes, stat1_info$sample, drop = FALSE]
  
  # 合并表达矩阵
  combined_expr <- cbind(c57_expr_sub, agb6_expr_sub, stat1_expr_sub)
  
  # 对行进行标准化
  expr_scaled <- t(scale(t(combined_expr)))
  expr_scaled[!is.finite(expr_scaled)] <- 0  # 处理无限值
  
  # 创建样本注释
  sample_annotation <- data.frame(
    Strain = rep(c("Ccs-11", "BALB/c", "STS"), 
                 times = c(ncol(c57_expr_sub), ncol(agb6_expr_sub), ncol(stat1_expr_sub))),
    Condition = factor(c(
      as.character(c57_info$condition),
      as.character(agb6_info$condition),
      as.character(stat1_info$condition)
    ), levels = c("Control", "Infected"))
  )
  rownames(sample_annotation) <- colnames(combined_expr)
  
  # 设置注释颜色
  annotation_colors <- list(
    Strain = c("Ccs-11" = "#1F77B4", "BALB/c" = "#FF7F0E", "STS" = "#2CA02C"),
    Condition = c("Control" = "#708090", "Infected" = "#A52A2A")
  )
  
  # 绘制热图
  heatmap_obj <- pheatmap(expr_scaled,
                          color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                          scale = "row",
                          annotation_col = sample_annotation,
                          annotation_colors = annotation_colors,
                          cluster_rows = TRUE,
                          cluster_cols = FALSE,  # 不聚类列
                          show_rownames = TRUE,
                          show_colnames = TRUE,
                          fontsize_row = 8,
                          fontsize_col = 8,
                          fontsize = 10,
                          border_color = NA,
                          main = title,
                          gaps_col = c(ncol(c57_expr_sub), ncol(c57_expr_sub) + ncol(agb6_expr_sub)),
                          annotation_names_col = TRUE,
                          annotation_legend = TRUE)
  
  return(list(plot = heatmap_obj, data = expr_scaled, annotation = sample_annotation, genes = all_genes))
}

# 10. 创建热图 -------------------------------------------------
cat("\n=== 创建热图 ===\n")
if(length(common_genes) >= 5) {
  cat("使用共同基因创建热图...\n")
  heatmap_result <- create_heatmap_for_image(
    c57_expr = c57_expr,
    agb6_expr = agb6_expr,
    stat1_expr = stat1_expr,
    c57_info = c57_info,
    agb6_info = agb6_info,
    stat1_info = stat1_info,
    gene_list = common_genes,
    title = "Gene Expression Heatmap: MPXV Infected vs Control\n(Genes from provided image)"
  )
} else {
  cat("共同基因不足，使用图片中的所有基因...\n")
  heatmap_result <- create_heatmap_for_image(
    c57_expr = c57_expr,
    agb6_expr = agb6_expr,
    stat1_expr = stat1_expr,
    c57_info = c57_info,
    agb6_info = agb6_info,
    stat1_info = stat1_info,
    gene_list = genes_from_image,
    title = "Gene Expression Heatmap: MPXV Infected vs Control\n(All genes from provided image)"
  )
}

# 保存热图
if(!is.null(heatmap_result)) {
  pdf("MPXV_Heatmap_From_Image.pdf", width = 12, height = 14)
  print(heatmap_result$plot)
  dev.off()
  cat("热图已保存为 'MPXV_Heatmap_From_Image.pdf'\n")
  
  # 保存热图数据
  write.csv(heatmap_result$data, "Heatmap_Expression_Data.csv")
  write.csv(heatmap_result$annotation, "Heatmap_Sample_Annotation.csv")
  write.csv(heatmap_result$genes, "Heatmap_Genes.csv", row.names = FALSE)
  
  # 显示热图
  print(heatmap_result$plot)
} else {
  cat("无法创建热图\n")
}

# 11. 创建样本数量图 -------------------------------------------------
sample_summary <- data.frame(
  Strain = c("Ccs-11", "BALB/c", "STS"),
  Control_Samples = c(sum(c57_info$condition == "Control"),
                      sum(agb6_info$condition == "Control"),
                      sum(stat1_info$condition == "Control")),
  Infected_Samples = c(sum(c57_info$condition == "Infected"),
                       sum(agb6_info$condition == "Infected"),
                       sum(stat1_info$condition == "Infected"))
)

cat("\n=== 样本信息汇总 ===\n")
print(sample_summary)

# 转换为长格式
sample_summary_long <- sample_summary %>%
  pivot_longer(cols = c(Control_Samples, Infected_Samples), names_to = "Condition", values_to = "Count")





# 12. 创建每个品系单独的热图 -------------------------------------------------
cat("\n=== 创建各品系单独的热图 ===\n")

create_single_strain_heatmap <- function(expr_matrix, sample_info, strain_name, gene_list, n_genes = 30) {
  if(is.null(expr_matrix) || nrow(expr_matrix) == 0) {
    warning(paste("没有", strain_name, "的表达矩阵"))
    return(NULL)
  }
  
  # 过滤基因
  genes_in_data <- gene_list[gene_list %in% rownames(expr_matrix)]
  if(length(genes_in_data) < 5) {
    cat(strain_name, "中没有足够的图片基因，使用高表达基因...\n")
    # 使用平均表达最高的基因
    gene_means <- rowMeans(expr_matrix, na.rm = TRUE)
    top_genes <- names(sort(gene_means, decreasing = TRUE))[1:min(n_genes, length(gene_means))]
    genes_in_data <- top_genes
  }
  
  cat(strain_name, " - 用于热图的基因:", length(genes_in_data), "个\n")
  
  # 提取样本
  valid_samples <- intersect(colnames(expr_matrix), sample_info$sample)
  if(length(valid_samples) < 2) {
    warning(paste(strain_name, "中没有足够的样本"))
    return(NULL)
  }
  
  # 按条件排序样本
  sample_info_ordered <- sample_info[sample_info$sample %in% valid_samples, ]
  sample_info_ordered <- sample_info_ordered[order(sample_info_ordered$condition), ]
  valid_samples_ordered <- sample_info_ordered$sample
  
  # 提取数据
  expr_subset <- expr_matrix[genes_in_data, valid_samples_ordered, drop = FALSE]
  
  # 标准化
  expr_scaled <- t(scale(t(expr_subset)))
  expr_scaled[!is.finite(expr_scaled)] <- 0
  
  # 创建样本注释
  annotation_col <- data.frame(
    Condition = factor(sample_info_ordered$condition, levels = c("Control", "Infected"))
  )
  rownames(annotation_col) <- valid_samples_ordered
  
  # 修正后的代码
  ann_colors <- list(
    Condition = c("Control" = "#708090", "Infected" = "#A52A2A")
  )
    
    # 绘制热图
    p <- pheatmap(expr_scaled,
                  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                  scale = "row",
                  annotation_col = annotation_col,
                  annotation_colors = ann_colors,
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,  # 保持分组顺序
                  show_rownames = TRUE,
                  show_colnames = TRUE,
                  fontsize_row = 8,
                  fontsize_col = 8,
                  border_color = NA,
                  main = paste(strain_name, "Strain: MPXV Infected vs Control"),
                  gaps_col = sum(sample_info_ordered$condition == "Control"))
    
    return(p)
}

# 创建各品系单独的热图
c57_heatmap <- create_single_strain_heatmap(c57_expr, c57_info, "Ccs-11", genes_from_image)
if(!is.null(c57_heatmap)) {
  pdf("Ccs-11_Heatmap.pdf", width = 8, height = 10)
  print(c57_heatmap)
  dev.off()
  cat("Ccs-11 热图已保存为 'Ccs-11_Heatmap.pdf'\n")
}

agb6_heatmap <- create_single_strain_heatmap(agb6_expr, agb6_info, "BALB/c", genes_from_image)
if(!is.null(agb6_heatmap)) {
  pdf("BALB-c_Heatmap.pdf", width = 8, height = 10)
  print(agb6_heatmap)
  dev.off()
  cat("BALB/c 热图已保存为 'BALB-c_Heatmap.pdf'\n")
}

stat1_heatmap <- create_single_strain_heatmap(stat1_expr, stat1_info, "STS", genes_from_image)
if(!is.null(stat1_heatmap)) {
  pdf("STS_Heatmap.pdf", width = 8, height = 10)
  print(stat1_heatmap)
  dev.off()
  cat("STS 热图已保存为 'STS_Heatmap.pdf'\n")
}

# 13. 创建输出目录并保存结果 -------------------------------------------------
output_dir <- "MPXV_Image_Analysis_Results"
dir.create(output_dir, showWarnings = FALSE)

cat("\n=== 保存分析结果 ===\n")
cat("保存目录:", output_dir, "\n")
#


















