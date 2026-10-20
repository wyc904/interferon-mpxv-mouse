# 加载必要的包
library(tidyverse)
library(UpSetR)
library(openxlsx)
library(viridis)
library(RColorBrewer)

# 读取数据
data <- read.xlsx("绘制.xlsx", sheet = 1)

# 检查数据结构
cat("数据维度:", dim(data), "\n")
cat("列名:", colnames(data), "\n")
cat("time唯一值:", unique(data$time), "\n")
cat("mice唯一值:", unique(data$mice), "\n")

# 查看前几行
print(head(data))

# 1. 数据预处理 - 提取显著差异基因
# 定义显著性阈值
padj_threshold <- 0.05
log2fc_threshold <- 1

# 筛选显著差异基因
sig_data <- data %>%
  filter(!is.na(padj)) %>%
  mutate(
    direction = case_when(
      padj < padj_threshold & log2FoldChange > log2fc_threshold ~ "Up",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down",
      TRUE ~ "Not_sig"
    )
  ) %>%
  filter(direction != "Not_sig")  # 只保留显著差异基因

cat("\n显著差异基因总数:", nrow(sig_data), "\n")
cat("各条件显著基因数:\n")
print(table(sig_data$mice, sig_data$time, sig_data$direction))

# 2. 准备UpSet图数据
# 定义所有时间点
time_points <- c("3d", "5d", "7d", "10d")
mouse_strains <- c("c57", "8bc")

# 创建所有条件组合
conditions <- c()
for(mice_type in mouse_strains) {
  for(time_point in time_points) {
    conditions <- c(conditions, 
                    paste0(mice_type, "_", time_point, "_up"),
                    paste0(mice_type, "_", time_point, "_down"))
  }
}

cat("\n总条件数:", length(conditions), "\n")
cat("条件列表:", paste(conditions, collapse = ", "), "\n")

# 创建空的基因列表
gene_lists <- list()

# 填充各个条件的基因列表
for(cond in conditions) {
  parts <- str_split(cond, "_", simplify = TRUE)
  mice_type <- parts[1]
  time_point <- parts[2]
  direction_type <- parts[3]
  
  # 筛选符合条件的基因
  cond_genes <- sig_data %>%
    filter(mice == mice_type, 
           time == time_point, 
           direction == str_to_title(direction_type)) %>%
    pull(gene) %>%
    unique()
  
  gene_lists[[cond]] <- cond_genes
  cat(paste(cond, "基因数:", length(cond_genes), "\n"))
}

# 3. 创建二元矩阵用于UpSet图
# 获取所有显著差异基因
all_sig_genes <- unique(sig_data$gene)
cat("\n所有显著差异基因总数:", length(all_sig_genes), "\n")

# 创建二元矩阵
binary_matrix <- data.frame(Gene = all_sig_genes)

for(cond in conditions) {
  binary_matrix[[cond]] <- as.integer(binary_matrix$Gene %in% gene_lists[[cond]])
}

# 查看二元矩阵结构
cat("\n二元矩阵维度:", dim(binary_matrix), "\n")
cat("前5行前5列:\n")
print(binary_matrix[1:5, 1:6])

# 4. 使用UpSetR绘制基础UpSet图
# 4.1 上调基因的UpSet图
up_conditions <- conditions[str_detect(conditions, "_up$")]

cat("\n上调基因条件数:", length(up_conditions), "\n")
cat("上调基因条件:", paste(up_conditions, collapse = ", "), "\n")

# 准备上调基因的二元矩阵
up_binary <- binary_matrix[, c("Gene", up_conditions)]

# 计算每个基因所在的条件数
up_binary$n_sets <- rowSums(up_binary[, up_conditions])

# 只保留至少在1个条件中上调的基因
up_binary <- up_binary[up_binary$n_sets > 0, ]

cat("上调基因矩阵维度:", dim(up_binary), "\n")

# 绘制上调基因UpSet图
pdf("UpSet_Upregulated_Genes_4timepoints.pdf", width = 14, height = 10)
upset(
  up_binary[, up_conditions],  # 只包含二元列
  nsets = length(up_conditions),
  nintersects = 30,  # 显示前30个交集
  mb.ratio = c(0.6, 0.4),  # 矩阵和条形图的比例
  order.by = "freq",  # 按频率排序
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.0,
  line.size = 0.8,
  mainbar.y.label = "基因数量",
  sets.x.label = "各条件基因数",
  text.scale = c(1.2, 1.2, 1, 1, 1.4, 1),  # 调整文本大小
  set_size.show = TRUE,  # 显示集合大小
  set_size.scale_max = max(colSums(up_binary[, up_conditions])) * 1.1
)
dev.off()

# 4.2 下调基因的UpSet图
down_conditions <- conditions[str_detect(conditions, "_down$")]

cat("\n下调基因条件数:", length(down_conditions), "\n")
cat("下调基因条件:", paste(down_conditions, collapse = ", "), "\n")

# 准备下调基因的二元矩阵
down_binary <- binary_matrix[, c("Gene", down_conditions)]

# 计算每个基因所在的条件数
down_binary$n_sets <- rowSums(down_binary[, down_conditions])

# 只保留至少在1个条件中下调的基因
down_binary <- down_binary[down_binary$n_sets > 0, ]

cat("下调基因矩阵维度:", dim(down_binary), "\n")

# 绘制下调基因UpSet图
pdf("UpSet_Downregulated_Genes_4timepoints.pdf", width = 14, height = 10)
upset(
  down_binary[, down_conditions],
  nsets = length(down_conditions),
  nintersects = 30,
  mb.ratio = c(0.6, 0.4),
  order.by = "freq",
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.0,
  line.size = 0.8,
  mainbar.y.label = "基因数量",
  sets.x.label = "各条件基因数",
  text.scale = c(1.2, 1.2, 1, 1, 1.4, 1),
  set_size.show = TRUE,
  set_size.scale_max = max(colSums(down_binary[, down_conditions])) * 1.1
)
dev.off()

# 5. 分别按品系绘制UpSet图
# 5.1 C57品系
c57_conditions <- conditions[str_detect(conditions, "^c57")]
c57_binary <- binary_matrix[, c("Gene", c57_conditions)]
c57_binary <- c57_binary[rowSums(c57_binary[, c57_conditions]) > 0, ]

cat("\nC57条件数:", length(c57_conditions), "\n")
cat("C57基因矩阵维度:", dim(c57_binary), "\n")

pdf("UpSet_C57_4timepoints.pdf", width = 12, height = 8)
upset(
  c57_binary[, c57_conditions],
  nsets = length(c57_conditions),
  nintersects = 20,
  mb.ratio = c(0.6, 0.4),
  order.by = "freq",
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.0,
  line.size = 0.8,
  mainbar.y.label = "基因数量",
  sets.x.label = "C57各条件基因数",
  text.scale = c(1.2, 1.2, 1, 1, 1.4, 1),
  set_size.show = TRUE,
  set_size.scale_max = max(colSums(c57_binary[, c57_conditions])) * 1.1
)
dev.off()

# 5.2 8BC品系
bc_conditions <- conditions[str_detect(conditions, "^8bc")]
bc_binary <- binary_matrix[, c("Gene", bc_conditions)]
bc_binary <- bc_binary[rowSums(bc_binary[, bc_conditions]) > 0, ]

cat("\n8BC条件数:", length(bc_conditions), "\n")
cat("8BC基因矩阵维度:", dim(bc_binary), "\n")

pdf("UpSet_8BC_4timepoints.pdf", width = 12, height = 8)
upset(
  bc_binary[, bc_conditions],
  nsets = length(bc_conditions),
  nintersects = 20,
  mb.ratio = c(0.6, 0.4),
  order.by = "freq",
  decreasing = TRUE,
  number.angles = 30,
  point.size = 3.0,
  line.size = 0.8,
  mainbar.y.label = "基因数量",
  sets.x.label = "8BC各条件基因数",
  text.scale = c(1.2, 1.2, 1, 1, 1.4, 1),
  set_size.show = TRUE,
  set_size.scale_max = max(colSums(bc_binary[, bc_conditions])) * 1.1
)
dev.off()

# 6. 分别按时间点绘制UpSet图
for(time_point in time_points) {
  time_conditions <- conditions[str_detect(conditions, paste0("_", time_point, "_"))]
  time_binary <- binary_matrix[, c("Gene", time_conditions)]
  time_binary <- time_binary[rowSums(time_binary[, time_conditions]) > 0, ]
  
  cat(paste0("\n", time_point, "时间点条件数:", length(time_conditions), "\n"))
  cat(paste0(time_point, "基因矩阵维度:", dim(time_binary), "\n"))
  
  pdf(paste0("UpSet_", time_point, "_timepoint.pdf"), width = 10, height = 7)
  upset(
    time_binary[, time_conditions],
    nsets = length(time_conditions),
    nintersects = 15,
    mb.ratio = c(0.6, 0.4),
    order.by = "freq",
    decreasing = TRUE,
    number.angles = 30,
    point.size = 2.8,
    line.size = 0.7,
    mainbar.y.label = "基因数量",
    sets.x.label = paste0(time_point, "各条件基因数"),
    text.scale = c(1.2, 1.2, 1, 1, 1.3, 1),
    set_size.show = TRUE,
    set_size.scale_max = max(colSums(time_binary[, time_conditions])) * 1.1
  )
  dev.off()
}

# 7. 绘制品系内时间动态比较
# 7.1 C57品系内时间点比较
c57_time_comparison <- list()
for(direction in c("up", "down")) {
  for(mice_type in mouse_strains) {
    # 获取该品系该方向的所有条件
    strain_dir_conditions <- conditions[str_detect(conditions, paste0("^", mice_type, "_.*_", direction, "$"))]
    
    if(length(strain_dir_conditions) > 0) {
      # 创建基因列表
      strain_dir_genes <- list()
      for(cond in strain_dir_conditions) {
        strain_dir_genes[[cond]] <- gene_lists[[cond]]
      }
      
      # 创建二元矩阵
      all_genes <- unique(unlist(strain_dir_genes))
      strain_binary <- data.frame(Gene = all_genes)
      for(cond in strain_dir_conditions) {
        strain_binary[[cond]] <- as.integer(strain_binary$Gene %in% strain_dir_genes[[cond]])
      }
      
      # 绘制UpSet图
      pdf(paste0("UpSet_", mice_type, "_", direction, "_time_comparison.pdf"), width = 10, height = 7)
      upset(
        strain_binary[, -1],  # 移除Gene列
        nsets = length(strain_dir_conditions),
        nintersects = 10,
        mb.ratio = c(0.6, 0.4),
        order.by = "freq",
        decreasing = TRUE,
        number.angles = 30,
        point.size = 2.8,
        line.size = 0.7,
        mainbar.y.label = "基因数量",
        sets.x.label = paste0(str_to_upper(mice_type), " ", str_to_title(direction), " 基因数"),
        text.scale = c(1.2, 1.2, 1, 1, 1.3, 1),
        set_size.show = TRUE,
        set_size.scale_max = max(colSums(strain_binary[, -1])) * 1.1
      )
      dev.off()
    }
  }
}

# 8. 统计分析函数
# 计算各种交集的统计信息
calculate_intersection_stats <- function(gene_lists) {
  # 获取所有条件
  all_conditions <- names(gene_lists)
  
  # 初始化结果数据框
  stats_df <- data.frame(
    Intersection = character(),
    Genes = integer(),
    Gene_List = character()
  )
  
  # 计算所有可能的交集
  for(i in 1:length(all_conditions)) {
    # 单个条件的组合
    combs <- combn(all_conditions, i, simplify = FALSE)
    
    for(comb in combs) {
      # 计算交集
      if(length(comb) == 1) {
        intersect_genes <- gene_lists[[comb]]
      } else {
        intersect_genes <- Reduce(intersect, gene_lists[comb])
      }
      
      # 保存结果
      intersection_name <- paste(comb, collapse = " ∩ ")
      stats_df <- rbind(stats_df, data.frame(
        Intersection = intersection_name,
        Genes = length(intersect_genes),
        Gene_List = paste(intersect_genes, collapse = ";")
      ))
    }
  }
  
  return(stats_df)
}

# 计算上调基因的交集统计
up_stats <- calculate_intersection_stats(gene_lists[str_detect(names(gene_lists), "_up$")])
cat("\n上调基因交集统计（前10个）:\n")
print(head(up_stats[order(-up_stats$Genes), ], 10))

# 计算下调基因的交集统计
down_stats <- calculate_intersection_stats(gene_lists[str_detect(names(gene_lists), "_down$")])
cat("\n下调基因交集统计（前10个）:\n")
print(head(down_stats[order(-down_stats$Genes), ], 10))

# 9. 计算特定交集的基因
# 9.1 计算所有时间点都上调的基因（每个品系内）
for(mice_type in mouse_strains) {
  up_time_conditions <- conditions[str_detect(conditions, paste0("^", mice_type, "_.*_up$"))]
  if(length(up_time_conditions) > 0) {
    all_time_up <- Reduce(intersect, gene_lists[up_time_conditions])
    cat(paste0("\n", str_to_upper(mice_type), "所有时间点共同上调基因数: ", length(all_time_up), "\n"))
  }
}

# 9.2 计算品系间共同上调的基因（每个时间点）
for(time_point in time_points) {
  time_up_conditions <- conditions[str_detect(conditions, paste0("_.*_", time_point, "_up$"))]
  if(length(time_up_conditions) >= 2) {
    strain_common_up <- Reduce(intersect, gene_lists[time_up_conditions])
    cat(paste0(time_point, "时间点品系间共同上调基因数: ", length(strain_common_up), "\n"))
  }
}

# 10. 保存统计结果
write.csv(up_stats, "Upregulated_Intersection_Stats_4timepoints.csv", row.names = FALSE)
write.csv(down_stats, "Downregulated_Intersection_Stats_4timepoints.csv", row.names = FALSE)
write.csv(binary_matrix, "Binary_Matrix_All_Genes_4timepoints.csv", row.names = FALSE)

# 11. 创建总结报告
cat("\n", strrep("=", 70), "\n", sep = "")
cat("UpSet分析总结报告 (4个时间点)\n")
cat(strrep("=", 70), "\n", sep = "")
cat("分析条件: C57和8BC小鼠，3d、5d、7d、10d四个时间点\n")
cat("显著性阈值: padj <", padj_threshold, ", |log2FC| >", log2fc_threshold, "\n")
cat("总显著差异基因数:", nrow(sig_data), "\n")
cat("上调基因总数:", sum(sig_data$direction == "Up"), "\n")
cat("下调基因总数:", sum(sig_data$direction == "Down"), "\n")

cat("\n各条件基因数统计:\n")
# 按品系和时间点分组统计
for(mice_type in mouse_strains) {
  cat(paste0("\n", str_to_upper(mice_type), ":\n"))
  for(time_point in time_points) {
    up_genes <- gene_lists[[paste0(mice_type, "_", time_point, "_up")]]
    down_genes <- gene_lists[[paste0(mice_type, "_", time_point, "_down")]]
    cat(paste0("  ", time_point, ": 上调 ", length(up_genes), 
               ", 下调 ", length(down_genes), "\n"))
  }
}

cat("\n关键交集发现:\n")
# 所有时间点都上调的基因
for(mice_type in mouse_strains) {
  up_time_conditions <- conditions[str_detect(conditions, paste0("^", mice_type, "_.*_up$"))]
  if(length(up_time_conditions) > 0) {
    all_time_up <- Reduce(intersect, gene_lists[up_time_conditions])
    cat(paste0(str_to_upper(mice_type), "所有时间点共同上调基因: ", length(all_time_up), "\n"))
  }
}

# 品系间共同上调基因（每个时间点）
cat("\n品系间共同上调基因:\n")
for(time_point in time_points) {
  c57_up <- gene_lists[[paste0("c57_", time_point, "_up")]]
  bc_up <- gene_lists[[paste0("8bc_", time_point, "_up")]]
  common_up <- intersect(c57_up, bc_up)
  cat(paste0(time_point, ": C57和8BC共同上调基因: ", length(common_up), "\n"))
}

# 品系特异性基因
cat("\n品系特异性上调基因（在C57中上调但8BC中不上调）:\n")
for(time_point in time_points) {
  c57_up <- gene_lists[[paste0("c57_", time_point, "_up")]]
  bc_up <- gene_lists[[paste0("8bc_", time_point, "_up")]]
  c57_specific <- setdiff(c57_up, bc_up)
  cat(paste0(time_point, ": C57特有上调基因: ", length(c57_specific), "\n"))
}

cat("\n生成的文件:\n")
cat("1. UpSet_Upregulated_Genes_4timepoints.pdf - 所有上调基因UpSet图\n")
cat("2. UpSet_Downregulated_Genes_4timepoints.pdf - 所有下调基因UpSet图\n")
cat("3. UpSet_C57_4timepoints.pdf - C57品系所有条件UpSet图\n")
cat("4. UpSet_8BC_4timepoints.pdf - 8BC品系所有条件UpSet图\n")
cat("5. UpSet_3d_timepoint.pdf - 3d时间点UpSet图\n")
cat("6. UpSet_5d_timepoint.pdf - 5d时间点UpSet图\n")
cat("7. UpSet_7d_timepoint.pdf - 7d时间点UpSet图\n")
cat("8. UpSet_10d_timepoint.pdf - 10d时间点UpSet图\n")
cat("9. UpSet_c57_up_time_comparison.pdf - C57上调基因时间比较\n")
cat("10. UpSet_c57_down_time_comparison.pdf - C57下调基因时间比较\n")
cat("11. UpSet_8bc_up_time_comparison.pdf - 8BC上调基因时间比较\n")
cat("12. UpSet_8bc_down_time_comparison.pdf - 8BC下调基因时间比较\n")
cat("13. Upregulated_Intersection_Stats_4timepoints.csv - 上调基因交集统计\n")
cat("14. Downregulated_Intersection_Stats_4timepoints.csv - 下调基因交集统计\n")
cat("15. Binary_Matrix_All_Genes_4timepoints.csv - 二元矩阵数据\n")

cat(strrep("=", 70), "\n", sep = "")