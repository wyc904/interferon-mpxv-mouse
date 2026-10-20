############################################################
# 40 个免疫轨迹基因热图 - 适配DESeq2分析流程
# 使用DESeq2的vst标准化数据
############################################################

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(grid)

## 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop/Project/10-干扰素小鼠/1-源/8bc_counts_gene", 
                     header = TRUE, sep = "\t", 
                     comment.char = "#", row.names = 1)

cat("原始数据维度:", dim(counts), "\n")

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:18]

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("8bc_10d_1_P", 
                            "8bc_3d_1_P", "8bc_3d_2_P", "8bc_3d_3_P",
                            "8bc_5d_1_P", "8bc_5d_2_P", "8bc_5d_3_P",
                            "8bc_7d_1_P", "8bc_7d_2_P", "8bc_7d_3_P",
                            "8bc_mock_1_P", "8bc_mock_2_P", "8bc_mock_3_P")

print("重命名后的列名:")
print(colnames(count_matrix))

## 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  timepoint = factor(c("10d", rep("3d", 3), rep("5d", 3), 
                       rep("7d", 3), rep("mock", 3)),
                     levels = c("mock", "3d", "5d", "7d", "10d")),
  condition = factor(c(rep("MPXV_lb", 10), rep("mock", 3)),
                     levels = c("mock", "MPXV_lb")),
  group_label = c(rep("MPXV_lb", 10), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

print("样本信息:")
print(sample_info)

## 3. 创建DESeq2对象并处理数据 --------------------------------------------
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

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 运行DESeq2分析
dds <- DESeq(dds)

# 获取vst标准化数据
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)

## 4. 定义40个免疫轨迹基因 ------------------------------------------------
immune_genes <- c("Oas1a","Oas1d","Oas1g","Oas1h","Ifit1bl1","Ifit1bl2","Ifit2","Ifit3","Ifi207","Ifi213","Usp18","Slfn2","Stat1","Iigp1","Tgtp1","Gbp6","Cxcl9","Ccl22","Il27","Il15ra","Socs1","S100a3","Aif1","Chit1","Ms4a2","Serping1","Tap1","Tapbpl","Ccl21a","Il7r","Cd4","Cd22","Igkc","Iglc1","Iglc3","Iglj1","Igkv4-83","Pax5","Bach2","Mmp8","Mmp27","Adamts14","Col5a1","Col4a3","Lgals1","Pdgfb","Aldh1a1","Cep55","Mcm7","Cdt1")

# 检查哪些基因在数据中
available_genes <- immune_genes[immune_genes %in% rownames(vst_matrix)]
missing_genes <- setdiff(immune_genes, rownames(vst_matrix))

cat("找到的免疫轨迹基因数量:", length(available_genes), "\n")
if (length(missing_genes) > 0) {
  cat("缺失的基因:\n")
  print(missing_genes)
}

if (length(available_genes) == 0) {
  stop("没有找到任何指定的免疫轨迹基因")
}

## 5. 提取并处理40个基因的表达数据 --------------------------------------
# 提取目标基因的vst标准化表达矩阵
expr_40 <- vst_matrix[available_genes, ]

# 对表达矩阵进行Z-score标准化（按行）
expr_40_scaled <- t(scale(t(expr_40)))

# 确保基因按照指定顺序排列
expr_40_ordered <- expr_40_scaled[available_genes, ]

## 6. 创建样本信息（适配横向热图布局）----------------------------------
# 生成新的样本标签
sample_info$new_label <- paste0(sample_info$timepoint, "-", 
                                gsub(".*_(\\d)_P", "\\1", sample_info$sample))

# 按时间点顺序对样本进行排序
sample_info_ordered <- sample_info %>%
  arrange(timepoint, sample)

# 重新排列表达矩阵的列
expr_40_final <- expr_40_ordered[, sample_info_ordered$sample]
colnames(expr_40_final) <- sample_info_ordered$new_label

# 转置矩阵，实现横向热图布局（样本在行，基因在列）
mat_transposed <- t(expr_40_final)

## 7. 创建注释信息 --------------------------------------------------------
# 行注释（样本注释）
annotation_row <- data.frame(
  Timepoint = factor(sample_info_ordered$timepoint, 
                     levels = c("mock", "3d", "5d", "7d", "10d")),
  Condition = factor(sample_info_ordered$group_label,
                     levels = c("MOCK", "MPXV_lb")),
  row.names = sample_info_ordered$new_label
)

# 创建颜色映射
timepoint_colors <- c("mock" = "#708090", "3d" = "#1B9E77", 
                      "5d" = "#D95F02", "7d" = "#7570B3", "10d" = "#E7298A")
condition_colors <- c("MOCK" = "#708090", "MPXV_lb" = "#A52A2A")

annotation_colors <- list(
  Timepoint = timepoint_colors,
  Condition = condition_colors
)

## 8. 绘制横向热图（不聚类基因，按指定顺序）-------------------------------
# 设置图形参数
color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 绘制热图
p_heat <- pheatmap(mat_transposed,
                   color = color_palette,
                   annotation_row = annotation_row,
                   annotation_colors = annotation_colors,
                   cluster_rows = FALSE,                 # 不聚类行（样本），保持时间顺序
                   cluster_cols = FALSE,                 # 不聚类列（基因），保持指定顺序
                   show_colnames = TRUE,                 # 显示基因名
                   show_rownames = TRUE,                 # 显示样本名
                   border_color = NA,
                   fontsize_row = 9,                    # 样本名字体大小
                   fontsize_col = 8,                    # 基因名字体大小
                   fontsize = 10,
                   main = "",
                   angle_col = 45,                      # 基因名倾斜45度
                   treeheight_col = 0,                  # 不显示基因聚类树
                   treeheight_row = 0)                  # 不显示样本聚类树





p_heat_vertical <- pheatmap(expr_40_final,
                            color = color_palette,
                            annotation_col = annotation_row,
                            annotation_colors = annotation_colors,
                            cluster_rows = FALSE,                 # 不聚类基因
                            cluster_cols = FALSE,                 # 不聚类样本
                            show_colnames = TRUE,                 # 显示样本名
                            show_rownames = TRUE,                 # 显示基因名
                            border_color = NA,
                            fontsize_row = 8,                     # 基因名字体大小
                            fontsize_col = 9,                     # 样本名字体大小
                            fontsize = 10,
                            main = "",
                            angle_col = 90,                      # 样本名倾斜45度
                            treeheight_col = 0,                   # 不显示样本聚类树
                            treeheight_row = 0)

p_heat_vertical 








## 9. 保存热图 -----------------------------------------------------------
# 保存PDF
# pdf("DESeq2_40Immune_Genes_Heatmap_Horizontal.pdf", width = 16, height = 10)
# grid::grid.newpage()
# grid::grid.draw(p_heat$gtable)
# dev.off()
# 
# # 保存PNG
# png("DESeq2_40Immune_Genes_Heatmap_Horizontal.png", width = 1600, height = 1000, res = 150)
# grid::grid.newpage()
# grid::grid.draw(p_heat$gtable)
# dev.off()

cat("✅ 横向热图已保存为：\n")
cat("   - DESeq2_40Immune_Genes_Heatmap_Horizontal.pdf\n")
cat("   - DESeq2_40Immune_Genes_Heatmap_Horizontal.png\n")

## 10. 保存基因列表和表达数据 --------------------------------------------
# 保存基因列表信息
gene_info <- data.frame(
  Gene = available_genes,
  Order = 1:length(available_genes)
)




# write.csv(gene_info, "DESeq2_40_Immune_Genes_List.csv", row.names = FALSE)
# 保存标准化后的表达矩阵
expr_output <- as.data.frame(expr_40_final) %>%
  tibble::rownames_to_column("Gene")
# write.csv(expr_output, "DESeq2_40_Immune_Genes_Expression.csv", row.names = FALSE)




cat("✅ 附加文件已保存：\n")
cat("   - DESeq2_40_Immune_Genes_List.csv\n")
cat("   - DESeq2_40_Immune_Genes_Expression.csv\n")

## 11. 统计信息输出 ------------------------------------------------------
cat("\n=== 热图生成统计信息 ===\n")
cat("总基因数量:", length(available_genes), "\n")
cat("样本数量:", ncol(expr_40_final), "\n")
cat("时间点分布:\n")
print(table(sample_info_ordered$timepoint))
cat("基因顺序（前10个）:\n")
print(head(available_genes, 10))
cat("缺失基因:", length(missing_genes), "\n")
if (length(missing_genes) > 0) {
  print(missing_genes)
}

## 12. 可选：绘制纵向热图（传统布局）------------------------------------
# 如果不转置矩阵，绘制传统的纵向热图
p_heat_vertical <- pheatmap(expr_40_final,
                            color = color_palette,
                            annotation_col = annotation_row,
                            annotation_colors = annotation_colors,
                            cluster_rows = FALSE,                 # 不聚类基因
                            cluster_cols = FALSE,                 # 不聚类样本
                            show_colnames = TRUE,                 # 显示样本名
                            show_rownames = TRUE,                 # 显示基因名
                            border_color = NA,
                            fontsize_row = 8,                     # 基因名字体大小
                            fontsize_col = 9,                     # 样本名字体大小
                            fontsize = 10,
                            main = "",
                            angle_col = 90,                      # 样本名倾斜45度
                            treeheight_col = 0,                   # 不显示样本聚类树
                            treeheight_row = 0)                   # 不显示基因聚类树



# 保存纵向热图
# pdf("DESeq2_40Immune_Genes_Heatmap_Vertical.pdf", width = 12, height = 10)
grid::grid.newpage()
grid::grid.draw(p_heat_vertical$gtable)
dev.off()
cat("✅ 纵向热图已保存为：DESeq2_40Immune_Genes_Heatmap_Vertical.pdf\n")






## 13. 与DESeq2差异分析结果整合 ------------------------------------------
# 进行各时间点vs mock的差异分析，并提取40个基因的结果
timepoints <- c("3d", "5d", "7d", "10d")
immune_genes_results <- list()

for(tp in timepoints) {
  cat("分析时间点:", tp, "\n")
  
  # 选择当前时间点和mock的样本
  tp_samples <- sample_info$sample[sample_info$timepoint %in% c(tp, "mock")]
  dds_sub <- dds[, tp_samples]
  
  # 确保condition因子水平正确（mock作为参考）
  dds_sub$condition <- relevel(dds_sub$condition, ref = "mock")
  
  # 运行DESeq2
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "MPXV_lb", "mock"), alpha = 0.05)
  
  # 提取40个免疫基因的结果
  immune_res <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(gene %in% available_genes) %>%
    arrange(match(gene, available_genes))
  
  if(nrow(immune_res) > 0) {
    immune_genes_results[[tp]] <- immune_res
    write.csv(immune_res, 
              paste0("DESeq2_Immune_Genes_", tp, "_vs_mock.csv"), 
              row.names = FALSE)
  }
}

cat("✅ 各时间点免疫基因差异分析结果已保存\n")

## 14. 创建汇总表格 ------------------------------------------------------
# 创建40个基因在各时间点的表达变化汇总表
summary_table <- data.frame(Gene = available_genes)

for(tp in timepoints) {
  if(tp %in% names(immune_genes_results)) {
    tp_results <- immune_genes_results[[tp]]
    # 添加log2FoldChange和padj列
    summary_table[[paste0(tp, "_log2FC")]] <- 
      tp_results$log2FoldChange[match(summary_table$Gene, tp_results$gene)]
    summary_table[[paste0(tp, "_padj")]] <- 
      tp_results$padj[match(summary_table$Gene, tp_results$gene)]
  }
}

write.csv(summary_table, "DESeq2_40_Immune_Genes_Summary.csv", row.names = FALSE)
cat("✅ 基因表达变化汇总表已保存: DESeq2_40_Immune_Genes_Summary.csv\n")

cat("\n=== 分析完成 ===\n")