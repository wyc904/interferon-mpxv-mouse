# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)    # 用于富集分析
library(org.Mm.eg.db)       # 小鼠基因注释数据库
library(DOSE)               # 疾病本体分析
library(openxlsx)           # 新增：用于Excel输出
library(viridis)            # 新增：用于颜色方案

# 设置输出目录
output_dir <- "D:/Desktop/C57_MPXV_analysis_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}






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
  condition = factor(c(rep("MPXV_lb", 12), rep("mock", 3)),  # 修正：10d只有1个样本
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
timepoints <- c("3d", "5d", "7d")  # 移除10d，因为只有一个重复
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
      change = factor(change, levels = c("Down-regulated", "Not Significant", "Up-regulated")),
      timepoint = tp  # 添加时间点信息
    ) %>%
    arrange(padj)
  
  deg_list[[tp]] <- res_df
  print(paste(tp, "vs mock: 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# ===============================================
# 7. 保存差异表达分析结果到文件
# ===============================================
message("保存差异表达分析结果...")

# 创建Excel工作簿
wb_deg <- createWorkbook()

# 保存每个时间点的详细差异表达结果
for(tp in timepoints) {
  if(!is.null(deg_list[[tp]])) {
    sheet_name <- paste0("DEGs_", tp)
    addWorksheet(wb_deg, sheet_name)
    writeData(wb_deg, sheet = sheet_name, deg_list[[tp]], rowNames = FALSE)
  }
}

# 合并所有时间点的DEGs并保存
all_degs <- bind_rows(deg_list)
addWorksheet(wb_deg, "All_DEGs_Combined")
writeData(wb_deg, sheet = "All_DEGs_Combined", all_degs, rowNames = FALSE)

# 保存显著DEGs汇总
sig_degs_summary <- all_degs %>%
  filter(change != "Not Significant") %>%
  group_by(timepoint, change) %>%
  summarise(Count = n(), .groups = "drop")

addWorksheet(wb_deg, "Significant_DEGs_Summary")
writeData(wb_deg, sheet = "Significant_DEGs_Summary", sig_degs_summary, rowNames = FALSE)

# 保存DEGs Excel文件
deg_excel_file <- file.path(output_dir, "8bc_MPXV_DEGs_Analysis.xlsx")
saveWorkbook(wb_deg, deg_excel_file, overwrite = TRUE)
message(paste("差异表达结果已保存到:", deg_excel_file))

# 同时保存CSV格式
write.csv(all_degs, file.path(output_dir, "8bc_MPXV_All_DEGs.csv"), row.names = FALSE)
message(paste("所有DEGs已保存为CSV:", file.path(output_dir, "8bc_MPXV_All_DEGs.csv")))

# ===============================================
# 8. 提取差异表达基因（DEGs）列表
# ===============================================
up_gene_list <- list()
down_gene_list <- list()
logFC_cutoff <- 1  # 定义 logFC 阈值

for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  if(nrow(res_df) > 0) {
    # 提取上调基因 (log2FoldChange > 1)
    up_genes <- res_df$gene[res_df$log2FoldChange > logFC_cutoff & res_df$padj < P_cutoff]
    if(length(up_genes) > 0) {
      up_gene_list[[tp]] <- up_genes
    }
    
    # 提取下调基因 (log2FoldChange < -1)
    down_genes <- res_df$gene[res_df$log2FoldChange < -logFC_cutoff & res_df$padj < P_cutoff]
    if(length(down_genes) > 0) {
      down_gene_list[[tp]] <- down_genes
    }
  }
}

# 检查基因列表是否为空
if (length(up_gene_list) == 0) {
  warning("在任何时间点均未找到显著上调基因。")
}
if (length(down_gene_list) == 0) {
  warning("在任何时间点均未找到显著下调基因。")
}

# ===============================================
# 9. 运行 GO 富集分析
# ===============================================
message("进行GO富集分析...")

# 创建GO结果工作簿
wb_go <- createWorkbook()

# 上调基因GO分析
if (length(up_gene_list) > 0) {
  message("正在对上调基因进行GO富集分析...")
  ck_up <- compareCluster(geneCluster = up_gene_list,
                          fun = "enrichGO",
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1)
  
  if (!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
    up_go_df <- as.data.frame(ck_up)
    addWorksheet(wb_go, "上调基因_GO富集")
    writeData(wb_go, sheet = "上调基因_GO富集", up_go_df, rowNames = FALSE)
    message("上调基因GO富集分析完成")
  }
}

# 下调基因GO分析
if (length(down_gene_list) > 0) {
  message("正在对下调基因进行GO富集分析...")
  ck_down <- compareCluster(geneCluster = down_gene_list,
                            fun = "enrichGO",
                            OrgDb = org.Mm.eg.db,
                            keyType = "SYMBOL",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.1)
  
  if (!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
    down_go_df <- as.data.frame(ck_down)
    addWorksheet(wb_go, "下调基因_GO富集")
    writeData(wb_go, sheet = "下调基因_GO富集", down_go_df, rowNames = FALSE)
    message("下调基因GO富集分析完成")
  }
}

# 保存GO结果
if (length(wb_go$sheet_names) > 0) {
  go_excel_file <- file.path(output_dir, "8bc_MPXV_GO_Analysis.xlsx")
  saveWorkbook(wb_go, go_excel_file, overwrite = TRUE)
  message(paste("GO分析结果已保存到:", go_excel_file))
} else {
  message("GO富集分析未产生结果，可能没有显著富集的通路")
}

# ===============================================
# 10. 可视化结果
# ===============================================
message("生成可视化图表...")

# 10.1 上调基因GO点图
if (exists("ck_up") && !is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
  p_up <- dotplot(ck_up, 
                  showCategory = 10, 
                  title = "8bc MPXV感染 - 上调基因GO生物学过程富集") + 
    scale_color_viridis_c(name = "调整P值") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存图片
  ggsave(file.path(output_dir, "8bc_上调基因_GO富集点图.pdf"), p_up, width = 12, height = 8)
  ggsave(file.path(output_dir, "8bc_上调基因_GO富集点图.png"), p_up, width = 12, height = 8, dpi = 300)
  message("上调基因GO点图已保存")
}

# 10.2 下调基因GO点图
if (exists("ck_down") && !is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
  p_down <- dotplot(ck_down, 
                    showCategory = 10, 
                    title = "8bc MPXV感染 - 下调基因GO生物学过程富集") +
    scale_color_viridis_c(name = "调整P值") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存图片
  ggsave(file.path(output_dir, "8bc_下调基因_GO富集点图.pdf"), p_down, width = 12, height = 8)
  ggsave(file.path(output_dir, "8bc_下调基因_GO富集点图.png"), p_down, width = 12, height = 8, dpi = 300)
  message("下调基因GO点图已保存")
}

# 10.3 差异基因数量柱状图
deg_counts <- all_degs %>%
  filter(change != "Not Significant") %>%
  group_by(timepoint, change) %>%
  summarise(count = n(), .groups = "drop")

if(nrow(deg_counts) > 0) {
  p_deg_count <- ggplot(deg_counts, aes(x = timepoint, y = count, fill = change)) +
    geom_col(position = "dodge") +
    labs(title = "8bc MPXV感染 - 差异基因数量", 
         x = "时间点", y = "基因数量") +
    scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
    theme_minimal()
  
  ggsave(file.path(output_dir, "8bc_差异基因数量.pdf"), p_deg_count, width = 8, height = 6)
  ggsave(file.path(output_dir, "8bc_差异基因数量.png"), p_deg_count, width = 8, height = 6, dpi = 300)
}

# 10.4 火山图（示例：7天时间点）
if("7d" %in% timepoints && !is.null(deg_list[["7d"]])) {
  volvano_data <- deg_list[["7d"]] %>%
    mutate(significant = ifelse(change != "Not Significant", "Significant", "Not Significant"))
  
  p_volcano <- ggplot(volvano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = "8bc MPXV感染 7天 vs Mock - 火山图",
         x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "8bc_7天_火山图.pdf"), p_volcano, width = 8, height = 6)
  ggsave(file.path(output_dir, "8bc_7天_火山图.png"), p_volcano, width = 8, height = 6, dpi = 300)
}

# ===============================================
# 11. 保存基因列表
# ===============================================
message("保存基因列表...")

# 创建基因列表工作簿
wb_genes <- createWorkbook()

# 保存上调基因列表
if(length(up_gene_list) > 0) {
  up_genes_df <- data.frame(
    Timepoint = rep(names(up_gene_list), sapply(up_gene_list, length)),
    Gene = unlist(up_gene_list)
  )
  addWorksheet(wb_genes, "上调基因列表")
  writeData(wb_genes, sheet = "上调基因列表", up_genes_df, rowNames = FALSE)
}

# 保存下调基因列表
if(length(down_gene_list) > 0) {
  down_genes_df <- data.frame(
    Timepoint = rep(names(down_gene_list), sapply(down_gene_list, length)),
    Gene = unlist(down_gene_list)
  )
  addWorksheet(wb_genes, "下调基因列表")
  writeData(wb_genes, sheet = "下调基因列表", down_genes_df, rowNames = FALSE)
}

# 保存基因列表Excel
genes_excel_file <- file.path(output_dir, "8bc_MPXV_Gene_Lists.xlsx")
if(length(wb_genes$sheet_names) > 0) {
  saveWorkbook(wb_genes, genes_excel_file, overwrite = TRUE)
  message(paste("基因列表已保存到:", genes_excel_file))
}

# ===============================================
# 12. 保存分析报告
# ===============================================
message("生成分析报告...")

# 创建分析摘要
analysis_summary <- data.frame(
  项目 = c("分析样本", "总基因数", "过滤后基因数", "分析时间点", 
         "3天差异基因", "5天差异基因", "7天差异基因"),
  结果 = c(ncol(dds), nrow(counts), nrow(dds), 
         paste(timepoints, collapse = ", "),
         ifelse("3d" %in% timepoints, 
                paste(sum(deg_list[["3d"]]$change != "Not Significant"), "个"), "N/A"),
         ifelse("5d" %in% timepoints, 
                paste(sum(deg_list[["5d"]]$change != "Not Significant"), "个"), "N/A"),
         ifelse("7d" %in% timepoints, 
                paste(sum(deg_list[["7d"]]$change != "Not Significant"), "个"), "N/A"))
)

write.csv(analysis_summary, file.path(output_dir, "8bc_分析摘要.csv"), row.names = FALSE)

# ===============================================
# 13. 保存工作空间
# ===============================================
save.image(file.path(output_dir, "8bc_MPXV_analysis_workspace.RData"))
message(paste("工作空间已保存到:", file.path(output_dir, "8bc_MPXV_analysis_workspace.RData")))

# ===============================================
# 14. 完成信息
# ===============================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("8bc MPXV感染分析完成！\n")
cat("输出目录:", output_dir, "\n")
cat("生成的文件:\n")
cat("1. 8bc_MPXV_DEGs_Analysis.xlsx - 差异表达分析结果\n")
cat("2. 8bc_MPXV_GO_Analysis.xlsx - GO富集分析结果\n")
cat("3. 8bc_MPXV_Gene_Lists.xlsx - 上下调基因列表\n")
cat("4. 各种可视化图表(PDF/PNG格式)\n")
cat("5. 8bc_MPXV_analysis_workspace.RData - 工作空间文件\n")
cat(paste(rep("=", 50), collapse = ""), "\n")