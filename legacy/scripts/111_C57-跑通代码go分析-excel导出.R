# 确保已安装并加载所需的库
library(clusterProfiler)
library(org.Mm.eg.db) # 假设您正在使用小鼠 (Mus musculus) 数据
library(ggplot2)
library(DOSE)
# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(openxlsx) # 新增：用于导出Excel文件

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

# ===============================================
# 1. 提取差异表达基因（DEGs）列表
# ===============================================
# 假设 `deg_list` 已从 DESeq2 结果生成，提取上调和下调基因
up_gene_list <- list()
down_gene_list <- list()
logFC_cutoff <- 1  # 定义 logFC 阈值

for(tp in timepoints) {
  res_df <- deg_list[[tp]]  # 从DESeq2结果中获取每个时间点的DEG数据
  
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

# 检查分离后的列表是否为空
if (length(up_gene_list) == 0) {
  stop("错误: 在任何时间点均未找到显著上调基因。请检查您的 logFC 阈值或 deg_list 数据。")
}
if (length(down_gene_list) == 0) {
  stop("错误: 在任何时间点均未找到显著下调基因。请检查您的 logFC 阈值或 deg_list 数据。")
}

# ===============================================
# 2. 运行 GO 富集比较分析 (上调基因)
# ===============================================
message("正在对上调基因进行 GO 富集比较分析...")
ck_up <- compareCluster(geneCluster = up_gene_list,
                        fun = "enrichGO",
                        OrgDb = org.Mm.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP", # 仅分析生物学过程 (BP)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1)

# ===============================================
# 3. 运行 GO 富集比较分析 (下调基因)
# ===============================================
message("正在对下调基因进行 GO 富集比较分析...")
ck_down <- compareCluster(geneCluster = down_gene_list,
                          fun = "enrichGO",
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1)

# ===============================================
# 4. 导出GO富集结果到Excel
# ===============================================

# 创建Excel工作簿
wb <- createWorkbook()

# 添加上调基因GO结果
if (!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
  up_go_df <- as.data.frame(ck_up)
  addWorksheet(wb, "上调基因_GO富集")
  writeData(wb, sheet = "上调基因_GO富集", up_go_df, rowNames = FALSE)
  message("上调基因GO富集结果已导出到Excel")
} else {
  message("注意: 上调基因的GO富集结果为空")
}

# 添加下调基因GO结果
if (!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
  down_go_df <- as.data.frame(ck_down)
  addWorksheet(wb, "下调基因_GO富集")
  writeData(wb, sheet = "下调基因_GO富集", down_go_df, rowNames = FALSE)
  message("下调基因GO富集结果已导出到Excel")
} else {
  message("注意: 下调基因的GO富集结果为空")
}

# 添加差异基因列表
addWorksheet(wb, "差异基因列表")
deg_summary <- data.frame()
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  if(nrow(res_df) > 0) {
    sig_genes <- res_df[res_df$change != "Not Significant", ]
    if(nrow(sig_genes) > 0) {
      sig_genes$Timepoint <- tp
      deg_summary <- rbind(deg_summary, sig_genes)
    }
  }
}
if(nrow(deg_summary) > 0) {
  writeData(wb, sheet = "差异基因列表", deg_summary, rowNames = FALSE)
  message("差异基因列表已导出到Excel")
}

# 保存Excel文件
excel_file <- "MPXV_GO_Analysis_Results.xlsx"
saveWorkbook(wb, excel_file, overwrite = TRUE)
message(paste("GO分析结果已保存到:", excel_file))

# ===============================================
# 5. 可视化：生成两张独立的点图
# ===============================================

# 5.1 上调基因的富集点图
if (!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
  message("生成上调基因 GO BP Dotplot...")
  p_up <- dotplot(ck_up, 
                  showCategory = 10, 
                  title = "上调基因 GO 生物学过程 (BP) 富集比较") + 
    scale_color_viridis_c(name = "调整P值") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_up)
  # 保存图片
  ggsave("上调基因_GO富集点图.pdf", p_up, width = 12, height = 8)
  ggsave("上调基因_GO富集点图.png", p_up, width = 12, height = 8, dpi = 300)
} else {
  message("注意: 上调基因的 GO 富集结果为空，跳过绘图。")
}

# 5.2 下调基因的富集点图
if (!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
  message("生成下调基因 GO BP Dotplot...")
  p_down <- dotplot(ck_down, 
                    showCategory = 10, 
                    title = "下调基因 GO 生物学过程 (BP) 富集比较") +
    scale_color_viridis_c(name = "调整P值") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_down)
  # 保存图片
  ggsave("下调基因_GO富集点图.pdf", p_down, width = 12, height = 8)
  ggsave("下调基因_GO富集点图.png", p_down, width = 12, height = 8, dpi = 300)
} else {
  message("注意: 下调基因的 GO 富集结果为空，跳过绘图。")
}

# ===============================================
# 6. 添加富集结果的详细分析表格
# ===============================================

# 创建详细分析工作簿
wb_detailed <- createWorkbook()

# 添加上调基因详细结果（按时间点分开）
if (!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
  up_go_df <- as.data.frame(ck_up)
  for(tp in timepoints) {
    sheet_name <- paste0("上调_", tp)
    tp_data <- up_go_df[up_go_df$Cluster == tp, ]
    if(nrow(tp_data) > 0) {
      addWorksheet(wb_detailed, sheet_name)
      writeData(wb_detailed, sheet = sheet_name, tp_data, rowNames = FALSE)
    }
  }
}

# 添加下调基因详细结果（按时间点分开）
if (!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
  down_go_df <- as.data.frame(ck_down)
  for(tp in timepoints) {
    sheet_name <- paste0("下调_", tp)
    tp_data <- down_go_df[down_go_df$Cluster == tp, ]
    if(nrow(tp_data) > 0) {
      addWorksheet(wb_detailed, sheet_name)
      writeData(wb_detailed, sheet = sheet_name, tp_data, rowNames = FALSE)
    }
  }
}

# 添加上调/下调基因列表
if(length(up_gene_list) > 0) {
  addWorksheet(wb_detailed, "上调基因列表")
  up_genes_df <- data.frame(
    Timepoint = rep(names(up_gene_list), sapply(up_gene_list, length)),
    Gene = unlist(up_gene_list)
  )
  writeData(wb_detailed, sheet = "上调基因列表", up_genes_df, rowNames = FALSE)
}

if(length(down_gene_list) > 0) {
  addWorksheet(wb_detailed, "下调基因列表")
  down_genes_df <- data.frame(
    Timepoint = rep(names(down_gene_list), sapply(down_gene_list, length)),
    Gene = unlist(down_gene_list)
  )
  writeData(wb_detailed, sheet = "下调基因列表", down_genes_df, rowNames = FALSE)
}

# 保存详细Excel文件
detailed_excel_file <- "MPXV_GO_Analysis_Detailed.xlsx"
saveWorkbook(wb_detailed, detailed_excel_file, overwrite = TRUE)
message(paste("详细GO分析结果已保存到:", detailed_excel_file))

message("\n分析完成！已生成以下文件：")
message("1. MPXV_GO_Analysis_Results.xlsx - 主要GO富集结果")
message("2. MPXV_GO_Analysis_Detailed.xlsx - 详细分时间点结果")
message("3. 上调基因_GO富集点图.pdf/png - 上调基因可视化")
message("4. 下调基因_GO富集点图.pdf/png - 下调基因可视化")