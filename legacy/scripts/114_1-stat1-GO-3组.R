



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





counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# print(dim(counts))
# counts2 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/AGB6_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# # print(dim(counts))
# counts3 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# # print(dim(counts))
# 
# 




colnames(counts)






# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 12:20]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix))

# 创建样本信息表
sample_info <- data.frame(
  sample = simplified_names,
  group = sapply(strsplit(simplified_names, "_"), function(x) paste(x[2:3], collapse = "_")),
  replicate = sapply(strsplit(simplified_names, "_"), function(x) x[4]),
  stringsAsFactors = FALSE
)

cat("\n=== 样本信息 ===\n")
print(sample_info)

# 1.1 构建样本信息表 (3组 vs 4组 vs 5组)
sample_info_raw <- data.frame(
  sample = c(
    # MPXV-Ilb 组 (Group 3)
    "JZ25055841.3.1.P.3.1.P_combined",
    "JZ25055842.3.2.P.3.2.P_combined", 
    "JZ25055843.3.3.P.3.3.P_combined",
    # MPXV-lb 组 (Group 4)
    "JZ25055844.4.1.P.4.1.P_combined",
    "JZ25055845.4.2.P.4.2.P_combined",
    "JZ25055846.4.3.P.4.3.P_combined",
    # Mock 组 (Group 5)
    "JZ25055847.5.1.P.5.1.P_combined",
    "JZ25055848.5.2.P.5.2.P_combined",
    "JZ25055849.5.3.P.5.3.P_combined"
  ),
  group_clear = c(
    "MPXV_Ilb", "MPXV_Ilb", "MPXV_Ilb",
    "MPXV_lb", "MPXV_lb", "MPXV_lb",
    "Mock", "Mock", "Mock"
  ),
  replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3)
)






# 创建新的样本信息表，匹配count_matrix中的样本名称
# 注意：这里的样本名称需要根据您实际的count_matrix列名进行调整
# 以下是假设的样本名称，您需要根据实际情况修改
sample_info <- data.frame(
  sample = colnames(count_matrix),  # 使用count_matrix的实际列名
  group = c(
    rep("MPXV_Ilb", 3),  # 假设前3个样本是MPXV_Ilb
    rep("MPXV_lb", 3),   # 中间3个样本是MPXV_lb
    rep("Mock", 3)       # 最后3个样本是Mock
  ),
  replicate = rep(1:3, 3),  # 每个组3个重复
  condition = c(
    rep("MPXV_Ilb", 3),
    rep("MPXV_lb", 3),
    rep("Mock", 3)
  ),
  stringsAsFactors = FALSE
)

# 设置行名为样本名
rownames(sample_info) <- sample_info$sample

# 查看样本信息
print("样本信息:")
print(sample_info)

# 3. 创建DESeq2对象 ----------------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 确保样本信息行名与count矩阵列名匹配
print("检查样本匹配:")
print(paste("count矩阵列数:", ncol(count_matrix_numeric)))
print(paste("样本信息行数:", nrow(sample_info)))
print("样本是否匹配:")
print(all(colnames(count_matrix_numeric) == rownames(sample_info)))

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ condition
)

print("DESeq2对象创建成功！")
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 设置参考水平为Mock
dds$condition <- relevel(dds$condition, ref = "Mock")

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行DESeq2分析 ----------------------------------------------------
dds <- DESeq(dds)

# 6. 进行差异分析（MPXV_lb vs Mock 和 MPXV_Ilb vs Mock）---------------
comparisons <- list(
  MPXV_lb_vs_Mock = c("condition", "MPXV_lb", "Mock"),
  MPXV_Ilb_vs_Mock = c("condition", "MPXV_Ilb", "Mock"),
  MPXV_lb_vs_MPXV_Ilb = c("condition", "MPXV_lb", "MPXV_Ilb")
)

results_list <- list()
deg_list <- list()

# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05

for(comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  print(paste("分析比较:", comp_name))
  
  # 提取结果
  res <- results(dds, contrast = comp, alpha = 0.05)
  
  # 保存完整结果
  results_list[[comp_name]] <- res
  
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
      comparison = comp_name  # 添加比较组信息
    ) %>%
    arrange(padj)
  
  deg_list[[comp_name]] <- res_df
  print(paste(comp_name, ": 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# ===============================================
# 7. 保存差异表达分析结果到文件
# ===============================================
message("保存差异表达分析结果...")

# 创建Excel工作簿
wb_deg <- createWorkbook()

# 保存每个比较的详细差异表达结果
for(comp_name in names(comparisons)) {
  if(!is.null(deg_list[[comp_name]])) {
    sheet_name <- gsub("_vs_", "_vs_", comp_name)
    addWorksheet(wb_deg, sheet_name)
    writeData(wb_deg, sheet = sheet_name, deg_list[[comp_name]], rowNames = FALSE)
  }
}

# 合并所有比较的DEGs并保存
all_degs <- bind_rows(deg_list)
addWorksheet(wb_deg, "All_DEGs_Combined")
writeData(wb_deg, sheet = "All_DEGs_Combined", all_degs, rowNames = FALSE)

# 保存显著DEGs汇总
sig_degs_summary <- all_degs %>%
  filter(change != "Not Significant") %>%
  group_by(comparison, change) %>%
  summarise(Count = n(), .groups = "drop")

addWorksheet(wb_deg, "Significant_DEGs_Summary")
writeData(wb_deg, sheet = "Significant_DEGs_Summary", sig_degs_summary, rowNames = FALSE)

# 保存DEGs Excel文件
deg_excel_file <- file.path(output_dir, "MPXV_DEGs_Analysis.xlsx")
saveWorkbook(wb_deg, deg_excel_file, overwrite = TRUE)
message(paste("差异表达结果已保存到:", deg_excel_file))

# 同时保存CSV格式
write.csv(all_degs, file.path(output_dir, "MPXV_All_DEGs.csv"), row.names = FALSE)
message(paste("所有DEGs已保存为CSV:", file.path(output_dir, "MPXV_All_DEGs.csv")))

# ===============================================
# 8. 提取差异表达基因（DEGs）列表
# ===============================================
# 分别提取每个比较的上下调基因
up_gene_list <- list()
down_gene_list <- list()
logFC_cutoff <- 1
P_cutoff <- 0.05

# 这里我们只对MPXV_lb_vs_Mock和MPXV_Ilb_vs_Mock进行GO分析
for(comp_name in c("MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock")) {
  if(comp_name %in% names(deg_list)) {
    res_df <- deg_list[[comp_name]]
    
    if(nrow(res_df) > 0) {
      # 提取上调基因
      up_genes <- res_df$gene[res_df$log2FoldChange > logFC_cutoff & 
                                res_df$padj < P_cutoff & 
                                !is.na(res_df$padj)]
      
      # 提取下调基因
      down_genes <- res_df$gene[res_df$log2FoldChange < -logFC_cutoff & 
                                  res_df$padj < P_cutoff & 
                                  !is.na(res_df$padj)]
      
      if(length(up_genes) > 0) {
        up_gene_list[[comp_name]] <- up_genes
        print(paste(comp_name, "上调基因数:", length(up_genes)))
      }
      
      if(length(down_genes) > 0) {
        down_gene_list[[comp_name]] <- down_genes
        print(paste(comp_name, "下调基因数:", length(down_genes)))
      }
    }
  }
}

# 检查基因列表是否为空
if (length(up_gene_list) == 0) {
  warning("未找到显著上调基因。")
} else {
  print(paste("总上调基因列表数:", length(up_gene_list)))
}

if (length(down_gene_list) == 0) {
  warning("未找到显著下调基因。")
} else {
  print(paste("总下调基因列表数:", length(down_gene_list)))
}





if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db")
}



# ===============================================
# 9. 运行 GO 富集分析
# ===============================================
# ===============================================
# 9. 运行 GO 富集分析
# ===============================================


# ===============================================
# 8. 提取差异表达基因（DEGs）列表
# ===============================================
# 分别提取每个比较的上下调基因
up_gene_list <- list()
down_gene_list <- list()
logFC_cutoff <- 1
P_cutoff <- 0.05

# 对三个比较都进行提取
for(comp_name in c("MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock", "MPXV_lb_vs_MPXV_Ilb")) {
  if(comp_name %in% names(deg_list)) {
    res_df <- deg_list[[comp_name]]
    
    if(nrow(res_df) > 0) {
      # 提取上调基因
      up_genes <- res_df$gene[res_df$log2FoldChange > logFC_cutoff & 
                                res_df$padj < P_cutoff & 
                                !is.na(res_df$padj)]
      
      # 提取下调基因
      down_genes <- res_df$gene[res_df$log2FoldChange < -logFC_cutoff & 
                                  res_df$padj < P_cutoff & 
                                  !is.na(res_df$padj)]
      
      if(length(up_genes) > 0) {
        up_gene_list[[comp_name]] <- up_genes
        print(paste(comp_name, "上调基因数:", length(up_genes)))
      }
      
      if(length(down_genes) > 0) {
        down_gene_list[[comp_name]] <- down_genes
        print(paste(comp_name, "下调基因数:", length(down_genes)))
      }
    }
  }
}

# 检查基因列表是否为空
if (length(up_gene_list) == 0) {
  warning("未找到显著上调基因。")
} else {
  print(paste("总上调基因列表数:", length(up_gene_list)))
  print("上调基因统计:")
  for(name in names(up_gene_list)) {
    print(paste(name, ":", length(up_gene_list[[name]]), "个基因"))
  }
}

if (length(down_gene_list) == 0) {
  warning("未找到显著下调基因。")
} else {
  print(paste("总下调基因列表数:", length(down_gene_list)))
  print("下调基因统计:")
  for(name in names(down_gene_list)) {
    print(paste(name, ":", length(down_gene_list[[name]]), "个基因"))
  }
}

# ===============================================
# 9. 运行 GO 富集分析
# ===============================================
message("进行GO富集分析...")

# 创建GO结果工作簿
wb_go <- createWorkbook()

# 检查gene list是否包含有效的基因
print("检查上调基因列表:")
for(name in names(up_gene_list)) {
  print(paste(name, ":", length(up_gene_list[[name]]), "个基因"))
}

print("检查下调基因列表:")
for(name in names(down_gene_list)) {
  print(paste(name, ":", length(down_gene_list[[name]]), "个基因"))
}

# 1) 上调基因GO分析
if (length(up_gene_list) > 0) {
  message("正在对上调基因进行GO富集分析...")
  
  # 检查是否有足够的基因进行GO分析
  total_up_genes <- sum(sapply(up_gene_list, length))
  if(total_up_genes < 5) {
    warning(paste("上调基因总数太少(", total_up_genes, ")，可能无法进行有效的GO富集分析"))
  } else {
    # 使用compareCluster进行富集分析
    tryCatch({
      ck_up <- compareCluster(geneCluster = up_gene_list,
                              fun = "enrichGO",
                              OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL",
                              ont = "BP",           # 生物学过程
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1,
                              pAdjustMethod = "BH",
                              minGSSize = 10,       # 最小基因集大小
                              maxGSSize = 500,      # 最大基因集大小
                              readable = TRUE)
      
      if (!is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
        up_go_df <- as.data.frame(ck_up)
        addWorksheet(wb_go, "上调基因_GO富集")
        writeData(wb_go, sheet = "上调基因_GO富集", up_go_df, rowNames = FALSE)
        
        # 保存上调基因的GO富集结果为单独的文件
        write.csv(up_go_df, 
                  file.path(output_dir, "上调基因_GO富集分析结果.csv"),
                  row.names = FALSE)
        
        message(paste("上调基因GO富集分析完成，找到", nrow(up_go_df), "个富集条目"))
        
        # 保存简化版本的GO结果
        up_go_simplified <- up_go_df %>%
          group_by(Cluster) %>%
          slice_head(n = 20) %>%  # 每个组取前20个富集结果
          arrange(pvalue)
        
        write.csv(up_go_simplified, 
                  file.path(output_dir, "上调基因_GO富集简化版.csv"),
                  row.names = FALSE)
        
      } else {
        warning("上调基因GO富集分析未产生结果")
      }
    }, error = function(e) {
      warning(paste("上调基因GO分析出错:", e$message))
    })
  }
}

# 2) 下调基因GO分析
if (length(down_gene_list) > 0) {
  message("正在对下调基因进行GO富集分析...")
  
  # 检查是否有足够的基因进行GO分析
  total_down_genes <- sum(sapply(down_gene_list, length))
  if(total_down_genes < 5) {
    warning(paste("下调基因总数太少(", total_down_genes, ")，可能无法进行有效的GO富集分析"))
  } else {
    tryCatch({
      ck_down <- compareCluster(geneCluster = down_gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Mm.eg.db,
                                keyType = "SYMBOL",
                                ont = "BP",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                pAdjustMethod = "BH",
                                minGSSize = 10,
                                maxGSSize = 500,
                                readable = TRUE)
      
      if (!is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
        down_go_df <- as.data.frame(ck_down)
        addWorksheet(wb_go, "下调基因_GO富集")
        writeData(wb_go, sheet = "下调基因_GO富集", down_go_df, rowNames = FALSE)
        
        # 保存下调基因的GO富集结果为单独的文件
        write.csv(down_go_df, 
                  file.path(output_dir, "下调基因_GO富集分析结果.csv"),
                  row.names = FALSE)
        
        message(paste("下调基因GO富集分析完成，找到", nrow(down_go_df), "个富集条目"))
        
        # 保存简化版本的GO结果
        down_go_simplified <- down_go_df %>%
          group_by(Cluster) %>%
          slice_head(n = 20) %>%  # 每个组取前20个富集结果
          arrange(pvalue)
        
        write.csv(down_go_simplified, 
                  file.path(output_dir, "下调基因_GO富集简化版.csv"),
                  row.names = FALSE)
        
      } else {
        warning("下调基因GO富集分析未产生结果")
      }
    }, error = function(e) {
      warning(paste("下调基因GO分析出错:", e$message))
    })
  }
}

# 保存GO结果到Excel
if (length(wb_go$sheet_names) > 0) {
  go_excel_file <- file.path(output_dir, "MPXV_GO_Analysis.xlsx")
  saveWorkbook(wb_go, go_excel_file, overwrite = TRUE)
  message(paste("GO分析结果已保存到:", go_excel_file))
} else {
  message("GO富集分析未产生结果，可能没有显著富集的通路")
}

# ===============================================
# 10. 运行KEGG富集分析
# ===============================================
message("进行KEGG富集分析...")

# 创建KEGG结果工作簿
wb_kegg <- createWorkbook()

# 1) 上调基因KEGG分析
if (length(up_gene_list) > 0) {
  message("正在对上调基因进行KEGG富集分析...")
  
  tryCatch({
    ck_up_kegg <- compareCluster(geneCluster = up_gene_list,
                                 fun = "enrichKEGG",
                                 organism = "mmu",  # 小鼠
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.1)
    
    if (!is.null(ck_up_kegg) && nrow(as.data.frame(ck_up_kegg)) > 0) {
      up_kegg_df <- as.data.frame(ck_up_kegg)
      addWorksheet(wb_kegg, "上调基因_KEGG富集")
      writeData(wb_kegg, sheet = "上调基因_KEGG富集", up_kegg_df, rowNames = FALSE)
      
      write.csv(up_kegg_df, 
                file.path(output_dir, "上调基因_KEGG富集分析结果.csv"),
                row.names = FALSE)
      
      message(paste("上调基因KEGG富集分析完成，找到", nrow(up_kegg_df), "个富集通路"))
    }
  }, error = function(e) {
    warning(paste("上调基因KEGG分析出错:", e$message))
  })
}

# 2) 下调基因KEGG分析
if (length(down_gene_list) > 0) {
  message("正在对下调基因进行KEGG富集分析...")
  
  tryCatch({
    ck_down_kegg <- compareCluster(geneCluster = down_gene_list,
                                   fun = "enrichKEGG",
                                   organism = "mmu",  # 小鼠
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1)
    
    if (!is.null(ck_down_kegg) && nrow(as.data.frame(ck_down_kegg)) > 0) {
      down_kegg_df <- as.data.frame(ck_down_kegg)
      addWorksheet(wb_kegg, "下调基因_KEGG富集")
      writeData(wb_kegg, sheet = "下调基因_KEGG富集", down_kegg_df, rowNames = FALSE)
      
      write.csv(down_kegg_df, 
                file.path(output_dir, "下调基因_KEGG富集分析结果.csv"),
                row.names = FALSE)
      
      message(paste("下调基因KEGG富集分析完成，找到", nrow(down_kegg_df), "个富集通路"))
    }
  }, error = function(e) {
    warning(paste("下调基因KEGG分析出错:", e$message))
  })
}

# 保存KEGG结果
if (length(wb_kegg$sheet_names) > 0) {
  kegg_excel_file <- file.path(output_dir, "MPXV_KEGG_Analysis.xlsx")
  saveWorkbook(wb_kegg, kegg_excel_file, overwrite = TRUE)
  message(paste("KEGG分析结果已保存到:", kegg_excel_file))
}











# ===============================================
# 10. 可视化结果
# ===============================================
message("生成可视化图表...")

# 10.1 上调基因GO点图
if (exists("ck_up") && !is.null(ck_up) && nrow(as.data.frame(ck_up)) > 0) {
  p_up <- dotplot(ck_up, 
                  showCategory = 10, 
                  title = "MPXV感染 - 上调基因GO生物学过程富集") + 
    scale_color_viridis_c(name = "调整P值") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存图片
  ggsave(file.path(output_dir, "上调基因_GO富集点图.pdf"), p_up, width = 12, height = 8)
  ggsave(file.path(output_dir, "上调基因_GO富集点图.png"), p_up, width = 12, height = 8, dpi = 300)
  message("上调基因GO点图已保存")
}

# 10.2 下调基因GO点图
if (exists("ck_down") && !is.null(ck_down) && nrow(as.data.frame(ck_down)) > 0) {
  p_down <- dotplot(ck_down, 
                    showCategory = 10, 
                    title = "MPXV感染 - 下调基因GO生物学过程富集") +
    scale_color_viridis_c(name = "调整P值") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存图片
  ggsave(file.path(output_dir, "下调基因_GO富集点图.pdf"), p_down, width = 12, height = 8)
  ggsave(file.path(output_dir, "下调基因_GO富集点图.png"), p_down, width = 12, height = 8, dpi = 300)
  message("下调基因GO点图已保存")
}

# 10.3 差异基因数量柱状图
if(exists("sig_degs_summary") && nrow(sig_degs_summary) > 0) {
  p_deg_count <- ggplot(sig_degs_summary, aes(x = comparison, y = Count, fill = change)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3) +
    labs(title = "MPXV感染 - 差异基因数量", 
         x = "比较组", y = "基因数量", fill = "表达变化") +
    scale_fill_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "差异基因数量.pdf"), p_deg_count, width = 8, height = 6)
  ggsave(file.path(output_dir, "差异基因数量.png"), p_deg_count, width = 8, height = 6, dpi = 300)
  message("差异基因数量柱状图已保存")
}

# 10.4 火山图
for(comp_name in c("MPXV_lb_vs_Mock", "MPXV_Ilb_vs_Mock")) {
  if(comp_name %in% names(deg_list)) {
    volcano_data <- deg_list[[comp_name]] %>%
      mutate(
        significant = ifelse(change != "Not Significant", as.character(change), "Not Significant"),
        significant = factor(significant, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
      )
    
    # 计算显著基因数量
    up_count <- sum(volcano_data$change == "Up-regulated")
    down_count <- sum(volcano_data$change == "Down-regulated")
    
    p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
      scale_color_manual(
        values = c("Up-regulated" = "#E41A1C", 
                   "Down-regulated" = "#377EB8", 
                   "Not Significant" = "grey70"),
        name = "表达变化"
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
      labs(
        title = paste0(gsub("_", " ", comp_name), " - 火山图"),
        x = "log2 Fold Change", 
        y = "-log10(adjusted p-value)",
        caption = paste0("上调: ", up_count, " | 下调: ", down_count)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
      )
    
    # 保存图片
    file_name <- gsub("_vs_", "_vs_", comp_name)
    ggsave(file.path(output_dir, paste0(file_name, "_火山图.pdf")), p_volcano, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(file_name, "_火山图.png")), p_volcano, width = 8, height = 6, dpi = 300)
    message(paste(comp_name, "火山图已保存"))
  }
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
    Comparison = rep(names(up_gene_list), sapply(up_gene_list, length)),
    Gene = unlist(up_gene_list)
  )
  addWorksheet(wb_genes, "上调基因列表")
  writeData(wb_genes, sheet = "上调基因列表", up_genes_df, rowNames = FALSE)
  
  # 保存为CSV
  write.csv(up_genes_df, 
            file.path(output_dir, "上调基因列表.csv"),
            row.names = FALSE)
}

# 保存下调基因列表
if(length(down_gene_list) > 0) {
  down_genes_df <- data.frame(
    Comparison = rep(names(down_gene_list), sapply(down_gene_list, length)),
    Gene = unlist(down_gene_list)
  )
  addWorksheet(wb_genes, "下调基因列表")
  writeData(wb_genes, sheet = "下调基因列表", down_genes_df, rowNames = FALSE)
  
  # 保存为CSV
  write.csv(down_genes_df, 
            file.path(output_dir, "下调基因列表.csv"),
            row.names = FALSE)
}

# 保存基因列表Excel
genes_excel_file <- file.path(output_dir, "MPXV_Gene_Lists.xlsx")
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
  项目 = c("分析样本数", "总基因数", "过滤后基因数", 
         "比较组数", "上调基因总数", "下调基因总数"),
  结果 = c(
    ncol(dds),
    nrow(count_matrix_numeric),
    nrow(dds),
    length(comparisons),
    ifelse(exists("up_genes_df"), nrow(up_genes_df), 0),
    ifelse(exists("down_genes_df"), nrow(down_genes_df), 0)
  )
)

write.csv(analysis_summary, file.path(output_dir, "分析摘要.csv"), row.names = FALSE)

# 保存详细的统计信息
deg_stats <- all_degs %>%
  filter(change != "Not Significant") %>%
  group_by(comparison, change) %>%
  summarise(基因数 = n(), .groups = "drop")

write.csv(deg_stats, file.path(output_dir, "差异基因统计.csv"), row.names = FALSE)

# ===============================================
# 13. 保存工作空间
# ===============================================
save.image(file.path(output_dir, "MPXV_analysis_workspace.RData"))
message(paste("工作空间已保存到:", file.path(output_dir, "MPXV_analysis_workspace.RData")))

# ===============================================
# 14. 完成信息
# ===============================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("MPXV感染差异表达分析完成！\n")
cat("输出目录:", output_dir, "\n")
cat("生成的文件:\n")
cat("1. MPXV_DEGs_Analysis.xlsx - 差异表达分析结果\n")
cat("2. MPXV_GO_Analysis.xlsx - GO富集分析结果\n")
cat("3. MPXV_Gene_Lists.xlsx - 上下调基因列表\n")
cat("4. 上调基因_GO富集分析结果.csv - 上调基因GO结果\n")
cat("5. 下调基因_GO富集分析结果.csv - 下调基因GO结果\n")
cat("6. 上调基因列表.csv - 上调基因列表\n")
cat("7. 下调基因列表.csv - 下调基因列表\n")
cat("8. 各种可视化图表(PDF/PNG格式)\n")
cat("9. MPXV_analysis_workspace.RData - 工作空间文件\n")
cat(paste(rep("=", 50), collapse = ""), "\n")