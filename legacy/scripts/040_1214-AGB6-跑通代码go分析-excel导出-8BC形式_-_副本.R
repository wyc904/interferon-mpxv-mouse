# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)    # 用于富集分析
library(org.Mm.eg.db)       # 小鼠基因注释数据库
library(openxlsx)           # 用于Excel输出
library(viridis)            # 用于颜色方案

# 设置输出目录
output_dir <- "D:/Desktop/AGB6_MPXV_analysis_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
cat("输出目录:", output_dir, "\n")

# 1. 读取并整理计数数据 -------------------------------------------------
cat("正在读取计数数据...\n")
counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/AGB6_counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 检查列名
cat("原始列名:\n")
print(colnames(counts))

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:19]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的列名:\n")
print(colnames(count_matrix))

# 2. 构建样本信息表 -------------------------------------------------
cat("构建样本信息表...\n")

# 定义我们要分析的6个样本
target_samples <- c(
  "JZ25070166.AGB6.MOCK.1.P.AGB6.MOCK.1.P_combined",              
  "JZ25070167.AGB6.MOCK.2.P.AGB6.MOCK.2.P_combined",              
  "JZ25070168.AGB6.MOCK.3.P.AGB6.MOCK.3.P_combined",              
  "JZ25070169.AGB6.MPXV_lb.1.P.AGB6.MPXV_lb.1.P_combined",        
  "JZ25070170.AGB6.MPXV_lb.2.P.AGB6.MPXV_lb.2.P_combined",        
  "JZ25070171.AGB6.MPXV_lb.3.P.AGB6.MPXV_lb.3.P_combined"
)

# 检查这些样本是否在count矩阵中
missing_samples <- setdiff(target_samples, colnames(count_matrix))
if (length(missing_samples) > 0) {
  stop("以下样本在count矩阵中找不到: ", paste(missing_samples, collapse = ", "))
}

# 提取这6个样本的数据
count_matrix_subset <- count_matrix[, target_samples]

# 创建样本信息
sample_info <- data.frame(
  sample = target_samples,
  condition = c(rep("AGB6.MOCK", 3), rep("AGB6.MPXV_lb", 3)),
  group_label = c(rep("Control", 3), rep("Infection", 3)),
  replicate = c(1, 2, 3, 1, 2, 3),
  stringsAsFactors = FALSE
)

# 将condition设为因子，确保Control为参考水平
sample_info$condition <- factor(sample_info$condition, levels = c("AGB6.MOCK", "AGB6.MPXV_lb"))
sample_info$group_label <- factor(sample_info$group_label, levels = c("Control", "Infection"))

rownames(sample_info) <- sample_info$sample

cat("样本信息表:\n")
print(sample_info)

# 检查样本数量
cat("Control样本数(AGB6.MOCK):", sum(sample_info$condition == "AGB6.MOCK"), "\n")
cat("Infection样本数(AGB6.MPXV_lb):", sum(sample_info$condition == "AGB6.MPXV_lb"), "\n")

# 3. 创建DESeq2对象 ----------------------------------------------------
cat("创建DESeq2对象...\n")
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix_subset)
mode(count_matrix_numeric) <- "numeric"

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ condition
)

cat("DESeq2对象创建成功！\n")
cat("基因数:", nrow(dds), "\n")
cat("样本数:", ncol(dds), "\n")

# 4. 过滤低表达基因 ----------------------------------------------------
cat("过滤低表达基因...\n")
# 因为只有6个样本，我们调整过滤条件
keep <- rowSums(counts(dds) >= 10) >= 2  # 至少在2个样本中表达量≥10
dds <- dds[keep, ]
cat("过滤后保留", nrow(dds), "个基因\n")

# 5. 运行DESeq2分析 ----------------------------------------------------
cat("运行DESeq2分析...\n")
dds <- DESeq(dds)

# 6. 提取差异表达结果 ----------------------------------------------------
cat("提取差异表达结果...\n")
# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05

# 获取结果 (AGB6.MPXV_lb vs AGB6.MOCK)
res <- results(dds, contrast = c("condition", "AGB6.MPXV_lb", "AGB6.MOCK"), alpha = 0.05)

# 提取显著差异基因并标记上下调
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(
    change = case_when(
      !is.na(padj) & padj < P_cutoff & log2FoldChange > FC_cutoff ~ "Up-regulated",
      !is.na(padj) & padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down-regulated",
      TRUE ~ "Not Significant"
    ),
    change = factor(change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
  ) %>%
  arrange(padj)

cat("差异基因统计:\n")
cat("上调基因数:", sum(res_df$change == "Up-regulated", na.rm = TRUE), "\n")
cat("下调基因数:", sum(res_df$change == "Down-regulated", na.rm = TRUE), "\n")
cat("总差异基因数:", sum(res_df$change != "Not Significant", na.rm = TRUE), "\n")

# 7. 保存差异表达分析结果到文件 ----------------------------------------
cat("保存差异表达分析结果...\n")

# 创建Excel工作簿
wb_deg <- createWorkbook()

# 保存详细差异表达结果
addWorksheet(wb_deg, "DEGs_Detailed")
writeData(wb_deg, sheet = "DEGs_Detailed", res_df, rowNames = FALSE)

# 保存显著DEGs
sig_degs <- res_df %>%
  filter(change != "Not Significant") %>%
  arrange(padj)

if (nrow(sig_degs) > 0) {
  addWorksheet(wb_deg, "Significant_DEGs")
  writeData(wb_deg, sheet = "Significant_DEGs", sig_degs, rowNames = FALSE)
}

# 保存统计汇总
deg_stats <- data.frame(
  Category = c("上调基因", "下调基因", "总差异基因"),
  Count = c(
    sum(res_df$change == "Up-regulated", na.rm = TRUE),
    sum(res_df$change == "Down-regulated", na.rm = TRUE),
    sum(res_df$change != "Not Significant", na.rm = TRUE)
  )
)

addWorksheet(wb_deg, "DEGs_Statistics")
writeData(wb_deg, sheet = "DEGs_Statistics", deg_stats, rowNames = FALSE)

# 保存样本信息
addWorksheet(wb_deg, "Sample_Information")
writeData(wb_deg, sheet = "Sample_Information", sample_info, rowNames = TRUE)

# 保存分析参数
params <- data.frame(
  Parameter = c("Log2FC cutoff", "P-value cutoff", "Low expression filter", "Comparison"),
  Value = c(FC_cutoff, P_cutoff, "Count >= 10 in >= 2 samples", "AGB6.MPXV_lb vs AGB6.MOCK")
)
addWorksheet(wb_deg, "Analysis_Parameters")
writeData(wb_deg, sheet = "Analysis_Parameters", params, rowNames = FALSE)

# 保存DEGs Excel文件
deg_excel_file <- file.path(output_dir, "AGB6_MPXV_DEGs_Analysis.xlsx")
saveWorkbook(wb_deg, deg_excel_file, overwrite = TRUE)
cat("差异表达结果已保存到:", deg_excel_file, "\n")

# 同时保存CSV格式
write.csv(res_df, file.path(output_dir, "AGB6_MPXV_All_Genes.csv"), row.names = FALSE)
if (nrow(sig_degs) > 0) {
  write.csv(sig_degs, file.path(output_dir, "AGB6_MPXV_Significant_DEGs.csv"), row.names = FALSE)
}
cat("CSV文件已保存\n")

# 8. 提取差异表达基因（DEGs）列表 -------------------------------------
cat("提取差异表达基因列表...\n")

# 提取上调基因
up_genes <- res_df$gene[res_df$change == "Up-regulated"]
down_genes <- res_df$gene[res_df$change == "Down-regulated"]

cat("上调基因数:", length(up_genes), "\n")
cat("下调基因数:", length(down_genes), "\n")

# 保存基因列表
if (length(up_genes) > 0) {
  write.table(up_genes, file.path(output_dir, "AGB6_Upregulated_Genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("上调基因列表已保存\n")
}

if (length(down_genes) > 0) {
  write.table(down_genes, file.path(output_dir, "AGB6_Downregulated_Genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("下调基因列表已保存\n")
}













# 9

cat("进行GO富集分析（Up/Down，输出Excel）...\n")

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(openxlsx)

# 参数
pvalue_cutoff <- 0.1
qvalue_cutoff <- 0.2
minGS <- 5
maxGS <- 1000
padj_method <- "BH"
ont_list <- c("BP","CC","MF")   # 只要BP就改成 c("BP")

# --- 自动猜ID类型（SYMBOL 或 ENSEMBL） ---
guess_id_type <- function(x){
  x <- na.omit(unique(x))
  if(length(x)==0) return("SYMBOL")
  frac_ens <- mean(grepl("^ENS", x))   # ENSMUSG/ENSG...
  if(frac_ens > 0.5) "ENSEMBL" else "SYMBOL"
}

# --- 转ENTREZ（失败则换类型重试）---
to_entrez <- function(ids, orgdb){
  ids <- na.omit(unique(ids))
  if(length(ids)==0) return(data.frame())
  ft <- guess_id_type(ids)
  
  out <- tryCatch({
    bitr(ids, fromType = ft, toType = "ENTREZID", OrgDb = orgdb)
  }, error = function(e) NULL)
  
  # 如果猜错了，再试另一个
  if(is.null(out) || nrow(out)==0){
    alt <- if(ft=="SYMBOL") "ENSEMBL" else "SYMBOL"
    out <- tryCatch({
      bitr(ids, fromType = alt, toType = "ENTREZID", OrgDb = orgdb)
    }, error = function(e) NULL)
  }
  
  if(is.null(out)) data.frame() else out
}

# 背景基因（建议用dds里过滤后的所有基因）
bg_ids <- rownames(dds)

# 转换 Up / Down / Universe 到 ENTREZ
up_map   <- to_entrez(up_genes, org.Mm.eg.db)
down_map <- to_entrez(down_genes, org.Mm.eg.db)
bg_map   <- to_entrez(bg_ids, org.Mm.eg.db)

cat("Up 转ENTREZ成功数:", nrow(up_map), "\n")
cat("Down 转ENTREZ成功数:", nrow(down_map), "\n")
cat("Universe 转ENTREZ成功数:", nrow(bg_map), "\n")

up_entrez   <- unique(up_map$ENTREZID)
down_entrez <- unique(down_map$ENTREZID)
bg_entrez   <- unique(bg_map$ENTREZID)

# enrichGO -> df，并计算 EnrichmentRatio（不筛掉）
enrich_to_df <- function(ego){
  df <- as.data.frame(ego)
  if(nrow(df)==0) return(df)
  df %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      BgRatio_num   = sapply(strsplit(BgRatio, "/"),   function(x) as.numeric(x[1]) / as.numeric(x[2])),
      EnrichmentRatio = GeneRatio_num / BgRatio_num
    ) %>% arrange(p.adjust)
}

run_enrich <- function(gene_entrez, universe_entrez, ont){
  if(length(gene_entrez) == 0) return(NULL)
  tryCatch({
    enrichGO(
      gene          = gene_entrez,
      universe      = universe_entrez,   # 这里类型已统一成ENTREZ
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = padj_method,
      pvalueCutoff  = pvalue_cutoff,
      qvalueCutoff  = qvalue_cutoff,
      minGSSize     = minGS,
      maxGSSize     = maxGS,
      readable      = TRUE
    )
  }, error = function(e){
    cat("enrichGO失败:", ont, " | ", e$message, "\n")
    NULL
  })
}

# 写Excel
wb <- createWorkbook()

# 写基因列表
addWorksheet(wb, "Up_Genes")
writeData(wb, "Up_Genes", data.frame(Gene=up_genes))
addWorksheet(wb, "Down_Genes")
writeData(wb, "Down_Genes", data.frame(Gene=down_genes))

# 逐个ont跑 Up/Down
for(ont in ont_list){
  # Up
  ego_up <- run_enrich(up_entrez, bg_entrez, ont)
  sheet_up <- substr(paste0("Up_GO_", ont), 1, 31)
  addWorksheet(wb, sheet_up)
  if(!is.null(ego_up) && nrow(as.data.frame(ego_up))>0){
    writeData(wb, sheet_up, enrich_to_df(ego_up))
  } else {
    writeData(wb, sheet_up, data.frame(Note="No enriched terms (Up)"))
  }
  
  # Down
  ego_down <- run_enrich(down_entrez, bg_entrez, ont)
  sheet_down <- substr(paste0("Down_GO_", ont), 1, 31)
  addWorksheet(wb, sheet_down)
  if(!is.null(ego_down) && nrow(as.data.frame(ego_down))>0){
    writeData(wb, sheet_down, enrich_to_df(ego_down))
  } else {
    writeData(wb, sheet_down, data.frame(Note="No enriched terms (Down)"))
  }
}

# 参数表
addWorksheet(wb, "GO_Params")
writeData(wb, "GO_Params", data.frame(
  Parameter=c("pvalue_cutoff","qvalue_cutoff","minGSSize","maxGSSize","pAdjustMethod","Ontology"),
  Value=c(pvalue_cutoff,qvalue_cutoff,minGS,maxGS,padj_method,paste(ont_list,collapse=","))
))

go_excel_file <- file.path(output_dir, "AGB6_MPXV_GO_UpDown.xlsx")
saveWorkbook(wb, go_excel_file, overwrite=TRUE)
cat("✅ GO Up/Down 结果已保存：", go_excel_file, "\n")
