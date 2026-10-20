# 确保已安装并加载所需的库
# install.packages(c("clusterProfiler", "org.Mm.eg.db", "ggplot2", "DOSE"))
library(clusterProfiler)
library(org.Mm.eg.db) # 假设您正在使用小鼠 (Mus musculus) 数据
library(ggplot2)
library(DOSE)

# ===============================================
# 1. 模拟数据 (请替换为您的真实 deg_list)
# ===============================================
timepoints <- c("3d", "5d", "7d", "10d")
deg_list <- list()

# 模拟不同时间点的差异表达基因 (DEGs)
set.seed(42)
mock_genes <- sample(keys(org.Mm.eg.db, keytype = "SYMBOL"), 500)

for(tp in timepoints) {
  # 随机选择基因
  n_degs <- sample(200:300, 1)
  genes_subset <- sample(mock_genes, n_degs)
  
  # 随机生成 log2FoldChange
  logFC <- rnorm(n_degs, mean = 0, sd = 2)
  logFC[logFC > 1.5] <- logFC[logFC > 1.5] * 2 # 增加一些显著上调
  logFC[logFC < -1.5] <- logFC[logFC < -1.5] * 2 # 增加一些显著下调
  
  deg_list[[tp]] <- data.frame(
    gene = genes_subset,
    log2FoldChange = logFC,
    stringsAsFactors = FALSE
  )
}

# ===============================================
# 2. 准备基因列表：分离上调和下调基因
# ===============================================
up_gene_list <- list()
down_gene_list <- list()
logFC_cutoff <- 1 # 定义 logFC 阈值

for(tp in timepoints) {
  deg_df <- deg_list[[tp]]
  
  if(nrow(deg_df) > 0) {
    # 提取上调基因 (log2FoldChange > 1)
    up_genes <- deg_df$gene[deg_df$log2FoldChange > logFC_cutoff]
    if(length(up_genes) > 0) {
      up_gene_list[[paste0(tp)]] <- up_genes
    }
    
    # 提取下调基因 (log2FoldChange < -1)
    down_genes <- deg_df$gene[deg_df$log2FoldChange < -logFC_cutoff]
    if(length(down_genes) > 0) {
      down_gene_list[[paste0(tp)]] <- down_genes
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
# 3. 运行 GO 富集比较分析 (上调基因)
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
# 4. 运行 GO 富集比较分析 (下调基因)
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
} else {
  message("注意: 下调基因的 GO 富集结果为空，跳过绘图。")
}

message("\n分析完成。已生成两张独立的点图，分别展示上调和下调基因的富集模式。")