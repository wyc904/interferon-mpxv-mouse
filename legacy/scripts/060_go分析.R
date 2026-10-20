# 对每个时间点进行独立的GO富集分析
timepoints <- c("3d", "5d", "7d", "10d")
go_results <- list()

for(tp in timepoints) {
  # 获取该时间点的差异基因
  degs <- deg_list[[tp]]$gene  # 假设deg_list已包含各时间点DEGs
  
  if(length(degs) > 0) {
    # GO富集分析
    go_enrich <- enrichGO(gene = degs,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont = "ALL",  # 同时分析BP, MF, CC
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1,
                          readable = TRUE)
    
    go_results[[tp]] <- go_enrich
    
    # 保存结果
    write.csv(as.data.frame(go_enrich), 
              paste0("GO_analysis_", tp, "_vs_mock.csv"))
  }
}






###############






# 比较不同时间点的GO富集模式
library(clusterProfiler)

# 准备基因列表（各时间点的上调和下调基因）
gene_list <- list()
for(tp in timepoints) {
  deg_df <- deg_list[[tp]]
  if(nrow(deg_df) > 0) {
    # 分别获取上调和下调基因
    up_genes <- deg_df$gene[deg_df$log2FoldChange > 1]
    down_genes <- deg_df$gene[deg_df$log2FoldChange < -1]
    
    gene_list[[paste0(tp, "_UP")]] <- up_genes
    gene_list[[paste0(tp, "_DOWN")]] <- down_genes
  }
}

# 比较聚类分析
ck <- compareCluster(geneCluster = gene_list, 
                     fun = "enrichGO",
                     OrgDb = org.Mm.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pvalueCutoff = 0.05)

# 可视化比较结果
dotplot(ck, showCategory = 5, 
        title = "GO Biological Process Enrichment Across Timepoints")








####################


# 分析重要GO条目随时间的富集程度变化
library(ggplot2)
library(dplyr)

# 提取关键GO条目的时间动态
key_go_terms <- c("type I interferon signaling pathway",
                  "response to virus",
                  "inflammatory response",
                  "T cell activation",
                  "apoptotic process")

go_trend_data <- data.frame()

for(tp in timepoints) {
  go_df <- as.data.frame(go_results[[tp]])
  if(nrow(go_df) > 0) {
    for(term in key_go_terms) {
      term_match <- go_df[grep(term, go_df$Description, ignore.case = TRUE), ]
      if(nrow(term_match) > 0) {
        go_trend_data <- rbind(go_trend_data,
                               data.frame(Timepoint = tp,
                                          GOTerm = term,
                                          Count = term_match$Count[1],
                                          GeneRatio = term_match$GeneRatio[1],
                                          pvalue = term_match$pvalue[1],
                                          p.adjust = term_match$p.adjust[1]))
      }
    }
  }
}

# 绘制GO条目时间趋势图
ggplot(go_trend_data, aes(x = Timepoint, y = -log10(p.adjust), 
                          color = GOTerm, group = GOTerm)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "关键GO条目富集程度的时间动态",
       x = "时间点", y = "-log10(调整后p值)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




##################




# 完整的干扰素GO分析代码（修复版）
library(dplyr)
library(tidyr)

# 提取干扰素相关GO条目
ifn_related_terms <- c("interferon", "ISG", "antiviral", "JAK-STAT", "IRF", "response to virus")

ifn_go_list <- list()

for(tp in timepoints) {
  if(tp %in% names(go_results) && !is.null(go_results[[tp]])) {
    go_df <- as.data.frame(go_results[[tp]])
    if(nrow(go_df) > 0) {
      # 筛选干扰素相关条目
      ifn_terms <- go_df[grepl(paste(ifn_related_terms, collapse = "|"), 
                               go_df$Description, ignore.case = TRUE), ]
      if(nrow(ifn_terms) > 0) {
        ifn_terms$Timepoint <- tp
        ifn_go_list[[tp]] <- ifn_terms
      }
    }
  }
}

# 合并所有时间点的结果
ifn_go_results <- do.call(rbind, ifn_go_list)

if(!is.null(ifn_go_results) && nrow(ifn_go_results) > 0) {
  # 转换为宽格式（修复版）
  ifn_go_wide <- ifn_go_results %>%
    dplyr::select(Timepoint, Description, p.adjust) %>%
    tidyr::pivot_wider(
      names_from = Timepoint, 
      values_from = p.adjust,
      values_fill = list(p.adjust = 1)  # 用1填充缺失值（表示不显著）
    )
  
  # 检查结果
  print("转换后的宽格式数据:")
  print(head(ifn_go_wide))
  
  # 转换为矩阵用于热图
  ifn_go_matrix <- as.matrix(ifn_go_wide[, -1])
  rownames(ifn_go_matrix) <- ifn_go_wide$Description
  
  # 绘制热图
  pheatmap::pheatmap(-log10(ifn_go_matrix),
                     main = "干扰素相关GO条目富集热图",
                     color = colorRampPalette(c("white", "red"))(100),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     na_col = "gray90",
                     fontsize_row = 8)
} else {
  print("未找到显著的干扰素相关GO条目")
}










































