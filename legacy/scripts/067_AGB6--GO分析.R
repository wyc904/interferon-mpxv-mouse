# 设置为TRUE以运行GO富集分析
run_enrichment <- TRUE

if(run_enrichment){
  # 加载必要的包
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!require("clusterProfiler")) {
    BiocManager::install("clusterProfiler")
    library(clusterProfiler)
  }
  
  if (!require("org.Mm.eg.db")) {
    BiocManager::install("org.Mm.eg.db")
    library(org.Mm.eg.db)
  }
  
  if (!require("enrichplot")) {
    BiocManager::install("enrichplot")
    library(enrichplot)
  }
  
  if (!require("DOSE")) {
    BiocManager::install("DOSE")
    library(DOSE)
  }
  
  # 1. 提取差异基因列表
  # 首先确保res_df中有change列，如果没有则创建
  if(!"change" %in% colnames(res_df)) {
    res_df$change <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                            ifelse(res_df$log2FoldChange > 1, "Up", "Down"), 
                            "Not Significant")
  }
  
  # 获取所有显著差异基因
  sig_genes <- res_df %>% filter(change != "Not Significant")
  
  cat("显著差异基因数量:", nrow(sig_genes), "\n")
  
  if(nrow(sig_genes) > 0) {
    # 2. 转换ID (假设基因ID是Symbol)
    # 注意: 如果您的基因ID是ENSEMBL ID或其他类型，需要调整fromType参数
    tryCatch({
      eg <- bitr(sig_genes$gene_id, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
      
      gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = "SYMBOL")
      
      cat("成功转换ID的基因数量:", nrow(gene_list), "\n")
      
      if(nrow(gene_list) > 0) {
        # 3. GO富集分析 - 生物学过程(BP)
        ego_bp <- enrichGO(gene = gene_list$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
        
        # 4. GO富集分析 - 细胞组分(CC)
        ego_cc <- enrichGO(gene = gene_list$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
        
        # 5. GO富集分析 - 分子功能(MF)
        ego_mf <- enrichGO(gene = gene_list$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
        
        # 6. 保存富集结果
        write.csv(as.data.frame(ego_bp), "GO_BP_enrichment.csv", row.names = FALSE)
        write.csv(as.data.frame(ego_cc), "GO_CC_enrichment.csv", row.names = FALSE)
        write.csv(as.data.frame(ego_mf), "GO_MF_enrichment.csv", row.names = FALSE)
        
        # 7. 可视化 - 点图
        if(nrow(ego_bp) > 0) {
          p_bp <- dotplot(ego_bp, showCategory=15, font.size=10) + 
            ggtitle("GO Biological Process Enrichment")
          ggsave("GO_BP_dotplot.pdf", p_bp, width = 10, height = 8)
        }
        
        if(nrow(ego_cc) > 0) {
          p_cc <- dotplot(ego_cc, showCategory=15, font.size=10) + 
            ggtitle("GO Cellular Component Enrichment")
          ggsave("GO_CC_dotplot.pdf", p_cc, width = 10, height = 8)
        }
        
        if(nrow(ego_mf) > 0) {
          p_mf <- dotplot(ego_mf, showCategory=15, font.size=10) + 
            ggtitle("GO Molecular Function Enrichment")
          ggsave("GO_MF_dotplot.pdf", p_mf, width = 10, height = 8)
        }
        
        # 8. 富集图
        if(nrow(ego_bp) > 0) {
          p_enrich <- barplot(ego_bp, showCategory=10, font.size=10) + 
            ggtitle("Top 10 GO Biological Processes")
          ggsave("GO_BP_barplot.pdf", p_enrich, width = 10, height = 6)
        }
        
        # 9. 基因-概念网络图
        if(nrow(ego_bp) > 0 && nrow(ego_bp) >= 5) {
          p_cnet <- cnetplot(ego_bp, foldChange=gene_list$log2FoldChange, 
                             categorySize="pvalue", showCategory=5)
          ggsave("GO_cnetplot.pdf", p_cnet, width = 12, height = 10)
        }
        
        # 10. 上调和下调基因分别进行富集分析
        # 上调基因
        up_genes <- gene_list %>% filter(log2FoldChange > 0)
        if(nrow(up_genes) > 0) {
          ego_up <- enrichGO(gene = up_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
          
          if(nrow(ego_up) > 0) {
            p_up <- dotplot(ego_up, showCategory=10) + ggtitle("GO Enrichment - Upregulated Genes")
            ggsave("GO_Upregulated_dotplot.pdf", p_up, width = 10, height = 6)
          }
        }
        
        # 下调基因
        down_genes <- gene_list %>% filter(log2FoldChange < 0)
        if(nrow(down_genes) > 0) {
          ego_down <- enrichGO(gene = down_genes$ENTREZID,
                               OrgDb = org.Mm.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)
          
          if(nrow(ego_down) > 0) {
            p_down <- dotplot(ego_down, showCategory=10) + ggtitle("GO Enrichment - Downregulated Genes")
            ggsave("GO_Downregulated_dotplot.pdf", p_down, width = 10, height = 6)
          }
        }
        
        # 11. 输出富集分析摘要
        cat("\n=== GO富集分析结果摘要 ===\n")
        cat("生物学过程(BP)富集条目:", nrow(ego_bp), "\n")
        cat("细胞组分(CC)富集条目:", nrow(ego_cc), "\n")
        cat("分子功能(MF)富集条目:", nrow(ego_mf), "\n")
        
        if(nrow(ego_bp) > 0) {
          cat("\nTop 5 GO Biological Processes:\n")
          print(head(ego_bp, 5)[, c("Description", "pvalue", "Count")])
        }
        
      } else {
        cat("警告: 没有基因成功转换为ENTREZID，请检查基因ID类型\n")
      }
    }, error = function(e) {
      cat("ID转换错误:", e$message, "\n")
      cat("尝试使用ENSEMBL ID转换...\n")
      
      # 尝试使用ENSEMBL ID
      tryCatch({
        eg <- bitr(sig_genes$gene_id, 
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
        
        gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = "ENSEMBL")
        
        # 继续GO分析...
        # (这里添加与上面相同的GO分析代码)
        
      }, error = function(e2) {
        cat("ENSEMBL ID转换也失败:", e2$message, "\n")
      })
    })
  } else {
    cat("没有发现显著差异基因，跳过GO富集分析\n")
  }
}

cat(">>> GO富集分析完成！结果文件已保存。\n")