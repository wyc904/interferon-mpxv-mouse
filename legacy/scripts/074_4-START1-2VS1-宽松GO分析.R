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
  
  # 1. 诊断问题
  cat("=== GO富集分析问题诊断 ===\n")
  
  # 检查显著差异基因
  if(!"change" %in% colnames(res_df)) {
    res_df$change <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                            ifelse(res_df$log2FoldChange > 1, "Up", "Down"), 
                            "Not Significant")
  }
  
  sig_genes <- res_df %>% filter(change != "Not Significant")
  cat("显著差异基因数量:", nrow(sig_genes), "\n")
  
  if(nrow(sig_genes) == 0) {
    cat("没有显著差异基因，尝试放宽阈值...\n")
    # 尝试使用更宽松的阈值
    sig_genes <- res_df %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5)
    cat("放宽阈值后的显著基因数量:", nrow(sig_genes), "\n")
  }
  
  if(nrow(sig_genes) > 0) {
    # 2. 检查基因ID类型
    cat("\n检查基因ID类型...\n")
    
    # 测试前几个基因的ID类型
    test_genes <- head(sig_genes$gene_id, 5)
    cat("前5个基因ID:", paste(test_genes, collapse = ", "), "\n")
    
    # 尝试多种ID类型转换
    id_types <- c("SYMBOL", "ENSEMBL", "ENTREZID", "ALIAS")
    successful_conversion <- FALSE
    
    for(id_type in id_types) {
      cat("尝试使用", id_type, "进行转换...\n")
      
      tryCatch({
        # 尝试转换
        eg <- bitr(sig_genes$gene_id, 
                   fromType = id_type, 
                   toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)
        
        if(nrow(eg) > 0) {
          cat("成功转换", nrow(eg), "个基因，使用", id_type, "类型\n")
          gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = id_type)
          successful_conversion <- TRUE
          used_id_type <- id_type
          break
        }
      }, error = function(e) {
        cat(id_type, "转换失败:", e$message, "\n")
      })
    }
    
    if(successful_conversion) {
      cat("成功转换ID的基因数量:", nrow(gene_list), "\n")
      
      # 3. 调整GO富集分析参数（放宽条件）
      cat("\n进行GO富集分析（使用宽松参数）...\n")
      
      # 使用更宽松的参数
      ego_bp <- enrichGO(gene = gene_list$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,    # 放宽p值阈值
                         qvalueCutoff = 0.2,
                         minGSSize = 5,        # 减小最小基因集大小
                         maxGSSize = 500,      # 限制最大基因集大小
                         readable = TRUE)
      
      ego_cc <- enrichGO(gene = gene_list$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,
                         qvalueCutoff = 0.2,
                         minGSSize = 5,
                         maxGSSize = 500,
                         readable = TRUE)
      
      ego_mf <- enrichGO(gene = gene_list$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,
                         qvalueCutoff = 0.2,
                         minGSSize = 5,
                         maxGSSize = 500,
                         readable = TRUE)
      
      # 4. 检查富集结果
      cat("\n=== GO富集分析结果 ===\n")
      cat("生物学过程(BP)富集条目:", nrow(ego_bp), "\n")
      cat("细胞组分(CC)富集条目:", nrow(ego_cc), "\n")
      cat("分子功能(MF)富集条目:", nrow(ego_mf), "\n")
      
      # 如果仍然没有结果，尝试其他方法
      if(nrow(ego_bp) == 0 & nrow(ego_cc) == 0 & nrow(ego_mf) == 0) {
        cat("\nGO富集无结果，尝试替代方法...\n")
        
        # 方法1: 使用更宽松的阈值
        ego_bp_loose <- enrichGO(gene = gene_list$ENTREZID,
                                 OrgDb = org.Mm.eg.db,
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.5,  # 非常宽松的p值
                                 qvalueCutoff = 0.5,
                                 minGSSize = 3,       # 更小的基因集
                                 readable = TRUE)
        
        cat("宽松阈值BP富集条目:", nrow(ego_bp_loose), "\n")
        
        # 方法2: 尝试KEGG富集分析
        if(require("org.Mm.eg.db")) {
          kegg <- enrichKEGG(gene = gene_list$ENTREZID,
                             organism = "mmu",  # 小鼠
                             pvalueCutoff = 0.1,
                             qvalueCutoff = 0.2)
          cat("KEGG通路富集条目:", nrow(kegg), "\n")
          
          if(nrow(kegg) > 0) {
            write.csv(as.data.frame(kegg), "KEGG_enrichment.csv", row.names = FALSE)
            p_kegg <- dotplot(kegg, showCategory=15) + ggtitle("KEGG Pathway Enrichment")
            ggsave("KEGG_dotplot.pdf", p_kegg, width = 10, height = 8)
          }
        }
        
        # 方法3: 尝试DO（疾病本体）富集
        if(require("DOSE")) {
          do <- enrichDO(gene = gene_list$ENTREZID,
                         pvalueCutoff = 0.1,
                         qvalueCutoff = 0.2)
          cat("疾病本体富集条目:", nrow(do), "\n")
        }
        
      } else {
        # 如果有GO富集结果，保存和可视化
        write.csv(as.data.frame(ego_bp), "GO_BP_enrichment.csv", row.names = FALSE)
        write.csv(as.data.frame(ego_cc), "GO_CC_enrichment.csv", row.names = FALSE)
        write.csv(as.data.frame(ego_mf), "GO_MF_enrichment.csv", row.names = FALSE)
        
        # 可视化
        if(nrow(ego_bp) > 0) {
          p_bp <- dotplot(ego_bp, showCategory=15) + ggtitle("GO Biological Process")
          ggsave("GO_BP_dotplot.pdf", p_bp, width = 10, height = 6)
        }
        if(nrow(ego_cc) > 0) {
          p_cc <- dotplot(ego_cc, showCategory=15) + ggtitle("GO Cellular Component")
          ggsave("GO_CC_dotplot.pdf", p_cc, width = 10, height = 6)
        }
        if(nrow(ego_mf) > 0) {
          p_mf <- dotplot(ego_mf, showCategory=15) + ggtitle("GO Molecular Function")
          ggsave("GO_MF_dotplot.pdf", p_mf, width = 10, height = 6)
        }
      }
      
    } else {
      cat("所有ID转换方法都失败了。请检查基因ID格式。\n")
      cat("常见的基因ID格式示例:\n")
      cat("SYMBOL: Actb, Gapdh, Tnf\n")
      cat("ENSEMBL: ENSMUSG00000059552, ENSMUSG00000031710\n")
      cat("ENTREZID: 11461, 14433, 21953\n")
    }
    
  } else {
    cat("没有显著差异基因可用于GO富集分析。\n")
    cat("建议:\n")
    cat("1. 检查差异分析阈值是否合适\n")
    cat("2. 考虑使用更宽松的阈值(padj < 0.1, |log2FC| > 0.5)\n")
    cat("3. 检查数据质量\n")
  }
}

cat(">>> GO富集分析完成！\n")