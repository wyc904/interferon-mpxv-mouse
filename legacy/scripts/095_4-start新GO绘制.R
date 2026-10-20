# --------------------------------------------------------------------------
# 模块 A: 设置、加载包与核心绘图函数定义 (针对MPXV分组比较)
# --------------------------------------------------------------------------

# 设置输出目录 (请确保该目录存在或R有权限创建)
output_dir <- "D:/Desktop/1" 
dir.create(output_dir, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 收集列表初始化 (请在运行任何模块前运行此代码)
# --------------------------------------------------------------------------
# 用于收集所有比较组上调基因的GO富集结果
all_ego_up_list <- list() 
# 用于收集所有比较组下调基因的GO富集结果
all_ego_down_list <- list()

# 设置为TRUE以运行GO富集分析
run_enrichment <- TRUE

if(run_enrichment){
  
  # I. 包加载 (确保所有必要的包都已加载)
  # ------------------------------------------------------------
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  
  required_packages_go <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE")
  for(pkg in required_packages_go) {
    if (!require(pkg, character.only = TRUE)) {
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # II. 核心绘图函数定义 (针对MPXV分组比较)
  # ------------------------------------------------------------
  
  #' 定义核心绘图函数：用于生成 Dotplot, Barplot, Cnetplot 三种图
  plot_go_results_multifaceted <- function(ego_result, ont_name, current_gene_list, comparison) {
    
    if (is.null(ego_result) || nrow(ego_result) == 0) return(invisible(NULL))
    
    # 1. Dotplot (点图)
    p_dot <- dotplot(ego_result, showCategory=15, font.size=10) + 
      ggtitle(paste0(comparison, " - GO ", ont_name, " Dotplot"))
    ggsave(file.path(output_dir, paste0("GO_", ont_name, "_dotplot_", comparison, ".pdf")), 
           p_dot, width = 10, height = 8)
    
    # 2. Barplot (条形图)
    p_bar <- barplot(ego_result, showCategory=10, font.size=10) + 
      ggtitle(paste0(comparison, " - Top 10 GO ", ont_name, " Barplot"))
    ggsave(file.path(output_dir, paste0("GO_", ont_name, "_barplot_", comparison, ".pdf")), 
           p_bar, width = 10, height = 6)
    
    # 3. Cnetplot (基因-概念网络图) - 依赖于基因的log2FoldChange
    if (nrow(ego_result) >= 5 && "ENTREZID" %in% colnames(current_gene_list)) {
      fold_changes <- current_gene_list$log2FoldChange
      names(fold_changes) <- current_gene_list$ENTREZID
      
      p_cnet <- cnetplot(ego_result, foldChange=fold_changes, 
                         categorySize="pvalue", showCategory=5, 
                         layout = "kk")
      ggsave(file.path(output_dir, paste0("GO_", ont_name, "_cnetplot_", comparison, ".pdf")), 
             p_cnet, width = 12, height = 10)
    }
  }
  
  # III. 主要 GO 分析流程函数定义 (针对MPXV分组比较)
  # ------------------------------------------------------------
  
  #' 定义主要的GO分析和绘图流程：执行富集、保存和调用绘图函数
  run_go_analysis_and_plot <- function(current_gene_list, comparison) {
    
    # 3. GO富集分析 (BP, CC, MF)
    ego_bp <- enrichGO(gene = current_gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego_cc <- enrichGO(gene = current_gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "CC",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego_mf <- enrichGO(gene = current_gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    
    # 4. 保存富集结果
    write.csv(as.data.frame(ego_bp), file.path(output_dir, paste0("GO_BP_enrichment_", comparison, ".csv")), row.names = FALSE)
    write.csv(as.data.frame(ego_cc), file.path(output_dir, paste0("GO_CC_enrichment_", comparison, ".csv")), row.names = FALSE)
    write.csv(as.data.frame(ego_mf), file.path(output_dir, paste0("GO_MF_enrichment_", comparison, ".csv")), row.names = FALSE)
    
    # 5/7/8/9. 可视化 (BP, CC, MF)
    plot_go_results_multifaceted(ego_bp, "BP", current_gene_list, comparison)
    plot_go_results_multifaceted(ego_cc, "CC", current_gene_list, comparison)
    plot_go_results_multifaceted(ego_mf, "MF", current_gene_list, comparison)
    
    # 10. 上调和下调基因分别进行富集分析 (只针对BP) 并收集结果
    # ------------------------------------------------------------
    up_genes <- current_gene_list %>% filter(log2FoldChange > 0)
    down_genes <- current_gene_list %>% filter(log2FoldChange < 0)
    
    # 上调基因 GO-BP
    if(nrow(up_genes) > 0) {
      ego_up <- enrichGO(gene = up_genes$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
      
      if(nrow(ego_up) > 0) {
        # 将结果保存到全局列表
        all_ego_up_list[[comparison]] <<- ego_up 
        # 同时保留单点图的生成
        p_up <- dotplot(ego_up, showCategory=10) + ggtitle(paste0(comparison, " - GO Enrichment - Upregulated BP"))
        ggsave(file.path(output_dir, paste0("GO_Upregulated_BP_dotplot_", comparison, ".pdf")), p_up, width = 10, height = 6)
      }
    }
    
    # 下调基因 GO-BP
    if(nrow(down_genes) > 0) {
      ego_down <- enrichGO(gene = down_genes$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
      
      if(nrow(ego_down) > 0) {
        # 将结果保存到全局列表
        all_ego_down_list[[comparison]] <<- ego_down
        # 同时保留单点图的生成
        p_down <- dotplot(ego_down, showCategory=10) + ggtitle(paste0(comparison, " - GO Enrichment - Downregulated BP"))
        ggsave(file.path(output_dir, paste0("GO_Downregulated_BP_dotplot_", comparison, ".pdf")), p_down, width = 10, height = 6)
      }
    }
    
    # 11. 输出富集分析摘要
    cat(paste("\n=== ", comparison, " GO富集分析结果摘要 ===\n"))
    cat("生物学过程(BP)富集条目:", ifelse(exists("ego_bp"), nrow(ego_bp), 0), "\n")
    cat("细胞组分(CC)富集条目:", ifelse(exists("ego_cc"), nrow(ego_cc), 0), "\n")
    cat("分子功能(MF)富集条目:", ifelse(exists("ego_mf"), nrow(ego_mf), 0), "\n")
    
    if(exists("ego_bp") && nrow(ego_bp) > 0) {
      cat(paste("\n", comparison, " Top 5 GO Biological Processes:\n"))
      print(head(ego_bp, 5)[, c("Description", "pvalue", "Count")])
    }
  }
} else {
  cat("GO富集分析已设置为不运行 (run_enrichment = FALSE)\n")
}

# --------------------------------------------------------------------------
# 模块 B: 比较组循环与 ID 转换前的准备 (针对MPXV分组比较)
# --------------------------------------------------------------------------

if(run_enrichment){
  # 定义要分析的比较组
  comparisons <- c("MPXV_IIb_vs_MOCK", "MPXV_Ib_vs_MOCK", "MPXV_IIb_vs_MPXV_Ib")
  
  for(comp in comparisons) {
    cat("\n======================================================\n")
    cat(paste(">>> 正在对比较组:", comp, "进行GO富集分析...\n"))
    cat("======================================================\n")
    
    # 0. 获取当前比较组的差异分析结果
    if(!comp %in% names(deg_list) || is.null(deg_list[[comp]])) {
      cat(paste("警告: 未找到比较组", comp, "的差异表达基因列表 (deg_list[[", comp, "]])，跳过。\n"))
      next
    }
    
    res_df_current <- deg_list[[comp]]
    gene_list <- data.frame() # 初始化基因列表
    
    # 1. 提取显著差异基因
    sig_genes <- res_df_current %>% filter(change != "Not Significant")
    sig_genes <- sig_genes %>% rename(gene_id = gene) 
    
    cat(paste(comp, "显著差异基因数量:", nrow(sig_genes), "\n"))
    
    if(nrow(sig_genes) > 0) {
      
      # 2. 尝试 ID 转换 (SYMBOL -> ENTREZID)
      tryCatch({
        eg <- bitr(sig_genes$gene_id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
        gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = "SYMBOL")
        
        cat(paste(comp, "成功转换ID的基因数量 (SYMBOL):", nrow(gene_list), "\n"))
        if(nrow(gene_list) > 0) {
          # 执行GO分析
          run_go_analysis_and_plot(gene_list, comp)
        } else {
          stop("SYMBOL ID 转换成功但结果为空，尝试 ENSEMBL ID")
        }
      }, error = function(e) {
        cat(paste(comp, " ID转换错误 (SYMBOL):", e$message, "\n"))
        
        # 尝试 ENSEMBL ID 转换
        tryCatch({
          eg <- bitr(sig_genes$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
          gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = "ENSEMBL")
          
          cat(paste(comp, "成功转换ID的基因数量 (ENSEMBL):", nrow(gene_list), "\n"))
          if(nrow(gene_list) > 0) {
            # 执行GO分析
            run_go_analysis_and_plot(gene_list, comp)
          } else {
            cat(paste(comp, "警告: ENSEMBL ID 转换成功但结果为空，跳过GO分析。\n"))
          }
        }, error = function(e2) {
          cat(paste(comp, "警告: ENSEMBL ID 转换也失败:", e2$message, "，跳过GO分析。\n"))
        })
      })
      
    } else {
      cat(paste(comp, "没有发现显著差异基因，跳过GO富集分析\n"))
    }
  } # 循环结束
}

# --------------------------------------------------------------------------
# 模块 D: 总结 (针对MPXV分组比较)
# --------------------------------------------------------------------------

cat("\n>>> MPXV分组比较GO富集分析完成！所有比较组结果文件和图表已保存到:", output_dir, "目录中。\n")

# --------------------------------------------------------------------------
# 模块 E: 合并所有比较组的GO富集结果可视化 (可选)
# --------------------------------------------------------------------------

# 如果有多个比较组的GO结果，可以进行比较分析
if(length(all_ego_up_list) > 1) {
  cat("\n>>> 生成比较组GO富集结果比较图...\n")
  
  # 合并上调基因的GO结果
  ego_up_combined <- merge_result(all_ego_up_list)
  if(nrow(ego_up_combined) > 0) {
    p_compare_up <- dotplot(ego_up_combined, showCategory=10) + 
      ggtitle("Comparison of GO BP Enrichment - Upregulated Genes")
    ggsave(file.path(output_dir, "GO_BP_Upregulated_Comparison.pdf"), p_compare_up, width = 12, height = 8)
  }
}

if(length(all_ego_down_list) > 1) {
  # 合并下调基因的GO结果
  ego_down_combined <- merge_result(all_ego_down_list)
  if(nrow(ego_down_combined) > 0) {
    p_compare_down <- dotplot(ego_down_combined, showCategory=10) + 
      ggtitle("Comparison of GO BP Enrichment - Downregulated Genes")
    ggsave(file.path(output_dir, "GO_BP_Downregulated_Comparison.pdf"), p_compare_down, width = 12, height = 8)
  }
}