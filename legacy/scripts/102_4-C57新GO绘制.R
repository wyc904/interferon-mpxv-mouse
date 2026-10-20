# --------------------------------------------------------------------------
# 模块 A: 设置、加载包与核心绘图函数定义
# --------------------------------------------------------------------------

# 设置输出目录 (请确保该目录存在或R有权限创建)
output_dir <- "D:/Desktop/1" 
dir.create(output_dir, showWarnings = FALSE)






# --------------------------------------------------------------------------
# 收集列表初始化 (请在运行任何模块前运行此代码)
# --------------------------------------------------------------------------
# 用于收集所有时间点上调基因的GO富集结果
all_ego_up_list <- list() 
# 用于收集所有时间点下调基因的GO富集结果
all_ego_down_list <- list()








# 设置为TRUE以运行GO富集分析
run_enrichment <- TRUE

if(run_enrichment){
  
  # I. 包加载 (确保所有必要的包都已加载)
  # ------------------------------------------------------------
  # 这一步仅在交互式环境中运行一次，确保依赖包到位
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  
  required_packages_go <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE")
  for(pkg in required_packages_go) {
    if (!require(pkg, character.only = TRUE)) {
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # II. 核心绘图函数定义
  # ------------------------------------------------------------
  
  #' 定义核心绘图函数：用于生成 Dotplot, Barplot, Cnetplot 三种图
  plot_go_results_multifaceted <- function(ego_result, ont_name, current_gene_list, timepoint) {
    
    if (is.null(ego_result) || nrow(ego_result) == 0) return(invisible(NULL))
    
    # 1. Dotplot (点图)
    p_dot <- dotplot(ego_result, showCategory=15, font.size=10) + 
      ggtitle(paste0(timepoint, " - GO ", ont_name, " Dotplot"))
    ggsave(file.path(output_dir, paste0("GO_", ont_name, "_dotplot_", timepoint, ".pdf")), 
           p_dot, width = 10, height = 8)
    
    # 2. Barplot (条形图)
    p_bar <- barplot(ego_result, showCategory=10, font.size=10) + 
      ggtitle(paste0(timepoint, " - Top 10 GO ", ont_name, " Barplot"))
    ggsave(file.path(output_dir, paste0("GO_", ont_name, "_barplot_", timepoint, ".pdf")), 
           p_bar, width = 10, height = 6)
    
    # 3. Cnetplot (基因-概念网络图) - 依赖于基因的log2FoldChange
    if (nrow(ego_result) >= 5 && "ENTREZID" %in% colnames(current_gene_list)) {
      fold_changes <- current_gene_list$log2FoldChange
      names(fold_changes) <- current_gene_list$ENTREZID
      
      p_cnet <- cnetplot(ego_result, foldChange=fold_changes, 
                         categorySize="pvalue", showCategory=5, 
                         layout = "kk")
      ggsave(file.path(output_dir, paste0("GO_", ont_name, "_cnetplot_", timepoint, ".pdf")), 
             p_cnet, width = 12, height = 10)
    }
  }
  
  # III. 主要 GO 分析流程函数定义
  # ------------------------------------------------------------
  
  #' 定义主要的GO分析和绘图流程：执行富集、保存和调用绘图函数
  run_go_analysis_and_plot <- function(current_gene_list, timepoint) {
    
    # 3. GO富集分析 (BP, CC, MF)
    ego_bp <- enrichGO(gene = current_gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego_cc <- enrichGO(gene = current_gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "CC",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego_mf <- enrichGO(gene = current_gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    
    # 4. 保存富集结果
    write.csv(as.data.frame(ego_bp), file.path(output_dir, paste0("GO_BP_enrichment_", timepoint, ".csv")), row.names = FALSE)
    write.csv(as.data.frame(ego_cc), file.path(output_dir, paste0("GO_CC_enrichment_", timepoint, ".csv")), row.names = FALSE)
    write.csv(as.data.frame(ego_mf), file.path(output_dir, paste0("GO_MF_enrichment_", timepoint, ".csv")), row.names = FALSE)
    
    # 5/7/8/9. 可视化 (BP, CC, MF)
    plot_go_results_multifaceted(ego_bp, "BP", current_gene_list, timepoint)
    plot_go_results_multifaceted(ego_cc, "CC", current_gene_list, timepoint)
    plot_go_results_multifaceted(ego_mf, "MF", current_gene_list, timepoint)
    
    # 10. 上调和下调基因分别进行富集分析 (只针对BP)
    # 10. 上调和下调基因分别进行富集分析 (只针对BP) 并收集结果
    # ------------------------------------------------------------
    up_genes <- current_gene_list %>% filter(log2FoldChange > 0)
    down_genes <- current_gene_list %>% filter(log2FoldChange < 0)
    
    # 上调基因 GO-BP
    if(nrow(up_genes) > 0) {
      ego_up <- enrichGO(gene = up_genes$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
      
      if(nrow(ego_up) > 0) {
        # 将结果保存到全局列表 (使用 <<- 跨作用域赋值)
        all_ego_up_list[[timepoint]] <<- ego_up 
        # 同时保留单点图的生成 (如果需要)
        p_up <- dotplot(ego_up, showCategory=10) + ggtitle(paste0(timepoint, " - GO Enrichment - Upregulated BP"))
        ggsave(file.path(output_dir, paste0("GO_Upregulated_BP_dotplot_", timepoint, ".pdf")), p_up, width = 10, height = 6)
      }
    }
    
    # 下调基因 GO-BP
    if(nrow(down_genes) > 0) {
      ego_down <- enrichGO(gene = down_genes$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
      
      if(nrow(ego_down) > 0) {
        # 将结果保存到全局列表 (使用 <<- 跨作用域赋值)
        all_ego_down_list[[timepoint]] <<- ego_down
        # 同时保留单点图的生成 (如果需要)
        p_down <- dotplot(ego_down, showCategory=10) + ggtitle(paste0(timepoint, " - GO Enrichment - Downregulated BP"))
        ggsave(file.path(output_dir, paste0("GO_Downregulated_BP_dotplot_", timepoint, ".pdf")), p_down, width = 10, height = 6)
      }
    }
    # 11. 输出富集分析摘要
    cat(paste("\n=== ", timepoint, " GO富集分析结果摘要 ===\n"))
    cat("生物学过程(BP)富集条目:", ifelse(exists("ego_bp"), nrow(ego_bp), 0), "\n")
    cat("细胞组分(CC)富集条目:", ifelse(exists("ego_cc"), nrow(ego_cc), 0), "\n")
    cat("分子功能(MF)富集条目:", ifelse(exists("ego_mf"), nrow(ego_mf), 0), "\n")
    
    if(exists("ego_bp") && nrow(ego_bp) > 0) {
      cat(paste("\n", timepoint, " Top 5 GO Biological Processes:\n"))
      print(head(ego_bp, 5)[, c("Description", "pvalue", "Count")])
    }
  }
} else {
  cat("GO富集分析已设置为不运行 (run_enrichment = FALSE)\n")
}






# --------------------------------------------------------------------------
# 模块 B: 时间点循环与 ID 转换前的准备
# --------------------------------------------------------------------------

if(run_enrichment){
  for(tp in timepoints) {
    cat("\n======================================================\n")
    cat(paste(">>> 正在对时间点:", tp, "进行GO富集分析...\n"))
    cat("======================================================\n")
    
    # 0. 获取当前时间点的差异分析结果
    if(!tp %in% names(deg_list) || is.null(deg_list[[tp]])) {
      cat(paste("警告: 未找到时间点", tp, "的差异表达基因列表 (deg_list[[", tp, "]])，跳过。\n"))
      next
    }
    
    res_df_current <- deg_list[[tp]]
    gene_list <- data.frame() # 初始化基因列表
    
    # 1. 提取显著差异基因
    sig_genes <- res_df_current %>% filter(change != "Not Significant")
    sig_genes <- sig_genes %>% rename(gene_id = gene) 
    
    cat(paste(tp, "显著差异基因数量:", nrow(sig_genes), "\n"))
    
    if(nrow(sig_genes) > 0) {
      
      # 2. 尝试 ID 转换 (SYMBOL -> ENTREZID)
      tryCatch({
        eg <- bitr(sig_genes$gene_id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
        gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = "SYMBOL")
        
        cat(paste(tp, "成功转换ID的基因数量 (SYMBOL):", nrow(gene_list), "\n"))
        if(nrow(gene_list) > 0) {
          # 跳转到 模块 C: 执行分析
          run_go_analysis_and_plot(gene_list, tp)
        } else {
          stop("SYMBOL ID 转换成功但结果为空，尝试 ENSEMBL ID")
        }
      }, error = function(e) {
        cat(paste(tp, " ID转换错误 (SYMBOL):", e$message, "\n"))
        
        # 模块 C（仅在 ID 转换失败时执行的 ENSEMBL 尝试）
        # 尝试 ENSEMBL ID 转换
        tryCatch({
          eg <- bitr(sig_genes$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
          gene_list <- merge(sig_genes, eg, by.x = "gene_id", by.y = "ENSEMBL")
          
          cat(paste(tp, "成功转换ID的基因数量 (ENSEMBL):", nrow(gene_list), "\n"))
          if(nrow(gene_list) > 0) {
            # 跳转到 模块 C: 执行分析
            run_go_analysis_and_plot(gene_list, tp)
          } else {
            cat(paste(tp, "警告: ENSEMBL ID 转换成功但结果为空，跳过GO分析。\n"))
          }
        }, error = function(e2) {
          cat(paste(tp, "警告: ENSEMBL ID 转换也失败:", e2$message, "，跳过GO分析。\n"))
        })
      })
      
    } else {
      cat(paste(tp, "没有发现显著差异基因，跳过GO富集分析\n"))
    }
  } # 循环结束
}





# --------------------------------------------------------------------------
# 模块 D: 总结
# --------------------------------------------------------------------------

cat("\n>>> GO富集分析完成！所有时间点结果文件和图表已保存到:", output_dir, "目录中。\n")








