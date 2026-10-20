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
  
  if (!require("openxlsx")) {
    install.packages("openxlsx")
    library(openxlsx)
  }
  
  # 创建输出目录
  go_output_dir <- "GO_Enrichment_Analysis"
  if (!dir.exists(go_output_dir)) {
    dir.create(go_output_dir, recursive = TRUE)
  }
  
  cat("GO分析结果将保存到目录:", go_output_dir, "\n")
  
  # 1. 检查res_df是否存在
  if(!exists("res_df")) {
    stop("错误：未找到res_df数据框。请确保已运行差异表达分析并生成res_df。")
  }
  
  cat("输入数据维度:", dim(res_df), "\n")
  cat("列名:", colnames(res_df), "\n")
  
  # 2. 创建change列（使用更宽松的阈值）
  cat("创建change列...\n")
  
  # 查看数据统计
  cat("\n=== 数据统计 ===\n")
  cat("总基因数:", nrow(res_df), "\n")
  cat("padj < 0.05 的基因数:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("padj为NA的基因数:", sum(is.na(res_df$padj)), "\n")
  cat("log2FoldChange > 1 的基因数:", sum(res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
  cat("log2FoldChange < -1 的基因数:", sum(res_df$log2FoldChange < -1, na.rm = TRUE), "\n")
  
  # 使用更宽松的阈值：padj < 0.1 或 pvalue < 0.05
  if(!"change" %in% colnames(res_df)) {
    # 方法1：使用pvalue而不是padj（更宽松）
    res_df$change <- ifelse(
      (!is.na(res_df$pvalue) & res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > 0.5) |
        (!is.na(res_df$padj) & res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 0.5),
      ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
      "Not Significant"
    )
  }
  
  # 统计change列
  cat("\n=== Change列统计 ===\n")
  cat("上调基因数(Up):", sum(res_df$change == "Up", na.rm = TRUE), "\n")
  cat("下调基因数(Down):", sum(res_df$change == "Down", na.rm = TRUE), "\n")
  cat("不显著基因数:", sum(res_df$change == "Not Significant", na.rm = TRUE), "\n")
  
  # 获取所有显著差异基因
  sig_genes <- res_df %>% filter(change != "Not Significant")
  cat("\n显著差异基因数量:", nrow(sig_genes), "\n")
  
  if(nrow(sig_genes) > 0) {
    # 3. 使用基因符号直接进行GO分析（跳过ID转换）
    cat("\n使用基因符号直接进行GO分析...\n")
    
    # 创建Excel工作簿
    wb <- createWorkbook()
    
    # 初始化统计变量
    up_go_count <- 0
    down_go_count <- 0
    
    # 上调基因
    up_genes <- sig_genes %>% filter(change == "Up")
    if(nrow(up_genes) >= 3) {  # 降低到至少3个基因
      cat("上调基因数量:", nrow(up_genes), "\n")
      
      # 查看前几个上调基因
      cat("前5个上调基因:", head(up_genes$gene_id, 5), "\n")
      
      # 上调基因GO分析 - 使用基因符号
      ego_up <- enrichGO(gene = up_genes$gene_id,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",  # 使用SYMBOL类型
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,   # 放宽阈值
                         qvalueCutoff = 0.3,   # 放宽阈值
                         minGSSize = 3,        # 降低最小基因集
                         maxGSSize = 1000,
                         readable = TRUE)
      
      # 保存上调基因GO结果到Excel
      if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
        ego_up_df <- as.data.frame(ego_up)
        # 添加标注列
        ego_up_df$Regulation <- "上调通路"
        ego_up_df$Regulation_English <- "Up-regulated"
        up_go_count <- nrow(ego_up_df)
        
        addWorksheet(wb, "上调通路_GO分析")
        writeData(wb, sheet = "上调通路_GO分析", ego_up_df, rowNames = FALSE)
        cat("上调基因GO分析结果已保存，条目数:", up_go_count, "\n")
        
        # 绘制上调基因GO点图
        p_up <- dotplot(ego_up, showCategory=15, font.size=10) + 
          ggtitle("GO Enrichment - 上调通路")
        ggsave(file.path(go_output_dir, "上调通路_GO点图.pdf"), p_up, width = 10, height = 8)
        
        # 绘制条形图
        p_up_bar <- barplot(ego_up, showCategory=10, font.size=10) + 
          ggtitle("Top 10 GO Terms - 上调通路")
        ggsave(file.path(go_output_dir, "上调通路_GO条形图.pdf"), p_up_bar, width = 10, height = 6)
        
        # 输出top结果
        cat("\nTop 5 上调通路:\n")
        top_up <- head(ego_up_df, 5)
        for(i in 1:nrow(top_up)) {
          cat(i, ". ", top_up$Description[i], 
              " (p=", format(top_up$pvalue[i], scientific=TRUE, digits=3), 
              ", padj=", format(top_up$p.adjust[i], scientific=TRUE, digits=3), 
              ")\n", sep="")
        }
      } else {
        cat("上调基因GO分析无显著结果\n")
      }
    } else {
      cat("上调基因数量不足3个，跳过上调基因GO分析\n")
    }
    
    # 下调基因
    down_genes <- sig_genes %>% filter(change == "Down")
    if(nrow(down_genes) >= 3) {  # 降低到至少3个基因
      cat("\n下调基因数量:", nrow(down_genes), "\n")
      
      # 查看前几个下调基因
      cat("前5个下调基因:", head(down_genes$gene_id, 5), "\n")
      
      # 下调基因GO分析 - 使用基因符号
      ego_down <- enrichGO(gene = down_genes$gene_id,
                           OrgDb = org.Mm.eg.db,
                           keyType = "SYMBOL",  # 使用SYMBOL类型
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.1,   # 放宽阈值
                           qvalueCutoff = 0.3,   # 放宽阈值
                           minGSSize = 3,        # 降低最小基因集
                           maxGSSize = 1000,
                           readable = TRUE)
      
      # 保存下调基因GO结果到Excel
      if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
        ego_down_df <- as.data.frame(ego_down)
        # 添加标注列
        ego_down_df$Regulation <- "下调通路"
        ego_down_df$Regulation_English <- "Down-regulated"
        down_go_count <- nrow(ego_down_df)
        
        addWorksheet(wb, "下调通路_GO分析")
        writeData(wb, sheet = "下调通路_GO分析", ego_down_df, rowNames = FALSE)
        cat("下调基因GO分析结果已保存，条目数:", down_go_count, "\n")
        
        # 绘制下调基因GO点图
        p_down <- dotplot(ego_down, showCategory=15, font.size=10) + 
          ggtitle("GO Enrichment - 下调通路")
        ggsave(file.path(go_output_dir, "下调通路_GO点图.pdf"), p_down, width = 10, height = 8)
        
        # 绘制条形图
        p_down_bar <- barplot(ego_down, showCategory=10, font.size=10) + 
          ggtitle("Top 10 GO Terms - 下调通路")
        ggsave(file.path(go_output_dir, "下调通路_GO条形图.pdf"), p_down_bar, width = 10, height = 6)
        
        # 输出top结果
        cat("\nTop 5 下调通路:\n")
        top_down <- head(ego_down_df, 5)
        for(i in 1:nrow(top_down)) {
          cat(i, ". ", top_down$Description[i], 
              " (p=", format(top_down$pvalue[i], scientific=TRUE, digits=3), 
              ", padj=", format(top_down$p.adjust[i], scientific=TRUE, digits=3), 
              ")\n", sep="")
        }
      } else {
        cat("下调基因GO分析无显著结果\n")
      }
    } else {
      cat("下调基因数量不足3个，跳过下调基因GO分析\n")
    }
    
    # 4. 如果上下调都没有结果，尝试使用所有差异基因
    if(up_go_count == 0 && down_go_count == 0 && nrow(sig_genes) >= 5) {
      cat("\n上下调基因均无GO结果，尝试使用所有差异基因进行分析...\n")
      
      ego_all <- enrichGO(gene = sig_genes$gene_id,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.3,
                          minGSSize = 3,
                          maxGSSize = 1000,
                          readable = TRUE)
      
      if (!is.null(ego_all) && nrow(as.data.frame(ego_all)) > 0) {
        ego_all_df <- as.data.frame(ego_all)
        ego_all_df$Regulation <- "所有差异基因"
        ego_all_df$Regulation_English <- "All DEGs"
        
        addWorksheet(wb, "所有差异基因_GO分析")
        writeData(wb, sheet = "所有差异基因_GO分析", ego_all_df, rowNames = FALSE)
        cat("所有差异基因GO分析结果已保存，条目数:", nrow(ego_all_df), "\n")
        
        # 绘制点图
        p_all <- dotplot(ego_all, showCategory=15, font.size=10) + 
          ggtitle("GO Enrichment - 所有差异基因")
        ggsave(file.path(go_output_dir, "所有差异基因_GO点图.pdf"), p_all, width = 10, height = 8)
      }
    }
    
    # 5. 创建汇总表
    cat("\n创建通路类型汇总表...\n")
    go_summary <- data.frame(
      分析类型 = character(),
      基因数 = integer(),
      GO通路数 = integer(),
      stringsAsFactors = FALSE
    )
    
    if (exists("up_genes") && nrow(up_genes) > 0) {
      go_summary <- rbind(go_summary, 
                          data.frame(分析类型 = "上调通路", 
                                     基因数 = nrow(up_genes), 
                                     GO通路数 = up_go_count))
    }
    
    if (exists("down_genes") && nrow(down_genes) > 0) {
      go_summary <- rbind(go_summary, 
                          data.frame(分析类型 = "下调通路", 
                                     基因数 = nrow(down_genes), 
                                     GO通路数 = down_go_count))
    }
    
    if (nrow(go_summary) > 0) {
      addWorksheet(wb, "GO_分析汇总")
      writeData(wb, sheet = "GO_分析汇总", go_summary, rowNames = FALSE)
    }
    
    # 6. 添加差异基因列表工作表
    addWorksheet(wb, "差异基因列表")
    writeData(wb, sheet = "差异基因列表", sig_genes, rowNames = FALSE)
    
    # 添加上调基因列表工作表
    if (exists("up_genes") && nrow(up_genes) > 0) {
      addWorksheet(wb, "上调基因列表")
      writeData(wb, sheet = "上调基因列表", up_genes, rowNames = FALSE)
    }
    
    # 添加下调基因列表工作表
    if (exists("down_genes") && nrow(down_genes) > 0) {
      addWorksheet(wb, "下调基因列表")
      writeData(wb, sheet = "下调基因列表", down_genes, rowNames = FALSE)
    }
    
    # 7. 添加分析参数工作表
    params_df <- data.frame(
      参数 = c("P值阈值", "Q值阈值", "p值校正方法", 
             "GO类别", "最小基因集大小", "最大基因集大小", 
             "显著性阈值", "倍数变化阈值", "基因ID类型"),
      值 = c("0.1", "0.3", "BH", "ALL (BP+CC+MF)", "3", "1000", 
            "pvalue < 0.05 或 padj < 0.1", "|log2FC| > 0.5", "SYMBOL")
    )
    addWorksheet(wb, "分析参数")
    writeData(wb, sheet = "分析参数", params_df, rowNames = FALSE)
    
    # 8. 保存Excel文件
    excel_file <- file.path(go_output_dir, "GO_富集分析_结果.xlsx")
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    cat("\nGO富集分析结果已保存到Excel文件:", excel_file, "\n")
    
    # 9. 同时保存CSV文件
    if (exists("ego_up_df") && nrow(ego_up_df) > 0) {
      write.csv(ego_up_df, 
                file.path(go_output_dir, "上调通路_GO富集分析.csv"), 
                row.names = FALSE)
    }
    if (exists("ego_down_df") && nrow(ego_down_df) > 0) {
      write.csv(ego_down_df, 
                file.path(go_output_dir, "下调通路_GO富集分析.csv"), 
                row.names = FALSE)
    }
    write.csv(sig_genes, 
              file.path(go_output_dir, "差异基因列表.csv"), 
              row.names = FALSE)
    
    # 10. 输出最终摘要
    cat("\n" + strrep("=", 60) + "\n")
    cat("GO富集分析最终摘要\n")
    cat(strrep("=", 60) + "\n")
    cat("总差异基因数:", nrow(sig_genes), "\n")
    cat("上调基因数:", ifelse(exists("up_genes"), nrow(up_genes), 0), "\n")
    cat("下调基因数:", ifelse(exists("down_genes"), nrow(down_genes), 0), "\n")
    cat("上调通路数:", up_go_count, "\n")
    cat("下调通路数:", down_go_count, "\n")
    
    # 输出文件列表
    cat("\n生成的文件列表:\n")
    cat("Excel文件:\n")
    cat("  -", excel_file, "\n")
    
    cat("\nPDF图形文件:\n")
    pdf_files <- list.files(go_output_dir, pattern = "\\.pdf$", full.names = FALSE)
    for (file in pdf_files) {
      cat("  -", file, "\n")
    }
    
    cat("\nCSV文件:\n")
    csv_files <- list.files(go_output_dir, pattern = "\\.csv$", full.names = FALSE)
    for (file in csv_files) {
      cat("  -", file, "\n")
    }
    
    cat(strrep("=", 60) + "\n")
    
  } else {
    cat("没有发现显著差异基因，跳过GO富集分析\n")
  }
}

cat(">>> GO富集分析完成！结果文件已保存。\n")