# 下调基因
down_genes_tbl <- gene_list %>% dplyr::filter(log2FoldChange < 0)
down_ids <- unique(na.omit(down_genes_tbl$ENTREZID))
cat("Down genes (ENTREZ) 数量:", length(down_ids), "\n")

# 输出文件（放到output_dir）
if(!exists("output_dir")) output_dir <- getwd()
excel_file <- file.path(output_dir, "GO_Downregulated_p_down_data.xlsx")

wb <- createWorkbook()

# 1) 写入参与分析的Down基因（无论是否富集成功都写）
addWorksheet(wb, "Down_Genes_Used")
writeData(wb, "Down_Genes_Used",
          down_genes_tbl %>% dplyr::select(gene, ENTREZID, log2FoldChange, padj))

# 如果Down基因太少，直接输出并结束
if (length(down_ids) < 3) {
  addWorksheet(wb, "Down_GO_BP_NOTE")
  writeData(wb, "Down_GO_BP_NOTE",
            data.frame(Note = paste0("Down基因ENTREZ数 < 3 (", length(down_ids), "), 无法进行富集")))
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat("✅ 已输出Excel(无富集):", excel_file, "\n")
} else {
  
  # 2) 尝试不同阈值，尽量让 down 能有结果
  cutoffs <- c(0.05, 0.1, 0.2)
  ego_down <- NULL
  used_cutoff <- NA
  
  for (pcut in cutoffs) {
    tmp <- tryCatch({
      enrichGO(
        gene          = down_ids,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = pcut,
        qvalueCutoff  = 0.2,
        minGSSize     = 3,      # 放宽一点，更容易出结果
        maxGSSize     = 1000,
        readable      = TRUE
      )
    }, error = function(e) NULL)
    
    if (!is.null(tmp) && nrow(as.data.frame(tmp)) > 0) {
      ego_down <- tmp
      used_cutoff <- pcut
      break
    }
  }
  
  # 3) 无论是否显著，都导出“底层表”
  addWorksheet(wb, "Params")
  writeData(wb, "Params", data.frame(
    Parameter = c("ont","pAdjustMethod","qvalueCutoff","minGSSize","maxGSSize","pvalueCutoff_used"),
    Value     = c("BP","BH",0.2,3,1000, ifelse(is.na(used_cutoff), "NONE(sig=0)", used_cutoff))
  ))
  
  if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
    
    # 有显著结果：导出表 + 画p_down
    down_df <- as.data.frame(ego_down) %>%
      dplyr::mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        BgRatio_num   = sapply(strsplit(BgRatio,   "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        EnrichmentRatio = GeneRatio_num / BgRatio_num
      ) %>% dplyr::arrange(p.adjust)
    
    addWorksheet(wb, "Down_GO_BP_sig")
    writeData(wb, "Down_GO_BP_sig", down_df)
    
    # p_down绘图（showCategory=10）
    p_down <- dotplot(ego_down, showCategory = 10) +
      ggtitle(paste0("GO Enrichment - Downregulated Genes (BP), pCut=", used_cutoff))
    ggsave(file.path(output_dir, "GO_Downregulated_dotplot.pdf"),
           p_down, width = 10, height = 6)
    
    cat("✅ Down BP 有结果，已绘制p_down并写入Excel\n")
    
  } else {
    # 仍然无显著：输出“原始富集表”（pvalueCutoff=1），用于排查/后续自己筛
    ego_down_raw <- enrichGO(
      gene          = down_ids,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      minGSSize     = 3,
      maxGSSize     = 1000,
      readable      = TRUE
    )
    
    raw_df <- as.data.frame(ego_down_raw)
    addWorksheet(wb, "Down_GO_BP_raw_all")
    if (nrow(raw_df) > 0) {
      raw_df <- raw_df %>%
        dplyr::mutate(
          GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
          BgRatio_num   = sapply(strsplit(BgRatio,   "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
          EnrichmentRatio = GeneRatio_num / BgRatio_num
        ) %>% dplyr::arrange(p.adjust)
      writeData(wb, "Down_GO_BP_raw_all", raw_df)
    } else {
      writeData(wb, "Down_GO_BP_raw_all",
                data.frame(Note="即使不做显著性筛选(pCut=1,qCut=1)，也没有任何GO条目返回；请检查ENTREZID是否属于Mm.eg.db或是否为空/异常"))
    }
    
    cat("⚠️ Down BP 在0.05/0.1/0.2都无显著，已导出raw_all用于排查\n")
  }
  
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat("✅ 已输出Excel:", excel_file, "\n")
}
