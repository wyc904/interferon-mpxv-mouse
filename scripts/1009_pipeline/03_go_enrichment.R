# ============================================================================
# stage3_go.R — GO BP 富集分析（compareCluster）+ 气泡图（上调红系 / 下调蓝系）
# 输入: _cache/<id>.rds
# 输出: 4-结果/<dataset>/05-GO富集/GO_Bubble_{Up,Down}_<id>.pdf|png
#       4-结果/<dataset>/00-表格/GO富集_<id>.xlsx
#       _cache/plots_<id>_go.rds（供 stage5 合成图）
# ============================================================================
source("scripts/1009_pipeline/utils/common.R", encoding = "UTF-8")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(stringr)
})

open_log("stage3_go")

RED_COLORS  <- c("#FFF0F0", "#FFBBBB", "#FF7777", "#EE3333", "#BB0000", "#880000", "#550000")
BLUE_COLORS <- c("#F0F4FF", "#C5D5F0", "#8BAFDE", "#4B80CA", "#2255AA", "#113388", "#081D5C")

tryCatch({

  run_compare_cluster <- function(gene_list) {
    if (length(gene_list) == 0) return(NULL)
    comp <- tryCatch(
      compareCluster(geneClusters = gene_list, fun = "enrichGO", OrgDb = org.Mm.eg.db,
                     keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2),
      error = function(e) { cat("  [GO] 错误:", conditionMessage(e), "\n"); NULL })
    if (is.null(comp)) return(NULL)
    df <- as.data.frame(comp)
    if (nrow(df) == 0) return(NULL)
    df <- df %>%
      mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        BgRatio_num   = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        RichFactor     = GeneRatio_num,
        FoldEnrichment = GeneRatio_num / BgRatio_num,
        LogP           = -log10(p.adjust))
    df
  }

  make_bubble <- function(df, colors, title) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    sig <- df %>% filter(p.adjust < 0.05)
    if (nrow(sig) == 0) return(NULL)
    top_terms <- sig %>%
      group_by(Cluster) %>% arrange(p.adjust) %>%
      slice_head(n = TOP_GO_PER_CLUSTER) %>% ungroup() %>%
      pull(Description) %>% unique()
    plot_data <- df %>% filter(Description %in% top_terms, p.adjust < 0.05)
    plot_data$Description <- str_trunc(plot_data$Description, width = 58)
    keep_lv <- unname(COMP_LABELS[COMP_ORDER])
    keep_lv <- keep_lv[keep_lv %in% unique(as.character(plot_data$Cluster))]
    plot_data$Cluster <- factor(as.character(plot_data$Cluster), levels = keep_lv)
    ord <- plot_data %>% group_by(Description) %>%
      summarise(first_cl = min(as.integer(Cluster)), best_p = min(p.adjust), .groups = "drop") %>%
      arrange(first_cl, best_p)
    plot_data$Description <- factor(plot_data$Description, levels = ord$Description)

    ggplot(plot_data, aes(x = Cluster, y = Description)) +
      geom_point(aes(size = Count, color = LogP), alpha = 0.85) +
      scale_color_gradientn(colors = colors, name = "-Log10(P.adj)",
        guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = unit(3.5, "cm"), barheight = unit(0.4, "cm"),
                               direction = "horizontal", order = 1)) +
      scale_size_continuous(range = c(2.5, 8), name = "Gene Count",
        guide = guide_legend(title.position = "top", title.hjust = 0.5,
                             direction = "horizontal", nrow = 1, order = 2)) +
      labs(title = title, x = NULL, y = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
            axis.text.y = element_text(size = 9, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.title.position = "plot",
            legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center",
            legend.title = element_text(face = "bold", size = 10),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(12, 12, 6, 12))
  }

  go_plots <- list()

  for (ds in DATASETS) {
    cache <- readRDS(cache_path(ds$id))
    cat("\n", strrep("=", 72), "\n## GO:", ds$label, " [", ds$id, "]\n", strrep("=", 72), "\n", sep = "")

    gene_up <- list(); gene_dn <- list()
    for (cid in COMP_ORDER) {
      r <- cache$results[[cid]]
      up <- unique(r$Geneid[r$change == "Up"]); dn <- unique(r$Geneid[r$change == "Down"])
      if (length(up) > 0) gene_up[[COMP_LABELS[[cid]]]] <- up
      if (length(dn) > 0) gene_dn[[COMP_LABELS[[cid]]]] <- dn
    }
    cat("上调基因集:", length(gene_up), "组 [", paste(sapply(gene_up, length), collapse = ","), "]\n")
    cat("下调基因集:", length(gene_dn), "组 [", paste(sapply(gene_dn, length), collapse = ","), "]\n")

    df_up <- run_compare_cluster(gene_up)
    df_dn <- run_compare_cluster(gene_dn)
    if (!is.null(df_up)) cat("上调 GO terms:", nrow(df_up), " (p.adjust<0.05:", sum(df_up$p.adjust < 0.05), ")\n")
    if (!is.null(df_dn)) cat("下调 GO terms:", nrow(df_dn), " (p.adjust<0.05:", sum(df_dn$p.adjust < 0.05), ")\n")

    p_up <- make_bubble(df_up, RED_COLORS,  paste0("GO BP Enrichment (Up) — ", cache$label))
    p_dn <- make_bubble(df_dn, BLUE_COLORS, paste0("GO BP Enrichment (Down) — ", cache$label))

    d5 <- file.path(ds_res_dir(ds), "05-GO富集"); dir.create(d5, showWarnings = FALSE, recursive = TRUE)
    if (!is.null(p_up)) {
      h <- max(7, length(levels(p_up$data$Description)) * 0.26 + 2.5)
      save_fig(p_up, file.path(d5, paste0("GO_Bubble_Up_", ds$id)), 11.5, h)
      cat("上调气泡图完成 (", length(levels(p_up$data$Description)), " terms )\n")
    }
    if (!is.null(p_dn)) {
      h <- max(7, length(levels(p_dn$data$Description)) * 0.26 + 2.5)
      save_fig(p_dn, file.path(d5, paste0("GO_Bubble_Down_", ds$id)), 11.5, h)
      cat("下调气泡图完成 (", length(levels(p_dn$data$Description)), " terms )\n")
    }

    # ---- Excel ----
    f <- file.path(ds_res_dir(ds), "00-表格", paste0("GO富集_", ds$id, ".xlsx"))
    wb <- createWorkbook()
    addWorksheet(wb, "README")
    desc <- data.frame(Item = c("数据集", "分析", "参数", "对比组", "筛选", "Sheet 说明"),
                       Value = c(ds$label,
                                 "clusterProfiler::compareCluster + enrichGO (BP)",
                                 "OrgDb=org.Mm.eg.db; keyType=SYMBOL; ont=BP; BH; pvalueCutoff=0.05; qvalueCutoff=0.2; 不做 simplify",
                                 paste(unname(COMP_LABELS[COMP_ORDER]), collapse = "; "),
                                 sprintf("每个比较组取 p.adjust 最小 Top %d 进入绘图；表中为全部结果", TOP_GO_PER_CLUSTER),
                                 "GO_Up_All/GO_Down_All=完整结果; UP_*/DN_*=按比较组; Plot_Up/Plot_Down=绘图数据; Summary=统计"))
    writeData(wb, "README", desc); setColWidths(wb, "README", cols = 1:2, widths = c(20, 80))

    if (!is.null(df_up)) {
      addWorksheet(wb, "GO_Up_All")
      writeData(wb, "GO_Up_All", df_up %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
      freezePane(wb, "GO_Up_All", firstRow = TRUE)
      for (cl in unique(df_up$Cluster)) {
        d <- df_up %>% filter(Cluster == cl) %>% arrange(p.adjust)
        sn <- substr(paste0("UP_", gsub(" ", "", cl)), 1, 31)
        addWorksheet(wb, sn); writeData(wb, sn, d %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
      }
    }
    if (!is.null(df_dn)) {
      addWorksheet(wb, "GO_Down_All")
      writeData(wb, "GO_Down_All", df_dn %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
      freezePane(wb, "GO_Down_All", firstRow = TRUE)
      for (cl in unique(df_dn$Cluster)) {
        d <- df_dn %>% filter(Cluster == cl) %>% arrange(p.adjust)
        sn <- substr(paste0("DN_", gsub(" ", "", cl)), 1, 31)
        addWorksheet(wb, sn); writeData(wb, sn, d %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
      }
    }
    if (!is.null(p_up)) { addWorksheet(wb, "Plot_Up");   writeData(wb, "Plot_Up",   p_up$data %>% mutate(across(where(is.numeric), ~ round(.x, 6)))) }
    if (!is.null(p_dn)) { addWorksheet(wb, "Plot_Down"); writeData(wb, "Plot_Down", p_dn$data %>% mutate(across(where(is.numeric), ~ round(.x, 6)))) }
    sm <- bind_rows(
      if (!is.null(df_up)) df_up %>% group_by(Cluster) %>% summarise(GO_Terms = sum(p.adjust < 0.05), .groups = "drop") %>% mutate(Direction = "Up") else NULL,
      if (!is.null(df_dn)) df_dn %>% group_by(Cluster) %>% summarise(GO_Terms = sum(p.adjust < 0.05), .groups = "drop") %>% mutate(Direction = "Down") else NULL)
    if (nrow(sm) > 0) { addWorksheet(wb, "Summary"); writeData(wb, "Summary", sm %>% dplyr::select(Direction, Cluster, GO_Terms)) }
    saveWorkbook(wb, f, overwrite = TRUE)
    cat("Excel 已保存:", f, "\n")

    go_plots[[ds$id]] <- list(up = p_up, dn = p_dn)
  }

  saveRDS(go_plots, file.path(CACHE, "plots_go_all.rds"))
  cat("\nSTAGE3 OK\n")

}, error = function(e) {
  cat("\n!! ERROR:", conditionMessage(e), "\n")
})

close_log()
