# ============================================================================
# stage1_deseq2.R — 按数据集(品系 × 组织)运行 DESeq2，提取 6 个比较结果
# 输出:
#   _cache/<id>.rds                    每个数据集的 vst 矩阵 / 结果 / 样本信息
#   4-结果/<dataset>/00-表格/DEG_全部比较_<id>.xlsx
#   4-结果/00-汇总/汇总_DEG数量_四套.xlsx
# ============================================================================
source("scripts/1009_pipeline/utils/common.R", encoding = "UTF-8")

open_log("stage1_deseq2")

tryCatch({

  cnt <- read_counts_matrix()
  ss  <- read_sample_sheet()
  gene_meta_all <- attr(cnt, "gene_meta")
  cat("计数矩阵:", nrow(cnt), "基因 ×", ncol(cnt), "样本\n")

  # 校验样本表与矩阵列一一对应
  stopifnot(setequal(ss$col, colnames(cnt)))
  stopifnot(nrow(ss) == ncol(cnt))

  summary_rows <- list()

  for (ds in DATASETS) {
    cat("\n", strrep("=", 72), "\n## 数据集: ", ds$label, "  [", ds$id, "]\n", strrep("=", 72), "\n", sep = "")

    si <- ss[ss$strain_c == ds$strain & ss$tissue == ds$tissue, ]
    si$group <- factor(si$group, levels = GROUP_LEVELS)
    si <- si[order(si$group, si$rep), ]
    rownames(si) <- si$col
    cat("样本数:", nrow(si), " | 各组:", paste(names(table(si$group)), table(si$group), sep = "=", collapse = ", "), "\n")

    sc <- cnt[, si$col, drop = FALSE]

    # 低表达过滤: count >= MIN_COUNT 的样本数 >= 最小分组样本数
    min_grp <- min(table(si$group))
    keep <- rowSums(sc >= MIN_COUNT) >= min_grp
    sc_f <- sc[keep, , drop = FALSE]
    cat("过滤: 保留", sum(keep), "/", nrow(sc), "基因 (最小分组 n =", min_grp, ")\n")

    dds <- DESeqDataSetFromMatrix(countData = sc_f, colData = si, design = ~ group)
    dds <- DESeq(dds)
    vsd <- vst(dds, blind = FALSE)

    res_l <- list()
    for (cmp in COMPARISONS) {
      r  <- results(dds, contrast = c("group", cmp$trt, cmp$ref), alpha = PADJ_CUT)
      rdf <- as.data.frame(r) %>%
        rownames_to_column("Geneid") %>%
        mutate(comparison = cmp$id, comparison_label = COMP_LABELS[[cmp$id]]) %>%
        classify_change() %>%
        arrange(padj) %>%
        select(Geneid, comparison, comparison_label, change, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
      res_l[[cmp$id]] <- rdf
      n_up <- sum(rdf$change == "Up"); n_dn <- sum(rdf$change == "Down")
      cat(sprintf("  %-8s  Up: %5d   Down: %5d   合计: %5d\n", cmp$id, n_up, n_dn, n_up + n_dn))
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        dataset = ds$id, dataset_label = ds$label, comparison = cmp$id,
        comparison_label = COMP_LABELS[[cmp$id]], Up = n_up, Down = n_dn, Total = n_up + n_dn)
    }

    # ---- 缓存 ----
    cache <- list(
      id = ds$id, label = ds$label, dir = ds$dir,
      vst = assay(vsd),
      gene_meta = gene_meta_all[keep, , drop = FALSE],
      results = res_l,
      sample_info = si[, c("col", "label_short", "group", "rep", "run_id", "seq_name")],
      filter = list(min_group = min_grp, kept = sum(keep), total = nrow(sc))
    )
    saveRDS(cache, cache_path(ds$id))
    cat("缓存已保存:", cache_path(ds$id), "\n")

    # ---- 每个数据集的 DEG 明细 Excel ----
    out_dir <- file.path(ds_res_dir(ds), "00-表格")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    f <- file.path(out_dir, paste0("DEG_全部比较_", ds$id, ".xlsx"))
    wb <- createWorkbook()

    addWorksheet(wb, "README")
    desc <- data.frame(Item = c("数据集", "品系", "组织", "样本数", "分组(样本数)",
                                "DESeq2 设计", "低表达过滤", "差异阈值", "比较组数", "生成时间",
                                "Sheet 说明"),
                       Value = c(ds$label, ds$strain, ds$tissue, nrow(si),
                                 paste(names(table(si$group)), table(si$group), sep = "=", collapse = ", "),
                                 "~ group (G1 为参照水平)",
                                 sprintf("count >= %d 的样本数 >= %d; 保留 %d/%d 基因", MIN_COUNT, min_grp, sum(keep), nrow(sc)),
                                 sprintf("padj < %g 且 |log2FC| > %g", PADJ_CUT, LFC_CUT),
                                 length(COMPARISONS), format(Sys.time(), "%Y-%m-%d %H:%M"),
                                 "All_Sig_DEGs=全部显著DEG; G2vsG1 等=按比较组; DEG_Summary=汇总; VST_*=归一化表达"))
    writeData(wb, "README", desc)
    setColWidths(wb, "README", cols = 1:2, widths = c(24, 70))

    all_sig <- bind_rows(lapply(res_l, function(r) r[r$change != "NS", ])) %>%
      mutate(across(where(is.numeric), ~ round(.x, 6)))
    addWorksheet(wb, "All_Sig_DEGs"); writeData(wb, "All_Sig_DEGs", all_sig)
    freezePane(wb, "All_Sig_DEGs", firstRow = TRUE)

    for (cid in COMP_ORDER) {
      d <- res_l[[cid]] %>% mutate(across(where(is.numeric), ~ round(.x, 6)))
      addWorksheet(wb, cid); writeData(wb, cid, d); freezePane(wb, cid, firstRow = TRUE)
    }

    sm <- bind_rows(lapply(res_l, function(r) data.frame(
      comparison = r$comparison[1], comparison_label = r$comparison_label[1],
      Up = sum(r$change == "Up"), Down = sum(r$change == "Down"),
      Total_Sig = sum(r$change != "NS"),
      Total_Tested = nrow(r),
      NotSig = sum(r$change == "NS"),
      NA_padj = sum(is.na(r$padj)))))
    addWorksheet(wb, "DEG_Summary"); writeData(wb, "DEG_Summary", sm)

    vst_df <- as.data.frame(assay(vsd)) %>%
      rownames_to_column("Geneid") %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)))
    addWorksheet(wb, "VST_Expression_AllGenes"); writeData(wb, "VST_Expression_AllGenes", vst_df)
    freezePane(wb, "VST_Expression_AllGenes", firstRow = TRUE)

    deg_genes <- unique(all_sig$Geneid)
    vst_deg <- vst_df[vst_df$Geneid %in% deg_genes, ]
    addWorksheet(wb, "VST_Expression_DEGs"); writeData(wb, "VST_Expression_DEGs", vst_deg)
    freezePane(wb, "VST_Expression_DEGs", firstRow = TRUE)

    saveWorkbook(wb, f, overwrite = TRUE)
    cat("Excel 已保存:", f, "\n")
  }

  # ---- 四套汇总 ----
  sm_all <- bind_rows(summary_rows)
  f2 <- file.path(SUB_DIR, "汇总_DEG数量_四套.xlsx")
  wb2 <- createWorkbook()
  addWorksheet(wb2, "DEG数量")
  writeData(wb2, "DEG数量", sm_all)
  freezePane(wb2, "DEG数量", firstRow = TRUE)
  addWorksheet(wb2, "分组结构")
  gs <- do.call(rbind, lapply(names(DATASETS), function(id) {
    cache <- readRDS(cache_path(id))
    tb <- table(cache$sample_info$group)
    data.frame(dataset = id, label = cache$label,
               G1 = as.integer(tb["G1"]), G2 = as.integer(tb["G2"]),
               G3 = as.integer(tb["G3"]), G4 = as.integer(tb["G4"]),
               kept_genes = cache$filter$kept)
  }))
  writeData(wb2, "分组结构", gs)
  saveWorkbook(wb2, f2, overwrite = TRUE)
  cat("\n汇总 Excel:", f2, "\n")
  cat("\nSTAGE1 OK\n")

}, error = function(e) {
  cat("\n!! ERROR:", conditionMessage(e), "\n")
})

close_log()
