# ============================================================================
# stage4b_heatmap_top50.R — Top DEGs 热图（每比较组 Top 8，每图合计 ≤48 个基因）
# 输入: _cache/<id>.rds
# 输出: 4-结果/<dataset>/06-热图/TopDEGs_Heatmap_<id>.pdf|png
#       4-结果/<dataset>/00-表格/热图数据_<id>.xlsx
# 说明: 原 stage4_heatmap.R 为每比较组 Top50 并集（148–188 个基因，21–26 in 高）；
#       本脚本改为每比较组 Top 8（去重后 ≤48），其余流程（VST → 行 Z-score ±3、
#       RdBu 色阶、按比较组分块、列注释仅 Group）与 stage4 一致；
#       另：画布尺寸改为按 gtable 自然尺寸自适应（原公式为长图设计，短图会留白）。
# ============================================================================
source("scripts/1009_pipeline/utils/common.R", encoding = "UTF-8")
suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

TOP_HEAT_PER_COMP <- 8   # 每比较组 padj 最小 Top N（6 组合计、去重后每图 ≤48）

RUN_IDS <- commandArgs(trailingOnly = TRUE)
if (length(RUN_IDS) == 0) RUN_IDS <- names(DATASETS)
RUN_DS <- DATASETS[RUN_IDS]
stopifnot(length(RUN_DS) > 0)

open_log("stage4b_heatmap_top50")

tryCatch({

  for (ds in RUN_DS) {
    cache <- readRDS(cache_path(ds$id))
    cat("\n", strrep("=", 72), "\n## 热图:", ds$label, " [", ds$id, "]\n", strrep("=", 72), "\n", sep = "")

    # ---- 每比较组 Top N 显著基因（并集）----
    sig <- bind_rows(lapply(COMP_ORDER, function(cid) {
      cache$results[[cid]] %>%
        filter(change != "NS") %>% arrange(padj) %>% head(TOP_HEAT_PER_COMP) %>%
        mutate(src = cid)
    }))
    sel_genes <- unique(sig$Geneid)

    # 每个基因的“来源比较组”：显著性最优的比较组
    best <- bind_rows(lapply(COMP_ORDER, function(cid) cache$results[[cid]] %>% filter(change != "NS"))) %>%
      arrange(padj) %>% distinct(Geneid, .keep_all = TRUE) %>%
      select(Geneid, comparison, baseMean, log2FoldChange, pvalue, padj)
    info <- best[match(sel_genes, best$Geneid), ]

    # 行排序：按来源比较组顺序 + padj
    ord <- order(match(info$comparison, COMP_ORDER), info$padj)
    sel_genes <- sel_genes[ord]; info <- info[ord, ]
    cat("热图基因数:", length(sel_genes), " (每比较组 Top", TOP_HEAT_PER_COMP, " 并集)\n", sep = "")
    print(table(factor(info$comparison, levels = COMP_ORDER)))

    # ---- 表达矩阵 ----
    expr <- cache$vst[sel_genes, , drop = FALSE]
    si <- cache$sample_info
    expr <- expr[, si$col, drop = FALSE]
    colnames(expr) <- si$label_short
    expr_z <- t(scale(t(expr)))
    expr_z[!is.finite(expr_z)] <- 0
    expr_z[expr_z >  3] <-  3
    expr_z[expr_z < -3] <- -3

    # ---- 注释 ----
    ann_col <- data.frame(Group = si$group, row.names = si$label_short)
    ann_colors <- list(Group = GROUP_COLORS[GROUP_LEVELS])
    sz_col <- as.integer(table(si$group)); sz_col <- sz_col[sz_col > 0]
    gcol <- cumsum(sz_col)[-length(sz_col)]
    sz_row <- as.integer(table(factor(info$comparison, levels = COMP_ORDER))); sz_row <- sz_row[sz_row > 0]
    grow <- cumsum(sz_row)[-length(sz_row)]

    ph <- pheatmap(
      expr_z,
      color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
      annotation_col = ann_col, annotation_colors = ann_colors,
      cluster_rows = FALSE, cluster_cols = FALSE,
      show_rownames = TRUE, show_colnames = TRUE,
      fontsize_row = if (nrow(expr_z) > 220) 4.2 else if (nrow(expr_z) > 140) 5 else 6.5,
      fontsize_col = 8, fontsize = 10,
      border_color = NA,
      main = sprintf("Top %d DEGs — %s", nrow(expr_z), ds$label),
      gaps_col = if (length(gcol) > 0) gcol else NULL,
      gaps_row = if (length(grow) > 0) grow else NULL,
      cellwidth = 13, cellheight = 8, treeheight_row = 0, treeheight_col = 0,
      angle_col = 45,
      legend_breaks = seq(-3, 3, 1), legend_labels = as.character(seq(-3, 3, 1)),
      silent = TRUE)

    # ---- 保存 ----
    d6 <- file.path(ds_res_dir(ds), "06-热图"); dir.create(d6, showWarnings = FALSE, recursive = TRUE)
    # 画布按内容紧凑自适应（stage4 原公式为 150+ 行设计，行数骤减后会留大片空白）
    w <- max(4.5, grid::convertWidth(gtable::gtable_width(ph$gtable), "in", TRUE) + 0.25)
    h <- max(3.5, grid::convertHeight(gtable::gtable_height(ph$gtable), "in", TRUE) + 0.25)
    base <- file.path(d6, paste0("TopDEGs_Heatmap_", ds$id))
    pdf(paste0(base, ".pdf"), width = w, height = h, useDingbats = FALSE)
    grid.newpage(); grid.draw(ph$gtable); dev.off()
    ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 300, background = "white")
    grid.newpage(); grid.draw(ph$gtable); dev.off()
    cat(sprintf("热图已保存: %.1f x %.1f in\n", w, h))

    # ---- Excel ----
    f <- file.path(ds_res_dir(ds), "00-表格", paste0("热图数据_", ds$id, ".xlsx"))
    wb <- createWorkbook()
    addWorksheet(wb, "README")
    desc <- data.frame(Item = c("数据集", "选基因", "标准化", "色阶", "Sheet 说明"),
                       Value = c(ds$label,
                                 sprintf("每个比较组 padj 最小 Top %d 的显著 DEGs 并集，共 %d 个", TOP_HEAT_PER_COMP, length(sel_genes)),
                                 "DESeq2 VST -> 行 Z-score（截断 ±3）",
                                 "RdBu 反向（蓝-白-红）",
                                 "Zscore=绘图矩阵; VST=VST表达值; DEGinfo=每个基因最优比较组的统计量"))
    writeData(wb, "README", desc); setColWidths(wb, "README", cols = 1:2, widths = c(16, 80))

    z_df <- as.data.frame(expr_z) %>% rownames_to_column("Geneid") %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
    addWorksheet(wb, "Zscore"); writeData(wb, "Zscore", z_df); freezePane(wb, "Zscore", firstRow = TRUE)
    raw_df <- as.data.frame(expr) %>% rownames_to_column("Geneid") %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
    addWorksheet(wb, "VST"); writeData(wb, "VST", raw_df); freezePane(wb, "VST", firstRow = TRUE)
    info_df <- info %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
    addWorksheet(wb, "DEGinfo"); writeData(wb, "DEGinfo", info_df); freezePane(wb, "DEGinfo", firstRow = TRUE)
    saveWorkbook(wb, f, overwrite = TRUE)
    cat("Excel 已保存:", f, "\n")
  }

  cat("\nSTAGE4B OK\n")

}, error = function(e) {
  cat("\n!! ERROR:", conditionMessage(e), "\n")
})

close_log()
