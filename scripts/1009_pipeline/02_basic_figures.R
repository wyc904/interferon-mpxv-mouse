# ============================================================================
# stage2_figures_basic.R — 四套基础图：PCA / DEG数量柱状图 / 韦恩图 / 火山图
# 输入: _cache/<id>.rds
# 输出: 4-结果/<dataset>/{01-PCA,02-DEG数量,03-韦恩图,04-火山图}/*.pdf|png
#       _cache/plots_<id>.rds  (供 stage5 合成图使用)
# ============================================================================
source("scripts/1009_pipeline/utils/common.R", encoding = "UTF-8")
suppressPackageStartupMessages({
  library(ggrepel)
  library(patchwork)
  library(ggVennDiagram)
})

open_log("stage2_figures_basic")

tryCatch({

  pca_data <- function(cache) {
    vst <- cache$vst
    v   <- apply(vst, 1, var)
    top <- names(sort(v, decreasing = TRUE))[1:min(500, length(v))]
    pc  <- prcomp(t(vst[top, , drop = FALSE]))
    pct <- round(100 * (pc$sdev^2 / sum(pc$sdev^2))[1:2], 1)
    si  <- cache$sample_info
    list(df = data.frame(sample = si$label_short, group = si$group,
                         PC1 = round(pc$x[, 1], 3), PC2 = round(pc$x[, 2], 3)),
         pct = pct, ntop = length(top))
  }

  pca_plot <- function(cache) {
    pd <- pca_data(cache)
    df <- pd$df
    pct <- pd$pct
    ggplot(df, aes(PC1, PC2, color = group, shape = group)) +
      geom_point(size = 4, alpha = 0.9) +
      geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, max.overlaps = 30, seed = 1) +
      scale_color_manual(values = GROUP_COLORS) +
      scale_shape_manual(values = c(G1 = 16, G2 = 17, G3 = 15, G4 = 18)) +
      labs(x = paste0("PC1: ", pct[1], "% variance"),
           y = paste0("PC2: ", pct[2], "% variance"),
           title = paste0("PCA — ", cache$label),
           subtitle = "G1 uninfected · G2 infected · G3/G4 treatment",
           color = "Group", shape = "Group") +
      theme_pub(13)
  }

  bar_data <- function(cache) {
    sm <- bind_rows(lapply(COMP_ORDER, function(cid) {
      r <- cache$results[[cid]]
      data.frame(comparison = cid, comparison_label = COMP_LABELS[[cid]],
                 direction = c("Up", "Down"),
                 n = c(sum(r$change == "Up"), sum(r$change == "Down")))
    }))
    sm$comparison_label <- factor(sm$comparison_label, levels = unname(COMP_LABELS[COMP_ORDER]))
    sm$direction <- factor(sm$direction, levels = c("Up", "Down"))
    sm
  }

  bar_plot <- function(cache) {
    sm <- bar_data(cache)
    ggplot(sm, aes(x = comparison_label, y = n, fill = direction)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      geom_text(aes(label = n), position = position_dodge(width = 0.85), vjust = -0.3, size = 3.2) +
      scale_fill_manual(values = c(Up = UP_COL, Down = DOWN_COL)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(x = NULL, y = "Number of Differentially Expressed Genes", fill = NULL,
           title = paste0("DEG Counts — ", cache$label),
           subtitle = sprintf("padj < %g & |log2FC| > %g", PADJ_CUT, LFC_CUT)) +
      theme_pub(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
            legend.position = "top")
  }

  venn_sets <- list(
    G3 = list(cids = c("G2vsG1", "G3vsG1", "G3vsG2"),
              labels = c("G2 vs G1", "G3 vs G1", "G3 vs G2")),
    G4 = list(cids = c("G2vsG1", "G4vsG1", "G4vsG2"),
              labels = c("G2 vs G1", "G4 vs G1", "G4 vs G2"))
  )

  make_venn <- function(cache, trt, mode = c("All", "Up", "Down"), fill_high = "#4C72B0") {
    mode <- match.arg(mode)
    s <- venn_sets[[trt]]
    lst <- setNames(lapply(s$cids, function(cid) {
      r <- cache$results[[cid]]
      g <- if (mode == "All") r$Geneid[r$change != "NS"] else r$Geneid[r$change == mode]
      unique(g)
    }), s$labels)
    if (sum(lengths(lst)) == 0) return(NULL)
    ggVennDiagram(lst, label = "count", label_alpha = 0, edge_size = 0.8) +
      scale_fill_gradient(low = "#F7F7F7", high = fill_high) +
      scale_color_manual(values = c("#A52A2A", "#1B9E77", "#7570B3")) +
      labs(title = paste0("DEG overlap (", mode, ") — ", cache$label),
           subtitle = paste(s$labels, collapse = "  |  ")) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
            legend.position = "none",
            plot.margin = margin(10, 10, 10, 10))
  }

  venn_region_df <- function(cache, trt, mode = c("All", "Up", "Down")) {
    mode <- match.arg(mode)
    s <- venn_sets[[trt]]
    lst <- setNames(lapply(s$cids, function(cid) {
      r <- cache$results[[cid]]
      g <- if (mode == "All") r$Geneid[r$change != "NS"] else r$Geneid[r$change == mode]
      unique(g)
    }), s$labels)
    allg <- unique(unlist(lst))
    if (length(allg) == 0) return(NULL)
    member <- sapply(lst, function(g) allg %in% g)
    key <- apply(member, 1, function(x) paste(which(x), collapse = "&"))
    df <- data.frame(Geneid = allg, key = key, stringsAsFactors = FALSE)
    agg <- df %>% group_by(key) %>%
      summarise(Count = n(), Genes = paste(Geneid, collapse = ", "), .groups = "drop")
    agg$Region <- sapply(strsplit(agg$key, "&"), function(ix) paste(s$labels[as.integer(ix)], collapse = " & "))
    agg %>% arrange(desc(Count)) %>% select(Region, Count, Genes)
  }

  volcano_plot <- function(cache, cid) {
    r <- cache$results[[cid]]
    r$neglogp <- -log10(r$padj)
    ycap <- quantile(r$neglogp, 0.995, na.rm = TRUE) * 1.05
    r$neglogp_c <- pmin(r$neglogp, ycap)
    top <- r %>% filter(change != "NS") %>% arrange(padj) %>% head(15) %>% pull(Geneid)
    ggplot(r, aes(log2FoldChange, neglogp_c)) +
      geom_point(aes(color = change), alpha = 0.65, size = 1.1) +
      scale_color_manual(values = c(Up = UP_COL, Down = DOWN_COL, NS = NS_COL)) +
      geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", color = "grey45", linewidth = 0.4) +
      geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", color = "grey45", linewidth = 0.4) +
      geom_text_repel(data = r[r$Geneid %in% top, ], aes(label = Geneid, color = change),
                      size = 2.6, show.legend = FALSE, max.overlaps = 20, seed = 2,
                      min.segment.length = 0, segment.size = 0.3) +
      labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~Adjusted~P-value),
           title = COMP_LABELS[[cid]],
           subtitle = sprintf("Up: %d · Down: %d", sum(r$change == "Up"), sum(r$change == "Down"))) +
      theme_pub(10) +
      theme(legend.position = "none")
  }

  bar_all <- list()

  for (ds in DATASETS) {
    cache <- readRDS(cache_path(ds$id))
    outdir <- ds_res_dir(ds)
    cat("\n## ", ds$label, " [", ds$id, "]\n", sep = "")

    # ---- PCA ----
    d1 <- file.path(outdir, "01-PCA"); dir.create(d1, showWarnings = FALSE, recursive = TRUE)
    p_pca <- pca_plot(cache)
    save_fig(p_pca, file.path(d1, paste0("PCA_", ds$id)), 6.5, 5.5)
    cat("  PCA 完成\n")

    # ---- DEG 数量 ----
    d2 <- file.path(outdir, "02-DEG数量"); dir.create(d2, showWarnings = FALSE, recursive = TRUE)
    p_bar <- bar_plot(cache)
    save_fig(p_bar, file.path(d2, paste0("DEGcount_", ds$id)), 7.5, 5.5)
    bar_all[[ds$id]] <- bar_data(cache) %>% mutate(dataset = ds$id, label = ds$label)
    cat("  DEG 数量柱状图完成\n")

    # ---- 韦恩图 ----
    d3 <- file.path(outdir, "03-韦恩图"); dir.create(d3, showWarnings = FALSE, recursive = TRUE)
    venn_plots <- list()
    for (trt in c("G3", "G4")) {
      for (mode in c("All", "Up", "Down")) {
        fill_high <- switch(mode, All = "#4C72B0", Up = "#C0392B", Down = "#21618C")
        pv <- make_venn(cache, trt, mode, fill_high)
        if (!is.null(pv)) {
          save_fig(pv, file.path(d3, paste0("Venn_", trt, "_", mode, "_", ds$id)), 6, 6)
          if (mode == "All") venn_plots[[trt]] <- pv
        }
      }
    }
    cat("  韦恩图完成 (", length(venn_plots), " 个 All 版本 )\n", sep = "")

    # ---- 火山图 ----
    d4 <- file.path(outdir, "04-火山图"); dir.create(d4, showWarnings = FALSE, recursive = TRUE)
    pl <- lapply(COMP_ORDER, function(cid) volcano_plot(cache, cid))
    for (i in seq_along(COMP_ORDER)) {
      ggsave(file.path(d4, paste0("Volcano_", COMP_ORDER[i], "_", ds$id, ".pdf")), pl[[i]],
             width = 4.3, height = 4.5, device = "pdf", useDingbats = FALSE)
    }
    volc_grid <- (pl[[1]] | pl[[2]] | pl[[3]]) / (pl[[4]] | pl[[5]] | pl[[6]]) +
      plot_annotation(title = paste0("Volcano — ", cache$label),
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15)))
    save_fig(volc_grid, file.path(d4, paste0("Volcano_Grid_", ds$id)), 12.5, 8.6)
    cat("  火山图完成 (6 单图 + 1 网格)\n")

    # ---- 供 stage5 合成图用 ----
    saveRDS(list(pca = p_pca, bar = p_bar, venn_g3 = venn_plots[["G3"]], venn_g4 = venn_plots[["G4"]],
                 volc_grid = volc_grid),
            file.path(CACHE, paste0("plots_", ds$id, ".rds")))

    # ---- 绘图依据 Excel（PCA / DEG数量 / 韦恩） ----
    f <- file.path(outdir, "00-表格", paste0("绘图依据_PCA-DEG数量-韦恩_", ds$id, ".xlsx"))
    wb <- createWorkbook()
    addWorksheet(wb, "README")
    desc <- data.frame(
      Item = c("说明", "PCA坐标", "DEG数量", "韦恩", "火山图", "其他"),
      Value = c(paste0(ds$label, " — 各图绘图依据数据"),
                "VST 方差最大 500 基因的 PC1/PC2 坐标（VST 表达值见 DEG_全部比较_*.xlsx 的 VST sheet）",
                "6 个比较组的上/下调显著基因数（padj<0.05 且 |log2FC|>1）",
                "每个治疗组 × All/Up/Down：三个集合重叠区域（分区）的计数与基因列表",
                "见 DEG_全部比较_*.xlsx 各比较组 sheet 与 All_Sig_DEGs（即为火山图数据）",
                "热图见 热图数据_*.xlsx；GO 见 GO富集_*.xlsx 的 Plot_Up/Plot_Down"))
    writeData(wb, "README", desc)
    setColWidths(wb, "README", cols = 1:2, widths = c(14, 95))

    pd <- pca_data(cache)
    addWorksheet(wb, "PCA坐标")
    writeData(wb, "PCA坐标", pd$df)
    writeData(wb, "PCA坐标", t(data.frame(说明 = sprintf("PC1 方差占比 %.1f%% | PC2 %.1f%% | 使用基因数 %d",
                                                            pd$pct[1], pd$pct[2], pd$ntop))),
              startRow = nrow(pd$df) + 3, colNames = FALSE)

    bd <- bind_rows(lapply(COMP_ORDER, function(cid) {
      r <- cache$results[[cid]]
      data.frame(comparison = COMP_LABELS[[cid]], Up = sum(r$change == "Up"),
                 Down = sum(r$change == "Down"), Total_Sig = sum(r$change != "NS"))
    }))
    addWorksheet(wb, "DEG数量")
    writeData(wb, "DEG数量", bd)
    freezePane(wb, "DEG数量", firstRow = TRUE)

    for (trt in c("G3", "G4")) {
      for (mode in c("All", "Up", "Down")) {
        vd <- venn_region_df(cache, trt, mode)
        sn <- paste0("Venn_", trt, "_", mode)
        addWorksheet(wb, sn)
        if (is.null(vd)) vd <- data.frame(Region = "（该方向无显著基因）", Count = 0, Genes = "")
        writeData(wb, sn, vd)
        setColWidths(wb, sn, cols = 1:3, widths = c(40, 8, 120))
      }
    }
    saveWorkbook(wb, f, overwrite = TRUE)
    cat("绘图依据 Excel 已保存:", f, "\n")
  }

  # ---- 四套汇总 DEG 数量图 ----
  bar_df <- bind_rows(bar_all)
  bar_df$label <- factor(bar_df$label, levels = sapply(DATASETS, `[[`, "label"))
  p_all <- ggplot(bar_df, aes(x = comparison_label, y = n, fill = direction)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_text(aes(label = n), position = position_dodge(width = 0.85), vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = c(Up = UP_COL, Down = DOWN_COL)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    facet_wrap(~ label, scales = "free_y", nrow = 2) +
    labs(x = NULL, y = "Number of DEGs", fill = NULL,
         title = "Differentially Expressed Genes — all four datasets",
         subtitle = sprintf("DESeq2, padj < %g & |log2FC| > %g", PADJ_CUT, LFC_CUT)) +
    theme_pub(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
          legend.position = "top")
  save_fig(p_all, file.path(SUB_DIR, "DEG数量_四套汇总"), 14, 9)
  cat("\n四套汇总 DEG 数量图完成\nSTAGE2 OK\n")

}, error = function(e) {
  cat("\n!! ERROR:", conditionMessage(e), "\n")
})

close_log()
