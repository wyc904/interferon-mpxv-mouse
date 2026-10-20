# ============================================================================
# stage3b_go_top30.R — GO BP 富集气泡图（每比较组 Top 5，每图合计 ≤30 条通路）
# 输入: _cache/<id>.rds（复用 stage1 的 DESeq2 结果，不重跑差异分析）
# 输出: 4-结果/<dataset>/05-GO富集/GO_Bubble_{Up,Down}_<id>.pdf|png
#       4-结果/00-汇总/GO富集_Top30_汇总.xlsx
# 用法: Rscript stage3b_go_top30.R [数据集id...]   不带参数=全部四套
# ============================================================================
source("scripts/1009_pipeline/utils/common.R", encoding = "UTF-8")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(stringr)
  library(writexl)   # 汇总表用 writexl 写：openxlsx 输出含悬空 drawing 关系，openpyxl 读不了
})

TOP_N <- 5    # 每个比较组取 Top N 条通路（6 组合计、去重后每图 ≤ 30 条）

RED_COLORS  <- c("#FFF0F0", "#FFBBBB", "#FF7777", "#EE3333", "#BB0000", "#880000", "#550000")
BLUE_COLORS <- c("#F0F4FF", "#C5D5F0", "#8BAFDE", "#4B80CA", "#2255AA", "#113388", "#081D5C")

RUN_IDS <- commandArgs(trailingOnly = TRUE)
if (length(RUN_IDS) == 0) RUN_IDS <- names(DATASETS)
RUN_DS <- DATASETS[RUN_IDS]
stopifnot(length(RUN_DS) > 0)

open_log("stage3b_go_top30")

cat("== stage3b: GO 气泡图 每组 Top", TOP_N, "（每图合计 ≤", TOP_N * length(COMP_ORDER), "条）==\n")
cat("数据集:", paste(sapply(RUN_DS, `[[`, "id"), collapse = ", "), "\n\n")

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
    # 每个 Cluster 按 p.adjust 取前 TOP_N 条
    sel <- sig %>%
      group_by(Cluster) %>% arrange(p.adjust) %>%
      slice_head(n = TOP_N) %>% ungroup()
    # Y 轴顺序: 首次出现的比较组顺序 + 最优 p 值
    ord <- data.frame(Description = unique(sel$Description), stringsAsFactors = FALSE)
    ord$first_cl <- sapply(ord$Description, function(t) {
      cls <- unique(as.character(sel$Cluster[sel$Description == t]))
      min(vapply(cls, function(c) match(c, unname(COMP_LABELS[COMP_ORDER])), numeric(1)))
    })
    ord$best_p <- sapply(ord$Description, function(t) min(sel$p.adjust[sel$Description == t]))
    ord <- ord %>% arrange(first_cl, best_p)
    # 绘图数据（完整 Description 优先用于 Excel，截断仅用于图上的因子标签）
    plot_data <- sel
    plot_data$Description_tr <- str_trunc(plot_data$Description, width = 58)
    keep_lv <- unname(COMP_LABELS[COMP_ORDER])
    keep_lv <- keep_lv[keep_lv %in% unique(as.character(plot_data$Cluster))]
    plot_data$Cluster <- factor(as.character(plot_data$Cluster), levels = keep_lv)
    plot_data$Description <- factor(plot_data$Description_tr,
                                    levels = str_trunc(ord$Description, width = 58))

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

  sheets <- list()
  sheets[["README"]] <- data.frame(
    Item = c("说明", "分析", "参数", "筛选", "用法"),
    Value = c(paste("GO BP 富集气泡图（每比较组 Top", TOP_N, "，每图合计 ≤30 条）"),
              "clusterProfiler::compareCluster + enrichGO (BP)",
              "OrgDb=org.Mm.eg.db; keyType=SYMBOL; ont=BP; BH; pvalueCutoff=0.05; qvalueCutoff=0.2; 不做 simplify",
              sprintf("每个比较组按 p.adjust 升序取前 %d 条；各组合并去重后每图 ≤30 条进入绘图", TOP_N),
              "每数据集一个 sheet；Direction=Up/Down 两个区块，区块内按比较组列出进入绘图的通路"),
    stringsAsFactors = FALSE)

  for (ds in RUN_DS) {
    cache <- readRDS(cache_path(ds$id))
    cat("\n", strrep("=", 72), "\n## GO Top", TOP_N, ":", ds$label, " [", ds$id, "]\n", strrep("=", 72), "\n", sep = "")

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

    p_up <- make_bubble(df_up, RED_COLORS,  paste0("GO BP Enrichment (Up) — ", cache$label, " (Top ", TOP_N, " per Group)"))
    p_dn <- make_bubble(df_dn, BLUE_COLORS, paste0("GO BP Enrichment (Down) — ", cache$label, " (Top ", TOP_N, " per Group)"))

    d5 <- file.path(ds_res_dir(ds), "05-GO富集"); dir.create(d5, showWarnings = FALSE, recursive = TRUE)
    if (!is.null(p_up)) {
      n <- length(levels(p_up$data$Description))
      h <- max(9, n * 0.26 + 3)
      save_fig(p_up, file.path(d5, paste0("GO_Bubble_Up_", ds$id)), 11.5, h)
      cat("上调气泡图完成 (", n, " terms, 高", round(h, 1), "in )\n")
    } else cat("上调：无显著通路，跳过\n")
    if (!is.null(p_dn)) {
      n <- length(levels(p_dn$data$Description))
      h <- max(9, n * 0.26 + 3)
      save_fig(p_dn, file.path(d5, paste0("GO_Bubble_Down_", ds$id)), 11.5, h)
      cat("下调气泡图完成 (", n, " terms, 高", round(h, 1), "in )\n")
    } else cat("下调：无显著通路，跳过\n")

    # ---- 汇总表：每数据集一个 sheet，Up/Down 两区块 ----
    parts <- list()
    for (dir_tag in c("Up", "Down")) {
      d <- if (dir_tag == "Up") p_up else p_dn
      if (is.null(d)) {
        parts[[dir_tag]] <- data.frame(Direction = dir_tag, Cluster = "无显著通路", stringsAsFactors = FALSE)
        next
      }
      pd <- d$data %>%
        mutate(Cluster = as.character(Cluster)) %>%
        dplyr::select(Cluster, Description, ID, GeneRatio, BgRatio,
                      RichFactor, FoldEnrichment, p.adjust, LogP, Count) %>%
        mutate(Direction = dir_tag,
               Description = as.character(Description),
               across(where(is.numeric), ~ round(.x, 6))) %>%
        arrange(Cluster, p.adjust) %>%
        relocate(Direction, .before = 1)
      parts[[dir_tag]] <- pd
    }
    all_cols <- union(colnames(parts[[1]]), colnames(parts[[2]]))
    pad_df <- function(x) {
      miss <- setdiff(all_cols, colnames(x))
      if (length(miss) > 0) for (m in miss) x[[m]] <- NA
      x[, all_cols, drop = FALSE]
    }
    pad_row <- as.data.frame(setNames(replicate(length(all_cols), NA, simplify = FALSE), all_cols))
    sheets[[ds$id]] <- rbind(pad_df(parts[["Up"]]), pad_row, pad_df(parts[["Down"]]))
    cat("汇总 sheet [", ds$id, "]:", nrow(sheets[[ds$id]]), "行 (含分隔行)\n", sep = "")
  }

  f <- file.path(SUB_DIR, "GO富集_Top30_汇总.xlsx")
  write_xlsx(sheets, f)
  cat("\n汇总 Excel:", f, "\n")
  cat("\nSTAGE3B OK\n")

}, error = function(e) {
  cat("\n!! ERROR:", conditionMessage(e), "\n")
})

close_log()
