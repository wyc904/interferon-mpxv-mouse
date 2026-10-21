# ============================================================================
# stage5_composite.R — 每套数据集的合成图（模板式 a–f：PCA/DEG数量/韦恩×2/火山网格/GO×2）
# 输入: _cache/plots_<id>.rds + _cache/plots_go_all.rds
# 输出: 4-结果/<dataset>/07-合成图/Composite_<id>.pdf|png
# ============================================================================
source("scripts/1009_pipeline/utils/common.R", encoding = "UTF-8")
suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(grid)
  library(ggVennDiagram)
})

open_log("stage5_composite")

go_h <- function(p) {
  if (is.null(p)) return(0)
  n <- length(levels(p$data$Description))
  if (is.na(n) || n == 0) n <- length(unique(p$data$Description))
  max(8, n * 0.26 + 2.5)
}

tryCatch({

  go_all <- readRDS(file.path(CACHE, "plots_go_all.rds"))

  for (ds in DATASETS) {
    cat("\n## 合成图:", ds$label, " [", ds$id, "]\n", sep = "")
    pl <- readRDS(file.path(CACHE, paste0("plots_", ds$id, ".rds")))
    go <- go_all[[ds$id]]

    # ---- 第一行: a PCA | b DEG 数量 | c 韦恩G3 | d 韦恩G4 ----
    cache <- readRDS(cache_path(ds$id))
    build_venn_small <- function(trt) {
      s <- list(G3 = c("G2vsG1", "G3vsG1", "G3vsG2"),
                G4 = c("G2vsG1", "G4vsG1", "G4vsG2"))[[trt]]
      labs <- gsub("vs", " vs ", s)
      lst <- setNames(lapply(s, function(cid) {
        r <- cache$results[[cid]]
        unique(r$Geneid[r$change != "NS"])
      }), labs)
      ggVennDiagram(lst, label = "count", label_alpha = 0, edge_size = 0.8,
                    set_size = 2.5, label_size = 3) +
        scale_fill_gradient(low = "#F7F7F7", high = "#4C72B0") +
        scale_color_manual(values = c("#A52A2A", "#1B9E77", "#7570B3")) +
        labs(title = paste0("DEG overlap (All) — ", cache$label),
             subtitle = paste(labs, collapse = "  |  ")) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
              plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
              legend.position = "none",
              aspect.ratio = 1,
              plot.margin = margin(4, 4, 4, 4))
    }
    r1 <- list(pl$pca, pl$bar)
    if (!is.null(pl$venn_g3)) r1 <- c(r1, list(build_venn_small("G3")))
    if (!is.null(pl$venn_g4)) r1 <- c(r1, list(build_venn_small("G4")))
    row1 <- plot_grid(plotlist = r1, nrow = 1, labels = letters[seq_along(r1)], label_size = 17,
                      rel_widths = c(1.0, 1.1, 0.95, 0.95)[seq_along(r1)])

    # ---- 第二行: e 火山网格（整体一个字母）----
    e_letter <- letters[length(r1) + 1]
    row2 <- ggdraw(patchworkGrob(pl$volc_grid)) +
      draw_label(e_letter, x = 0.008, y = 0.992, hjust = 0, vjust = 1, size = 17, fontface = "bold")

    rows <- list(row1, row2)
    relh <- c(6.0, 10.0)
    h_total <- 1.0 + 6.0 + 10.0

    # ---- 第三行: f GO 上调/下调 ----
    go_pair <- NULL
    if (!is.null(go$up) && !is.null(go$dn)) go_pair <- (go$up | go$dn)
    else if (!is.null(go$up)) go_pair <- go$up
    else if (!is.null(go$dn)) go_pair <- go$dn
    if (!is.null(go_pair)) {
      f_letter <- letters[length(r1) + 2]
      rg <- ggdraw(patchworkGrob(go_pair)) +
        draw_label(f_letter, x = 0.008, y = 0.992, hjust = 0, vjust = 1, size = 17, fontface = "bold")
      rows <- c(rows, list(rg))
      h_go <- max(go_h(go$up), go_h(go$dn))
      relh <- c(relh, h_go)
      h_total <- h_total + h_go
    }

    title <- ggdraw() + draw_label(
      paste0("MPXV infection — ", ds$label, "    (G1 uninfected | G2 infected | G3 Tecovirimat | G4 DMSO vehicle)"),
      fontface = "bold", size = 16)

    final <- plot_grid(title, plot_grid(plotlist = rows, ncol = 1, rel_heights = relh),
                       ncol = 1, rel_heights = c(0.9, h_total))

    d7 <- file.path(ds_res_dir(ds), "07-合成图"); dir.create(d7, showWarnings = FALSE, recursive = TRUE)
    base <- file.path(d7, paste0("Composite_", ds$id))
    W <- 17
    H <- 1.2 + h_total
    ggsave(paste0(base, ".pdf"), final, width = W, height = H, device = "pdf", useDingbats = FALSE)
    ggsave(paste0(base, ".png"), final, width = W, height = H, dpi = 150, device = ragg::agg_png, bg = "white")
    cat(sprintf("合成图完成: %.1f x %.1f in\n", W, H))
  }

  cat("\nSTAGE5 OK\n")

}, error = function(e) {
  cat("\n!! ERROR:", conditionMessage(e), "\n")
})

close_log()
