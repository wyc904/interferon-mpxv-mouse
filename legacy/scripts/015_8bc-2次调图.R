library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)

FC_cutoff <- 1
P_cutoff  <- 0.05

color_up   <- "#C75D69"
color_down <- "#4682B4"
color_ns   <- "#D3D3D3"

# 关键：保证是整数计数（DESeq2推荐）
count_matrix_numeric <- as.matrix(count_matrix)
storage.mode(count_matrix_numeric) <- "integer"

# 全部样本的 dds（只建一次）
sample_info$group_tp <- factor(as.character(sample_info$timepoint),
                               levels = c("mock","3d","5d","7d","10d"))

dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData   = sample_info,
  design    = ~ group_tp
)

# 可选：全局过滤（如果你想保留“只在某个时间点出现”的基因，建议先别全局过滤）
# keep_all <- rowSums(counts(dds_all) >= 10) >= 3
# dds_all  <- dds_all[keep_all, ]

dds_all <- DESeq(dds_all)

get_res_df <- function(tp, pairwise = TRUE) {
  if (pairwise) {
    # 子集：只取 tp + mock 样本，重新跑（更可比）
    keep_samp <- colData(dds_all)$group_tp %in% c("mock", tp)
    dds_sub <- dds_all[, keep_samp]
    dds_sub$group_tp <- droplevels(dds_sub$group_tp)
    design(dds_sub) <- ~ group_tp
    
    # 子集过滤（按该对比的样本数更合理）
    keep_gene <- rowSums(counts(dds_sub) >= 10) >= 3
    dds_sub <- dds_sub[keep_gene, ]
    
    dds_sub <- DESeq(dds_sub)
    res <- results(dds_sub, contrast = c("group_tp", tp, "mock"), alpha = P_cutoff)
  } else {
    # 全数据一起跑的结果
    res <- results(dds_all, contrast = c("group_tp", tp, "mock"), alpha = P_cutoff)
  }
  
  # 不要直接 filter 掉 NA：否则点会“消失”，让你误以为缺基因
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      padj_plot = ifelse(is.na(padj), 1, padj),
      lfc_plot  = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
      change = case_when(
        padj_plot < P_cutoff & lfc_plot >=  FC_cutoff  ~ "Up-regulated",
        padj_plot < P_cutoff & lfc_plot <= -FC_cutoff  ~ "Down-regulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # ✅检查：显著基因有没有被漏标
  sig_raw <- rownames(res)[!is.na(res$padj) &
                             res$padj < P_cutoff &
                             !is.na(res$log2FoldChange) &
                             abs(res$log2FoldChange) >= FC_cutoff]
  sig_marked <- res_df$gene[res_df$change != "Not Significant"]
  
  missing <- setdiff(sig_raw, sig_marked)
  if (length(missing) > 0) {
    message(tp, " vs mock：发现漏标显著基因数量 = ", length(missing))
    print(head(missing, 20))
  } else {
    message(tp, " vs mock：未发现漏标显著基因 ✅")
  }
  
  res_df
}

plot_volcano <- function(res_df, tp) {
  top_genes <- res_df %>%
    filter(change != "Not Significant") %>%
    arrange(padj_plot) %>%
    head(15) %>%
    pull(gene)
  
  ggplot(res_df, aes(x = lfc_plot, y = -log10(padj_plot))) +
    geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down-regulated" = color_down,
                                  "Not Significant" = color_ns,
                                  "Up-regulated" = color_up)) +
    geom_vline(xintercept = c(-FC_cutoff, FC_cutoff),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(P_cutoff),
               linetype = "dashed", color = "grey50") +
    geom_text_repel(data = subset(res_df, gene %in% top_genes),
                    aes(label = gene, color = change),
                    size = 2.5, show.legend = FALSE, max.overlaps = 15) +
    labs(title = paste("Volcano:", tp, "vs Mock"),
         x = expression(Log[2]~Fold~Change),
         y = expression(-Log[10]~Adjusted~P-value)) +
    theme_bw() +
    theme(legend.position = "bottom")
}

timepoints <- c("3d", "5d", "7d", "10d")

# 选择模式：TRUE=每个时间点子集重跑（推荐用来对齐你另一张图）
pairwise_mode <- TRUE

deg_list <- list()
for (tp in timepoints) {
  res_df <- get_res_df(tp, pairwise = pairwise_mode)
  deg_list[[tp]] <- res_df
  
  p <- plot_volcano(res_df, tp)
  print(p)
  
  # 保存：PDF 是矢量图，最高清；PNG 可用 ragg 更清晰
  ggsave(file.path(output_dir, paste0(tp, "_vs_Mock.pdf")),
         p, width = 4, height = 4.5, device = cairo_pdf)
}
