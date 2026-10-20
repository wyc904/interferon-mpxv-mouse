suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(grid)
})

# ==============================================================================
# 0. 前置检查
# ==============================================================================
if (!exists("count_matrix")) stop("错误：请先加载 count_matrix！（行=gene，列=sample）")

# 输出目录
out_dir <- "D:/Desktop/1/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 如果 count_matrix 不是整数，DESeq2 会报错，这里强制转整数（更稳）
count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"

cat(">>> count_matrix 维度：", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# ==============================================================================
# 1. 构建样本信息表：Tecovirimat / DMSO / mock（严格按你给的样本名）
# ==============================================================================
tec_samples <- c(
  "JZ25055835.1.1.P.1.1.P_combined",
  "JZ25055836.1.2.P.1.2.P_combined",
  "JZ25055837.1.3.P.1.3.P_combined"
)

dmso_samples <- c(
  "JZ25055838.2.1.P.2.1.P_combined",
  "JZ25055839.2.2.P.2.2.P_combined",
  "JZ25055840.2.3.P.2.3.P_combined"
)

mock_samples <- c(
  "JZ25055847.5.1.P.5.1.P_combined",
  "JZ25055848.5.2.P.5.2.P_combined",
  "JZ25055849.5.3.P.5.3.P_combined"
)

sample_info_raw <- data.frame(
  sample = c(tec_samples, dmso_samples, mock_samples),
  group_clear = c(
    rep("Tecovirimat", length(tec_samples)),
    rep("DMSO",        length(dmso_samples)),
    rep("mock",        length(mock_samples))
  ),
  replicate = c(
    1:length(tec_samples),
    1:length(dmso_samples),
    1:length(mock_samples)
  ),
  stringsAsFactors = FALSE
)

# 检查样本是否都在 count_matrix 列中
missing <- setdiff(sample_info_raw$sample, colnames(count_matrix))
if (length(missing) > 0) {
  stop(paste0("错误：以下样本在 count_matrix 中未找到：\n", paste(missing, collapse = "\n")))
}

# 构建 coldata（DESeq2 用）
coldata <- sample_info_raw %>%
  mutate(
    # 设定因子顺序：DMSO 为参考组（你需要谁当参考就把它放第一个）
    condition = factor(group_clear, levels = c("DMSO", "Tecovirimat", "mock"))
  ) %>%
  column_to_rownames("sample")

# 取子矩阵并保证顺序一致
counts_subset <- count_matrix[, rownames(coldata), drop = FALSE]

cat(">>> 数据准备完成。组别：", paste(levels(coldata$condition), collapse = ", "), "\n")
cat(">>> 样本数：", ncol(counts_subset), "\n")

# ==============================================================================
# 2. DESeq2 主流程（一次建模，多次 contrast）
# ==============================================================================
dds <- DESeqDataSetFromMatrix(
  countData = counts_subset,
  colData   = coldata,
  design    = ~ condition
)

# 基因过滤：至少 2 个样本 counts >= 10（你原逻辑保留）
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat(">>> 过滤后保留基因数:", nrow(dds), "\n")

dds <- DESeq(dds)

# VST（后面热图可能用到）
vsd <- vst(dds, blind = FALSE)

# ==============================================================================
# 3. 可视化配色
# ==============================================================================
color_up   <- "#C75D69"
color_down <- "#4682B4"
color_ns   <- "#D3D3D3"

# ==============================================================================
# 4. pheatmap 保存辅助函数（可选）
# ==============================================================================
save_pheatmap_pdf <- function(x, filename, width = 7, height = 9) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# ==============================================================================
# 5. 封装：对任意两组做差异分析 + Volcano +（可选Heatmap）+ 保存
# ==============================================================================
run_pairwise <- function(dds, vsd, groupA, groupB,
                         FC_cutoff = 1.0, P_cutoff = 0.05,
                         top_label_n = 15, heat_top_n = 30,
                         out_dir = "D:/Desktop/1/") {
  
  cat("\n============================================================\n")
  cat(">>> 组间比较：", groupA, "vs", groupB, "\n")
  
  # contrast = (condition, groupA, groupB) => log2FC >0 表示 groupA 上调
  res <- results(dds, contrast = c("condition", groupA, groupB))
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)
  
  # 标记显著性（保持你原逻辑）
  res_df$change <- "Not Significant"
  res_df$change[!is.na(res_df$padj) & res_df$padj < P_cutoff & res_df$log2FoldChange >  FC_cutoff] <- "Up-regulated"
  res_df$change[!is.na(res_df$padj) & res_df$padj < P_cutoff & res_df$log2FoldChange < -FC_cutoff] <- "Down-regulated"
  res_df$change <- factor(res_df$change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
  
  # 统计
  sig_summary <- table(res_df$change)
  cat(">>> 差异基因统计:\n")
  cat("上调基因:", ifelse("Up-regulated" %in% names(sig_summary), sig_summary[["Up-regulated"]], 0), "\n")
  cat("下调基因:", ifelse("Down-regulated" %in% names(sig_summary), sig_summary[["Down-regulated"]], 0), "\n")
  cat("不显著基因:", ifelse("Not Significant" %in% names(sig_summary), sig_summary[["Not Significant"]], 0), "\n")
  
  # ---------- 保存结果 ----------
  prefix <- paste0("DESeq2_", groupA, "_vs_", groupB)
  write.csv(res_df, file.path(out_dir, paste0(prefix, "_all.csv")), row.names = FALSE)
  
  sig_genes <- res_df %>% filter(change != "Not Significant")
  write.csv(sig_genes, file.path(out_dir, paste0(prefix, "_sig.csv")), row.names = FALSE)
  write.csv(sig_genes %>% filter(change == "Up-regulated"),
            file.path(out_dir, paste0(prefix, "_up.csv")), row.names = FALSE)
  write.csv(sig_genes %>% filter(change == "Down-regulated"),
            file.path(out_dir, paste0(prefix, "_down.csv")), row.names = FALSE)
  
  # ---------- Volcano（关键修复：padj==0 -> 极小正数，避免Inf点被丢） ----------
  # 这里不再 filter(!is.na(padj))，而是让 NA 的 padj 落在 y=0；0的padj转成有限值
  plot_df <- res_df %>%
    mutate(
      padj_plot = case_when(
        is.na(padj) ~ 1,
        padj == 0 ~ {
          min_pos <- suppressWarnings(min(padj[padj > 0], na.rm = TRUE))
          if (is.finite(min_pos)) min_pos * 0.1 else 1e-300
        },
        TRUE ~ padj
      ),
      minuslog10padj = -log10(padj_plot),
      log2FoldChange_plot = ifelse(is.na(log2FoldChange), 0, log2FoldChange)
    )
  
  # 仅标注显著且最小 padj 的 top
  top_genes <- plot_df %>%
    filter(change != "Not Significant") %>%
    arrange(padj) %>%
    head(top_label_n) %>%
    pull(gene_id)
  
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange_plot, y = minuslog10padj)) +
    geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Down-regulated"   = color_down,
      "Not Significant"  = color_ns,
      "Up-regulated"     = color_up
    )) +
    geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "grey50") +
    geom_text_repel(
      data = subset(plot_df, gene_id %in% top_genes),
      aes(label = gene_id),
      size = 3, show.legend = FALSE, max.overlaps = 50
    ) +
    labs(
      title = paste0(groupA, " vs ", groupB),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~Adjusted~P-value)
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  print(p_vol)
  
  # ✅ 正确保存：用 groupA/groupB 拼文件名，不再用 tp
  ggsave(
    filename = file.path(out_dir, paste0(prefix, "_volcano.pdf")),
    plot     = p_vol,
    width    = 4,
    height   = 4.5,
    dpi      = 800
  )
  
  # （可选）返回结果方便你后面继续做图
  return(list(res_df = res_df, plot_df = plot_df, p_vol = p_vol))
}

# ==============================================================================
# 6. 运行：你想比较哪两组，就调用一次
# ==============================================================================
# 示例1：Tecovirimat vs DMSO
r1 <- run_pairwise(dds, vsd, groupA = "Tecovirimat", groupB = "DMSO", out_dir = out_dir)

# 示例2：Tecovirimat vs mock
r2 <- run_pairwise(dds, vsd, groupA = "Tecovirimat", groupB = "mock", out_dir = out_dir)

# 示例3：mock vs DMSO（看你需不需要）
r3 <- run_pairwise(dds, vsd, groupA = "mock", groupB = "DMSO", out_dir = out_dir)
