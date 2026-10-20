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
if(!exists("count_matrix")) stop("错误：请先加载 count_matrix！（行=gene，列=sample）")

# 如果 count_matrix 不是整数，DESeq2 会报错，这里强制转整数（更稳）
count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"

cat(">>> count_matrix 维度：", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# ==============================================================================
# 1. 构建样本信息表：Tecovirimat / DMSO / mock
#    （严格按你给的样本名）
# ==============================================================================
colnames(count_matrix) <- trimws(colnames(count_matrix))

tec_samples <- c(
  "JZ25070172.AGB6.Tecovirimat.1.P.AGB6.Tecovirimat.1.P_combined",
  "JZ25070173.AGB6.Tecovirimat.2.P.AGB6.Tecovirimat.2.P_combined",
  "JZ25070174.AGB6.Tecovirimat.3.P.AGB6.Tecovirimat.3.P_combined"
)

dmso_samples <- c(
  "JZ25070169.AGB6.MPXV_lb.1.P.AGB6.MPXV_lb.1.P_combined",
  "JZ25070170.AGB6.MPXV_lb.2.P.AGB6.MPXV_lb.2.P_combined",
  "JZ25070171.AGB6.MPXV_lb.3.P.AGB6.MPXV_lb.3.P_combined"
)

mock_samples <- c(
  "JZ25070166.AGB6.MOCK.1.P.AGB6.MOCK.1.P_combined",
  "JZ25070167.AGB6.MOCK.2.P.AGB6.MOCK.2.P_combined",
  "JZ25070168.AGB6.MOCK.3.P.AGB6.MOCK.3.P_combined"
)

sample_info_raw <- data.frame(
  sample = c(tec_samples, dmso_samples, mock_samples),
  group_clear = c(
    rep("Tecovirimat", length(tec_samples)),
    rep("DMSO", length(dmso_samples)),
    rep("mock", length(mock_samples))
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
if(length(missing) > 0){
  stop(paste0("错误：以下样本在 count_matrix 中未找到：\n", paste(missing, collapse = "\n")))
}

# 构建 coldata（DESeq2 用）
coldata <- sample_info_raw %>%
  mutate(
    condition = factor(group_clear, levels = c("mock", "DMSO", "Tecovirimat")) # 设定基准/顺序
  ) %>%
  column_to_rownames("sample")

# 取子矩阵并保证顺序一致
counts_subset <- count_matrix[, rownames(coldata)]

cat(">>> 数据准备完成。组别：", paste(levels(coldata$condition), collapse = ", "), "\n")
cat(">>> 样本数：", ncol(counts_subset), "\n")

# ==============================================================================
# 2. DESeq2 主流程（一次建模，多次 contrast）
# ==============================================================================
dds <- DESeqDataSetFromMatrix(
  countData = counts_subset,
  colData = coldata,
  design = ~ condition
)

# 基因过滤：至少 2 个样本 counts >= 10
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat(">>> 过滤后保留基因数:", nrow(dds), "\n")

# 运行DESeq分析
dds <- DESeq(dds)

# ==============================================================================
# 3. 可视化配色（按三组）
# ==============================================================================
color_group <- c(
  "Tecovirimat" = "#8B4513",  # 深棕
  "DMSO" = "#708090",         # 灰
  "mock" = "#2E8B57"          # 绿
)

color_up <- "#C75D69"
color_down <- "#4682B4"
color_ns <- "#D3D3D3"

# ==============================================================================
# 4. 先做一次 PCA（全样本）
# ==============================================================================
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 按我们表里的 replicate 生成标签
pcaData$replicate <- coldata[pcaData$name, "replicate"]
pcaData$new_label <- paste0(pcaData$condition, "_", pcaData$replicate)

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = color_group) +
  geom_text_repel(aes(label = new_label), size = 3, show.legend = FALSE, nudge_y = 0.5) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    title = "PCA: Tecovirimat vs DMSO vs mock",
    color = "Condition", shape = "Condition"
  ) +
  theme_bw()

print(p_pca)
ggsave("Figure_PCA_AllGroups.pdf", p_pca, width = 6.5, height = 5.5)

# ==============================================================================
# 5. 封装：对任意两组做差异分析 + Volcano + Heatmap + 保存
# ==============================================================================
save_pheatmap_pdf <- function(x, filename, width = 7, height = 9) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

run_pairwise <- function(dds, vsd, groupA, groupB,
                         FC_cutoff = 1.0, P_cutoff = 0.05,
                         top_label_n = 15, heat_top_n = 30) {
  
  cat("\n============================================================\n")
  cat(">>> 组间比较：", groupA, "vs", groupB, "\n")
  
  # contrast = (condition, groupA, groupB) => log2FC >0 表示 groupA 上调
  res <- results(dds, contrast = c("condition", groupA, groupB), alpha = P_cutoff)
  
  res_df <- as.data.frame(res) %>%
    mutate(gene_id = rownames(.)) %>%
    arrange(padj)
  
  # 标记显著性
  res_df$change <- "Not Significant"
  res_df$change[!is.na(res_df$padj) & res_df$padj < P_cutoff & res_df$log2FoldChange > FC_cutoff] <- "Up-regulated"
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
  write.csv(res_df, paste0(prefix, "_all.csv"), row.names = FALSE)
  
  sig_genes <- res_df %>% filter(change != "Not Significant")
  if (nrow(sig_genes) > 0) {
    write.csv(sig_genes, paste0(prefix, "_sig.csv"), row.names = FALSE)
    
    write.csv(sig_genes %>% filter(change == "Up-regulated"),
              paste0(prefix, "_up.csv"), row.names = FALSE)
    write.csv(sig_genes %>% filter(change == "Down-regulated"),
              paste0(prefix, "_down.csv"), row.names = FALSE)
  } else {
    cat(">>> 没有显著差异基因\n")
  }
  
  # ---------- Volcano ----------
  plot_df <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(minuslog10padj = -log10(padj))
  
  # 仅标注显著且最小 padj 的 top
  top_genes <- plot_df %>%
    filter(change != "Not Significant") %>%
    arrange(padj) %>%
    head(top_label_n) %>%
    pull(gene_id)
  
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange, y = minuslog10padj)) +
    geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Down-regulated" = color_down,
      "Not Significant" = color_ns,
      "Up-regulated" = color_up
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
  ggsave(paste0("Figure_Volcano_", groupA, "_vs_", groupB, ".pdf"),
         p_vol, width = 6, height = 5.5, dpi = 300)
  
  # ---------- 热图 ----------
  if (nrow(sig_genes) > 0) {
    # 选择差异最显著的基因
    sig_genes_ordered <- res_df %>%
      filter(change != "Not Significant") %>%
      arrange(padj) %>%
      head(min(heat_top_n, nrow(sig_genes)))
    
    if (nrow(sig_genes_ordered) > 1) {
      # 获取vsd矩阵
      vsd_mat <- assay(vsd)
      
      # 只取两个比较组的样本
      sample_names <- rownames(coldata)[coldata$condition %in% c(groupA, groupB)]
      heatmap_mat <- vsd_mat[sig_genes_ordered$gene_id, sample_names, drop = FALSE]
      
      # 样本注释
      annotation_col <- data.frame(
        Condition = factor(coldata[sample_names, "condition"])
      )
      rownames(annotation_col) <- sample_names
      
      # 颜色映射
      ann_colors <- list(Condition = color_group[levels(annotation_col$Condition)])
      
      # 绘制热图
      pheat <- pheatmap(heatmap_mat,
                        scale = "row",
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        show_rownames = ifelse(nrow(heatmap_mat) <= 50, TRUE, FALSE),
                        show_colnames = TRUE,
                        annotation_col = annotation_col,
                        annotation_colors = ann_colors,
                        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                        fontsize_row = 8,
                        main = paste0("Heatmap: ", groupA, " vs ", groupB, "\nTop ", nrow(sig_genes_ordered), " DE genes"),
                        silent = TRUE
      )
      
      # 保存热图
      print(pheat)
      save_pheatmap_pdf(pheat, paste0("Figure_Heatmap_", groupA, "_vs_", groupB, ".pdf"))
    } else {
      cat(">>> 差异基因数量不足，跳过热图绘制\n")
    }
  } else {
    cat(">>> 没有显著差异基因，跳过热图绘制\n")
  }
  
  return(res_df)
}

# ==============================================================================
# 6. 运行三组两两比较
# ==============================================================================
cat("\n>>> 开始进行两两比较分析...\n")

comparisons <- list(
  c("Tecovirimat", "DMSO"),
  c("Tecovirimat", "mock"),
  c("DMSO", "mock")
)

# 存储所有结果
all_results <- list()

for (cmp in comparisons) {
  groupA <- cmp[1]
  groupB <- cmp[2]
  cat("\n============================================================\n")
  cat(">>> 正在分析：", groupA, "vs", groupB, "\n")
  
  result <- run_pairwise(
    dds = dds,
    vsd = vsd,
    groupA = groupA,
    groupB = groupB,
    FC_cutoff = 1.0,
    P_cutoff = 0.05,
    top_label_n = 15,
    heat_top_n = 30
  )
  
  all_results[[paste0(groupA, "_vs_", groupB)]] <- result
}

# ==============================================================================
# 7. 保存样本信息和分析总结
# ==============================================================================
write.csv(as.data.frame(coldata), "Sample_info_Tecovirimat_DMSO_mock.csv")

# 创建分析总结
summary_df <- data.frame(
  Comparison = names(all_results),
  Upregulated = sapply(all_results, function(x) sum(x$change == "Up-regulated", na.rm = TRUE)),
  Downregulated = sapply(all_results, function(x) sum(x$change == "Down-regulated", na.rm = TRUE)),
  Total_DE = sapply(all_results, function(x) sum(x$change != "Not Significant", na.rm = TRUE))
)

cat("\n============================================================\n")
cat(">>> 分析总结：\n")
print(summary_df)
write.csv(summary_df, "DESeq2_analysis_summary.csv", row.names = FALSE)

# 保存R工作空间以便后续分析
save(dds, vsd, coldata, all_results, summary_df, 
     file = "DESeq2_analysis_results.RData")

cat("\n>>> 全部分析完成！已输出以下文件：\n")
cat("1. PCA图: Figure_PCA_AllGroups.pdf\n")
cat("2. 各比较组火山图: Figure_Volcano_*.pdf\n")
cat("3. 各比较组热图: Figure_Heatmap_*.pdf\n")
cat("4. 差异基因结果: DESeq2_*_all/sig/up/down.csv\n")
cat("5. 样本信息: Sample_info_Tecovirimat_DMSO_mock.csv\n")
cat("6. 分析总结: DESeq2_analysis_summary.csv\n")
cat("7. R工作空间: DESeq2_analysis_results.RData\n")