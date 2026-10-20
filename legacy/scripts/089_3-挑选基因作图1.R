############################################################
# MPXV_lb 时间序列差异表达分析 & 免疫反应-恢复过程标志物热图（完整代码）
############################################################

# 0. 加载必要的包 ---------------------------------------------------------
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(grid)

# 1. 读取并整理计数数据 ---------------------------------------------------
counts <- read.delim(
  "D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt",
  header = TRUE, sep = "\t", comment.char = "#", row.names = 1
)

# 提取 count 矩阵（从第 6 列开始是真正的 count 数据）
count_matrix <- counts[, 6:20]

# 重命名列名为简洁的样本名（共 15 个样本）
colnames(count_matrix) <- c(
  "C57_10d_1_P", "C57_10d_2_P", "C57_10d_3_P",
  "C57_3d_1_P",  "C57_3d_2_P",  "C57_3d_3_P",
  "C57_5d_1_P",  "C57_5d_2_P",  "C57_5d_3_P",
  "C57_7d_1_P",  "C57_7d_2_P",  "C57_7d_3_P",
  "c57_mock_1_P","c57_mock_2_P","c57_mock_3_P"
)

cat("重命名后的列名:\n")
print(colnames(count_matrix))

# 2. 创建样本信息数据框 ---------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  timepoint = factor(
    c(rep("10d", 3), rep("3d", 3), rep("5d", 3), 
      rep("7d", 3),  rep("mock", 3)),
    levels = c("mock", "3d", "5d", "7d", "10d")  # 指定时间顺序
  ),
  condition = factor(
    c(rep("MPXV_lb", 12), rep("mock", 3)),
    levels = c("mock", "MPXV_lb")                # mock 为对照
  ),
  group_label = c(rep("MPXV_lb", 12), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

cat("样本信息:\n")
print(sample_info)

# 3. 创建 DESeq2 对象 -----------------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 创建 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData   = sample_info,
  design    = ~ condition
)

cat("DESeq2 对象创建成功！\n")
cat("基因数:", nrow(dds), "\n")
cat("样本数:", ncol(dds), "\n")

# 4. 过滤低表达基因 -------------------------------------------------------
# 至少在 3 个样本中表达量 ≥ 10
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("过滤后保留", nrow(dds), "个基因\n")

# 5. 运行 DESeq2 分析 -----------------------------------------------------
dds <- DESeq(dds)

# 6. 变换：VST（用于热图等可视化）---------------------------------------
vsd <- vst(dds, blind = FALSE)

# 7. 按时间点进行差异分析 ------------------------------------------------
timepoints   <- c("3d", "5d", "7d", "10d")
results_list <- list()
deg_list     <- list()

# 定义差异基因阈值
FC_cutoff <- 1.0   # |log2FoldChange| > 1 (即 2 倍差异)
P_cutoff  <- 0.05  # padj < 0.05

for (tp in timepoints) {
  cat("分析时间点:", tp, "\n")
  
  # 选择当前时间点和 mock 的样本
  tp_samples <- sample_info$sample[sample_info$timepoint %in% c(tp, "mock")]
  dds_sub    <- dds[, tp_samples]
  
  # 确保 condition 因子水平正确（mock 作为参考）
  dds_sub$condition <- relevel(dds_sub$condition, ref = "mock")
  
  # 运行 DESeq2（子集）
  dds_sub <- DESeq(dds_sub)
  res <- results(
    dds_sub,
    contrast = c("condition", "MPXV_lb", "mock"),
    alpha    = 0.05
  )
  
  # 保存完整结果
  results_list[[tp]] <- res
  
  # 提取显著差异基因并标记上下调
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      change = case_when(
        padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up-regulated",
        padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      change = factor(
        change,
        levels = c("Down-regulated", "Not Significant", "Up-regulated")
      )
    ) %>%
    arrange(padj)
  
  deg_list[[tp]] <- res_df
  
  cat(tp, "vs mock: 发现",
      sum(res_df$change != "Not Significant"),
      "个显著 DEGs\n")
}

# 8. 全局配色设置 ---------------------------------------------------------
color_mock <- "#708090"   # 灰色
color_mpxv <- "#A52A2A"   # 深棕色
color_up   <- "#C75D69"   # 上调基因颜色（浅红）
color_down <- "#4682B4"   # 下调基因颜色（蓝色）
color_ns   <- "#D3D3D3"   # 非显著（浅灰）

# 9. 合并所有时间点的热图（Top 30 DEGs，列按时间排序 + 聚类行）------------
# 提取所有时间点的显著差异基因
all_sig_genes <- unique(unlist(lapply(deg_list, function(x) {
  x %>%
    filter(change != "Not Significant") %>%
    pull(gene)
})))

if (length(all_sig_genes) > 10) {
  # 计算每个基因在不同时间点出现的最小 padj
  gene_min_padj <- sapply(all_sig_genes, function(gene) {
    min(sapply(deg_list, function(df) {
      idx <- which(df$gene == gene)
      if (length(idx) > 0) df$padj[idx] else 1
    }), na.rm = TRUE)
  })
  
  # 取前 30 个最显著的基因
  top_genes <- names(sort(gene_min_padj))[1:min(30, length(gene_min_padj))]
  
  # ☆ 按时间点 mock → 3d → 5d → 7d → 10d 排列样本列
  sample_info_ordered_all <- sample_info %>%
    arrange(timepoint)
  
  # 提取表达矩阵并标准化（按行 Z-score）
  mat <- assay(vsd)[top_genes, sample_info_ordered_all$sample]
  mat <- t(scale(t(mat)))
  
  # 创建样本注释（按时间点顺序）
  annotation_col <- data.frame(
    Timepoint = sample_info_ordered_all$timepoint,
    Condition = sample_info_ordered_all$group_label
  )
  rownames(annotation_col) <- sample_info_ordered_all$sample
  
  # 创建颜色映射
  annotation_colors <- list(
    Timepoint = c(
      "mock" = color_mock,
      "3d"   = "#1B9E77",
      "5d"   = "#D95F02",
      "7d"   = "#7570B3",
      "10d"  = "#E7298A"
    ),
    Condition = c(
      "MOCK"    = color_mock,
      "MPXV_lb" = color_mpxv
    )
  )
  
  # 绘制热图：列固定时间顺序，按基因聚类（体现梯度变化）
  p_heat <- pheatmap(
    mat,
    color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    annotation_col    = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows      = TRUE,   # 基因聚类
    cluster_cols      = FALSE,  # 列不打乱，保持 mock → 3d → 5d → 7d → 10d
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 8,
    fontsize_col      = 9,
    main              = "MPXV_lb 时间序列差异表达基因热图（Top 30，按时间排序）"
  )
  
  # 保存热图为 PDF
  pdf("Heatmap_TimeSeries_All_DEGs_Top30_TimeOrdered.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  dev.off()
}

# 10. 各时间点 DEGs 数量统计柱状图 --------------------------------------
deg_counts <- data.frame(
  timepoint = timepoints,
  up = sapply(timepoints, function(x) {
    if (x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$change == "Up-regulated")
    } else { 0 }
  }),
  down = sapply(timepoints, function(x) {
    if (x %in% names(deg_list) && nrow(deg_list[[x]]) > 0) {
      sum(deg_list[[x]]$change == "Down-regulated")
    } else { 0 }
  })
)

deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c(up, down),
               names_to = "direction", values_to = "count")

p_bar <- ggplot(deg_counts_long,
                aes(x = timepoint, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "各时间点差异基因数量统计",
    x     = "时间点",
    y     = "DEGs 数量"
  ) +
  scale_fill_manual(values = c("down" = color_down, "up" = color_up),
                    labels = c("down" = "Down-regulated", "up" = "Up-regulated")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_bar)
ggsave("DEGs_Count_By_Timepoint.pdf", p_bar, width = 8, height = 6)

# 11. 免疫反应 -> 恢复过程标志物基因热图（至少 40 个基因）----------------

# 11.1 定义“免疫反应 + 恢复/调节”大列表（约 48 个基因）
# 分 4 类：IFN/抗病毒、炎症趋化、适应性免疫、恢复/修复与调节
genes_ifn <- c(
  "Ddx58", "Ifih1", "Tlr3", "Tlr7", "Sting1",
  "Irf7", "Irf9", "Stat1",
  "Oas1a", "Oas1b", "Oas2", "Oas3",
  "Mx1", "Mx2",
  "Ifit1", "Ifit2", "Ifit3",
  "Isg15", "Rsad1", "Gbp2b"
)

genes_inflam <- c(
  "Ccl2", "Ccl3", "Ccl5", "Ccl21a",
  "Cxcl9", "Cxcl10", "Cxcl11",
  "Il6", "Il1b", "Tnf", "Ifng"
)

genes_adapt <- c(
  "Cd3e", "Cd4", "Cd8a",      # T 细胞
  "Cd19", "Cd79a",            # B 细胞
  "Nkg7", "Gzmb"              # NK / CTL 细胞毒
)

genes_recovery <- c(
  "Cd68", "Mrc1", "Havcr1", "Mmp9",  # 巨噬/修复相关
  "Il10", "Il10ra", "Il10rb",        # 抗炎/恢复
  "Tgfb1", "Socs1", "Socs3"          # 免疫负调控
)

immune_recovery_genes <- c(
  genes_ifn,
  genes_inflam,
  genes_adapt,
  genes_recovery
)

# 去除重复
immune_recovery_genes <- unique(immune_recovery_genes)

# 11.2 在表达矩阵中筛选真正存在的基因 -----------------------------------
genes_present <- immune_recovery_genes[immune_recovery_genes %in% rownames(vsd)]

cat("在表达矩阵中找到的免疫反应/恢复相关基因数量：", length(genes_present), "\n")
print(genes_present)

if (length(genes_present) > 0) {
  
  # 提取这些基因的表达量
  marker_expr <- assay(vsd)[genes_present, ]
  
  # 按时间顺序排列样本（mock → 3d → 5d → 7d → 10d）
  sample_info_ordered <- sample_info %>%
    arrange(timepoint)
  
  marker_expr <- marker_expr[, sample_info_ordered$sample]
  
  # 行标准化（每个基因自身 Z-score）
  marker_expr_scaled <- t(scale(t(marker_expr)))
  
  # 11.3 构建列注释：时间点 + 处理组 ------------------------------------
  annotation_col2 <- data.frame(
    Timepoint = sample_info_ordered$timepoint,
    Condition = sample_info_ordered$group_label
  )
  rownames(annotation_col2) <- sample_info_ordered$sample
  
  # 11.4 构建行注释：基因类别（IFN/炎症/适应性/恢复）----------------------
  gene_category <- sapply(rownames(marker_expr_scaled), function(g) {
    if (g %in% genes_ifn)        return("IFN/抗病毒")
    if (g %in% genes_inflam)     return("炎症/趋化因子")
    if (g %in% genes_adapt)      return("适应性免疫")
    if (g %in% genes_recovery)   return("恢复/调节")
    return("其他")
  })
  
  annotation_row <- data.frame(
    Category = factor(
      gene_category,
      levels = c("IFN/抗病毒", "炎症/趋化因子", "适应性免疫", "恢复/调节", "其他")
    )
  )
  rownames(annotation_row) <- rownames(marker_expr_scaled)
  
  # 11.5 颜色映射（列 + 行注释一起放在一个 list）---------------------------
  annotation_colors3 <- list(
    Timepoint = c(
      "mock" = color_mock,
      "3d"   = "#1B9E77",
      "5d"   = "#D95F02",
      "7d"   = "#7570B3",
      "10d"  = "#E7298A"
    ),
    Condition = c(
      "MOCK"    = color_mock,
      "MPXV_lb" = color_mpxv
    ),
    Category = c(
      "IFN/抗病毒"     = "#E41A1C",
      "炎症/趋化因子" = "#377EB8",
      "适应性免疫"     = "#4DAF4A",
      "恢复/调节"       = "#984EA3",
      "其他"           = "#999999"
    )
  )
  
  # 11.6 画免疫反应 -> 恢复过程标志物热图：列按时间，按基因聚类 ----------
  p_core <- pheatmap(
    marker_expr_scaled,
    color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    annotation_col    = annotation_col2,
    annotation_row    = annotation_row,
    annotation_colors = annotation_colors3,
    cluster_rows      = TRUE,   # 基因聚类
    cluster_cols      = FALSE,  # 列固定时间轴：mock → 3d → 5d → 7d → 10d
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 9,
    main              = "猴痘感染 C57 小鼠免疫反应到恢复过程标志物热图"
  )
  
  # 保存 PDF
  pdf("Heatmap_Immune_Response_Recovery_Markers_TimeSeries.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_core$gtable)
  dev.off()
  
} else {
  cat("警告：定义的免疫反应/恢复相关基因在 vsd 中都没有找到，请检查基因名与行名是否一致（大小写、物种等）。\n")
}

############################################################
# 代码结束
############################################################
