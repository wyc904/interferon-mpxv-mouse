############################################################
# 高水平 Venn 图（Cell 风格）三组比较分析代码
# - 读入计数
# - DESeq2 差异分析
# - 提取三组 DEGs
# - 画所有 DEGs / 上调 / 下调 的 Venn 图
# - 导出高质量 PDF & 600 dpi TIFF
############################################################

## 0. 环境准备 -----------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

# 加载必要 R 包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(VennDiagram)
library(ggVennDiagram)

## 0.1 定义统一的"Cell 风格"主题 ----------------------------------------
theme_cell <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border       = element_blank(),
      panel.grid         = element_blank(),
      axis.line          = element_line(linewidth = 0.5),
      axis.ticks         = element_line(linewidth = 0.4),
      plot.title         = element_text(
        hjust = 0.5, face = "bold", size = base_size + 2
      ),
      plot.subtitle      = element_text(
        hjust = 0.5, size = base_size
      ),
      legend.title       = element_blank(),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      legend.box.background = element_blank()
    )
}
theme_set(theme_cell())

## 0.2 输出目录 -----------------------------------------------------------
out_dir  <- "Venn_Figures_Three_Groups"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

############################################################
# 1. 读取并整理计数数据
############################################################
counts <- read.delim("D:/Desktop/Project/10-干扰素小鼠/1-源/STA1_counts_gene.txt", 
                     header = TRUE, sep = "\t", 
                     comment.char = "#", row.names = 1)

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:20]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix))

# 筛选出需要的9个样本（3组×3重复）
target_samples <- c(
  # mock 组 (Group 1)
  "JZ25055847.5.1.P.5.1.P_combined.Aligned.sortedByCoord.out.bam",
  "JZ25055848.5.2.P.5.2.P_combined.Aligned.sortedByCoord.out.bam", 
  "JZ25055849.5.3.P.5.3.P_combined.Aligned.sortedByCoord.out.bam",
  # MPXV-Ilb 组 (Group 3)
  "JZ25055841.3.1.P.3.1.P_combined.Aligned.sortedByCoord.out.bam"
,"JZ25055842.3.2.P.3.2.P_combined.Aligned.sortedByCoord.out.bam"
,"JZ25055843.3.3.P.3.3.P_combined.Aligned.sortedByCoord.out.bam",
  # MPXV-lb 组 (Group 4)
"JZ25055844.4.1.P.4.1.P_combined.Aligned.sortedByCoord.out.bam"
,"JZ25055845.4.2.P.4.2.P_combined.Aligned.sortedByCoord.out.bam"
,"JZ25055846.4.3.P.4.3.P_combined.Aligned.sortedByCoord.out.bam"
)

# 检查样本是否存在
missing_samples <- target_samples[!target_samples %in% colnames(count_matrix)]
if (length(missing_samples) > 0) {
  cat("警告：以下样本在计数矩阵中未找到:\n")
  print(missing_samples)
  # 使用实际存在的样本
  target_samples <- target_samples[target_samples %in% colnames(count_matrix)]
}

count_matrix_subset <- count_matrix[, target_samples]





# 创建样本信息表
sample_info <- data.frame(
  sample = target_samples,
  condition = c(
    rep("mock", 3),
    rep("MPXV_Ilb", 3),
    rep("MPXV_lb", 3)
  ),
  row.names = target_samples
)

# 确保condition为因子，设置mock为参考水平
sample_info$condition <- factor(sample_info$condition, 
                                levels = c("mock", "MPXV_Ilb", "MPXV_lb"))

cat("\n=== 样本信息 ===\n")
print(sample_info)

# 检查count数据的统计信息
cat("\n=== Count数据统计 ===\n")
cat("总基因数:", nrow(count_matrix_subset), "\n")
cat("总测序reads数:", sum(count_matrix_subset), "\n")
cat("每个样本的reads数:\n")
print(colSums(count_matrix_subset))

############################################################
# 2. 创建 DESeq2 对象 & 过滤低表达基因
############################################################
count_matrix_numeric <- as.matrix(count_matrix_subset)
mode(count_matrix_numeric) <- "numeric"

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData   = sample_info,
  design    = ~ condition
)

cat("\nDESeq2 对象创建成功：\n")
cat("  基因数：", nrow(dds), "\n")
cat("  样本数：", ncol(dds), "\n")

# 过滤低表达基因：至少 3 个样本中 read ≥ 10
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("过滤后保留基因数：", nrow(dds), "\n")

# 运行DESeq2
dds <- DESeq(dds)

############################################################
# 3. 执行两两比较，提取差异表达基因
############################################################
FC_cutoff <- 1.0   # |log2FC| > 1
P_cutoff  <- 0.05  # padj < 0.05

# 3.1 MPXV_Ilb vs mock
cat("\n===== 比较1: MPXV_Ilb vs mock =====\n")
res_Ilb_vs_mock <- results(dds, 
                           contrast = c("condition", "MPXV_Ilb", "mock"),
                           alpha = 0.05)

# 整理结果
res_Ilb_df <- as.data.frame(res_Ilb_vs_mock) %>%
  rownames_to_column("gene") %>%
  mutate(
    change = case_when(
      padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up",
      padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down",
      TRUE                                        ~ "NotSig"
    ),
    change = factor(change, levels = c("Down", "NotSig", "Up"))
  ) %>%
  arrange(padj)

cat("MPXV_Ilb vs mock：显著 DEGs 数 = ", 
    sum(res_Ilb_df$change != "NotSig", na.rm = TRUE), "\n")
cat("上调基因数：", sum(res_Ilb_df$change == "Up", na.rm = TRUE), "\n")
cat("下调基因数：", sum(res_Ilb_df$change == "Down", na.rm = TRUE), "\n")

# 3.2 MPXV_lb vs mock
cat("\n===== 比较2: MPXV_lb vs mock =====\n")
res_lb_vs_mock <- results(dds, 
                          contrast = c("condition", "MPXV_lb", "mock"),
                          alpha = 0.05)

# 整理结果
res_lb_df <- as.data.frame(res_lb_vs_mock) %>%
  rownames_to_column("gene") %>%
  mutate(
    change = case_when(
      padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up",
      padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down",
      TRUE                                        ~ "NotSig"
    ),
    change = factor(change, levels = c("Down", "NotSig", "Up"))
  ) %>%
  arrange(padj)

cat("MPXV_lb vs mock：显著 DEGs 数 = ", 
    sum(res_lb_df$change != "NotSig", na.rm = TRUE), "\n")
cat("上调基因数：", sum(res_lb_df$change == "Up", na.rm = TRUE), "\n")
cat("下调基因数：", sum(res_lb_df$change == "Down", na.rm = TRUE), "\n")

# 3.3 MPXV_lb vs MPXV_Ilb
cat("\n===== 比较3: MPXV_lb vs MPXV_Ilb =====\n")
res_lb_vs_Ilb <- results(dds, 
                         contrast = c("condition", "MPXV_lb", "MPXV_Ilb"),
                         alpha = 0.05)

# 整理结果
res_lb_vs_Ilb_df <- as.data.frame(res_lb_vs_Ilb) %>%
  rownames_to_column("gene") %>%
  mutate(
    change = case_when(
      padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up",
      padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down",
      TRUE                                        ~ "NotSig"
    ),
    change = factor(change, levels = c("Down", "NotSig", "Up"))
  ) %>%
  arrange(padj)

cat("MPXV_lb vs MPXV_Ilb：显著 DEGs 数 = ", 
    sum(res_lb_vs_Ilb_df$change != "NotSig", na.rm = TRUE), "\n")
cat("上调基因数：", sum(res_lb_vs_Ilb_df$change == "Up", na.rm = TRUE), "\n")
cat("下调基因数：", sum(res_lb_vs_Ilb_df$change == "Down", na.rm = TRUE), "\n")

############################################################
# 4. 构建 DEGs 列表用于韦恩图
############################################################

# 4.1 提取显著DEGs（全部）
Ilb_vs_mock_genes <- res_Ilb_df %>% 
  filter(change != "NotSig") %>% 
  pull(gene)

lb_vs_mock_genes <- res_lb_df %>% 
  filter(change != "NotSig") %>% 
  pull(gene)

lb_vs_Ilb_genes <- res_lb_vs_Ilb_df %>% 
  filter(change != "NotSig") %>% 
  pull(gene)

# 4.2 提取上调基因
Ilb_vs_mock_up <- res_Ilb_df %>% 
  filter(change == "Up") %>% 
  pull(gene)

lb_vs_mock_up <- res_lb_df %>% 
  filter(change == "Up") %>% 
  pull(gene)

lb_vs_Ilb_up <- res_lb_vs_Ilb_df %>% 
  filter(change == "Up") %>% 
  pull(gene)

# 4.3 提取下调基因
Ilb_vs_mock_down <- res_Ilb_df %>% 
  filter(change == "Down") %>% 
  pull(gene)

lb_vs_mock_down <- res_lb_df %>% 
  filter(change == "Down") %>% 
  pull(gene)

lb_vs_Ilb_down <- res_lb_vs_Ilb_df %>% 
  filter(change == "Down") %>% 
  pull(gene)

# 4.4 构建列表
venn_gene_lists <- list(
  "MPXV_Ilb_vs_mock" = Ilb_vs_mock_genes,
  "MPXV_lb_vs_mock" = lb_vs_mock_genes,
  "MPXV_lb_vs_MPXV_Ilb" = lb_vs_Ilb_genes
)

up_gene_lists <- list(
  "MPXV_Ilb_vs_mock" = Ilb_vs_mock_up,
  "MPXV_lb_vs_mock" = lb_vs_mock_up,
  "MPXV_lb_vs_MPXV_Ilb" = lb_vs_Ilb_up
)

down_gene_lists <- list(
  "MPXV_Ilb_vs_mock" = Ilb_vs_mock_down,
  "MPXV_lb_vs_mock" = lb_vs_mock_down,
  "MPXV_lb_vs_MPXV_Ilb" = lb_vs_Ilb_down
)

############################################################
# 5. 重叠基因模式统计
############################################################

extract_overlap_genes <- function(gene_lists, output_dir, prefix = "all_DEGs") {
  all_genes <- unique(unlist(gene_lists))
  if (length(all_genes) == 0) return(NULL)
  
  # 构建 presence 矩阵（基因 × 组）
  presence_df <- sapply(gene_lists, function(g) all_genes %in% g) %>%
    as.data.frame() %>%
    mutate(gene = all_genes, .before = 1)
  
  # 统计每种 TRUE/FALSE 组合对应的基因数
  overlap_counts <- presence_df %>%
    group_by(across(-gene)) %>%
    summarise(n = n(), .groups = "drop")
  
  # 保存重叠模式
  write.csv(
    overlap_counts,
    file.path(output_dir,
              paste0("overlap_pattern_counts_", prefix, ".csv")),
    row.names = FALSE
  )
  
  # 核心重叠基因：所有组都 TRUE
  core_genes <- presence_df %>%
    filter(if_all(-gene, ~ .x == TRUE)) %>%
    pull(gene)
  
  if (length(core_genes) > 0) {
    writeLines(
      core_genes,
      file.path(output_dir,
                paste0("core_DEGs_all_groups_", prefix, ".txt"))
    )
    cat("  发现", length(core_genes),
        "个在所有比较中均差异表达的核心基因（", prefix, "）\n", sep = "")
  } else {
    cat("  没有在所有比较中都差异表达的核心基因（", prefix, "）\n", sep = "")
  }
  
  list(
    overlap_counts = overlap_counts,
    core_genes     = core_genes
  )
}

# 对"所有显著 DEGs"做重叠分析
overlap_all <- extract_overlap_genes(venn_gene_lists, out_dir, "all")

############################################################
# 6. 绘制三组 Venn 图
############################################################

## 6.1 所有显著 DEGs 的 Venn 图 -----------------------------------------
cat("\n绘制三组比较的Venn图...\n")

# 设置颜色
group_colors <- c("#4E79A7", "#F28E2B", "#E15759")  # 蓝、橙、红

p_venn_all <- ggVennDiagram(
  venn_gene_lists,
  label       = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#F5F5F5", high = "#4C72B0") +
  scale_color_manual(values = group_colors) +
  labs(
    title    = "Venn Diagram of Differentially Expressed Genes",
    subtitle = "Three-way comparison between conditions"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    plot.margin   = margin(5, 5, 5, 5, "mm")
  )

# PDF（矢量） & 600 dpi TIFF
pdf(file.path(out_dir, "Venn_All_DEGs_Three_Groups.pdf"), width = 7, height = 6)
print(p_venn_all)
dev.off()

tiff(file.path(out_dir, "Venn_All_DEGs_Three_Groups_600dpi.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
print(p_venn_all)
dev.off()

cat("已输出：Venn_All_DEGs_Three_Groups.pdf / Venn_All_DEGs_Three_Groups_600dpi.tiff\n")

## 6.2 上调基因 Venn 图 --------------------------------------------------
if (length(unlist(up_gene_lists)) > 0) {
  overlap_up <- extract_overlap_genes(up_gene_lists, out_dir, "Up")
  
  p_venn_up <- ggVennDiagram(
    up_gene_lists,
    label       = "count",
    label_alpha = 0
  ) +
    scale_fill_gradient(low = "#FDEDEC", high = "#C0392B") +  # 红色系
    scale_color_manual(values = group_colors) +
    labs(title = "Venn Diagram of Up-regulated Genes",
         subtitle = "Three-way comparison between conditions") +
    theme_void() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  pdf(file.path(out_dir, "Venn_Up_DEGs_Three_Groups.pdf"), width = 7, height = 6)
  print(p_venn_up)
  dev.off()
  
  tiff(file.path(out_dir, "Venn_Up_DEGs_Three_Groups_600dpi.tiff"),
       width = 7, height = 6, units = "in", res = 600, compression = "lzw")
  print(p_venn_up)
  dev.off()
  
  cat("已输出：Venn_Up_DEGs_Three_Groups.pdf / Venn_Up_DEGs_Three_Groups_600dpi.tiff\n")
} else {
  cat("没有上调基因，跳过上调基因Venn图绘制\n")
}

## 6.3 下调基因 Venn 图 ---------------------------------------------------
if (length(unlist(down_gene_lists)) > 0) {
  overlap_down <- extract_overlap_genes(down_gene_lists, out_dir, "Down")
  
  p_venn_down <- ggVennDiagram(
    down_gene_lists,
    label       = "count",
    label_alpha = 0
  ) +
    scale_fill_gradient(low = "#EBF5FB", high = "#21618C") +  # 蓝色系
    scale_color_manual(values = group_colors) +
    labs(title = "Venn Diagram of Down-regulated Genes",
         subtitle = "Three-way comparison between conditions") +
    theme_void() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  pdf(file.path(out_dir, "Venn_Down_DEGs_Three_Groups.pdf"), width = 7, height = 6)
  print(p_venn_down)
  dev.off()
  
  tiff(file.path(out_dir, "Venn_Down_DEGs_Three_Groups_600dpi.tiff"),
       width = 7, height = 6, units = "in", res = 600, compression = "lzw")
  print(p_venn_down)
  dev.off()
  
  cat("已输出：Venn_Down_DEGs_Three_Groups.pdf / Venn_Down_DEGs_Three_Groups_600dpi.tiff\n")
} else {
  cat("没有下调基因，跳过下调基因Venn图绘制\n")
}










# 7. 重叠模式柱状图
############################################################

if (!is.null(overlap_all$overlap_counts)) {
  # 使用基础R方法创建pattern列
  overlap_summary <- overlap_all$overlap_counts
  
  # 提取除了n列以外的列（这些是逻辑值列）
  logic_cols <- overlap_summary[, !names(overlap_summary) %in% "n", drop = FALSE]
  
  # 为每一行创建模式标签
  pattern_labels <- apply(logic_cols, 1, function(x) {
    labs <- names(venn_gene_lists)[x]
    if (length(labs) == 0) "None" else paste(labs, collapse = " & ")
  })
  
  overlap_summary$pattern <- pattern_labels
  
  # 按n降序排列
  overlap_summary <- overlap_summary[order(overlap_summary$n, decreasing = TRUE), ]
  
  p_overlap_bar <- ggplot(overlap_summary,
                          aes(x = reorder(pattern, n), y = n)) +
    geom_col(width = 0.6, fill = "#2E86AB", alpha = 0.8) +
    geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(
      title = "Number of DEGs in Each Overlap Pattern",
      x     = "Overlap pattern",
      y     = "Gene count"
    ) +
    theme_bw() +
    theme(
      panel.border       = element_blank(),
      panel.grid         = element_blank(),
      axis.line          = element_line(linewidth = 0.5),
      axis.ticks         = element_line(linewidth = 0.4),
      plot.title         = element_text(
        hjust = 0.5, face = "bold", size = 16
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y = element_text(size = 10)
    )
  
  pdf(file.path(out_dir, "Overlap_Pattern_Barplot_Three_Groups.pdf"),
      width = 8, height = 5)
  print(p_overlap_bar)
  dev.off()
  
  tiff(file.path(out_dir, "Overlap_Pattern_Barplot_Three_Groups_600dpi.tiff"),
       width = 8, height = 5, units = "in", res = 600, compression = "lzw")
  print(p_overlap_bar)
  dev.off()
  
  cat("已输出：Overlap_Pattern_Barplot_Three_Groups.pdf / Overlap_Pattern_Barplot_Three_Groups_600dpi.tiff\n")
}



























############################################################
# 8. 保存详细结果
############################################################

# 保存每个比较的完整结果
write.csv(res_Ilb_df, file.path(out_dir, "MPXV_Ilb_vs_mock_DEGs.csv"), row.names = FALSE)
write.csv(res_lb_df, file.path(out_dir, "MPXV_lb_vs_mock_DEGs.csv"), row.names = FALSE)
write.csv(res_lb_vs_Ilb_df, file.path(out_dir, "MPXV_lb_vs_MPXV_Ilb_DEGs.csv"), row.names = FALSE)

# 保存DEGs列表
for (group in names(venn_gene_lists)) {
  genes <- venn_gene_lists[[group]]
  if (length(genes) > 0) {
    writeLines(genes, file.path(out_dir, paste0("DEGs_", group, ".txt")))
  }
}

# 保存session信息
sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n" + strrep("=", 60) + "\n")
cat("分析完成！\n")
cat("输出目录：", normalizePath(out_dir), "\n")
cat("\n分析汇总：\n")
cat("1. MPXV_Ilb vs mock: ", length(Ilb_vs_mock_genes), " 个 DEGs\n")
cat("2. MPXV_lb vs mock: ", length(lb_vs_mock_genes), " 个 DEGs\n")
cat("3. MPXV_lb vs MPXV_Ilb: ", length(lb_vs_Ilb_genes), " 个 DEGs\n")
cat(strrep("=", 60), "\n")