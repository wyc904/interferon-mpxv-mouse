# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(VennDiagram)
library(ggVennDiagram)

## 0.1 定义统一的"Cell风格"主题 ----------------------------------------
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
out_dir  <- "Venn_Figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

############################################################
# 1. 读取并整理计数数据
############################################################

# 读取counts文件
counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/8bc_counts_gene", 
                     header = TRUE, sep = "\t", 
                     comment.char = "#", row.names = 1)

cat("原始数据维度:", dim(counts), "\n")
cat("列名:\n")
print(colnames(counts))

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:18]

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("8bc_10d_1_P", 
                            "8bc_3d_1_P", "8bc_3d_2_P", "8bc_3d_3_P",
                            "8bc_5d_1_P", "8bc_5d_2_P", "8bc_5d_3_P",
                            "8bc_7d_1_P", "8bc_7d_2_P", "8bc_7d_3_P",
                            "8bc_mock_1_P", "8bc_mock_2_P", "8bc_mock_3_P")

print("重命名后的列名:")
print(colnames(count_matrix))

############################################################
# 2. 样本信息（根据实际数据修正）
############################################################

sample_info <- data.frame(
  sample     = colnames(count_matrix),
  timepoint  = factor(
    c("10d",                    # 8bc_10d_1_P
      rep("3d",  3),            # 8bc_3d_1_P, 8bc_3d_2_P, 8bc_3d_3_P
      rep("5d",  3),            # 8bc_5d_1_P, 8bc_5d_2_P, 8bc_5d_3_P
      rep("7d",  3),            # 8bc_7d_1_P, 8bc_7d_2_P, 8bc_7d_3_P
      rep("mock", 3)),          # 8bc_mock_1_P, 8bc_mock_2_P, 8bc_mock_3_P
    levels = c("mock", "3d", "5d", "7d", "10d")
  ),
  condition  = factor(
    c(rep("MPXV_lb", 10), rep("mock", 3)),  # 10个感染样本 + 3个mock
    levels = c("mock", "MPXV_lb")
  ),
  group_label = c(rep("MPXV_lb", 10), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

cat("\n样本信息概览：\n")
print(sample_info)

# 验证样本数量
cat("总样本数:", nrow(sample_info), "\n")
cat("MPXV_lb样本数:", sum(sample_info$condition == "MPXV_lb"), "\n")
cat("mock样本数:", sum(sample_info$condition == "mock"), "\n")

############################################################
# 3. 创建DESeq2对象 & 过滤低表达基因
############################################################

count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData   = sample_info,
  design    = ~ condition
)

cat("\nDESeq2对象创建成功：\n")
cat("  基因数：", nrow(dds), "\n")
cat("  样本数：", ncol(dds), "\n")

# 过滤低表达基因：至少3个样本中read ≥ 10
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("过滤后保留基因数：", nrow(dds), "\n")

############################################################
# 4. 分时间点做差异分析（MPXV_lb vs mock）
############################################################

timepoints   <- c("3d", "5d", "7d", "10d")
results_list <- list()
deg_list     <- list()

FC_cutoff <- 1.0   # |log2FC| > 1
P_cutoff  <- 0.05  # padj < 0.05

for (tp in timepoints) {
  cat("\n===== 分析时间点：", tp, " =====\n")
  
  # 当前时间点 + mock样本
  tp_samples <- sample_info$sample[sample_info$timepoint %in% c(tp, "mock")]
  dds_sub    <- dds[, tp_samples]
  
  # 检查样本数量
  cat("  用于比较的样本数：", length(tp_samples), "\n")
  cat("  MPXV_lb样本：", sum(dds_sub$condition == "MPXV_lb"), "\n")
  cat("  mock样本：", sum(dds_sub$condition == "mock"), "\n")
  
  # 设置对照组
  dds_sub$condition <- relevel(dds_sub$condition, ref = "mock")
  
  # 差异分析
  dds_sub <- DESeq(dds_sub)
  res     <- results(
    dds_sub,
    contrast = c("condition", "MPXV_lb", "mock"),
    alpha    = 0.05
  )
  
  results_list[[tp]] <- res
  
  # 整理为data.frame，并标记上下调
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      change = case_when(
        padj < P_cutoff & log2FoldChange >  FC_cutoff ~ "Up",
        padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down",
        TRUE                                        ~ "NotSig"
      ),
      change = factor(change,
                      levels = c("Down", "NotSig", "Up"))
    ) %>%
    arrange(padj)
  
  deg_list[[tp]] <- res_df
  
  # 统计结果
  n_sig <- sum(res_df$change != "NotSig", na.rm = TRUE)
  n_up <- sum(res_df$change == "Up", na.rm = TRUE)
  n_down <- sum(res_df$change == "Down", na.rm = TRUE)
  
  cat(tp, "vs mock：显著DEGs数 =", n_sig, 
      "(上调 =", n_up, "下调 =", n_down, ")\n")
}

############################################################
# 5. 构建各时间点DEGs列表（全部/上调/下调）
############################################################

# 5.1 所有显著DEGs
venn_gene_lists <- list()
for (tp in timepoints) {
  sig_genes <- deg_list[[tp]] %>%
    filter(change != "NotSig") %>%
    pull(gene)
  if (length(sig_genes) > 0) {
    venn_gene_lists[[tp]] <- sig_genes
  }
}

if (length(venn_gene_lists) < 2) {
  stop("显著DEGs的时间点少于2个，无法绘制Venn图。")
}

cat("\n用于Venn图的时间点：\n")
print(names(venn_gene_lists))

# 5.2 上调 & 下调DEGs列表
up_gene_lists   <- list()
down_gene_lists <- list()

for (tp in timepoints) {
  if (!tp %in% names(deg_list)) next
  
  up_genes <- deg_list[[tp]] %>%
    filter(change == "Up") %>%
    pull(gene)
  down_genes <- deg_list[[tp]] %>%
    filter(change == "Down") %>%
    pull(gene)
  
  if (length(up_genes) > 0)
    up_gene_lists[[tp]] <- up_genes
  
  if (length(down_genes) > 0)
    down_gene_lists[[tp]] <- down_genes
}

cat("上调基因列表长度：", length(up_gene_lists), "\n")
cat("下调基因列表长度：", length(down_gene_lists), "\n")

############################################################
# 6. 重叠基因模式统计
############################################################

extract_overlap_genes <- function(gene_lists, output_dir, prefix = "all_DEGs") {
  all_genes <- unique(unlist(gene_lists))
  if (length(all_genes) == 0) return(NULL)
  
  # 构建presence矩阵（基因×组）
  presence_df <- sapply(gene_lists, function(g) all_genes %in% g) %>%
    as.data.frame() %>%
    mutate(gene = all_genes, .before = 1)
  
  # 统计每种TRUE/FALSE组合对应的基因数
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
  
  # 核心重叠基因：所有组都TRUE
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
        "个在所有组中均差异表达的核心基因（", prefix, "）\n", sep = "")
  } else {
    cat("  没有在所有组都差异表达的核心基因（", prefix, "）\n", sep = "")
  }
  
  list(
    overlap_counts = overlap_counts,
    core_genes     = core_genes
  )
}

# 对"所有显著DEGs"做一次重叠分析
overlap_all <- extract_overlap_genes(venn_gene_lists, out_dir, "all")

############################################################
# 7. 绘制Venn图（Cell风格）
############################################################

## 7.1 所有显著DEGs的Venn图 -----------------------------------------
p_venn_all <- ggVennDiagram(
  venn_gene_lists,
  label       = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#FFF5EB", high = "#FD8D3C") +  # 浅橙到深橙
  scale_color_manual(values = rep("#D94701", length(venn_gene_lists)))  # 深红色 # 标准学术蓝 +
  labs(
    title    = "Differentially Expressed Genes at Each Time Point",
    subtitle = paste("MPXV_lb vs mock (",
                     paste(names(venn_gene_lists), collapse = ", "),
                     ")", sep = "")
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    plot.margin   = margin(5, 5, 5, 5, "mm")
  )





# PDF（矢量）& 600 dpi TIFF
pdf(file.path(out_dir, "Venn_All_DEGs.pdf"), width = 6, height = 5)
print(p_venn_all)
dev.off()

tiff(file.path(out_dir, "Venn_All_DEGs_600dpi.tiff"),
     width = 6, height = 5, units = "in", res = 600, compression = "lzw")
print(p_venn_all)
dev.off()

cat("\n已输出：Venn_All_DEGs.pdf / Venn_All_DEGs_600dpi.tiff\n")

## 7.2 上调基因Venn图 --------------------------------------------------
if (length(up_gene_lists) >= 2) {
  overlap_up <- extract_overlap_genes(up_gene_lists, out_dir, "Up")
  
  p_venn_up <- ggVennDiagram(
    up_gene_lists,
    label       = "count",
    label_alpha = 0
  ) +
    scale_fill_gradient(low = "#FDEDEC", high = "#C0392B") +  # 红色系
    scale_color_manual(values = rep("#C0392B", length(up_gene_lists))) +
    labs(title = "Up-regulated DEGs Overlap") +
    theme_void() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  pdf(file.path(out_dir, "Venn_Up_DEGs.pdf"), width = 6, height = 5)
  print(p_venn_up)
  dev.off()
  
  tiff(file.path(out_dir, "Venn_Up_DEGs_600dpi.tiff"),
       width = 6, height = 5, units = "in", res = 600, compression = "lzw")
  print(p_venn_up)
  dev.off()
  
  cat("已输出：Venn_Up_DEGs.pdf / Venn_Up_DEGs_600dpi.tiff\n")
}

## 7.3 下调基因Venn图 ---------------------------------------------------
if (length(down_gene_lists) >= 2) {
  overlap_down <- extract_overlap_genes(down_gene_lists, out_dir, "Down")
  
  p_venn_down <- ggVennDiagram(
    down_gene_lists,
    label       = "count",
    label_alpha = 0
  ) +
    scale_fill_gradient(low = "#EBF5FB", high = "#21618C") +  # 蓝色系
    scale_color_manual(values = rep("#21618C", length(down_gene_lists))) +
    labs(title = "Down-regulated DEGs Overlap") +
    theme_void() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  pdf(file.path(out_dir, "Venn_Down_DEGs.pdf"), width = 6, height = 5)
  print(p_venn_down)
  dev.off()
  
  tiff(file.path(out_dir, "Venn_Down_DEGs_600dpi.tiff"),
       width = 6, height = 5, units = "in", res = 600, compression = "lzw")
  print(p_venn_down)
  dev.off()
  
  cat("已输出：Venn_Down_DEGs.pdf / Venn_Down_DEGs_600dpi.tiff\n")
}

############################################################
# 8. 重叠模式柱状图（用于补充图）
############################################################

if (!is.null(overlap_all$overlap_counts)) {
  overlap_summary <- overlap_all$overlap_counts %>%
    mutate(
      pattern = apply(select(., -n), 1, function(x) {
        labs <- names(venn_gene_lists)[as.logical(x)]
        if (length(labs) == 0) "None" else paste(labs, collapse = " & ")
      })
    ) %>%
    arrange(desc(n))
  
  p_overlap_bar <- ggplot(overlap_summary,
                          aes(x = reorder(pattern, n), y = n)) +
    geom_col(width = 0.6, fill = "#2E86AB", alpha = 0.8) +
    coord_flip() +
    labs(
      title = "Number of DEGs in Each Overlap Pattern",
      x     = "Overlap pattern",
      y     = "Gene count"
    ) +
    theme_cell() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  pdf(file.path(out_dir, "Overlap_Pattern_Barplot.pdf"),
      width = 6, height = 4)
  print(p_overlap_bar)
  dev.off()
  
  tiff(file.path(out_dir, "Overlap_Pattern_Barplot_600dpi.tiff"),
       width = 6, height = 4, units = "in", res = 600, compression = "lzw")
  print(p_overlap_bar)
  dev.off()
  
  cat("已输出：Overlap_Pattern_Barplot.pdf / Overlap_Pattern_Barplot_600dpi.tiff\n")
}

############################################################
# 9. 保存分析结果
############################################################

# 保存DEGs列表
for(tp in timepoints) {
  if(tp %in% names(deg_list)) {
    write.csv(deg_list[[tp]], 
              file.path(out_dir, paste0("DEGs_", tp, "_vs_mock.csv")),
              row.names = FALSE)
  }
}

# 保存会话信息
save.image(file.path(out_dir, "Venn_Analysis_Workspace.RData"))

cat("\n=== 分析完成！ ===\n")
cat("输出目录：", normalizePath(out_dir), "\n")
cat("生成的文件：\n")
cat("  - Venn_All_DEGs.pdf / .tiff (所有DEGs韦恩图)\n")
if(length(up_gene_lists) >= 2) cat("  - Venn_Up_DEGs.pdf / .tiff (上调基因韦恩图)\n")
if(length(down_gene_lists) >= 2) cat("  - Venn_Down_DEGs.pdf / .tiff (下调基因韦恩图)\n")
cat("  - Overlap_Pattern_Barplot.pdf / .tiff (重叠模式柱状图)\n")
cat("  - 各时间点DEGs列表CSV文件\n")
cat("  - 重叠模式统计CSV文件\n")
cat("  - 分析工作空间RData文件\n")