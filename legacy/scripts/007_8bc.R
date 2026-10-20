




library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(grid)

## 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop/Project/10-干扰素小鼠/1-源/8bc_counts_gene", 
                     header = TRUE, sep = "\t", 
                     comment.char = "#", row.names = 1)

cat("原始数据维度:", dim(counts), "\n")

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

## 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  timepoint = factor(c("10d", rep("3d", 3), rep("5d", 3), 
                       rep("7d", 3), rep("mock", 3)),
                     levels = c("mock", "3d", "5d", "7d", "10d")),
  condition = factor(c(rep("MPXV_lb", 10), rep("mock", 3)),
                     levels = c("mock", "MPXV_lb")),
  group_label = c(rep("MPXV_lb", 10), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

print("样本信息:")
print(sample_info)

## 3. 创建DESeq2对象并处理数据 --------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ condition
)

print("DESeq2对象创建成功！")
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 运行DESeq2分析
dds <- DESeq(dds)

# 获取vst标准化数据
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)





## 4. 为“各时间点 vs mock”创建分组变量并重新跑一次DESeq2 -----------------

# 创建一个用于对比的组：mock / 3d / 5d / 7d / 10d
sample_info$group_tp <- factor(as.character(sample_info$timepoint),
                               levels = c("mock", "3d", "5d", "7d", "10d"))

# 用同一套count矩阵重新构建dds_tp（不影响你原来的dds对象）
dds_tp <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ group_tp
)

# 同样过滤低表达（保持你原逻辑：>=10 的count至少出现在3个样本）
keep_tp <- rowSums(counts(dds_tp) >= 10) >= 3
dds_tp <- dds_tp[keep_tp, ]

dds_tp <- DESeq(dds_tp)

## 5. 生成 deg_list/timepoints 和阈值（供火山图循环直接使用） -------------
timepoints <- c("3d", "5d", "7d", "10d")

# 阈值：如果你之前脚本已有就不要重复定义；这里给一个常用默认
FC_cutoff <- 1
P_cutoff <- 0.05

deg_list <- list()



for (tp in timepoints) {
  # tp vs mock
  res <- results(dds_tp, contrast = c("group_tp", tp, "mock"))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    # 去掉NA，避免-log10(padj) 出现Inf/NA
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      change = case_when(
        padj < P_cutoff & log2FoldChange >= FC_cutoff  ~ "Up-regulated",
        padj < P_cutoff & log2FoldChange <= -FC_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      )
    )
  
  deg_list[[tp]] <- res_df
}







# 7.2 各时间点火山图
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
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
    
    print(p_vol)
    # 修改ggsave调用：添加dpi=600和指定路径
    ggsave(paste0("D:/Desktop/1/", tp, "_vs_Mock.pdf"), p_vol, width = 4, height = 4.5, dpi = 800)
  }
}
# 7.2 各时间点火山图
for(tp in timepoints) {
  res_df <- deg_list[[tp]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
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
    
    print(p_vol)
    # 修改ggsave调用：添加dpi=600和指定路径
    ggsave(paste0("D:/Desktop/1/", tp, "_vs_Mock.pdf"), p_vol, width = 4, height = 4.5, dpi = 800)
  }
}
















