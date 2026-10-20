

# 加载包
library(edgeR)
library(limma)
library(ggplot2)
# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# 1. 读取并整理计数数据 -------------------------------------------------
counts <- read.delim("D:/Desktop//Project/10-干扰素小鼠/1-源/STA1_counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 查看原始列名，确认样本顺序
print("原始列名（第7-15列）:")
print(colnames(counts)[1:20])

# 提取count矩阵（第7-15列是我们要分析的样本）
count_matrix <- counts[, 12:20]

# 根据您提供的信息，正确重命名列名
# 第7-9列: MPXV_IIb样本
# 第10-12列: MPXV_Ib样本  
# 第13-15列: MOCK样本

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("MPXV_IIb_1", "MPXV_IIb_2", "MPXV_IIb_3",
                            "MPXV_Ib_1", "MPXV_Ib_2", "MPXV_Ib_3",
                            "MOCK_1", "MOCK_2", "MOCK_3")

print("重命名后的列名:")
print(colnames(count_matrix))

# 验证样本对应关系
sample_mapping <- data.frame(
  Original_Name = colnames(counts)[12:20],
  New_Name = colnames(count_matrix),
  Group = c(rep("MPXV_IIb", 3), rep("MPXV_Ib", 3), rep("MOCK", 3))
)

print("样本对应关系:")
print(sample_mapping)

# 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  group = factor(c(rep("MPXV_IIb", 3), rep("MPXV_Ib", 3), rep("MOCK", 3)),
                 levels = c("MOCK", "MPXV_Ib", "MPXV_IIb")),
  replicate = rep(1:3, 3)
)

rownames(sample_info) <- sample_info$sample

print("样本信息:")
print(sample_info)

# 3. 创建DESeq2对象 ----------------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ group
)

print("DESeq2对象创建成功！")
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行DESeq2分析 ----------------------------------------------------
dds <- DESeq(dds)
print("DESeq2分析完成！")

# 后续可以继续添加差异分析等代码...








# 6. 进行差异分析比较 ----------------------------------------------
comparisons <- list(
  c("MPXV_IIb", "MOCK"),
  c("MPXV_Ib", "MOCK"),
  c("MPXV_IIb", "MPXV_Ib")
)

results_list <- list()
deg_list <- list()

# 定义差异基因阈值
FC_cutoff <- 1.0  # log2FoldChange > 1 (即2倍差异)
P_cutoff <- 0.05  # padj < 0.05

for(comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  comp_name <- paste0(group1, "_vs_", group2)
  
  print(paste("分析比较:", comp_name))
  
  # 确保参考水平正确
  dds$group <- relevel(dds$group, ref = group2)
  
  # 重新运行DESeq2
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", group1, group2), alpha = 0.05)
  
  # 保存完整结果
  results_list[[comp_name]] <- res
  
  # 提取显著差异基因并标记上下调
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      change = case_when(
        padj < P_cutoff & log2FoldChange > FC_cutoff ~ "Up-regulated",
        padj < P_cutoff & log2FoldChange < -FC_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      change = factor(change, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
    ) %>%
    arrange(padj)
  
  deg_list[[comp_name]] <- res_df
  print(paste(comp_name, ": 发现", sum(res_df$change != "Not Significant"), "个显著DEGs"))
}

# 7. 可视化 --------------------------------------------------------

# 全局配色设置
color_mock <- "#708090"    # 灰色
color_IIb <- "#A52A2A"     # 深棕色 - MPXV IIb
color_Ib <- "#8B4513"      # 马鞍棕色 - MPXV Ib
color_up <- "#C75D69"      # 浅红
color_down <- "#4682B4"    # 蓝色
color_ns <- "#D3D3D3"      # 浅灰


# # 检查文件的前几行
# Y_A_B_counts<-readLines("D:/Desktop/Project/生信/3-差异基因分析/out/A-B-counts.txt", n = 10)
# # 读取 featureCounts 生成的文件
# counts <- Y_A_B_counts
# 
# head(Y_A_B_counts)
# 
# 
# 
# count_data <- counts[, 7:ncol(counts)] 
# 
# 
# D:\Desktop\10-干扰素小鼠\1-源
# 

counts <- read.delim("D:/Desktop/Project/10-干扰素小鼠/1-源/STA1_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# print(dim(counts))
# counts2 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/AGB6_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# # print(dim(counts))
# counts3 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# # print(dim(counts))
# 
# 




colnames(counts)






# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:20]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix))

# 创建样本信息表
sample_info <- data.frame(
  sample = simplified_names,
  group = sapply(strsplit(simplified_names, "_"), function(x) paste(x[2:3], collapse = "_")),
  replicate = sapply(strsplit(simplified_names, "_"), function(x) x[4]),
  stringsAsFactors = FALSE
)

cat("\n=== 样本信息 ===\n")
print(sample_info)

# 检查count数据的统计信息
cat("\n=== Count数据统计 ===\n")
cat("总基因数:", nrow(count_matrix), "\n")
cat("总测序reads数:", sum(count_matrix), "\n")
cat("每个样本的reads数:\n")
print(colSums(count_matrix))

# 检查是否有零表达基因
zero_genes <- rowSums(count_matrix) == 0
cat("零表达基因数:", sum(zero_genes), "\n")
cat("占比:", round(mean(zero_genes) * 100, 2), "%\n")
# 7.2 各比较组火山图（拆分版本）

# 7.2.1 MPXV_IIb_vs_MOCK 火山图
if("MPXV_IIb_vs_MOCK" %in% names(deg_list)) {
  res_df <- deg_list[["MPXV_IIb_vs_MOCK"]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol_IIb_vs_MOCK <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
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
      labs(title = "Volcano: MPXV IIb vs MOCK", 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      xlim(-13, 13) +  # 添加x轴范围限制
      theme_bw() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_vol_IIb_vs_MOCK)
    ggsave("Volcano_MPXV_IIb_vs_MOCK.pdf", p_vol_IIb_vs_MOCK, width = 7, height = 6)
    cat("已生成 MPXV IIb vs MOCK 火山图\n")
  }
}






# 7.2.2 MPXV_Ib_vs_MOCK 火山图
if("MPXV_Ib_vs_MOCK" %in% names(deg_list)) {
  res_df <- deg_list[["MPXV_Ib_vs_MOCK"]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol_Ib_vs_MOCK <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
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
      labs(title = "Volcano: MPXV Ib vs MOCK", 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value)) +
      xlim(-15, 15) +  # 设置x轴范围为-15到15
      ylim(0, 200) +   # 设置y轴范围为0到200
      theme_bw() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_vol_Ib_vs_MOCK)
    ggsave("Volcano_MPXV_Ib_vs_MOCK.pdf", p_vol_Ib_vs_MOCK, width = 7, height = 6)
    cat("已生成 MPXV Ib vs MOCK 火山图\n")
  }
}







# 7.2.3 MPXV_IIb_vs_MPXV_Ib 火山图
if("MPXV_IIb_vs_MPXV_Ib" %in% names(deg_list)) {
  res_df <- deg_list[["MPXV_IIb_vs_MPXV_Ib"]]
  
  if(nrow(res_df) > 0) {
    # 选择最显著的基因进行标记
    top_genes <- res_df %>% 
      filter(change != "Not Significant") %>%
      arrange(padj) %>% 
      head(15) %>% 
      pull(gene)
    
    p_vol_IIb_vs_Ib <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
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
      labs(title = "Volcano: MPXV IIb vs MPXV Ib", 
           x = expression(Log[2]~Fold~Change), 
           y = expression(-Log[10]~Adjusted~P-value))+
      xlim(-13, 13) +  # 设置x轴范围为-15到15 +
      theme_bw() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_vol_IIb_vs_Ib)
    ggsave(paste0("D:/Desktop/1/", tp, "_vs_Mock.pdf"), p_vol, width = 4, height = 4.5, dpi = 800)
    cat("已生成 MPXV IIb vs MPXV Ib 火山图\n")
  }
}





# 检查是否所有比较组都成功生成了火山图
generated_plots <- c()
if(exists("p_vol_IIb_vs_MOCK")) generated_plots <- c(generated_plots, "MPXV_IIb_vs_MOCK")
if(exists("p_vol_Ib_vs_MOCK")) generated_plots <- c(generated_plots, "MPXV_Ib_vs_MOCK")
if(exists("p_vol_IIb_vs_Ib")) generated_plots <- c(generated_plots, "MPXV_IIb_vs_MPXV_Ib")

cat("成功生成的火山图:", paste(generated_plots, collapse = ", "), "\n")


















