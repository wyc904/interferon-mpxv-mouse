# 加载必要的包
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)    # 用于富集分析
library(org.Mm.eg.db)       # 小鼠基因注释数据库
library(DOSE)               # 疾病本体分析

counts <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/C57-counts_gene.txt", 
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix <- counts[, 6:20]

# 重命名列名为简洁的样本名
colnames(count_matrix) <- c("C57_10d_1_P", "C57_10d_2_P", "C57_10d_3_P",
                            "C57_3d_1_P", "C57_3d_2_P", "C57_3d_3_P",
                            "C57_5d_1_P", "C57_5d_2_P", "C57_5d_3_P",
                            "C57_7d_1_P", "C57_7d_2_P", "C57_7d_3_P",
                            "c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P")

print("重命名后的列名:")
print(colnames(count_matrix))

# 2. 创建样本信息数据框 -------------------------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  timepoint = factor(c(rep("10d", 3), rep("3d", 3), rep("5d", 3), 
                       rep("7d", 3), rep("mock", 3)),
                     levels = c("mock", "3d", "5d", "7d", "10d")),
  condition = factor(c(rep("MPXV_lb", 12), rep("mock", 3)),
                     levels = c("mock", "MPXV_lb")),
  group_label = c(rep("MPXV_lb", 12), rep("MOCK", 3))
)

rownames(sample_info) <- sample_info$sample

print("样本信息:")
print(sample_info)

# 3. 创建DESeq2对象 ----------------------------------------------------
# 确保计数矩阵只包含数值数据
count_matrix_numeric <- as.matrix(count_matrix)
mode(count_matrix_numeric) <- "numeric"

# 创建DESeq2数据集 - 使用timepoint作为主要因子
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_numeric,
  colData = sample_info,
  design = ~ timepoint  # 修改为以timepoint为主要因子
)

print("DESeq2对象创建成功！")
print(paste("基因数:", nrow(dds)))
print(paste("样本数:", ncol(dds)))

# 4. 过滤低表达基因 ----------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3  # 至少在3个样本中表达量≥10
dds <- dds[keep, ]
print(paste("过滤后保留", nrow(dds), "个基因"))

# 5. 运行LRT (似然比检验) 分析 -----------------------------------------
# 使用LRT来识别随时间点变化的基因
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)  # 简化模型为截距模型

# 提取LRT结果
lrt_results <- results(dds_lrt, alpha = 0.01)  # 使用更严格的阈值
print("LRT分析完成！")

# 筛选显著的时间依赖性基因 (padj < 0.01)
sig_genes <- rownames(lrt_results)[which(lrt_results$padj < 0.01 & !is.na(lrt_results$padj))]
print(paste("发现", length(sig_genes), "个随时间显著变化的基因 (padj < 0.01)"))


# 如果基因数太多，可以使用更严格的阈值
if(length(sig_genes) > 2000) {
  sig_genes_strict <- rownames(lrt_results)[which(lrt_results$padj < 0.001 & !is.na(lrt_results$padj))]
  print(paste("使用更严格阈值 (padj < 0.001) 得到", length(sig_genes_strict), "个基因"))
  sig_genes <- sig_genes_strict
}





# 6. 生物学功能过滤 ---------------------------------------------------
# 6.1 基因ID转换（假设当前基因ID是ENSEMBL ID）
# 如果基因ID是基因符号，可以跳过这一步
# 这里假设基因ID是ENSEMBL ID，需要转换为ENTREZID用于富集分析

# 首先尝试识别基因ID类型
gene_ids <- rownames(counts)[1:10]  # 查看前10个基因ID
print("前10个基因ID:")
print(gene_ids)

# 函数：尝试转换基因ID
convert_gene_ids <- function(genes) {
  # 尝试不同的ID转换方法
  tryCatch({
    # 如果是ENSEMBL ID（以ENSMUSG开头）
    if(any(grepl("^ENSMUSG", genes[1:10]))) {
      # 使用biomaRt进行转换
      if(!require(biomaRt)) {
        BiocManager::install("biomaRt")
        library(biomaRt)
      }
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      id_map <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "mgi_symbol"),
                      filters = "ensembl_gene_id", 
                      values = genes, 
                      mart = mart)
      return(id_map)
    } else {
      # 如果是基因符号，直接使用
      return(data.frame(mgi_symbol = genes, entrezgene_id = NA))
    }
  }, error = function(e) {
    print("基因ID转换失败，将使用原始基因符号进行富集分析")
    return(data.frame(mgi_symbol = genes, entrezgene_id = NA))
  })
}

# 转换显著基因的ID
sig_gene_symbols <- sig_genes  # 默认使用原始ID

# 6.2 GO/KEGG富集分析
print("开始GO/KEGG富集分析...")

# 目标免疫相关通路
target_go_bp <- c(
  "immune system process",
  "response to virus", 
  "cytokine-mediated signaling pathway",
  "leukocyte migration"
)

target_kegg <- c(
  "Cytokine-cytokine receptor interaction",
  "Chemokine signaling pathway", 
  "Toll-like receptor signaling pathway",
  "T cell receptor signaling pathway"
)

# GO富集分析
go_enrich <- enrichGO(gene          = sig_gene_symbols,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = "SYMBOL",  # 如果是ENSEMBL ID需要修改
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

# KEGG富集分析
kegg_enrich <- enrichKEGG(gene         = sig_gene_symbols,
                          organism     = 'mmu',
                          pvalueCutoff = 0.05,
                          keyType      = "kegg")

print("富集分析完成！")

# 6.3 筛选目标通路的基因
immune_genes <- c()

# 从GO结果中提取目标通路的基因
if(!is.null(go_enrich)) {
  go_results <- as.data.frame(go_enrich)
  
  # 筛选目标GO通路
  target_go_results <- go_results %>%
    filter(grepl(paste(target_go_bp, collapse = "|"), 
                 Description, ignore.case = TRUE))
  
  if(nrow(target_go_results) > 0) {
    # 提取这些通路中的基因
    for(i in 1:nrow(target_go_results)) {
      genes_in_pathway <- strsplit(target_go_results$geneID[i], "/")[[1]]
      immune_genes <- unique(c(immune_genes, genes_in_pathway))
    }
    print(paste("从GO通路中找到", length(immune_genes), "个免疫相关基因"))
  }
}

# 从KEGG结果中提取目标通路的基因
if(!is.null(kegg_enrich)) {
  kegg_results <- as.data.frame(kegg_enrich)
  
  # 筛选目标KEGG通路
  target_kegg_results <- kegg_results %>%
    filter(grepl(paste(target_kegg, collapse = "|"), 
                 Description, ignore.case = TRUE))
  
  if(nrow(target_kegg_results) > 0) {
    # 提取这些通路中的基因
    for(i in 1:nrow(target_kegg_results)) {
      genes_in_pathway <- strsplit(target_kegg_results$geneID[i], "/")[[1]]
      immune_genes <- unique(c(immune_genes, genes_in_pathway))
    }
    print(paste("从KEGG通路中找到", length(immune_genes), "个免疫相关基因"))
  }
}

# 如果通过富集分析找到的基因太少，使用已知的免疫基因列表
if(length(immune_genes) < 50) {
  print("通过富集分析找到的基因较少，使用预定义的免疫基因列表...")
  
  # 预定义的小鼠免疫相关基因列表
  predefined_immune_genes <- c(
    "Il6", "Tnf", "Ifnb1", "Ifna1", "Cxcl1", "Cxcl2", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
    "Irf1", "Irf3", "Irf7", "Stat1", "Stat2", "Stat3", "Jak1", "Jak2", "Tyk2",
    "Tlr2", "Tlr3", "Tlr4", "Tlr7", "Tlr8", "Tlr9", "Myd88", "Traf6", "Nfkb1", "Rela",
    "Cd4", "Cd8a", "Cd8b1", "Cd3e", "Cd3g", "Cd3d", "Cd19", "Cd20", "Cd40", "Cd40lg",
    "Il2", "Il4", "Il10", "Il12a", "Il12b", "Il17a", "Il21", "Il23a", "Ifng", "Tgfb1",
    "B2m", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "H2-DMb2"
  )
  
  # 取显著基因与预定义免疫基因的交集
  immune_genes <- intersect(sig_gene_symbols, predefined_immune_genes)
  print(paste("使用预定义列表找到", length(immune_genes), "个免疫相关基因"))
}

# 7. 保存结果 ---------------------------------------------------------
output_dir <- "LRT_TimeCourse_Analysis"
dir.create(output_dir, showWarnings = FALSE)

# 保存LRT结果
lrt_df <- as.data.frame(lrt_results) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write.csv(lrt_df, file.path(output_dir, "LRT_time_dependent_genes_full.csv"), row.names = FALSE)






# 保存显著基因
sig_genes_df <- data.frame(gene = sig_genes)
write.csv(sig_genes_df, file.path(output_dir, "Significant_time_dependent_genes.csv"), row.names = FALSE)

# 保存免疫相关基因
immune_genes_df <- data.frame(gene = immune_genes)
write.csv(immune_genes_df, file.path(output_dir, "Immune_related_time_dependent_genes.csv"), row.names = FALSE)

# 保存富集分析结果
if(!is.null(go_enrich)) {
  write.csv(as.data.frame(go_enrich), file.path(output_dir, "GO_enrichment_results.csv"), row.names = FALSE)
}

if(!is.null(kegg_enrich)) {
  write.csv(as.data.frame(kegg_enrich), file.path(output_dir, "KEGG_enrichment_results.csv"), row.names = FALSE)
}

# 8. 可视化结果 ------------------------------------------------------
# 8.1 LRT结果MA图
ma_data <- as.data.frame(lrt_results) %>%
  rownames_to_column("gene") %>%
  mutate(significant = padj < 0.01 & !is.na(padj),
         immune_related = gene %in% immune_genes)

p_ma <- ggplot(ma_data, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significant, alpha = significant), size = 1) +
  scale_color_manual(values = c("grey", "red")) +
  scale_alpha_manual(values = c(0.3, 0.7)) +
  scale_x_log10() +
  labs(title = "MA Plot: Time-dependent Genes (LRT)",
       x = "Mean Expression", y = "Log2 Fold Change") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

if(sum(ma_data$immune_related) > 0) {
  p_ma <- p_ma + 
    geom_point(data = subset(ma_data, immune_related), 
               color = "blue", size = 1.5, alpha = 0.8) +
    geom_text_repel(data = subset(ma_data, immune_related & padj < 0.001),
                    aes(label = gene), size = 3, max.overlaps = 20)
}

print(p_ma)
ggsave(file.path(output_dir, "MA_plot_LRT.pdf"), p_ma, width = 10, height = 8)

# 8.2 富集分析结果气泡图
if(!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
  p_go <- dotplot(go_enrich, showCategory = 15) +
    labs(title = "GO Biological Process Enrichment")
  print(p_go)
  ggsave(file.path(output_dir, "GO_enrichment_dotplot.pdf"), p_go, width = 10, height = 8)
}

if(!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
  p_kegg <- dotplot(kegg_enrich, showCategory = 15) +
    labs(title = "KEGG Pathway Enrichment")
  print(p_kegg)
  ggsave(file.path(output_dir, "KEGG_enrichment_dotplot.pdf"), p_kegg, width = 10, height = 8)
}

# 9. 结果总结 ---------------------------------------------------------
cat("=== LRT时间序列分析完成 ===\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总基因数:", nrow(dds), "\n")
cat("LRT显著基因 (padj < 0.01):", length(sig_genes), "\n")
cat("免疫相关基因:", length(immune_genes), "\n\n")

cat("目标免疫通路基因统计:\n")
cat("- 免疫系统过程相关基因:", sum(grepl("immune", immune_genes, ignore.case = TRUE)), "\n")
cat("- 细胞因子相关基因:", sum(grepl("il|ifn|cxcl|ccl", immune_genes, ignore.case = TRUE)), "\n")
cat("- 信号通路相关基因:", sum(grepl("stat|jak|tlr|irf", immune_genes, ignore.case = TRUE)), "\n")

# 保存工作空间
save.image(file.path(output_dir, "LRT_time_course_analysis.RData"))

cat("\n分析完成！结果保存在:", output_dir, "目录中\n")
cat("主要输出文件:\n")
cat("- LRT_time_dependent_genes_full.csv: 完整的LRT分析结果\n")
cat("- Significant_time_dependent_genes.csv: 显著时间依赖性基因列表\n")
cat("- Immune_related_time_dependent_genes.csv: 免疫相关基因列表（用于Mfuzz分析）\n")
cat("- 富集分析结果和可视化图表\n")






















