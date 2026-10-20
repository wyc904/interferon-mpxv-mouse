# ==============================================================================
# 4. 进阶分析环境加载
# ==============================================================================
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠注释库
library(enrichplot)
library(pathview)
library(ggplot2)
library(ggrepel)


# install.packages("pathview",repos="http://R-Forge.R-project.org")


# 创建结果保存目录
dir.create("Advanced_Analysis", showWarnings = FALSE)

# ==============================================================================
# 5. ID转换 (Symbol -> Entrez ID)
# ==============================================================================
# GSEA和KEGG需要Entrez ID。我们将转换后的ID合并到结果表中。
cat(">>> 正在进行基因ID转换...\n")

# 移除版本号(如果有)并转换
genes <- rownames(res_df)
entrez_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# 将Entrez ID合并回 res_df
# 注意：这可能会因为多对一映射导致行数略微变化，或者去除无法匹配的基因
res_df_anno <- merge(res_df, entrez_ids, by.x="gene_id", by.y="SYMBOL")

# 去除重复的 Entrez ID (保留 padj 最小的那个)
res_df_anno <- res_df_anno %>% 
  arrange(padj) %>% 
  distinct(ENTREZID, .keep_all = TRUE)

cat(">>> ID转换完成。用于功能分析的基因数:", nrow(res_df_anno), "\n")

# ==============================================================================
# 6. MA Plot (M-vs-A plot)
# ==============================================================================
# 这种图用于检查数据是否存在系统性偏差（通常在Y=0轴上下对称）
# X轴: 平均表达量 (log10 Mean), Y轴: 差异倍数 (Log2FC)

cat(">>> 绘制 MA Plot...\n")

p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  # 绘制所有点（灰色）
  geom_point(col = "grey80", size = 0.8, alpha = 0.5) +
  # 绘制显著差异点
  geom_point(data = subset(res_df, change != "Not Significant"), 
             aes(color = change), size = 1.2, alpha = 0.8) +
  # 添加基准线
  geom_hline(yintercept = 0, col = "black", size = 0.5) +
  # 标记 Top 10 显著基因
  geom_text_repel(data = subset(res_df, gene_id %in% top_genes[1:10]), 
                  aes(label = gene_id), size = 3, max.overlaps = 20) +
  # 颜色设置
  scale_color_manual(values = c("Down-regulated" = color_down, "Up-regulated" = color_up)) +
  labs(title = "MA Plot: sta1.MPXV_lb vs mock", 
       x = "Log10 Mean Expression", 
       y = "Log2 Fold Change") +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_ma)
# ggsave("Advanced_Analysis/Figure_MA_Plot.pdf", p_ma, width = 6, height = 5)

# ==============================================================================
# 7. Gene Boxplots / Dotplots (单基因表达量图)
# ==============================================================================
# 这里我们自动选取差异最显著的 Top 4 基因进行展示
# 你也可以手动指定基因，例如: target_genes <- c("Stat1", "Mx1", "Il6", "Tnf")

cat(">>> 绘制单基因表达量图...\n")

# 选取 padj 最小的4个基因
target_genes <- res_df %>% arrange(padj) %>% head(4) %>% pull(gene_id)
# 如果你想看特定基因，请取消注释下一行并修改基因名：
# target_genes <- c("Stat1", "Ifnb1", "Isg15", "Ccl2") 

# 提取这些基因的标准化 Count 数据
plot_data_list <- list()
for(g in target_genes){
  # 检查基因是否存在
  if(g %in% rownames(dds)){
    d <- plotCounts(dds, gene=g, intgroup="condition", returnData=TRUE)
    d$gene <- g
    plot_data_list[[g]] <- d
  }
}
combined_plot_data <- do.call(rbind, plot_data_list)

if(!is.null(combined_plot_data)){
  p_gene <- ggplot(combined_plot_data, aes(x = condition, y = count, fill = condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2) +
    facet_wrap(~gene, scales = "free_y") + # 分面展示
    scale_fill_manual(values = c("Control" = color_control, "Infection" = color_infect)) +
    labs(title = "Top Significant Genes Expression", y = "Normalized Counts") +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_rect(fill="white"))
  
  print(p_gene)
  ggsave("Advanced_Analysis/Figure_Gene_Boxplots.pdf", p_gene, width = 8, height = 6)
}

# ==============================================================================
# 8. GSEA (Gene Set Enrichment Analysis) - KEGG
# ==============================================================================
# GSEA 不需要阈值筛选，而是使用所有基因的排序列表
cat(">>> 运行 GSEA 分析 (基于 KEGG 通路)...\n")

# 8.1 准备排序列表 (Ranked List)
# 排序指标通常使用 Log2FC，降序排列
geneList <- res_df_anno$log2FoldChange
names(geneList) <- res_df_anno$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# 8.2 运行 GSEA
# organism = 'mmu' 代表小鼠 (mouse), 'hsa' 代表人类
# pvalueCutoff 放宽到 1 以确保能跑出结果，后续再筛选
gsea_kegg <- gseKEGG(geneList     = geneList,
                     organism     = 'mmu',
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     verbose      = FALSE)

# 8.3 可视化 GSEA 结果
if(nrow(gsea_kegg) > 0){
  # 仅保留显著的通路
  gsea_sig <- filter(gsea_kegg, p.adjust < 0.25) # GSEA 推荐阈值有时较宽，可设为 0.05
  
  if(nrow(gsea_sig) > 0) {
    # 绘制山峰图 (选取 NES 分数最高的通路)
    top_pathway_id <- gsea_sig@result$ID[1]
    top_pathway_desc <- gsea_sig@result$Description[1]
    
    p_gsea <- gseaplot2(gsea_kegg, geneSetID = top_pathway_id, 
                        title = paste("GSEA:", top_pathway_desc),
                        pvalue_table = TRUE)
    
    print(p_gsea)
    ggsave("Advanced_Analysis/Figure_GSEA_Top1.pdf", p_gsea, width = 8, height = 6)
    
    # 保存结果表
    write.csv(as.data.frame(gsea_kegg), "Advanced_Analysis/GSEA_KEGG_Results.csv")
  } else {
    cat("GSEA 未发现显著富集的通路 (padj < 0.25)。\n")
  }
} else {
  cat("GSEA 分析未返回结果，可能是基因ID匹配率过低。\n")
}












# ==============================================================================
# 9. KEGG Pathway Map (通路映射图)
# ==============================================================================
# 使用 pathview 将 Log2FC 映射到具体的 KEGG 通路图上
cat(">>> 绘制 KEGG Pathway Map...\n")

if(nrow(gsea_kegg) > 0){
  # 自动选择最显著的那条通路进行绘制
  # 或者你可以手动指定，例如 "mmu04620" (Toll-like receptor)
  target_pathway_id <- gsea_kegg@result$ID[1] 
  
  # 注意：pathview 会直接在当前工作目录生成 .png 文件
  # out.suffix 设置文件后缀
  tryCatch({
    pathview(gene.data  = geneList,
             pathway.id = target_pathway_id,
             species    = "mmu",
             limit      = list(gene=max(abs(geneList)), cpd=1),
             out.suffix = "pathview_output",
             kegg.dir   = ".") # 这里的 . 代表当前目录，确保有写入权限
    
    cat(paste0("Pathview 图已生成: ", target_pathway_id, ".pathview_output.png\n"))
    
    # 移动文件到结果目录
    file.rename(paste0(target_pathway_id, ".pathview_output.png"), 
                paste0("Advanced_Analysis/", target_pathway_id, "_Pathview.png"))
    
  }, error = function(e){
    cat("Pathview 绘图失败，可能是网络连接问题或通路ID无效。\n")
  })
}

cat(">>> 所有进阶分析完成！请检查 Advanced_Analysis 文件夹。\n")




















