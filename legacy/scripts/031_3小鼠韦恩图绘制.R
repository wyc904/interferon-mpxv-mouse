# 需要的包
library(edgeR)
library(limma)
library(ggplot2)

library(ggVennDiagram)



# 加载包
library(edgeR)
library(limma)
library(ggplot2)


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

counts3 <- read.delim("D:/Desktop/10-干扰素小鼠/1-源/STA1_counts_gene.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)
# print(dim(counts))







# 提取count矩阵（从第6列开始是真正的count数据）
count_matrix3 <- counts3[, 6:20]

# 简化列名（去掉复杂的路径信息）
sample_names <- colnames(count_matrix3)
simplified_names <- gsub("STAR_aligned\\.trimmomatic\\.|\\.Aligned\\.sortedByCoord\\.out\\.bam", "", sample_names)
colnames(count_matrix3) <- simplified_names

cat("简化后的样本名称:\n")
print(colnames(count_matrix3))








sta1_samples <- c(
  "JZ25055844.4.1.P.4.1.P_combined",
  "JZ25055845.4.2.P.4.2.P_combined",
  "JZ25055846.4.3.P.4.3.P_combined"
)

# =========================
# 2) 合并 count 矩阵（取三者共有基因，最稳）
# =========================
common_genes <- Reduce(intersect, list(rownames(count_matrix),
                                       rownames(count_matrix2),
                                       rownames(count_matrix3)))

count_all <- cbind(
  count_matrix[common_genes, , drop = FALSE],
  count_matrix2[common_genes, , drop = FALSE],
  count_matrix3[common_genes, , drop = FALSE]
)

# 按你想对比的顺序排列列
sample_order <- c(c57_samples, agb6_samples, sta1_samples)

# 检查样本名是否都在矩阵里
missing <- setdiff(sample_order, colnames(count_all))
if (length(missing) > 0) {
  stop("这些样本名在合并后的 count 矩阵中找不到：\n",
       paste(missing, collapse = "\n"),
       "\n\n请运行：print(colnames(count_all)) 检查列名是否完全一致。")
}

count_all <- count_all[, sample_order, drop = FALSE]

group <- factor(c(rep("C57", length(c57_samples)),
                  rep("AGB6", length(agb6_samples)),
                  rep("STA1", length(sta1_samples))),
                levels = c("C57", "AGB6", "STA1"))

# =========================
# 3) edgeR 预处理 + limma-voom 差异分析
# =========================
dge <- DGEList(counts = count_all, group = group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

cont <- makeContrasts(
  AGB6_vs_C57   = AGB6 - C57,
  STA1_vs_C57   = STA1 - C57,
  STA1_vs_AGB6  = STA1 - AGB6,
  levels = design
)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

# =========================
# 4) 提取上调/下调差异基因（阈值可改）
# =========================
get_deg_sets <- function(coef_name, padj = 0.05, lfc = 1) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  
  sig <- tt[!is.na(tt$adj.P.Val) &
              tt$adj.P.Val < padj &
              abs(tt$logFC) >= lfc, , drop = FALSE]
  
  list(
    up   = sig$gene[sig$logFC > 0],
    down = sig$gene[sig$logFC < 0],
    table = tt
  )
}

# 你可以在这里改阈值：padj / lfc
res_agb6_c57  <- get_deg_sets("AGB6_vs_C57",  padj = 0.05, lfc = 1)
res_sta1_c57  <- get_deg_sets("STA1_vs_C57",  padj = 0.05, lfc = 1)
res_sta1_agb6 <- get_deg_sets("STA1_vs_AGB6", padj = 0.05, lfc = 1)

up_sets <- list(
  AGB6_vs_C57  = res_agb6_c57$up,
  STA1_vs_C57  = res_sta1_c57$up,
  STA1_vs_AGB6 = res_sta1_agb6$up
)

down_sets <- list(
  AGB6_vs_C57  = res_agb6_c57$down,
  STA1_vs_C57  = res_sta1_c57$down,
  STA1_vs_AGB6 = res_sta1_agb6$down
)

# =========================
# 5) 画韦恩图（上调 / 下调各一张）
# =========================
# 设置颜色
group_colors <- c("#4E79A7", "#F28E2B", "#E15759")  # 蓝、橙、红



# 修改下面的Venn图，使其风格与上面保持一致
# 设置相同的颜色
group_colors_three <- c("#4E79A7", "#F28E2B", "#E15759")  # 蓝、橙、红

# 上调基因Venn图
p_up <- ggVennDiagram(
  up_sets,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#FDEDEC", high = "#C0392B") +  # 红色系
  scale_color_manual(values = group_colors_three) +
  labs(
    title = "",
    subtitle = ""
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# 下调基因Venn图
p_down <- ggVennDiagram(
  down_sets,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#EBF5FB", high = "#21618C") +  # 蓝色系
  scale_color_manual(values = group_colors_three) +
  labs(
    title = "",
    subtitle = ""
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# 输出图形
print(p_up)
print(p_down)

# # 可选：保存图形文件
# pdf("Upregulated_DEGs_Venn_Styled.pdf", width = 7, height = 6)
# print(p_up)
# dev.off()
# 
# pdf("Downregulated_DEGs_Venn_Styled.pdf", width = 7, height = 6)
# print(p_down)
# dev.off()
# 
# cat("已输出风格统一的Venn图：Upregulated_DEGs_Venn_Styled.pdf 和 Downregulated_DEGs_Venn_Styled.pdf\n")
print(p_up)
print(p_down)
# 




# # 可选：保存图片
# ggsave("Venn_Up_DEGs.png", p_up, width = 7, height = 6, dpi = 300)
# ggsave("Venn_Down_DEGs.png", p_down, width = 7, height = 6, dpi = 300)
# 
# # 可选：把每个对比的结果表导出
# write.csv(res_agb6_c57$table,  "DEG_AGB6_vs_C57.csv",  row.names = TRUE)
# write.csv(res_sta1_c57$table,  "DEG_STA1_vs_C57.csv",  row.names = TRUE)
# write.csv(res_sta1_agb6$table, "DEG_STA1_vs_AGB6.csv", row.names = TRUE)
