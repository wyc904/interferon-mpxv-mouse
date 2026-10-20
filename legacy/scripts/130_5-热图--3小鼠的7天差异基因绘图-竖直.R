# 加载必要的包
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(viridis)

# 设置工作目录
setwd("D:/Desktop")

# 1. 读取Excel文件
cat("正在读取Excel文件...\n")
data <- read.xlsx("绘制3体重热图.xlsx", sheet = 1)
cat("数据维度:", dim(data), "\n")

# 2. 检查列名
cat("列名:\n")
print(colnames(data))

# 3. 定义选定的基因列表（按指定顺序）
selected_genes_ordered <- c(
  "Irf8", "Jak2", "Ifnar2", "Oasl2", "Oasl1a", "Isg20", "Ifit1", "Bst2", "Mx1", "Dhx58", "Adar1", "Ifna", "Ifng","Cd3e", "Tnfsf13b", "Cxcl9", "Ccl2", "Cenpi", "Aif1", "Slc11a1", "B2m", "Vcam1", "Sele", "Lgals9","Serpinb9b", "Serpinb9c", "Mmp8", "Mmp13", "Mmp3", "Nlrp3", "Casp3", "Serpine1"
  ,"C3", "C6","Thbs1", "Myo1f", "Plek", "Rhoc", "Acta2", "Myl6",
  
  
  "Enpp2"  ,   # 外核苷酸焦磷酸酶/磷酸二酯酶2
 "Ccl27a",   # 趋化因子
  "Ccr3",     # 趋化因子受体
  "Cxcr5",    # 趋化因子受体
  "Il11ra1",  # IL-11受体
 "Il34",     # 白细胞介素34
  
  "Csf1r",    # 集落刺激因子1受体
  "Clec4f",   # C型凝集素受体
  "Clec9a",   # C型凝集素受体
 
  "Adora3",   # 腺苷A3受体
  "Tnfaip2",  # TNFα诱导蛋白

  "Tnfaip8l1"# TNFα诱导蛋白样1
 
)

cat("指定顺序的基因数量:", length(selected_genes_ordered), "\n")

# 4. 筛选数据集中存在的基因
cat("数据集中Gene列的前10个值:\n")
print(head(data$Gene))

# 尝试不区分大小写匹配
data_genes <- data$Gene
selected_genes_in_data <- c()

# 检查每个选定的基因是否在数据集中
for(gene in selected_genes_ordered) {
  # 精确匹配
  if(gene %in% data_genes) {
    selected_genes_in_data <- c(selected_genes_in_data, gene)
  } else {
    # 尝试大写匹配
    gene_upper <- toupper(gene)
    if(gene_upper %in% toupper(data_genes)) {
      matched_gene <- data_genes[toupper(data_genes) == gene_upper][1]
      selected_genes_in_data <- c(selected_genes_in_data, matched_gene)
      cat("注意: 基因", gene, "匹配为", matched_gene, "\n")
    } else {
      cat("警告: 基因", gene, "在数据集中未找到\n")
    }
  }
}

cat("在数据集中找到的基因数量:", length(selected_genes_in_data), "\n")
cat("找到的基因:", paste(selected_genes_in_data, collapse = ", "), "\n")

# 5. 提取选定基因的数据
selected_data <- data[data$Gene %in% selected_genes_in_data, ]

# 如果找不到任何基因，尝试使用数据集中前40个基因
if(nrow(selected_data) == 0) {
  cat("警告: 未找到任何选定基因，使用数据集中前40个基因\n")
  selected_data <- data[1:min(40, nrow(data)), ]
} else if(nrow(selected_data) < 20) {
  cat("警告: 只找到少量选定基因，将使用所有找到的基因\n")
}

cat("筛选后的数据维度:", dim(selected_data), "\n")

# 6. 检查列顺序并重新排序
# 定义期望的列顺序
desired_columns <- c("Gene", "C57_Control", "8BC_Control", "STAT1_Control", 
                     "C57_Infected", "8BC_Infected", "STAT1_Infected")

# 检查并重新排列列
if(all(desired_columns %in% colnames(selected_data))) {
  selected_data <- selected_data[, desired_columns]
  cat("列顺序已调整\n")
} else {
  # 如果列名不完全匹配，尝试查找相似的列
  cat("期望的列名不完全匹配，当前列名:\n")
  print(colnames(selected_data))
  
  # 尝试找到匹配的列
  actual_columns <- colnames(selected_data)
  c57_control_col <- actual_columns[grep("C57.*Control|control.*C57", actual_columns, ignore.case = TRUE)]
  c57_infected_col <- actual_columns[grep("C57.*Infected|infected.*C57", actual_columns, ignore.case = TRUE)]
  bc_control_col <- actual_columns[grep("8BC.*Control|control.*8BC|BC.*Control", actual_columns, ignore.case = TRUE)]
  bc_infected_col <- actual_columns[grep("8BC.*Infected|infected.*8BC|BC.*Infected", actual_columns, ignore.case = TRUE)]
  stat1_control_col <- actual_columns[grep("STAT1.*Control|control.*STAT1", actual_columns, ignore.case = TRUE)]
  stat1_infected_col <- actual_columns[grep("STAT1.*Infected|infected.*STAT1", actual_columns, ignore.case = TRUE)]
  
  # 重新组合数据
  selected_data <- selected_data[, c("Gene", c57_control_col, bc_control_col, stat1_control_col,
                                     c57_infected_col, bc_infected_col, stat1_infected_col)]
  colnames(selected_data) <- desired_columns
  cat("列已重命名并重新排序\n")
}

# 7. 准备热图数据
# 设置行名
rownames(selected_data) <- selected_data$Gene

# 提取表达矩阵
expr_matrix <- as.matrix(selected_data[, -1])  # 移除Gene列
cat("表达矩阵维度:", dim(expr_matrix), "\n")

# 检查数据
cat("表达矩阵前3行:\n")
print(expr_matrix[1:min(3, nrow(expr_matrix)), ])

# 8. 对行进行标准化（Z-score）
cat("对行进行标准化...\n")
expr_matrix_scaled <- t(scale(t(expr_matrix)))
# 处理NA/无限值
expr_matrix_scaled[is.na(expr_matrix_scaled)] <- 0
expr_matrix_scaled[!is.finite(expr_matrix_scaled)] <- 0

# 9. 创建列注释
col_annotation <- data.frame(
  Strain = factor(rep(c("C57", "8BC", "STAT1"), 2)),
  Condition = factor(c(rep("Control", 3), rep("Infected", 3)), 
                     levels = c("Control", "Infected"))
)
rownames(col_annotation) <- colnames(expr_matrix)

# 10. 按照指定顺序排列行
# 只保留在数据集中找到的基因
ordered_genes_found <- selected_genes_ordered[selected_genes_ordered %in% rownames(expr_matrix_scaled)]

# 如果找到的基因与指定顺序不完全匹配，尝试大小写不敏感的匹配
if(length(ordered_genes_found) < length(selected_genes_in_data)) {
  # 尝试大小写不敏感的匹配
  temp_ordered_genes <- c()
  for(gene in selected_genes_ordered) {
    gene_lower <- tolower(gene)
    data_genes_lower <- tolower(rownames(expr_matrix_scaled))
    if(gene_lower %in% data_genes_lower) {
      idx <- which(data_genes_lower == gene_lower)[1]
      temp_ordered_genes <- c(temp_ordered_genes, rownames(expr_matrix_scaled)[idx])
    }
  }
  ordered_genes_found <- temp_ordered_genes
}

cat("按照指定顺序找到的基因数量:", length(ordered_genes_found), "\n")
cat("按指定顺序排列的基因:", paste(head(ordered_genes_found, 10), collapse = ", "), "...\n")

# 重新排列表达矩阵，按照指定顺序
expr_matrix_scaled_ordered <- expr_matrix_scaled[ordered_genes_found, , drop = FALSE]

# 11. 设置颜色
# 列注释颜色
strain_colors <- c("C57" = "#1F77B4", "8BC" = "#FF7F0E", "STAT1" = "#2CA02C")
condition_colors <- c("Control" = "#708090", "Infected" = "#A52A2A")

annotation_colors <- list(
  Strain = strain_colors,
  Condition = condition_colors
)

# 12. 绘制热图
cat("开始绘制热图...\n")

# 设置颜色映射
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# 绘制热图，不聚类行，按照指定顺序排列
heatmap_plot <- pheatmap(
  expr_matrix_scaled_ordered,
  color = color_palette,
  scale = "row",  # 对行进行标准化
  annotation_col = col_annotation,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # 不对行进行聚类，保持指定顺序
  cluster_cols = FALSE,  # 不对列进行聚类，保持指定顺序
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  border_color = NA,
  main = "",
  gaps_col = 3,  # 在Control和Infected之间添加分割线
  display_numbers = FALSE,  # 不显示数值
  number_color = "black",
  number_format = "%.1f",
  cellwidth = 20,
  cellheight = 10
)
# 
# # 13. 保存热图
# output_pdf <- "猴痘病毒感染_差异基因热图_按指定顺序.pdf"
# pdf(output_pdf, width = 10, height = 12)
# print(heatmap_plot)
# dev.off()
# cat("热图已保存为:", output_pdf, "\n")
# 
# # 也保存为PNG格式
# output_png <- "猴痘病毒感染_差异基因热图_按指定顺序.png"
# png(output_png, width = 2000, height = 2400, res = 300)
# print(heatmap_plot)
# dev.off()
# cat("热图已保存为:", output_png, "\n")
# 
# # 14. 保存热图数据
# # 保存标准化后的数据
# scaled_data <- as.data.frame(expr_matrix_scaled_ordered)
# scaled_data$Gene <- rownames(expr_matrix_scaled_ordered)
# scaled_data <- scaled_data[, c("Gene", "C57_Control", "8BC_Control", "STAT1_Control", 
#      
# 
# cat("热图数据已保存为CSV文件\n")
