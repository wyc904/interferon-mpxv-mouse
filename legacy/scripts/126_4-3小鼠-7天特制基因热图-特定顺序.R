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
  # 炎症和抗病毒反应

    "Irf8",     # 干扰素调节因子8，抗病毒免疫关键转录因子
    "Jak2",     # JAK-STAT信号通路关键激酶
    "Cd3e",     # T细胞受体复合物，T细胞激活标志
    "Ptprc",    # CD45，白细胞共同抗原
    "B2m",      # β2微球蛋白，MHC I类分子组分
    "Cd44",     # 细胞表面粘附分子，免疫细胞迁移
    "Tnfsf13b", # B细胞激活因子
    "Il13ra1",  # IL-13受体，调节Th2免疫应答
    "Ifnar2",   # 干扰素α/β受体，I型干扰素信号
    "Oasl2",    # 寡腺苷酸合成酶，抗病毒效应分子
    "Oasl1a",   # 寡腺苷酸合成酶家族成员
    "Isg20",    # 干扰素刺激基因，核酸外切酶活性
    "Ifit1",    # 干扰素诱导蛋白，抑制病毒复制
    "Ifit2",    # 干扰素诱导蛋白，限制病毒传播
    "Bst2",     # 抗病毒限制因子（tetherin）
    "Mx1",      # 干扰素诱导的GTP酶，抗病毒蛋白
    "Tnf",      # 肿瘤坏死因子，促炎细胞因子
    "Trem1",    # 髓样细胞触发受体-1，炎症放大
    "Trem2",    # 髓样细胞触发受体-2，炎症调节
    "Serpinb10",# 丝氨酸蛋白酶抑制剂，调节炎症
    "Serpinb9b",# 丝氨酸蛋白酶抑制剂，调节细胞凋亡
    "Serpinb9c",# 丝氨酸蛋白酶抑制剂家族成员
    "Mmp8",     # 基质金属蛋白酶8，组织重塑
    "Mmp13",    # 基质金属蛋白酶13，细胞外基质降解
    "Mmp3",     # 基质金属蛋白酶3，炎症反应
    "Mmp14",    # 基质金属蛋白酶14，细胞迁移
    "Tlr9",     # Toll样受体9，识别病毒DNA
    "Ddx60",    # 病毒RNA感受器
    "Dhx58",    # 病毒RNA感受器
    "Nlrp3",    # 炎症小体组分，IL-1β激活
    "Aim2",     # AIM2炎症小体，识别病毒DNA
    "Casp4",    # 炎症性caspase
    "Casp3",    # 凋亡执行caspase
    "Itgax",    # CD11c，树突细胞标志
    "Cd300lg",  # 免疫调节受体
    "Clec7a",   # C型凝集素受体，模式识别受体
    "C3",       # 补体C3，补体系统激活
    "C6",       # 补体C6，膜攻击复合体组分
    "Serpine1", # 纤溶酶原激活物抑制剂
    "Bcl2",     # 抗凋亡蛋白，调节细胞生存
    "Fas"  ,"Ddx58",
    "Sting1",
    "Ifi204",
    "Ifi206",
    "Ifi208",
    "Ifi214",
    "Usp18",
    "Stat2",
    "Tyk2",
    "Tlr12",
    "Unc93b1",
    "Samd9l",
    "Ccl19",
    "Ccl21a",
    "Ccl27a",
    "Ccr3",
    "Cxcr5",
    "Il7r",
    "Il11ra1",
    "Il17ra",
    "Il17re",
    "Il18",
    "Il34",
    "Tap1",
    "Tap2",
    "Psmb8",
    "H2-Eb1",
    "H2-K2",
    "Csf1r",
    "Clec4f",
    "Clec9a",
    "Ms4a4b",
    "Ms4a4c",
    "Adora2a",
    "Adora3",
    "Tnfaip2",
    "Tnfaip8",
    "Tnfaip8l1",
    "Ncf1",
    "Ctsw"
   
)

cat("指定顺序的基因数量:", length(selected_genes_ordered), "\n")



# 创建功能类别数据框
gene_info <- data.frame(
  Gene = selected_genes_ordered,
  Category = gene_categories_ordered,
  stringsAsFactors = FALSE
)

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

# 10. 创建行注释（按指定顺序和功能分类）
# 创建一个数据框存储基因信息
gene_info_in_data <- data.frame(
  Gene = selected_genes_in_data,
  Category = NA
)

# 根据指定顺序分配功能类别
# 创建一个从基因到类别的映射
gene_to_category <- setNames(gene_info$Category, gene_info$Gene)

# 为每个在数据中存在的基因分配功能类别
for(i in 1:length(selected_genes_in_data)) {
  gene <- selected_genes_in_data[i]
  if(gene %in% names(gene_to_category)) {
    gene_info_in_data$Category[i] <- gene_to_category[gene]
  } else {
    # 如果基因名不匹配，尝试不区分大小写
    gene_lower <- tolower(gene)
    gene_names_lower <- tolower(names(gene_to_category))
    if(gene_lower %in% gene_names_lower) {
      idx <- which(gene_names_lower == gene_lower)[1]
      gene_info_in_data$Category[i] <- gene_to_category[idx]
    } else {
      # 如果仍然找不到，分配默认类别
      gene_info_in_data$Category[i] <- "其他功能基因"
    }
  }
}

# 创建行注释数据框
row_annotation <- data.frame(
  Category = factor(gene_info_in_data$Category, 
                    levels = unique(gene_info_in_data$Category))
)
rownames(row_annotation) <- gene_info_in_data$Gene

# 11. 设置颜色
# 列注释颜色
strain_colors <- c("C57" = "#1F77B4", "8BC" = "#FF7F0E", "STAT1" = "#2CA02C")
condition_colors <- c("Control" = "#708090", "Infected" = "#A52A2A")

# 行注释颜色（使用Set3调色板）
unique_categories <- unique(gene_info_in_data$Category)
category_colors <- colorRampPalette(brewer.pal(min(9, length(unique_categories)), "Set3"))(length(unique_categories))
names(category_colors) <- unique_categories

annotation_colors <- list(
  Strain = strain_colors,
  Condition = condition_colors,
  Category = category_colors
)

# 12. 按照指定顺序排列行
# 只保留在数据集中找到的基因
ordered_genes_found <- selected_genes_ordered[selected_genes_ordered %in% rownames(expr_matrix_scaled)]

# 如果找到的基因与指定顺序不完全匹配，尝试大小写不敏感的匹配
if(length(ordered_genes_found) < length(selected_genes_in_data)) {
  # 尝试大小写不敏感的匹配
  ordered_genes_found <- c()
  for(gene in selected_genes_ordered) {
    gene_lower <- tolower(gene)
    data_genes_lower <- tolower(rownames(expr_matrix_scaled))
    if(gene_lower %in% data_genes_lower) {
      idx <- which(data_genes_lower == gene_lower)[1]
      ordered_genes_found <- c(ordered_genes_found, rownames(expr_matrix_scaled)[idx])
    }
  }
}

cat("按照指定顺序找到的基因数量:", length(ordered_genes_found), "\n")
cat("按指定顺序排列的基因:", paste(ordered_genes_found, collapse = ", "), "\n")

# 重新排列表达矩阵，按照指定顺序
expr_matrix_scaled_ordered <- expr_matrix_scaled[ordered_genes_found, , drop = FALSE]

# 相应地重新排列行注释
row_annotation_ordered <- row_annotation[ordered_genes_found, , drop = FALSE]

# 13. 计算行分割位置
# 根据功能类别计算分割位置
if(nrow(row_annotation_ordered) > 0) {
  # 获取每个类别的基因数
  category_counts <- table(row_annotation_ordered$Category)
  # 计算累积和，排除最后一个
  row_gaps <- cumsum(category_counts)[-length(category_counts)]
} else {
  row_gaps <- NULL
}

# 14. 绘制热图
cat("开始绘制热图...\n")

# 设置颜色映射
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# 绘制热图，不聚类行，按照指定顺序排列
heatmap_plot <- pheatmap(
  expr_matrix_scaled_ordered,
  color = color_palette,
  scale = "row",  # 对行进行标准化
  annotation_col = col_annotation,
  annotation_row = row_annotation_ordered,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # 不对行进行聚类，保持指定顺序
  cluster_cols = FALSE,  # 不对列进行聚类，保持指定顺序
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  border_color = NA,
  main = "猴痘病毒感染小鼠差异表达基因热图\n(C57, 8BC, STAT1品系, 感染后7天 vs Mock)",
  gaps_col = 3,  # 在Control和Infected之间添加分割线
  gaps_row = row_gaps,  # 在功能类别之间添加分割线
  display_numbers = FALSE,  # 不显示数值
  number_color = "black",
  number_format = "%.1f",
  cellwidth = 20,
  cellheight = 10
)

# 15. 保存热图
output_pdf <- "猴痘病毒感染_差异基因热图_按指定顺序.pdf"
pdf(output_pdf, width = 10, height = 12)
print(heatmap_plot)
dev.off()
cat("热图已保存为:", output_pdf, "\n")

# 也保存为PNG格式
output_png <- "猴痘病毒感染_差异基因热图_按指定顺序.png"
png(output_png, width = 2000, height = 2400, res = 300)
print(heatmap_plot)
dev.off()
cat("热图已保存为:", output_png, "\n")

# 16. 保存热图数据
# 保存标准化后的数据
scaled_data <- as.data.frame(expr_matrix_scaled_ordered)
scaled_data$Gene <- rownames(expr_matrix_scaled_ordered)
scaled_data <- scaled_data[, c("Gene", "C57_Control", "8BC_Control", "STAT1_Control", 
                               "C57_Infected", "8BC_Infected", "STAT1_Infected")]

write.csv(scaled_data, "热图数据_标准化表达量_按指定顺序.csv", row.names = FALSE)
write.csv(selected_data, "热图数据_原始表达量.csv", row.names = FALSE)

# 保存注释信息
write.csv(col_annotation, "列注释信息.csv")
write.csv(row_annotation_ordered, "行注释信息_按指定顺序.csv")

cat("热图数据已保存为CSV文件\n")

# 17. 创建各品系对比的热图
cat("创建各品系对比热图...\n")

# 计算感染/对照的比值
if(nrow(selected_data) > 0) {
  # 计算log2FoldChange
  c57_fc <- selected_data$C57_Infected - selected_data$C57_Control
  bc_fc <- selected_data$`8BC_Infected` - selected_data$`8BC_Control`
  stat1_fc <- selected_data$STAT1_Infected - selected_data$STAT1_Control
  
  # 创建FoldChange矩阵
  fc_matrix <- cbind(c57_fc, bc_fc, stat1_fc)
  rownames(fc_matrix) <- selected_data$Gene
  colnames(fc_matrix) <- c("C57_FC", "8BC_FC", "STAT1_FC")
  
  # 按照指定顺序排列FoldChange矩阵
  fc_matrix_ordered <- fc_matrix[ordered_genes_found, , drop = FALSE]
  
  # 绘制FoldChange热图
  fc_heatmap <- pheatmap(
    fc_matrix_ordered,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    scale = "none",  # 不进行标准化
    cluster_rows = FALSE,  # 不聚类行
    cluster_cols = FALSE,  # 不聚类列
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = NA,
    main = "各品系感染/对照log2FoldChange热图\n(按指定基因顺序)",
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    cellwidth = 40,
    cellheight = 10
  )
  
  # 保存FoldChange热图
  pdf("各品系_log2FoldChange_热图_按指定顺序.pdf", width = 8, height = 10)
  print(fc_heatmap)
  dev.off()
  cat("FoldChange热图已保存为: 各品系_log2FoldChange_热图_按指定顺序.pdf\n")
  
  # 保存FoldChange数据
  fc_data <- as.data.frame(fc_matrix_ordered)
  fc_data$Gene <- rownames(fc_matrix_ordered)
  write.csv(fc_data, "各品系_log2FoldChange_数据_按指定顺序.csv", row.names = FALSE)
}

# 18. 创建汇总统计
cat("\n=== 汇总统计 ===\n")
cat("使用的基因数量:", nrow(expr_matrix_scaled_ordered), "\n")
cat("样本信息:\n")
cat("  Control组: C57, 8BC, STAT1\n")
cat("  Infected组: C57, 8BC, STAT1\n")

# 计算各组的平均表达量
if(nrow(selected_data) > 0) {
  # 计算log2FoldChange
  selected_data$C57_log2FC <- selected_data$C57_Infected - selected_data$C57_Control
  selected_data$`8BC_log2FC` <- selected_data$`8BC_Infected` - selected_data$`8BC_Control`
  selected_data$STAT1_log2FC <- selected_data$STAT1_Infected - selected_data$STAT1_Control
  
  # 计算每个品系中差异表达最显著的基因
  cat("\n各品系表达变化最大的基因:\n")
  
  # C57
  c57_sig_up <- selected_data[order(selected_data$C57_log2FC, decreasing = TRUE)[1:5], c("Gene", "C57_log2FC")]
  c57_sig_down <- selected_data[order(selected_data$C57_log2FC, decreasing = FALSE)[1:5], c("Gene", "C57_log2FC")]
  cat("\nC57表达上调前5的基因:\n")
  print(c57_sig_up)
  cat("\nC57表达下调前5的基因:\n")
  print(c57_sig_down)
  
  # 8BC
  bc_sig_up <- selected_data[order(selected_data$`8BC_log2FC`, decreasing = TRUE)[1:5], c("Gene", "8BC_log2FC")]
  bc_sig_down <- selected_data[order(selected_data$`8BC_log2FC`, decreasing = FALSE)[1:5], c("Gene", "8BC_log2FC")]
  cat("\n8BC表达上调前5的基因:\n")
  print(bc_sig_up)
  cat("\n8BC表达下调前5的基因:\n")
  print(bc_sig_down)
  
  # STAT1
  stat1_sig_up <- selected_data[order(selected_data$STAT1_log2FC, decreasing = TRUE)[1:5], c("Gene", "STAT1_log2FC")]
  stat1_sig_down <- selected_data[order(selected_data$STAT1_log2FC, decreasing = FALSE)[1:5], c("Gene", "STAT1_log2FC")]
  cat("\nSTAT1表达上调前5的基因:\n")
  print(stat1_sig_up)
  cat("\nSTAT1表达下调前5的基因:\n")
  print(stat1_sig_down)
  
  # 保存差异表达数据
  write.csv(selected_data[, c("Gene", "C57_Control", "C57_Infected", "C57_log2FC", 
                              "8BC_Control", "8BC_Infected", "8BC_log2FC",
                              "STAT1_Control", "STAT1_Infected", "STAT1_log2FC")], 
            "差异表达数据_按指定顺序.csv", row.names = FALSE)
}

cat("\n分析完成！\n")
cat("输出文件:\n")
cat("1. 猴痘病毒感染_差异基因热图_按指定顺序.pdf - 主热图（按指定基因顺序）\n")
cat("2. 猴痘病毒感染_差异基因热图_按指定顺序.png - 主热图(PNG格式)\n")
cat("3. 各品系_log2FoldChange_热图_按指定顺序.pdf - 各品系FoldChange热图\n")
cat("4. 热图数据_标准化表达量_按指定顺序.csv - 标准化后的表达矩阵\n")
cat("5. 热图数据_原始表达量.csv - 原始表达矩阵\n")
cat("6. 各品系_log2FoldChange_数据_按指定顺序.csv - FoldChange数据\n")
cat("7. 差异表达数据_按指定顺序.csv - 包含所有差异表达信息\n")
cat("8. 列注释信息.csv - 样本注释信息\n")
cat("9. 行注释信息_按指定顺序.csv - 基因注释信息（按指定顺序）\n")

# 19. 创建额外的热图，按照品系和感染状态分组
cat("\n创建额外的热图，按品系分组...\n")

# 重新排列列，按品系分组
expr_matrix_scaled_ordered_groups <- expr_matrix_scaled_ordered[, c("C57_Control", "C57_Infected", "8BC_Control", "8BC_Infected", "STAT1_Control", "STAT1_Infected")]

# 创建新的列注释
col_annotation_groups <- data.frame(
  Strain = factor(c("C57", "C57", "8BC", "8BC", "STAT1", "STAT1"),
                  levels = c("C57", "8BC", "STAT1")),
  Condition = factor(c("Control", "Infected", "Control", "Infected", "Control", "Infected"),
                     levels = c("Control", "Infected"))
)
rownames(col_annotation_groups) <- colnames(expr_matrix_scaled_ordered_groups)

# 绘制分组热图
heatmap_groups <- pheatmap(
  expr_matrix_scaled_ordered_groups,
  color = color_palette,
  scale = "row",  # 对行进行标准化
  annotation_col = col_annotation_groups,
  annotation_row = row_annotation_ordered,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # 不对行进行聚类
  cluster_cols = FALSE,  # 不对列进行聚类
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  border_color = NA,
  main = "猴痘病毒感染小鼠差异表达基因热图\n(按品系分组)",
  gaps_col = c(2, 4),  # 在C57和8BC、8BC和STAT1之间添加分割线
  gaps_row = row_gaps,  # 在功能类别之间添加分割线
  display_numbers = FALSE,
  number_color = "black",
  number_format = "%.1f",
  cellwidth = 20,
  cellheight = 10
)

# 保存分组热图
pdf("猴痘病毒感染_差异基因热图_按品系分组.pdf", width = 10, height = 12)
print(heatmap_groups)
dev.off()
cat("分组热图已保存为: 猴痘病毒感染_差异基因热图_按品系分组.pdf\n")

cat("\n所有热图绘制完成！\n")