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

# 3. 定义选定的基因列表
selected_genes <- c(
  # 炎症和抗病毒反应
  "Serpinb10",  # 丝氨酸蛋白酶抑制剂，参与炎症反应
  "Irf8",       # 干扰素调节因子8，干扰素信号通路关键因子
  "Mmp8",       # 基质金属蛋白酶8，炎症和组织重塑
  "S100a3",     # S100钙结合蛋白，炎症反应
  "Trem1",      # 髓细胞触发受体1，炎症反应
  "Trem2",      # 髓细胞触发受体2，炎症反应
  "Ccl19",      # 趋化因子，免疫细胞招募
  "Zfp280b",    # 锌指蛋白，可能参与基因调控
  "Angptl1",    # 血管生成素样蛋白1，可能与炎症相关
  
  # 干扰素和抗病毒通路
  "Irgc1",      # 免疫相关GTP酶，抗病毒免疫
  "Kpna1",      # 核导入蛋白，可能参与干扰素信号转导
  
  # 细胞应激和凋亡
  "Btg2",       # 抗增殖蛋白，细胞应激反应
  "Mtch1",      # 线粒体载体蛋白，细胞凋亡
  "Dido1",      # 细胞凋亡调节因子
  "Cep55",      # 中心体蛋白，细胞周期
  
  # 免疫细胞功能
  "Btnl4",      # 嗜乳脂蛋白样4，免疫调节
  "Timd2",      # T细胞免疫球蛋白和粘蛋白结构域2
  "Ms4a2",      # 肥大细胞表面抗原，过敏反应
  
  # 信号转导
  "Shc2",       # SHC转换蛋白2，信号转导
  "Prkce",      # 蛋白激酶C epsilon，信号转导
  "Srpk1",      # 丝氨酸/精氨酸蛋白激酶1
  "Frs3",       # 成纤维细胞生长因子受体底物3
  "Eif4g1",     # 真核翻译起始因子4G1
  "Wnt7b",      # Wnt信号通路
  
  # 代谢和转运
  "Cyp1a2",     # 细胞色素P450，药物代谢
  "Atp8b4",     # ATP酶，脂质转运
  "Slc2a12",    # 葡萄糖转运蛋白
  "Gapdh",      # 糖酵解关键酶，也参与非代谢功能
  
  # 转录调控
  "Hira",       # 组蛋白调控因子
  "Ing5",       # 抑制因子5，转录调控
  "Katnb1",     # 微管切割蛋白，转录调控
  "Arhgap19",   # Rho GTP酶激活蛋白
  "Cmip",       # c-Maf诱导蛋白
  
  # 细胞结构
  "Dnah1",      # 动力蛋白，纤毛运动
  "Pkp4",       # 桥粒蛋白，细胞连接
  "Epn3",       # 网格蛋白连接蛋白，内存作用
  
  # 其他功能基因
  "Trnp1",      # 调节核仁和增殖
  "Ep400",      # 组蛋白乙酰转移酶复合体
  "Cnr1",       # 大麻素受体1
  "Gabrr2",     # GABA受体
  "Olfr59",     # 嗅觉受体
  "Col6a6",     # 胶原蛋白VI
  "Sgca",       # 肌聚糖α
  "Tns4"        # 张力蛋白4，细胞粘附
)

cat("选定的基因数量:", length(selected_genes), "\n")

# 4. 筛选数据集中存在的基因
# 首先检查数据中的Gene列名称
cat("数据集中Gene列的前10个值:\n")
print(head(data$Gene))

# 尝试不区分大小写匹配
data_genes <- data$Gene
selected_genes_in_data <- c()

# 检查每个选定的基因是否在数据集中
for(gene in selected_genes) {
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

# 10. 创建行注释（基因功能分类）
# 定义基因功能分类
gene_categories <- data.frame(
  Gene = selected_genes_in_data[selected_genes_in_data %in% selected_data$Gene],
  Category = NA
)

# 为每个基因分配功能类别
assign_category <- function(gene_name) {
  gene_lower <- tolower(gene_name)
  
  # 根据基因名关键词分配类别
  if(grepl("serpin|mmp|s100|trem|ccl|angptl", gene_lower)) {
    return("炎症/抗病毒")
  } else if(grepl("irf|irgc|kpna", gene_lower)) {
    return("干扰素通路")
  } else if(grepl("btg|mtch|dido|cep", gene_lower)) {
    return("细胞应激/凋亡")
  } else if(grepl("btnl|timd|ms4a", gene_lower)) {
    return("免疫细胞功能")
  } else if(grepl("shc|prk|srpk|frs|eif|wnt", gene_lower)) {
    return("信号转导")
  } else if(grepl("cyp|atp|slc|gapdh", gene_lower)) {
    return("代谢/转运")
  } else if(grepl("hir|ing|katn|arhg|cmip", gene_lower)) {
    return("转录调控")
  } else if(grepl("dnah|pkp|epn", gene_lower)) {
    return("细胞结构")
  } else {
    return("其他功能")
  }
}

# 为每个基因分配类别
for(i in 1:nrow(selected_data)) {
  gene_name <- rownames(selected_data)[i]
  gene_categories$Category[gene_categories$Gene == gene_name] <- assign_category(gene_name)
}

# 创建行注释数据框
row_annotation <- data.frame(
  Category = factor(gene_categories$Category, 
                    levels = c("炎症/抗病毒", "干扰素通路", "细胞应激/凋亡", 
                               "免疫细胞功能", "信号转导", "代谢/转运", 
                               "转录调控", "细胞结构", "其他功能"))
)
rownames(row_annotation) <- gene_categories$Gene

# 11. 设置颜色
# 列注释颜色
strain_colors <- c("C57" = "#1F77B4", "8BC" = "#FF7F0E", "STAT1" = "#2CA02C")
condition_colors <- c("Control" = "#708090", "Infected" = "#A52A2A")

# 行注释颜色（使用Set3调色板）
category_colors <- colorRampPalette(brewer.pal(9, "Set3"))(length(unique(row_annotation$Category)))
names(category_colors) <- unique(row_annotation$Category)

annotation_colors <- list(
  Strain = strain_colors,
  Condition = condition_colors,
  Category = category_colors
)

# 12. 绘制热图
cat("开始绘制热图...\n")

# 设置颜色映射
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# 绘制热图
heatmap_plot <- pheatmap(
  expr_matrix_scaled,
  color = color_palette,
  scale = "row",  # 对行进行标准化
  annotation_col = col_annotation,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,  # 对行进行聚类
  cluster_cols = FALSE,  # 不对列进行聚类，保持指定顺序
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  border_color = NA,
  main = "猴痘病毒感染小鼠差异表达基因热图\n(C57, 8BC, STAT1品系, 感染后7天 vs Mock)",
  gaps_col = 3,  # 在Control和Infected之间添加分割线
  display_numbers = FALSE,  # 不显示数值
  number_color = "black",
  number_format = "%.1f",
  cellwidth = 20,
  cellheight = 10,
  treeheight_row = 50,
  treeheight_col = 20
)

# 13. 保存热图
output_pdf <- "猴痘病毒感染_差异基因热图.pdf"
pdf(output_pdf, width = 10, height = 12)
print(heatmap_plot)
dev.off()
cat("热图已保存为:", output_pdf, "\n")

# 也保存为PNG格式
output_png <- "猴痘病毒感染_差异基因热图.png"
png(output_png, width = 2000, height = 2400, res = 300)
print(heatmap_plot)
dev.off()
cat("热图已保存为:", output_png, "\n")

# 14. 保存热图数据
# 保存标准化后的数据
scaled_data <- as.data.frame(expr_matrix_scaled)
scaled_data$Gene <- rownames(scaled_data)
scaled_data <- scaled_data[, c("Gene", "C57_Control", "8BC_Control", "STAT1_Control", 
                               "C57_Infected", "8BC_Infected", "STAT1_Infected")]

write.csv(scaled_data, "热图数据_标准化表达量.csv", row.names = FALSE)
write.csv(selected_data, "热图数据_原始表达量.csv", row.names = FALSE)

# 保存注释信息
write.csv(col_annotation, "列注释信息.csv")
write.csv(row_annotation, "行注释信息.csv")

cat("热图数据已保存为CSV文件\n")

# 15. 创建简化的热图（不聚类行）
cat("创建简化版热图（不聚类行）...\n")

# 按照基因功能类别排序
if(!is.null(row_annotation$Category)) {
  # 按照类别排序
  gene_order <- order(row_annotation$Category)
  expr_matrix_sorted <- expr_matrix_scaled[gene_order, ]
  row_annotation_sorted <- row_annotation[gene_order, , drop = FALSE]
  
  # 计算行分割位置
  category_counts <- table(row_annotation_sorted$Category)
  row_gaps <- cumsum(category_counts)[-length(category_counts)]
} else {
  expr_matrix_sorted <- expr_matrix_scaled
  row_gaps <- NULL
}

# 绘制不聚类行的热图
heatmap_simple <- pheatmap(
  expr_matrix_sorted,
  color = color_palette,
  scale = "row",
  annotation_col = col_annotation,
  annotation_row = row_annotation_sorted,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  border_color = NA,
  main = "猴痘病毒感染小鼠差异表达基因热图\n(按功能分类排序)",
  gaps_col = 3,  # 在Control和Infected之间添加分割线
  gaps_row = row_gaps,  # 在基因类别之间添加分割线
  display_numbers = FALSE,
  cellwidth = 20,
  cellheight = 10
)

# 保存简化版热图
pdf("猴痘病毒感染_差异基因热图_按功能排序.pdf", width = 10, height = 12)
print(heatmap_simple)
dev.off()
cat("简化版热图已保存为: 猴痘病毒感染_差异基因热图_按功能排序.pdf\n")

# 16. 创建各品系对比的热图
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
  
  # 绘制FoldChange热图
  fc_heatmap <- pheatmap(
    fc_matrix,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    scale = "none",  # 不进行标准化
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = NA,
    main = "各品系感染/对照log2FoldChange热图",
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    cellwidth = 40,
    cellheight = 10
  )
  
  # 保存FoldChange热图
  pdf("各品系_log2FoldChange_热图.pdf", width = 8, height = 10)
  print(fc_heatmap)
  dev.off()
  cat("FoldChange热图已保存为: 各品系_log2FoldChange_热图.pdf\n")
  
  # 保存FoldChange数据
  fc_data <- as.data.frame(fc_matrix)
  fc_data$Gene <- rownames(fc_matrix)
  write.csv(fc_data, "各品系_log2FoldChange_数据.csv", row.names = FALSE)
}

# 17. 创建汇总统计
cat("\n=== 汇总统计 ===\n")
cat("使用的基因数量:", nrow(selected_data), "\n")
cat("样本信息:\n")
cat("  Control组: C57, 8BC, STAT1\n")
cat("  Infected组: C57, 8BC, STAT1\n")

# 计算各组的平均表达量
if(nrow(selected_data) > 0) {
  cat("\n各组的平均表达量:\n")
  cat("  C57 Control平均:", mean(selected_data$C57_Control, na.rm = TRUE), "\n")
  cat("  8BC Control平均:", mean(selected_data$`8BC_Control`, na.rm = TRUE), "\n")
  cat("  STAT1 Control平均:", mean(selected_data$STAT1_Control, na.rm = TRUE), "\n")
  cat("  C57 Infected平均:", mean(selected_data$C57_Infected, na.rm = TRUE), "\n")
  cat("  8BC Infected平均:", mean(selected_data$`8BC_Infected`, na.rm = TRUE), "\n")
  cat("  STAT1 Infected平均:", mean(selected_data$STAT1_Infected, na.rm = TRUE), "\n")
}

cat("\n分析完成！\n")
cat("输出文件:\n")
cat("1. 猴痘病毒感染_差异基因热图.pdf - 主热图\n")
cat("2. 猴痘病毒感染_差异基因热图.png - 主热图(PNG格式)\n")
cat("3. 猴痘病毒感染_差异基因热图_按功能排序.pdf - 按功能分类排序的热图\n")
cat("4. 各品系_log2FoldChange_热图.pdf - 各品系FoldChange热图\n")
cat("5. 热图数据_标准化表达量.csv - 标准化后的表达矩阵\n")
cat("6. 热图数据_原始表达量.csv - 原始表达矩阵\n")
cat("7. 各品系_log2FoldChange_数据.csv - FoldChange数据\n")
cat("8. 列注释信息.csv - 样本注释信息\n")
cat("9. 行注释信息.csv - 基因注释信息\n")