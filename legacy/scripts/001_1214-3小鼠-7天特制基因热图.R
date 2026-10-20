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
data <- read.xlsx("热图绘制.xlsx", sheet = 1)
cat("数据维度:", dim(data), "\n")

# 检查列名
print(colnames(data))
# 输出: [1] "Gene"        "MOCK1"       "MOCK2"       "MOCK3"       "MPXV.IIB-1"  "MPXV.IIB-2"  "MPXV.IIB-33" "MPXV.IB-1"   "MPXV.IB-2"   "MPXV.IB-3"  

# 3. 定义选定的基因列表
selected_genes <- c(
  # 炎症和抗病毒反应
  
  "Cxcr1","Nfkb1","Chit1","S100a3","Ly6i","Il1r2","Il1r1","Il1rl2","Tnfsf9","Mmp8","Mmp3","Mmp13","Mmp27","Adamts9","Itga8","Col1a2","Vim","Gbp6","Oas1f","Oas1b","Trim62","Trim7","Ifi214","Ifi204","Ifi208","Clec9a","Cd209a","Clec1a","Clec1b","Cd209b","Clec4f","H2-Eb1","H2-K2","Ccl19","Ccl21a","Pax5","Bach2","Cd72","Fcrl5","Slamf1","Iglc1"
  
  
)

# 
# 
# 
# "Trem2",      # 髓细胞触发受体2，炎症反应
# "Ccl19",      # 趋化因子，免疫细胞招募
# "Zfp280b",    # 锌指蛋白，可能参与基因调控
# "Angptl1",    # 血管生成素样蛋白1，可能与炎症相关
# "Dido1",      # 细胞凋亡调节因子
# "Btnl4",      # 嗜乳脂蛋白样4，免疫调节
# "Timd2",      # T细胞免疫球蛋白和粘蛋白结构域2
# "Prkce",      # 蛋白激酶C epsilon，信号转导
# "Frs3",       # 成纤维细胞生长因子受体底物3
# "Eif4g1",     # 真核翻译起始因子4G1
# "Wnt7b",      # Wnt信号通路
# 
# 
# 
# 
# 
# # 代谢和转运
# 
# "Slc2a12",    # 葡萄糖转运蛋白
# # 转录调控
# "Hira",       # 组蛋白调控因子
# "Ing5",       # 抑制因子5，转录调控
# "Katnb1",     # 微管切割蛋白，转录调
# "Cmip",       # c-Maf诱导蛋白
# # 细胞结构
# "Dnah1",      # 动力蛋白，纤毛运动
# "Pkp4",       # 桥粒蛋白，细胞连接
# "Epn3",       # 网格蛋白连接蛋白，内存作用
# 
# # 其他功能基因
# "Trnp1",      # 调节核子核仁和增殖
# "Ep400",      # 组蛋白乙酰转移酶复合体
# "Cnr1",       # 大麻素受体1
# "Gabrr2",     # GABA受体
# "Col6a6",     # 胶原蛋白VI
# "Sgca",       # 肌聚糖α
# "Tns4"     ,  
# 
# 
# "Cyp1a2",     # 细胞色素P450，药物代谢
# "Atp8b4",     # ATP酶，脂质转运
# "Olfr59",     # 嗅觉受体
# "Irf8",       # 干扰素调节因子8，干扰素信号通路关键因子
# "Mmp8",       # 基质金属蛋白酶8，炎症和组织重塑
# "S100a3",     # S100钙结合蛋白，炎症反应
# "Gapdh",      # 糖酵解关键酶，也参与非代谢功能
# "Srpk1",      # 丝氨酸/精氨酸蛋白激酶1
# "Irgc1",      # 免疫相关GTP酶，抗病毒免疫
# "Kpna1",      # 核导入蛋白，可能参与干扰素信号转导
# # 细胞应激和凋亡
# "Btg2",       # 抗增殖蛋白，细胞应激反应
# "Mtch1",      # 线粒体载体蛋白，细胞凋亡
# 
# "Cep55",      # 中心体蛋白，细胞周期
# 
# # 免疫细胞功能
# 
# "Ms4a2",      # 肥大细胞表面抗原，过敏反应
# 
# # 信号转导
# "Shc2"     # SHC转换蛋白2，信号转导
# 
# 







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

# 6. 按照selected_genes的顺序重新排列行
# 确保selected_genes_in_data中的基因按照selected_genes的顺序排列
ordered_genes <- c()
for(gene in selected_genes) {
  if(gene %in% selected_genes_in_data) {
    ordered_genes <- c(ordered_genes, gene)
  }
}

# 按照ordered_genes的顺序重新排列数据
selected_data <- selected_data[match(ordered_genes, selected_data$Gene), ]

cat("按照指定顺序重新排列后的基因顺序:\n")
print(selected_data$Gene)

# 7. 根据实际列名调整列顺序
# 原始列名: "Gene", "MOCK1", "MOCK2", "MOCK3", "MPXV.IIB-1", "MPXV.IIB-2", "MPXV.IIB-33", "MPXV.IB-1", "MPXV.IB-2", "MPXV.IB-3"

# 定义期望的列顺序（按照MOCK, MPXV.IIB, MPXV.IB分组）
desired_columns <- c("Gene", "MOCK1", "MOCK2", "MOCK3", 
                     "MPXV.IIB-1", "MPXV.IIB-2", "MPXV.IIB-33", 
                     "MPXV.IB-1", "MPXV.IB-2", "MPXV.IB-3")

# 检查并重新排列列
if(all(desired_columns %in% colnames(selected_data))) {
  selected_data <- selected_data[, desired_columns]
  cat("列顺序已调整\n")
} else {
  cat("期望的列名不完全匹配，当前列名:\n")
  print(colnames(selected_data))
  stop("列名不匹配，请检查列名并修改desired_columns")
}

# 8. 准备热图数据
# 设置行名
rownames(selected_data) <- selected_data$Gene

# 提取表达矩阵
expr_matrix <- as.matrix(selected_data[, -1])  # 移除Gene列
cat("表达矩阵维度:", dim(expr_matrix), "\n")

# 检查数据
cat("表达矩阵前3行:\n")
print(expr_matrix[1:min(3, nrow(expr_matrix)), ])

# 9. 对行进行标准化（Z-score）
cat("对行进行标准化...\n")
expr_matrix_scaled <- t(scale(t(expr_matrix)))
# 处理NA/无限值
expr_matrix_scaled[is.na(expr_matrix_scaled)] <- 0
expr_matrix_scaled[!is.finite(expr_matrix_scaled)] <- 0

# 10. 创建列注释
# 根据实际列名创建分组
col_annotation <- data.frame(
  Group = factor(
    c(rep("MOCK", 3), rep("MPXV.IIB", 3), rep("MPXV.IB", 3)),
    levels = c("MOCK", "MPXV.IIB", "MPXV.IB")
  )
)


rownames(col_annotation) <- colnames(expr_matrix)

# 11. 设置颜色
# 列注释颜色
group_colors <- c(
  "MOCK" = "#708090",     # 灰色
  "MPXV.IIB" = "#A52A2A", # 棕色
  "MPXV.IB" = "#8B0000"   # 深红色
)

annotation_colors <- list(
  Group = group_colors
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
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # 不对行进行聚类，保持指定顺序
  cluster_cols = FALSE,  # 不对列进行聚类，保持指定顺序
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  border_color = NA,
  main = "猴痘病毒感染小鼠差异表达基因热图\n(MOCK vs MPXV.IIB vs MPXV.IB)",
  gaps_col = c(3, 6),  # 在MOCK和MPXV.IIB之间，以及MPXV.IIB和MPXV.IB之间添加分割线
  display_numbers = FALSE,  # 不显示数值
  number_color = "black",
  number_format = "%.1f",
  cellwidth = 20,
  cellheight = 10,
  treeheight_row = 0,  # 设置为0，不显示行聚类树
  treeheight_col = 0   # 设置为0，不显示列聚类树
)
# 









rownames(col_annotation) <- colnames(expr_matrix)

# 11. 设置颜色
# 列注释颜色
group_colors <- c(
  "MOCK" = "#708090",     # 灰色
  "MPXV.IIB" = "#A52A2A", # 棕色
  "MPXV.IB" = "#8B0000"   # 深红色
)

annotation_colors <- list(
  Group = group_colors
)

# 12. 绘制横向热图
cat("开始绘制横向热图...\n")

# 设置颜色映射
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# 转置矩阵，使基因在列，样本在行
expr_matrix_transposed <- t(expr_matrix_scaled)

# 创建行注释（原来的列注释现在变为行注释）
row_annotation <- col_annotation
rownames(row_annotation) <- rownames(expr_matrix_transposed)

# 绘制横向热图
heatmap_plot <- pheatmap(
  expr_matrix_transposed,  # 使用转置后的矩阵
  color = color_palette,
  scale = "column",  # 对列（基因）进行标准化
  annotation_row = row_annotation,  # 行注释
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # 不对行（样本）进行聚类
  cluster_cols = FALSE,  # 不对列（基因）进行聚类
  show_rownames = TRUE,  # 显示行名（样本）
  show_colnames = TRUE,  # 显示列名（基因）
  fontsize_row = 10,     # 行名字体大小
  fontsize_col = 8,      # 列名字体大小
  fontsize = 10,
  border_color = NA,
  main = "",
  gaps_row = c(3, 6),  # 在样本组之间添加分割线
  display_numbers = FALSE,  # 不显示数值
  number_color = "black",
  number_format = "%.1f",
 
  treeheight_row = 0,  # 设置为0，不显示行聚类树
  treeheight_col = 0,  # 设置为0，不显示列聚类树
  angle_col = 45       # 列名（基因名）旋转45度
)







