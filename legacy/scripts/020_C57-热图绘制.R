############################################################
# 40 个免疫轨迹基因热图（横向布局版本）- 按指定顺序排列
# 去除phase标注，不进行基因聚类
############################################################

library(readxl)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(tidyr)

## 1. 读取表达矩阵 --------------------------------------------------------

# 假设 Excel 第一列是基因名，其余是样本
expr <- read_excel("D:/Desktop/c57All_Genes_Expression_Analysis.xlsx")

# 去掉列名 & 第一列（基因名）中的前后空格
colnames(expr) <- trimws(colnames(expr))
expr[[1]]      <- trimws(expr[[1]])

# 把第一列重命名为 "Gene"
colnames(expr)[1] <- "Gene"
cat("表达矩阵列名：\n")
print(colnames(expr))

## 2. 定义基因顺序（按照您指定的顺序）---------------------------------------

genes_order <- c("Ciita", "H2-Ab1", "H2-D1", "H2-Eb1", "H2-K1", "Psmb8", "Psmb9", "Tap1", "Tap2", "Tapbp", "Tapbpl", "Cd40", "Cd86", "Cd209a", "Clec9a", "Flt3", "Irf8", "Itgax", "Ly75", "Siglech", "Spib", "Tcf4", "Ccl3", "Ccl4", "Il12b", "Tnf", "Ccl19", "Ccl21a", "Cxcl10", "Cxcl9")


# 
# "Cxcl17", "Itln1", "Colec10", "Cd22", "Cd72", "Sit1", 
# "Bach2", "Slamf1", "Stab1", "Trem2", "Angptl1",
# "Slfn2", "Oas1a", "Ifi213", "Ifit1bl1", "Ifit1bl2",  
# "Cyp4f18", "Sgsh", "Usp36", 
# "Rtp1", "Slc26a11", "Cyp11b2", "Sppl3", "Cnot8",
# "Klrc3", "Ykt6", "Washc4", "Rps28", "Dlg3", "Srpk1",
# "Nlrx1", "Cnr2", "Serping1", "Tcaim", "Tapbp", "Itgb7", 
# "Col5a1", "Mmrn2", "Egln2", "Nfe2l1"

## 3. 只保留这 40 个基因的表达 --------------------------------------------

expr_40 <- expr %>%
  dplyr::filter(Gene %in% genes_order)

# 提示哪些基因在表达矩阵中缺失（如果有）
missing_genes <- setdiff(genes_order, expr_40$Gene)
if (length(missing_genes) > 0) {
  cat("⚠ 下列基因在表达矩阵中未找到，将不出现在热图中：\n")
  print(missing_genes)
}

# 按指定顺序排序基因行
available_genes <- genes_order[genes_order %in% expr_40$Gene]
expr_40 <- expr_40 %>%
  dplyr::arrange(match(Gene, available_genes))

## 4. 构造数值表达矩阵（关键：确保是 numeric） ----------------------------

# 把除 Gene 以外所有列转成 numeric
expr_40_num <- expr_40 %>%
  dplyr::mutate(across(-Gene, ~ as.numeric(as.character(.))))

# 提取数值矩阵
expr_40_mat <- as.matrix(expr_40_num[, -1])
rownames(expr_40_mat) <- expr_40_num$Gene

# 强制转换为数值矩阵
expr_40_mat <- apply(expr_40_mat, 2, as.numeric)
rownames(expr_40_mat) <- expr_40_num$Gene
colnames(expr_40_mat) <- colnames(expr_40_num)[-1]

# 检查数据类型
cat("矩阵数据类型:", class(expr_40_mat), "\n")
cat("元素数据类型:", typeof(expr_40_mat), "\n")
cat("行名示例:", head(rownames(expr_40_mat)), "\n")

## 5. 按时间顺序重排列样本列 -----------------------------------------------

target_order <- c(
  "c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P",
  "C57_3d_1_P",   "C57_3d_2_P",   "C57_3d_3_P",
  "C57_5d_1_P",   "C57_5d_2_P",   "C57_5d_3_P",
  "C57_7d_1_P",   "C57_7d_2_P",   "C57_7d_3_P",
  "C57_10d_1_P",  "C57_10d_2_P",  "C57_10d_3_P"
)

current_samples   <- colnames(expr_40_mat)
available_samples <- target_order[target_order %in% current_samples]
missing_samples   <- setdiff(target_order, current_samples)

if (length(missing_samples) > 0) {
  cat("⚠ 下列样本在表达矩阵中未找到，将不会出现在热图中：\n")
  print(missing_samples)
}

expr_40_mat <- expr_40_mat[, available_samples, drop = FALSE]

## 6. 构建样本信息（模拟原始代码的sample_info结构）-----------------------

# 创建样本信息数据框，模拟原始代码中的sample_info
sample_info <- data.frame(
  sample_name = available_samples,
  timepoint = sapply(strsplit(available_samples, "_"), function(x) {
    if(x[2] == "mock") "mock" else x[2]
  }),
  group_label = ifelse(sapply(strsplit(available_samples, "_"), function(x) x[2]) == "mock", "MOCK", "MPXV_lb"),
  row.names = available_samples
)

# 定义期望的排序顺序
ordered_timepoints <- c("mock", "3d", "5d", "7d", "10d")

# 生成新的样本标签 (例如: "mock-1", "mock-2", "3d-1", ...)
sample_info$new_label <- ave(
  seq_along(sample_info$timepoint), 
  sample_info$timepoint, 
  FUN = function(x) paste0(sample_info$timepoint[x], "-", seq_along(x))
)

# 对样本信息表进行排序
sample_info_ordered <- sample_info %>%
  # 将 timepoint 转换为因子，确保排序顺序正确
  mutate(timepoint_fct = factor(timepoint, levels = ordered_timepoints)) %>%
  # 按时间点因子和新标签排序，确保时间组内顺序稳定
  arrange(timepoint_fct, new_label)

# 提取新的列名和原始的样本名（用于匹配矩阵）
ordered_cols_new_labels <- sample_info_ordered$new_label
ordered_cols_original_names <- rownames(sample_info_ordered) 

## 7. 转置矩阵并重排序样本（实现横向热图）-------------------------------

# 根据排序后的样本名对矩阵进行重排序和重命名
expr_40_mat_ordered <- expr_40_mat[, ordered_cols_original_names] 
colnames(expr_40_mat_ordered) <- ordered_cols_new_labels # 设置新列名

# 对基因表达进行Z-score标准化（按行）
mat_scaled <- t(scale(t(expr_40_mat_ordered)))

# 确保基因按照指定顺序排列
mat_scaled <- mat_scaled[available_genes, ]

# 转置矩阵，实现横向热图布局（样本在行，基因在列）
mat_transposed <- t(mat_scaled)

## 8. 创建注释信息 --------------------------------------------------------

# 行注释（样本注释）
annotation_row <- data.frame(
  Timepoint = factor(sample_info_ordered$timepoint, levels = c("mock", "3d", "5d", "7d", "10d")),
  Condition = factor(sample_info_ordered$group_label, levels = c("MOCK", "MPXV_lb")),
  row.names = sample_info_ordered$new_label
)

# 创建颜色映射（只保留样本相关的颜色）
timepoint_colors <- c("mock" = "#708090", "3d" = "#1B9E77", 
                      "5d" = "#D95F02", "7d" = "#7570B3", "10d" = "#E7298A")
condition_colors <- c("MOCK" = "#708090", "MPXV_lb" = "#A52A2A")

annotation_colors <- list(
  Timepoint = timepoint_colors,
  Condition = condition_colors
)

## 9. 绘制横向热图（不聚类基因，按指定顺序）-------------------------------

# 设置图形参数
color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 绘制热图
# p_heat <- pheatmap(mat_transposed,
#                    color = color_palette,
#                    annotation_row = annotation_row,      # 行注释（样本）
#                    annotation_colors = annotation_colors,
#                    cluster_rows = FALSE,                 # 不聚类行（样本），保持时间顺序
#                    cluster_cols = FALSE,                 # 不聚类列（基因），保持指定顺序
#                    show_colnames = TRUE,                 # 显示基因名
#                    show_rownames = TRUE,                 # 显示样本名
#                    border_color = NA,
#                    fontsize_row = 9,                     # 样本名字体大小
#                    fontsize_col = 8,                     # 基因名字体大小
#                    fontsize = 10,
#                    main = "",
#                    angle_col = 45,                       # 基因名倾斜45度
#                    treeheight_col = 0,                   # 不显示基因聚类树
#                    treeheight_row = 0)                   # 不显示样本聚类树
# 
# p_heat 


# dev.off()

# 
p_heat_vertical <- pheatmap(mat_scaled,
                            color = color_palette,
                            annotation_col = annotation_row,
                            annotation_colors = annotation_colors,
                            cluster_rows = FALSE,                 # 不聚类基因
                            cluster_cols = FALSE,                 # 不聚类样本
                            show_colnames = TRUE,                 # 显示样本名
                            show_rownames = TRUE,                 # 显示基因名
                            border_color = NA,
                            fontsize_row = 8,                     # 基因名字体大小
                            fontsize_col = 9,                     # 样本名字体大小
                            fontsize = 10,
                            main = "",
                            angle_col = 90,                      # 样本名倾斜45度
                            treeheight_col = 0,                   # 不显示样本聚类树
                            treeheight_row = 0)

p_heat_vertical 

# 
# 
# ## 10. 保存热图 -----------------------------------------------------------
# 
# # 保存PDF（横向布局，宽度大于高度）
# # pdf("Heatmap_40Genes_Ordered_Horizontal.pdf", width = 16, height = 10)
# grid::grid.newpage()
# grid::grid.draw(p_heat$gtable)
# dev.off()
# 
# # 保存PNG
# # png("Heatmap_40Genes_Ordered_Horizontal.png", width = 1600, height = 1000, res = 150)
# grid::grid.newpage()
# grid::grid.draw(p_heat$gtable)
# dev.off()
# 
# cat("✅ 横向热图已保存为：\n")
# cat("   - Heatmap_40Genes_Ordered_Horizontal.pdf\n")
# cat("   - Heatmap_40Genes_Ordered_Horizontal.png\n")
# 
# ## 11. 保存基因列表和表达数据 ---------------------------------------------
# 
# # 保存基因列表信息
# gene_info <- data.frame(
#   Gene = available_genes,
#   Order = 1:length(available_genes)
# )
# 
# # write.csv(gene_info, "40_Genes_Ordered_List.csv", row.names = FALSE)
# 
# # 保存标准化后的表达矩阵
# expr_output <- as.data.frame(mat_scaled) %>%
#   tibble::rownames_to_column("Gene")
# # write.csv(expr_output, "40_Genes_Scaled_Expression.csv", row.names = FALSE)


## 12. 统计信息输出 -------------------------------------------------------
