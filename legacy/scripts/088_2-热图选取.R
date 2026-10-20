############################################################
# 基因表达时间序列热图绘制
# 数据来源：D:\Desktop\1.xlsx
# 时间顺序：mock → 3d → 5d → 7d → 10d（每个时间点3个重复）
############################################################

# 加载所需包
library(readxl)    # 读取Excel文件
library(dplyr)     # 数据处理
library(pheatmap)  # 绘制热图
library(RColorBrewer)  # 颜色调色板
library(grid)      # 图形输出

# 1. 读取Excel数据 --------------------------------------------------------
# 注意：Windows路径用双反斜杠"\\"或正斜杠"/"
file_path <- "D:/Desktop/1.xlsx"  # 数据文件路径
expr_data <- read_excel(file_path)  # 读取Excel文件

# 查看数据结构（可选，用于调试）
cat("数据维度（行×列）：", nrow(expr_data), "×", ncol(expr_data), "\n")
cat("列名：", paste(colnames(expr_data), collapse = ", "), "\n")


# 2. 数据预处理 ----------------------------------------------------------
# 2.1 清理列名和基因名（去除前后空格）
colnames(expr_data) <- trimws(colnames(expr_data))  # 清理列名空格
expr_data$Gene <- trimws(expr_data$Gene)            # 清理基因名空格

# 2.2 定义目标样本顺序（按时间：mock→3d→5d→7d→10d，每个时间点3个重复）
target_sample_order <- c(
  "c57_mock_1_P", "c57_mock_2_P", "c57_mock_3_P",  # mock组（对照）
  "C57_3d_1_P",   "C57_3d_2_P",   "C57_3d_3_P",    # 3d组（感染后3天）
  "C57_5d_1_P",   "C57_5d_2_P",   "C57_5d_3_P",    # 5d组（感染后5天）
  "C57_7d_1_P",   "C57_7d_2_P",   "C57_7d_3_P",    # 7d组（感染后7天）
  "C57_10d_1_P",  "C57_10d_2_P",  "C57_10d_3_P"     # 10d组（感染后10天）
)

# 2.3 筛选并排序样本列（仅保留数据中存在的样本）
existing_samples <- colnames(expr_data)[colnames(expr_data) %in% target_sample_order]
ordered_samples <- existing_samples[order(match(existing_samples, target_sample_order))]  # 按目标顺序排列

# 检查缺失样本（可选，用于调试）
missing_samples <- setdiff(target_sample_order, existing_samples)
if (length(missing_samples) > 0) {
  cat("⚠️ 警告：以下样本在数据中不存在，已跳过：\n", paste(missing_samples, collapse = ", "), "\n")
}

# 2.4 提取表达矩阵（行=基因，列=样本，值=表达量）
# 确保表达量为数值型（处理可能的字符型数据）
expr_matrix <- expr_data[, c("Gene", ordered_samples)]  # 保留基因列和目标样本列
expr_matrix[, ordered_samples] <- lapply(expr_matrix[, ordered_samples], function(x) as.numeric(as.character(x)))  # 转换为数值

# 转换为矩阵（行名=基因名，关键！确保行标签是基因名称）
expr_mat <- as.matrix(expr_matrix[, ordered_samples])  # 样本表达矩阵
rownames(expr_mat) <- expr_matrix$Gene  # 行名设为基因名称（核心：避免数字索引）


# 3. 构建注释信息（可选，标注时间点） --------------------------------------
# 3.1 列注释：标记每个样本的时间点（mock/3d/5d/7d/10d）
time_points <- sapply(strsplit(ordered_samples, "_"), function(x) {
  if (tolower(x[2]) == "mock") "mock" else x[2]  # 提取时间点（如"3d"）
})
annotation_col <- data.frame(
  Timepoint = factor(time_points, levels = c("mock", "3d", "5d", "7d", "10d")),  # 时间点顺序
  row.names = ordered_samples  # 行名=样本名
)

# 3.2 注释颜色（时间点用不同颜色区分）
annotation_colors <- list(
  Timepoint = c(
    mock = "#708090",  # 灰色（对照）
    `3d` = "#1B9E77",  # 绿色（3天）
    `5d` = "#D95F02",  # 橙色（5天）
    `7d` = "#7570B3",  # 紫色（7天）
    `10d` = "#E7298A"  # 粉色（10天）
  )
)


# 4. 绘制热图 ------------------------------------------------------------
# 4.1 热图参数设置
heatmap_params <- list(
  mat = expr_mat,  # 表达矩阵（行=基因，列=样本）
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # 蓝-白-红渐变（低→高表达）
  scale = "row",  # 行标准化（Z-score）：突出基因在时间上的变化趋势
  cluster_rows = TRUE,  # 基因聚类（按表达模式分组）
  cluster_cols = FALSE,  # 列不聚类（固定时间顺序）
  show_rownames = TRUE,  # 显示基因名称（核心：行标签为基因名）
  show_colnames = TRUE,  # 显示样本名称
  fontsize_row = 8,  # 基因名称字体大小（根据基因数量调整，避免重叠）
  fontsize_col = 8,  # 样本名称字体大小
  main = "MPXV感染后C57小鼠免疫基因表达时间序列热图",  # 标题
  annotation_col = annotation_col,  # 列注释（时间点）
  annotation_colors = annotation_colors,  # 注释颜色
  gaps_col = c(3, 6, 9, 12),  # 在时间点间添加间隙（mock后、3d后、5d后、7d后）
  na_col = "grey",  # 缺失值用灰色表示
  border_color = NA  # 去除单元格边框（更美观）
)

# 4.2 绘制热图并保存
set.seed(123)  # 保证聚类结果可重复
p_heatmap <- do.call(pheatmap, heatmap_params)

# 保存为PDF（矢量图，高清）
pdf("s基因表达时间序列热图.pdf", width = 12, height =25)
grid.newpage()
grid.draw(p_heatmap$gtable)
dev.off()

# 保存为PNG（位图，方便插入文档）
png("基因表达时间序列热图.png", width = 2000, height =5000, res = 150)
grid.newpage()
grid.draw(p_heatmap$gtable)
dev.off()


# 5. 输出结果说明 --------------------------------------------------------
cat("\n✅ 热图绘制完成！\n")
cat("📄 输出文件：\n",
    "- 矢量图：基因表达时间序列热图.pdf\n",
    "- 位图：基因表达时间序列热图.png\n")
cat("🔍 热图说明：\n",
    "- 行：基因名称（如Prss46、Slfn2，无数字索引）\n",
    "- 列：样本（按时间顺序：mock→3d→5d→7d→10d，每个时间点3个重复）\n",
    "- 颜色：蓝色=低表达，白色=平均表达，红色=高表达\n",
    "- 注释：顶部颜色条标注时间点（灰色=mock，绿色=3d，橙色=5d，紫色=7d，粉色=10d）\n")