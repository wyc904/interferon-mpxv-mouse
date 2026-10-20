############################################################
# 按时间顺序重排表达矩阵 & 输出 Excel 方便挑选时间序列差异基因
# 原始文件：D:/Desktop/1.xlsx
# 原始列名：
# Gene, C57_10d_1_P, C57_10d_2_P, C57_10d_3_P,
#       C57_3d_1_P,  C57_3d_2_P,  C57_3d_3_P,
#       C57_5d_1_P,  C57_5d_2_P,  C57_5d_3_P,
#       C57_7d_1_P,  C57_7d_2_P,  C57_7d_3_P,
#       c57_mock_1_P, c57_mock_2_P, c57_mock_3_P
############################################################

# 如未安装 writexl，请先运行：install.packages("writexl")
library(readxl)
library(writexl)

## 1. 读取原始数据 ---------------------------------------------------------

expr <- read_excel("D:/Desktop/1.xlsx")

# 保证第一列叫 Gene（如果本来就是，就不会有影响）
colnames(expr)[1] <- "Gene"

cat("原始列名：\n")
print(colnames(expr))

## 2. 定义各时间点的样本列名 ----------------------------------------------

cols_10d <- c("C57_10d_1_P", "C57_10d_2_P", "C57_10d_3_P")
cols_3d  <- c("C57_3d_1_P",  "C57_3d_2_P",  "C57_3d_3_P")
cols_5d  <- c("C57_5d_1_P",  "C57_5d_2_P",  "C57_5d_3_P")
cols_7d  <- c("C57_7d_1_P",  "C57_7d_2_P",  "C57_7d_3_P")
cols_mock<- c("c57_mock_1_P","c57_mock_2_P","c57_mock_3_P")

# 想要的时间顺序：Gene → mock → 3d → 5d → 7d → 10d
target_order <- c(
  "Gene",
  cols_mock,
  cols_3d,
  cols_5d,
  cols_7d,
  cols_10d
)

# 只保留在数据中真实存在的列名（防止名字有微小差异时报错）
target_order <- target_order[target_order %in% colnames(expr)]

cat("按时间重排后的列顺序：\n")
print(target_order)

## 3. 按时间顺序重排列 ----------------------------------------------------

expr_ordered <- expr[, target_order]

## 4. 计算各时间点的平均表达（基于重排后的矩阵） --------------------------

# 防止列名缺失，这里再和实际列名取交集
cols_mock_use <- cols_mock[cols_mock %in% colnames(expr_ordered)]
cols_3d_use   <- cols_3d[cols_3d %in% colnames(expr_ordered)]
cols_5d_use   <- cols_5d[cols_5d %in% colnames(expr_ordered)]
cols_7d_use   <- cols_7d[cols_7d %in% colnames(expr_ordered)]
cols_10d_use  <- cols_10d[cols_10d %in% colnames(expr_ordered)]

# 先复制一份，用于存放平均值
expr_means <- expr_ordered

# 使用基础 R 计算每一行在对应时间点的平均值
expr_means$mean_mock <- if (length(cols_mock_use) > 0) {
  rowMeans(as.matrix(expr_ordered[, cols_mock_use, drop = FALSE]), na.rm = TRUE)
} else NA_real_

expr_means$mean_3d <- if (length(cols_3d_use) > 0) {
  rowMeans(as.matrix(expr_ordered[, cols_3d_use, drop = FALSE]), na.rm = TRUE)
} else NA_real_

expr_means$mean_5d <- if (length(cols_5d_use) > 0) {
  rowMeans(as.matrix(expr_ordered[, cols_5d_use, drop = FALSE]), na.rm = TRUE)
} else NA_real_

expr_means$mean_7d <- if (length(cols_7d_use) > 0) {
  rowMeans(as.matrix(expr_ordered[, cols_7d_use, drop = FALSE]), na.rm = TRUE)
} else NA_real_

expr_means$mean_10d <- if (length(cols_10d_use) > 0) {
  rowMeans(as.matrix(expr_ordered[, cols_10d_use, drop = FALSE]), na.rm = TRUE)
} else NA_real_

# 只保留 Gene + 各时间点平均值列
expr_means_out <- expr_means[, c("Gene", "mean_mock", "mean_3d", "mean_5d", "mean_7d", "mean_10d")]

## 5. 写出为新的 Excel 文件 -----------------------------------------------

out_path <- "D:/Desktop/1_time_ordered_with_means.xlsx"

write_xlsx(
  list(
    "表达矩阵_按时间排序" = expr_ordered,
    "各时间点平均表达"   = expr_means_out
  ),
  path = out_path
)

cat("✅ 已导出到：", out_path, "\n")
