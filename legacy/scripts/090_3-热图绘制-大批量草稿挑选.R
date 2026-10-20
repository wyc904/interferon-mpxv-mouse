############################################################
# 40 个免疫轨迹基因热图：
# 天然免疫启动 → 适应性免疫激活 → 效应执行与免疫调节 → 稳态恢复
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
expr <- read_excel("D:/Desktop/1.xlsx")

# 去掉列名 & 第一列（基因名）中的前后空格
colnames(expr) <- trimws(colnames(expr))
expr[[1]]      <- trimws(expr[[1]])

# 把第一列重命名为 "Gene"
colnames(expr)[1] <- "Gene"
cat("表达矩阵列名：\n")
print(colnames(expr))

## 2. 定义目标基因顺序 ----------------------------------------------------

target_genes <- c("Cxcl17", "Itln1", "Colec10", "Cd22", "Cd72", "Sit1", 
                  "Bach2", "Slamf1", "Stab1", "Trem2", "Angptl1",
                  
                  
                  
                  "Slfn2", "Oas1a", "Ifi213", "Ifit1bl1", "Ifit1bl2",  
                  "Cyp4f18",  "Sgsh",  "Usp36", 
                  "Rtp1",  "Slc26a11", "Cyp11b2", "Sppl3", "Cnot8",
                  
                  
                  
                  "Klrc3", "Ykt6", "Washc4", "Rps28", "Dlg3", "Srpk1",
                  
                  
                  "Nlrx1", "Cnr2", "Serping1", "Tcaim", "Tapbp", "Itgb7", 
                  "Col5a1", "Mmrn2", "Egln2", "Nfe2l1")

# 去除重复基因，保持第一次出现的顺序
target_genes <- target_genes[!duplicated(target_genes)]

cat("目标基因数量（去重后）：", length(target_genes), "\n")

## 3. 筛选目标基因并按照指定顺序排列 --------------------------------------

# 筛选存在的基因
existing_genes <- target_genes[target_genes %in% expr$Gene]
missing_genes <- target_genes[!target_genes %in% expr$Gene]

cat("存在的基因数量：", length(existing_genes), "\n")
if (length(missing_genes) > 0) {
  cat("缺失的基因：", paste(missing_genes, collapse = ", "), "\n")
}

# 按照目标基因顺序筛选和排序
expr_filtered <- expr[expr$Gene %in% existing_genes, ]

# 创建基因顺序因子并排序
gene_order <- factor(expr_filtered$Gene, levels = existing_genes)
expr_sorted <- expr_filtered[order(gene_order), ]

## 4. 重新排列样本列（按时间点顺序） --------------------------------------

# 定义时间点顺序
timepoint_order <- c("mock", "3d", "5d", "7d", "10d")

# 获取当前列名（除了Gene列）
sample_cols <- colnames(expr_sorted)[-1]

# 按照时间点顺序对样本列进行排序
sorted_cols <- c()
for (timepoint in timepoint_order) {
  timepoint_cols <- sample_cols[grepl(timepoint, sample_cols, ignore.case = TRUE)]
  sorted_cols <- c(sorted_cols, timepoint_cols)
}

# 重新排列数据框列
expr_final <- expr_sorted[, c("Gene", sorted_cols)]

## 5. 准备热图数据 --------------------------------------------------------

# 提取数值矩阵（转置，使基因在行，样本在列）
heatmap_data <- as.matrix(expr_final[, -1])
rownames(heatmap_data) <- expr_final$Gene

# 对行（基因）进行Z-score标准化
heatmap_data_scaled <- t(scale(t(heatmap_data)))

## 6. 创建样本注释（按时间点） --------------------------------------------

# 从列名提取时间点信息
sample_annotation <- data.frame(
  Timepoint = factor(
    sapply(colnames(heatmap_data), function(x) {
      if (grepl("mock", x, ignore.case = TRUE)) return("Mock")
      if (grepl("3d", x)) return("3d")
      if (grepl("5d", x)) return("5d")
      if (grepl("7d", x)) return("7d")
      if (grepl("10d", x)) return("10d")
      return("Unknown")
    }),
    levels = c("Mock", "3d", "5d", "7d", "10d")
  )
)
rownames(sample_annotation) <- colnames(heatmap_data)

## 7. 设置颜色 ------------------------------------------------------------

# 时间点颜色
timepoint_colors <- list(
  Timepoint = c(Mock = "#E6F2FF", "3d" = "#B3D9FF", "5d" = "#80BFFF", 
                "7d" = "#4DA6FF", "10d" = "#1A8CFF")
)

# 热图颜色
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

## 8. 绘制热图 ------------------------------------------------------------

pheatmap(heatmap_data_scaled,
         color = heatmap_colors,
         scale = "row",  # 行标准化
         cluster_rows = FALSE,  # 不聚类行（保持指定顺序）
         cluster_cols = FALSE,  # 不聚类列（保持时间点顺序）
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = sample_annotation,
         annotation_colors = timepoint_colors,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "免疫轨迹基因表达热图 (Mock → 3d → 5d → 7d → 10d)",
         gaps_col = cumsum(table(sample_annotation$Timepoint)[-length(table(sample_annotation$Timepoint))])
)

cat("热图绘制完成！\n")
cat("基因顺序：按照指定的56个基因顺序（去重后）\n")
cat("时间点顺序：Mock → 3d → 5d → 7d → 10d\n")