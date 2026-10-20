# 加载必要的包
library(tidyverse)
library(UpSetR)
library(openxlsx)
library(scales)
library(viridis)
library(grid)
library(gridExtra)
library(readr)

# 设置科学计算显示格式
options(scipen = 999)

# ============================================================================
# 1. 数据加载和预处理
# ============================================================================

# 读取数据
cat("正在读取数据...\n")
raw_data <- read.xlsx("C57-8BC-绘制UPSET-总.xlsx", sheet = 1)

# 检查数据结构
cat("原始数据维度:", dim(raw_data), "\n")
cat("列名:", colnames(raw_data), "\n")

# 查看前几行
cat("\n前6行数据预览:\n")
print(head(raw_data))

# 统一的数据清洗函数
clean_numeric_data <- function(data) {
  # 定义数值列
  numeric_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  
  # 创建清理后的数据
  cleaned_data <- data
  
  # 确保基本列类型正确
  cleaned_data <- cleaned_data %>%
    mutate(
      time = as.character(time),
      mice = as.character(mice),
      gene = as.character(gene)
    )
  
  # 使用parse_number安全转换数值列
  for(col in numeric_cols) {
    if(col %in% colnames(cleaned_data)) {
      # 转换为字符，然后使用parse_number
      cleaned_data[[col]] <- parse_number(
        as.character(cleaned_data[[col]]), 
        na = c("", "NA", "#N/A", "#VALUE!", "Inf", "-Inf", "inf", "-inf")
      )
    }
  }
  
  return(cleaned_data)
}

# 应用数据清洗
cat("\n开始数据清洗...\n")
data <- clean_numeric_data(raw_data)

# ============================================================================
# 2. 提取显著差异基因
# ============================================================================

# 定义显著性阈值
padj_threshold <- 0.05
log2fc_threshold <- 1

# 提取显著差异基因
sig_data <- data %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    direction = case_when(
      padj < padj_threshold & log2FoldChange > log2fc_threshold ~ "Up",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down",
      TRUE ~ "Not_sig"
    )
  ) %>%
  filter(direction != "Not_sig") %>%
  distinct(time, mice, gene, .keep_all = TRUE)

cat("\n显著差异基因提取结果:\n")
cat("  显著差异基因总数:", nrow(sig_data), "\n")

# 检查每个条件的基因数
cat("\n各条件下差异基因数量:\n")
sig_counts <- sig_data %>%
  group_by(mice, time) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(mice, time)
print(sig_counts)

# ============================================================================
# 3. 准备UpSet图数据
# ============================================================================

# 定义所有时间点和小鼠品系
time_points <- c("3d", "5d", "7d", "10d")
mouse_strains <- c("c57", "8bc")

# 定义条件顺序 - 按照c57-3d, 8bc-3d, c57-5d, 8bc-5d...的顺序
# 这是您要求的顺序：每个时间点内，先c57后8bc
conditions_ordered <- c()
for(time_pt in time_points) {
  for(strain in mouse_strains) {
    cond_name <- paste(strain, time_pt, sep = "_")
    conditions_ordered <- c(conditions_ordered, cond_name)
  }
}

# 定义显示名称 - 按照相同的顺序
display_names <- c(
  "c57_3d" = "C57-3d",
  "8bc_3d" = "8BC-3d",  # 注意：这里改为了8bc_3d
  "c57_5d" = "C57-5d", 
  "8bc_5d" = "8BC-5d",  # 注意：这里改为了8bc_5d
  "c57_7d" = "C57-7d",
  "8bc_7d" = "8BC-7d",  # 注意：这里改为了8bc_7d
  "c57_10d" = "C57-10d",
  "8bc_10d" = "8BC-10d"  # 注意：这里改为了8bc_10d
)

total_gene_lists <- list()

# 填充各个条件的基因列表
for(cond in conditions_ordered) {
  parts <- str_split(cond, "_", simplify = TRUE)
  mice_type <- parts[1]
  time_point <- parts[2]
  
  # 筛选该条件下所有显著差异基因
  cond_genes <- sig_data %>%
    filter(mice == mice_type, time == time_point) %>%
    pull(gene) %>%
    unique()
  
  total_gene_lists[[cond]] <- cond_genes
  cat(paste(display_names[cond], "总差异基因数:", length(cond_genes), "\n"))
}

# 创建总的二元矩阵
all_sig_genes <- unique(sig_data$gene)
cat("\n所有显著差异基因总数（去重后）:", length(all_sig_genes), "\n")

total_binary <- data.frame(Gene = all_sig_genes)

# 按照指定顺序添加二进制列
for(cond in conditions_ordered) {
  total_binary[[cond]] <- as.integer(total_binary$Gene %in% total_gene_lists[[cond]])
}

cat("\n二元矩阵维度:", dim(total_binary), "\n")
cat("\n列顺序:", colnames(total_binary)[-1], "\n")  # 显示列顺序

# 检查是否有基因存在
if (length(all_sig_genes) == 0) {
  stop("没有找到显著差异基因，无法生成UpSet图。")
}

# ============================================================================
# 4. 生成主UpSet图
# ============================================================================

cat("\n正在生成高质量UpSet图...\n")

# 设置图形参数
pdf("UpSet_Total_DEGs_Journal_Quality.pdf", 
    width = 12, 
    height = 9, 
    useDingbats = FALSE)

# 创建基础UpSet图
# 检查数据，确保有足够的交集
intersect_counts <- colSums(total_binary[, -1])
if (sum(intersect_counts) > 0) {
  # 创建颜色方案 - 为C57和8BC使用不同的颜色
  c57_color <- "#2E86AB"  # 蓝色
  b8bc_color <- "#A23B72"  # 红色
  
  # 创建自定义颜色向量，按照条件顺序
  set_colors <- c()
  for(cond in conditions_ordered) {
    if(startsWith(cond, "c57_")) {
      set_colors <- c(set_colors, c57_color)
    } else {
      set_colors <- c(set_colors, b8bc_color)
    }
  }
  
  # 生成UpSet图
  upset_plot <- upset(
    total_binary[, -1],  # 移除Gene列
    nsets = length(conditions_ordered),
    sets = conditions_ordered,  # 明确指定集合顺序
    nintersects = 30,
    mb.ratio = c(0.55, 0.45),
    order.by = "freq",
    decreasing = TRUE,
    number.angles = 0,
    point.size = 3,
    line.size = 0.8,
    mainbar.y.label = "Number of Genes",
    sets.x.label = "DEGs per Condition",
    text.scale = c(1.4, 1.4, 1, 1.2, 1.6, 1.2),
    set_size.show = TRUE,
    set_size.scale_max = max(colSums(total_binary[, -1])) * 1.2,
    set_size.numbers_size = 4,
    main.bar.color = "#2E86AB",
    sets.bar.color = set_colors,  # 使用自定义颜色
    shade.color = "#F5F5F5",
    matrix.color = "#333333",
    mainbar.y.max = NULL,
    shade.alpha = 0.2,
    matrix.dot.alpha = 0.8,
    empty.intersections = NULL,
    show.numbers = "yes",
    keep.order = TRUE  # 保持集合顺序
  )
  
  # 添加自定义注释
  grid.text("Intersection Analysis of Differentially Expressed Genes", 
            x = 0.5, y = 0.97, 
            gp = gpar(fontsize = 16, fontface = "bold", col = "#333333"))
  
  grid.text("C57 and 8BC mice at 3, 5, 7, 10 days post-infection", 
            x = 0.5, y = 0.935, 
            gp = gpar(fontsize = 12, col = "#333333", fontface = "italic"))
  
  # 添加图例
  legend_x_start <- 0.05
  legend_y_start <- 0.05
  legend_spacing <- 0.08
  
  # C57图例
  grid.rect(x = legend_x_start, y = legend_y_start, 
            width = 0.015, height = 0.02,
            gp = gpar(fill = c57_color, col = NA))
  grid.text("C57", 
            x = legend_x_start + 0.02, y = legend_y_start,
            just = "left", 
            gp = gpar(fontsize = 10, col = "#333333"))
  
  # 8BC图例
  grid.rect(x = legend_x_start, y = legend_y_start - 0.03, 
            width = 0.015, height = 0.02,
            gp = gpar(fill = b8bc_color, col = NA))
  grid.text("8BC", 
            x = legend_x_start + 0.02, y = legend_y_start - 0.03,
            just = "left", 
            gp = gpar(fontsize = 10, col = "#333333"))
  
} else {
  # 如果没有交集，创建一个简单的消息图
  plot.new()
  grid.text("No significant intersections found", 
            x = 0.5, y = 0.5, 
            gp = gpar(fontsize = 16, fontface = "bold", col = "#333333"))
  grid.text("The data does not contain overlapping differentially expressed genes", 
            x = 0.5, y = 0.45, 
            gp = gpar(fontsize = 12, col = "#333333"))
}

dev.off()

cat("主图已保存: UpSet_Total_DEGs_Journal_Quality.pdf\n")

# ============================================================================
# 5. 保存交集结果
# ============================================================================

# 创建一个函数来提取指定交集的基因
get_intersection_genes <- function(conditions) {
  subset_data <- total_binary
  for(cond in conditions_ordered) {
    if(cond %in% conditions) {
      subset_data <- subset_data[subset_data[[cond]] == 1, ]
    } else {
      subset_data <- subset_data[subset_data[[cond]] == 0, ]
    }
  }
  return(subset_data$Gene)
}

# 保存所有显著差异基因列表
cat("\n正在保存基因列表...\n")

# 创建Excel工作簿
wb <- createWorkbook()

# 保存各条件的基因列表
for(cond in conditions_ordered) {
  genes <- total_binary[total_binary[[cond]] == 1, ]$Gene
  if(length(genes) > 0) {
    gene_data <- data.frame(
      Gene = genes,
      Condition = display_names[cond]
    )
    addWorksheet(wb, display_names[cond])
    writeData(wb, display_names[cond], gene_data)
  }
}

# 保存交集信息
# 找出前20个最大的交集
intersection_sizes <- apply(total_binary[, -1], 1, function(row) {
  sum(row)
})

# 统计各大小交集的数量
intersection_counts <- table(intersection_sizes)
intersection_summary <- data.frame(
  IntersectionSize = as.numeric(names(intersection_counts)),
  NumberOfIntersections = as.numeric(intersection_counts),
  NumberOfGenes = sapply(as.numeric(names(intersection_counts)), function(size) {
    sum(intersection_sizes == size)
  })
)

addWorksheet(wb, "Intersection_Summary")
writeData(wb, "Intersection_Summary", intersection_summary)

# 保存二元矩阵
addWorksheet(wb, "Binary_Matrix")
writeData(wb, "Binary_Matrix", total_binary)

# 保存Excel文件
saveWorkbook(wb, "DEGs_Intersection_Analysis.xlsx", overwrite = TRUE)
cat("Excel文件已保存: DEGs_Intersection_Analysis.xlsx\n")

cat("\n分析完成！\n")
cat("=============================================\n")
cat("输出文件:\n")
cat("1. UpSet_Total_DEGs_Journal_Quality.pdf - 高质量UpSet图\n")
cat("2. DEGs_Intersection_Analysis.xlsx - 基因列表和交集信息\n")
cat("=============================================\n")