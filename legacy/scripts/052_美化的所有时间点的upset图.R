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
raw_data <- read.xlsx("绘制.xlsx", sheet = 1)

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

# 创建条件标签
display_names <- c(
  "c57_3d" = "C57-3d",
  "c57_5d" = "C57-5d", 
  "c57_7d" = "C57-7d",
  "c57_10d" = "C57-10d",
  "8bc_3d" = "8BC-3d",
  "8bc_5d" = "8BC-5d",
  "8bc_7d" = "8BC-7d",
  "8bc_10d" = "8BC-10d"
)

total_conditions <- names(display_names)
total_gene_lists <- list()

# 填充各个条件的基因列表
for(cond in total_conditions) {
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

# 添加二进制列
for(cond in total_conditions) {
  total_binary[[cond]] <- as.integer(total_binary$Gene %in% total_gene_lists[[cond]])
}

cat("\n二元矩阵维度:", dim(total_binary), "\n")

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
  upset_plot <- upset(
    total_binary[, -1],  # 移除Gene列
    nsets = length(total_conditions),
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
    sets.bar.color = "#A23B72",
    shade.color = "#F5F5F5",
    matrix.color = "#333333",
    mainbar.y.max = NULL,
    shade.alpha = 0.2,
    matrix.dot.alpha = 0.8,
    empty.intersections = NULL,
    show.numbers = "yes",
    keep.order = FALSE
  )
  
  # 添加自定义注释
  grid.text("Intersection Analysis of Differentially Expressed Genes", 
            x = 0.5, y = 0.97, 
            gp = gpar(fontsize = 16, fontface = "bold", col = "#333333"))
  
  grid.text("C57 and 8BC mice at 3, 5, 7, 10 days post-infection", 
            x = 0.5, y = 0.935, 
            gp = gpar(fontsize = 12, col = "#333333", fontface = "italic"))
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
# 5. 创建分面版本
