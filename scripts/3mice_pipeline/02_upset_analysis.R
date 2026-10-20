# Run from repository root
source("scripts/3mice_pipeline/utils/config.R", encoding = "UTF-8")

# ============================================================================
# 干扰素小鼠 RNA-seq Upset 图分析
# 小鼠品系: BALB/c + C57BL/6
# 时间点: mock, 3d, 5d, 7d, 10d
# 数据来源: featureCounts 基因计数矩阵
# 输出: 矢量 PDF Upset 图 + Excel 基因列表
# ============================================================================

# ---- 0. 安装/加载包 ----
required_pkgs <- c("DESeq2", "tidyverse", "UpSetR", "openxlsx", "grid", "viridis")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# 避免科学计数法
options(scipen = 999)

# ---- 1. 读取 featureCounts 数据 ----
cat("========================================\n")
cat("  读取 count 矩阵 ...\n")
cat("========================================\n")

read_fc <- function(file_path, label) {
  cat(sprintf("  读取 %s: %s\n", label, basename(file_path)))
  raw <- read.table(file_path, header = TRUE, sep = "\t",
                    skip = 1, comment.char = "#",
                    check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("    维度: %d 基因 × %d 列\n", nrow(raw), ncol(raw)))
  return(raw)
}

# 路径设置
data_dir <- DATA_DIR_3MICE
balb_raw <- read_fc(file.path(data_dir, "BALB_F_counts_gene(1).txt"), "BALB/c")
c57_raw  <- read_fc(file.path(data_dir, "C57_F_counts_gene(1).txt"),  "C57BL/6")

# ---- 2. 整理 count 矩阵和样本元数据 ----
cat("\n--- 整理样本注释 ---\n")

# 解析样本名称 -> (strain, time, replicate)
parse_sample <- function(colname) {
  # 处理类似 "BALBC_10d_3_F/Aligned.sortedByCoord.out.bam" 的列名
  base <- basename(colname)          # 去路径
  base <- sub("\\.bam$", "", base)   # 去 .bam
  base <- sub("_Aligned.*", "", base) # 去 _Aligned...
  parts <- strsplit(base, "_")[[1]]
  # parts 大致为: [strain_prefix, time, rep, "F"] 或 ["c57", "mock", rep, "F"]
  # 需要处理不同命名风格
  n <- length(parts)
  if (n >= 3) {
    rep_str  <- parts[n-1]  # 倒数第二个是重复编号
    time_str <- parts[n-2]  # 倒数第三个是时间点
    strain_prefix <- paste(parts[1:(n-3)], collapse = "_")
    if (strain_prefix == "") strain_prefix <- parts[1]
    # 统一品系名
    strain <- ifelse(grepl("BALB", strain_prefix, ignore.case = TRUE),
                     "BALB", "C57")
    # 时间点标准化
    time_map <- c("mock" = "mock", "3d" = "3d", "5d" = "5d", "7d" = "7d", "10d" = "10d")
    time_std <- time_map[tolower(time_str)]
    if (is.na(time_std)) time_std <- tolower(time_str)
    return(c(strain = strain, time = time_std, rep = rep_str))
  }
  return(c(strain = "unknown", time = "unknown", rep = "unknown"))
}

# 从每个 count 矩阵中提取 metadata 列(前6列)和 count 列
extract_counts <- function(raw, strain_label) {
  meta_cols <- raw[, 1:6, drop = FALSE]
  colnames(meta_cols) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- raw[, -(1:6), drop = FALSE]

  # 清理列名: 去掉 /Aligned.sortedByCoord.out.bam 后缀, 保留样本ID
  clean_names <- sapply(colnames(count_cols), function(nm) {
    nm <- sub("/Aligned\\.sortedByCoord\\.out\\.bam$", "", nm)
    return(nm)
  })
  colnames(count_cols) <- clean_names

  list(meta = meta_cols, counts = count_cols)
}

balb <- extract_counts(balb_raw, "BALB")
c57  <- extract_counts(c57_raw, "C57")

# 合并两个 count 矩阵 (按基因名)
counts_merged <- merge(
  data.frame(Geneid = balb$meta$Geneid, balb$counts, check.names = FALSE),
  data.frame(Geneid = c57$meta$Geneid,  c57$counts,  check.names = FALSE),
  by = "Geneid", all = FALSE
)
rownames(counts_merged) <- counts_merged$Geneid
counts_merged$Geneid <- NULL

cat(sprintf("\n合并后计数矩阵: %d 基因 × %d 样本\n", nrow(counts_merged), ncol(counts_merged)))
cat("样本列名:", paste(colnames(counts_merged), collapse = ", "), "\n")

# ---- 3. 创建 DESeq2 所需的 colData ----
all_samples <- data.frame(
  sample_id  = colnames(counts_merged),
  stringsAsFactors = FALSE
)

# 解析每个样本
sample_info <- do.call(rbind, lapply(colnames(counts_merged), parse_sample))
colnames(sample_info) <- c("strain", "time", "rep")
sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)

colData <- cbind(all_samples, sample_info)
colData$group <- with(colData, paste(strain, time, sep = "_"))
colData$strain <- as.factor(colData$strain)
colData$time   <- factor(colData$time, levels = c("mock", "3d", "5d", "7d", "10d"))
colData$group  <- as.factor(colData$group)

cat("\n样本分组信息:\n")
print(table(colData$strain, colData$time))

# ---- 4. DESeq2 差异分析 (按品系分别进行) ----
cat("\n========================================\n")
cat("  运行 DESeq2 差异表达分析\n")
cat("========================================\n")

padj_threshold <- 0.05
log2fc_threshold <- 1

# 为每个品系单独运行 DESeq2
run_deseq2_strain <- function(strain_name, counts_mat, colData) {
  cat(sprintf("\n--- %s 品系 DESeq2 分析 ---\n", strain_name))

  # 筛选该品系的样本
  strain_samples <- colData[colData$strain == strain_name, ]
  strain_counts <- counts_mat[, strain_samples$sample_id, drop = FALSE]

  cat(sprintf("  样本数: %d\n", nrow(strain_samples)))
  print(table(strain_samples$time))

  # DESeq2 设计: ~ time (将 mock 作为 reference level)
  strain_samples$time <- relevel(strain_samples$time, ref = "mock")

  dds <- DESeqDataSetFromMatrix(
    countData = strain_counts,
    colData   = strain_samples,
    design    = ~ time
  )

  # 过滤低表达基因 (至少 N 个样本 count >= 10, N = 最小重复数)
  smallest_group <- min(table(strain_samples$time))
  keep <- rowSums(counts(dds) >= 10) >= smallest_group
  dds <- dds[keep, ]
  cat(sprintf("  过滤后保留基因数: %d\n", sum(keep)))

  # 运行 DESeq2
  dds <- DESeq(dds)

  # 提取各时间点 vs mock 的结果
  time_levels <- setdiff(levels(strain_samples$time), "mock")
  deg_list <- list()

  for (tp in time_levels) {
    cat(sprintf("    → %s_vs_mock ...\n", tp))
    res <- results(dds, contrast = c("time", tp, "mock"), alpha = padj_threshold)
    res_df <- as.data.frame(res)
    res_df$Geneid <- rownames(res_df)
    res_df$comparison <- sprintf("%s_%s_vs_mock", strain_name, tp)
    res_df$strain <- strain_name
    res_df$time   <- tp
    deg_list[[tp]] <- res_df
  }

  # 合并该品系所有比较
  all_degs <- do.call(rbind, deg_list)
  rownames(all_degs) <- NULL
  return(all_degs)
}

# 分别跑 BALB 和 C57
balb_degs <- run_deseq2_strain("BALB", counts_merged, colData)
c57_degs  <- run_deseq2_strain("C57",  counts_merged, colData)

all_results <- rbind(balb_degs, c57_degs)

cat(sprintf("\n所有比较基因总数: %d\n", nrow(all_results)))

# ---- 5. 提取显著差异基因 ----
cat("\n========================================\n")
cat("  提取显著差异基因\n")
cat("========================================\n")

sig_genes <- all_results %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    direction = case_when(
      padj < padj_threshold & log2FoldChange >  log2fc_threshold ~ "Up",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down",
      TRUE ~ "Not_sig"
    )
  ) %>%
  filter(direction != "Not_sig") %>%
  distinct(comparison, Geneid, .keep_all = TRUE)

cat(sprintf("显著差异基因条目数: %d\n", nrow(sig_genes)))

# ---- 6. 构建 Upset 图数据 ----
cat("\n--- 构建 Upset 矩阵 ---\n")

comparisons <- unique(sig_genes$comparison)
cat(sprintf("比较组数: %d\n", length(comparisons)))
cat("比较组:", paste(comparisons, collapse = ", "), "\n")

# 构建二元矩阵: 基因 × 比较组
gene_lists <- split(sig_genes$Geneid, sig_genes$comparison)
all_sig_gene_ids <- unique(sig_genes$Geneid)
cat(sprintf("去重后显著基因总数: %d\n", length(all_sig_gene_ids)))

# 确保 conditions_ordered 只包含实际存在的比较组
# 排序: BALB_3d, C57_3d, BALB_5d, C57_5d, ...
time_order <- c("3d", "5d", "7d", "10d")
strain_order <- c("BALB", "C57")

conditions_ordered <- c()
for (tp in time_order) {
  for (st in strain_order) {
    cond_name <- sprintf("%s_%s_vs_mock", st, tp)
    if (cond_name %in% comparisons) {
      conditions_ordered <- c(conditions_ordered, cond_name)
    }
  }
}
cat("可用条件:", paste(conditions_ordered, collapse = ", "), "\n")

# 显示名称
display_names <- setNames(
  gsub("_vs_mock", "", conditions_ordered),
  conditions_ordered
)
# 美化: BALB_3d -> BALB-3d
display_names <- setNames(
  gsub("_", "-", gsub("_vs_mock", "", conditions_ordered)),
  conditions_ordered
)

# 创建二元矩阵
binary_matrix <- data.frame(Geneid = all_sig_gene_ids)
for (cond in conditions_ordered) {
  binary_matrix[[cond]] <- as.integer(binary_matrix$Geneid %in% gene_lists[[cond]])
}

cat(sprintf("二元矩阵维度: %d × %d\n", nrow(binary_matrix), ncol(binary_matrix) - 1))

# 打印各条件下基因数
cat("\n各条件下显著差异基因数:\n")
for (cond in conditions_ordered) {
  n_genes <- length(gene_lists[[cond]])
  cat(sprintf("  %-25s: %d (上: %d, 下: %d)\n",
      display_names[cond], n_genes,
      sum(sig_genes$comparison == cond & sig_genes$direction == "Up"),
      sum(sig_genes$comparison == cond & sig_genes$direction == "Down")))
}

# ---- 7. 绘制 UpSet 图 (矢量 PDF) ----
cat("\n========================================\n")
cat("  生成矢量 Upset 图\n")
cat("========================================\n")

output_dir <- OUT_DIR
pdf_file <- file.path(output_dir, "Upset_DEGs_Interferon_Mice.pdf")

# 定义颜色方案
balb_color <- "#E64B35"   # 红色系 - BALB/c
c57_color  <- "#4DBBD5"   # 蓝色系 - C57BL/6

set_colors <- ifelse(grepl("BALB", conditions_ordered), balb_color, c57_color)

# 计算每个交集的大小
set_sizes <- sapply(conditions_ordered, function(cond) sum(binary_matrix[[cond]]))
max_set_size <- max(set_sizes)

# PDF 输出 (矢量)
pdf(pdf_file, width = 14, height = 10, useDingbats = FALSE)

# 准备 Upset 输入 (去掉 Geneid 列)
upset_input <- binary_matrix[, conditions_ordered, drop = FALSE]
colnames(upset_input) <- display_names[conditions_ordered]

upset(
  upset_input,
  nsets              = length(conditions_ordered),
  sets               = display_names[conditions_ordered],
  keep.order          = TRUE,
  nintersects         = 40,
  mb.ratio            = c(0.55, 0.45),
  order.by            = "freq",
  decreasing          = TRUE,
  number.angles       = 0,
  point.size          = 3.5,
  line.size           = 1,
  mainbar.y.label     = "Intersection Size\n(Number of DEGs)",
  sets.x.label        = "DEGs per Comparison",
  text.scale          = c(1.6, 1.5, 1.2, 1.2, 2, 1.3),
  set_size.show       = TRUE,
  set_size.scale_max  = max_set_size * 1.15,
  set_size.numbers_size = 5,
  main.bar.color      = "#3C5488",
  sets.bar.color      = set_colors,
  shade.color         = "#F0F0F0",
  matrix.color        = "#333333",
  shade.alpha         = 0.15,
  matrix.dot.alpha    = 0.8,
  show.numbers        = "yes"
)

# 添加标题
grid.text(
  "Upset Plot: DEGs Intersection Analysis",
  x = 0.5, y = 0.97,
  gp = gpar(fontsize = 17, fontface = "bold", col = "#2C3E50")
)
grid.text(
  "BALB/c & C57BL/6 Mice — Interferon Response (3d, 5d, 7d, 10d vs Mock)",
  x = 0.5, y = 0.935,
  gp = gpar(fontsize = 12, col = "#7F8C8D", fontface = "italic")
)

# 图例
legend_x <- 0.05
legend_y <- 0.06

# BALB 图例
grid.rect(x = legend_x, y = legend_y,
          width = 0.018, height = 0.025,
          gp = gpar(fill = balb_color, col = balb_color))
grid.text("BALB/c",
          x = legend_x + 0.028, y = legend_y,
          just = "left",
          gp = gpar(fontsize = 11, col = "#2C3E50"))

# C57 图例
grid.rect(x = legend_x, y = legend_y - 0.035,
          width = 0.018, height = 0.025,
          gp = gpar(fill = c57_color, col = c57_color))
grid.text("C57BL/6",
          x = legend_x + 0.028, y = legend_y - 0.035,
          just = "left",
          gp = gpar(fontsize = 11, col = "#2C3E50"))

# 显著性阈值注释
grid.text(
  sprintf("Threshold: |log2FC| > %g,  padj < %g", log2fc_threshold, padj_threshold),
  x = 0.95, y = 0.02,
  just = "right",
  gp = gpar(fontsize = 9, col = "#95A5A6")
)

dev.off()
cat(sprintf("\n✓ PDF 矢量图已保存: %s\n", pdf_file))

# ---- 8. 导出基因列表 (Excel) ----
cat("\n--- 导出 Excel 基因列表 ---\n")

wb <- createWorkbook()

# 8a. 各比较组基因列表
for (cond in conditions_ordered) {
  cond_genes <- sig_genes %>%
    filter(comparison == cond) %>%
    arrange(padj) %>%
    select(Geneid, baseMean, log2FoldChange, lfcSE, pvalue, padj, direction)

  if (nrow(cond_genes) > 0) {
    sheet_name <- substr(display_names[cond], 1, 31) # Excel sheet名最长31字符
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, cond_genes)
    # 冻结首行
    freezePane(wb, sheet_name, firstRow = TRUE)
  }
}

# 8b. 交集摘要
intersection_sizes <- rowSums(upset_input)
intersection_summary <- data.frame(
  IntersectionSize = sort(unique(intersection_sizes)),
  GeneCount = sapply(sort(unique(intersection_sizes)), function(s) sum(intersection_sizes == s))
)
addWorksheet(wb, "Intersection_Summary")
writeData(wb, "Intersection_Summary", intersection_summary)

# 8c. 完整二元矩阵
addWorksheet(wb, "Binary_Matrix")
binary_export <- binary_matrix
colnames(binary_export)[-1] <- display_names[colnames(binary_export)[-1]]
writeData(wb, "Binary_Matrix", binary_export)

# 8d. 所有显著结果 (长格式)
addWorksheet(wb, "All_Sig_Results")
writeData(wb, "All_Sig_Results", sig_genes)

excel_file <- file.path(output_dir, "DEGs_Upset_Results.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat(sprintf("✓ Excel 结果已保存: %s\n", excel_file))

# ---- 9. 汇总输出 ----
cat("\n========================================\n")
cat("  分析完成！\n")
cat("========================================\n")
cat("\n输出文件:\n")
cat(sprintf("  1. %s  — 矢量 PDF Upset 图\n", pdf_file))
cat(sprintf("  2. %s — 基因列表与交集分析\n", excel_file))
cat("\n分析参数:\n")
cat(sprintf("  padj 阈值:       %g\n", padj_threshold))
cat(sprintf("  |log2FC| 阈值:   %g\n", log2fc_threshold))
cat(sprintf("  总比较组数:      %d\n", length(conditions_ordered)))
cat(sprintf("  总显著基因(去重): %d\n", length(all_sig_gene_ids)))
cat("========================================\n")
