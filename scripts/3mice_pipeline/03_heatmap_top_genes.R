# Run from repository root
source("scripts/3mice_pipeline/utils/config.R", encoding = "UTF-8")

# ============================================================================
# 干扰素小鼠 Top DEGs 热图 (按品系 C57 & BALB)
# 数据: featureCounts 计数矩阵
# 方法: DESeq2 + vst标准化 + Z-score行缩放 + pheatmap
# 输出: 矢量PDF + Excel
# ============================================================================

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(viridis)
library(grid)

options(scipen = 999)

# ---- 0. 参数设置 ----
TOP_N_PER_COMPARISON <- 30   # 每个比较组选取的top基因数
PADJ_THRESHOLD <- 0.05       # 显著性阈值
LOG2FC_THRESHOLD <- 1        # |log2FC| 阈值

out_dir <- OUT_DIR

# ---- 1. 读取数据 ----
cat("========================================\n")
cat("  读取 count 矩阵\n")
cat("========================================\n")

read_fc <- function(fp, label) {
  cat(sprintf("  %s: %s\n", label, basename(fp)))
  raw <- read.table(fp, header = TRUE, sep = "\t",
                    skip = 1, comment.char = "#",
                    check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("    维度: %d × %d\n", nrow(raw), ncol(raw)))
  return(raw)
}

data_dir <- DATA_DIR_3MICE
balb_raw <- read_fc(file.path(data_dir, "BALB_F_counts_gene(1).txt"), "BALB/c")
c57_raw  <- read_fc(file.path(data_dir, "C57_F_counts_gene(1).txt"),  "C57BL/6")

# ---- 2. 整理矩阵 ----
extract_counts <- function(raw) {
  meta_cols  <- raw[, 1:6, drop = FALSE]
  colnames(meta_cols) <- c("Geneid","Chr","Start","End","Strand","Length")
  count_cols <- raw[, -(1:6), drop = FALSE]
  # 清理列名: 去掉 /Aligned.sortedByCoord.out.bam
  clean <- sapply(colnames(count_cols), function(nm) {
    sub("/Aligned\\.sortedByCoord\\.out\\.bam$", "", nm)
  })
  colnames(count_cols) <- clean
  list(meta = meta_cols, counts = count_cols)
}

balb <- extract_counts(balb_raw)
c57  <- extract_counts(c57_raw)

# 合并
counts_all <- merge(
  data.frame(Geneid = balb$meta$Geneid, balb$counts, check.names = FALSE),
  data.frame(Geneid = c57$meta$Geneid,  c57$counts,  check.names = FALSE),
  by = "Geneid", all = FALSE
)
rownames(counts_all) <- counts_all$Geneid
counts_all$Geneid <- NULL

cat(sprintf("合并矩阵: %d genes × %d samples\n", nrow(counts_all), ncol(counts_all)))

# ---- 3. 样本注释 ----
parse_sample <- function(colname) {
  parts <- strsplit(colname, "_")[[1]]
  n <- length(parts)
  rep_str  <- parts[n-1]
  time_str <- parts[n-2]
  strain_prefix <- paste(parts[1:(n-3)], collapse = "_")
  if (strain_prefix == "") strain_prefix <- parts[1]
  strain <- ifelse(grepl("BALB", strain_prefix, ignore.case = TRUE), "BALB", "C57")
  time_map <- c("mock"="mock","3d"="3d","5d"="5d","7d"="7d","10d"="10d")
  time_std <- time_map[tolower(time_str)]
  if (is.na(time_std)) time_std <- tolower(time_str)
  c(strain = strain, time = time_std, rep = rep_str)
}

colData <- data.frame(
  sample_id = colnames(counts_all),
  stringsAsFactors = FALSE
)
info <- do.call(rbind, lapply(colnames(counts_all), parse_sample))
colnames(info) <- c("strain","time","rep")
colData <- cbind(colData, as.data.frame(info, stringsAsFactors = FALSE))
colData$strain <- as.factor(colData$strain)
colData$time   <- factor(colData$time, levels = c("mock","3d","5d","7d","10d"))
colData$group  <- with(colData, factor(paste(strain, time, sep = "_")))

cat("\n样本分组:\n")
print(table(colData$strain, colData$time))

# ---- 4. DESeq2 按品系运行, 并获取 vst 归一化表达 ----
cat("\n========================================\n")
cat("  DESeq2 + VST 归一化\n")
cat("========================================\n")

run_deseq_and_vst <- function(strain_name, counts_mat, colData) {
  cat(sprintf("\n--- %s ---\n", strain_name))
  strain_samples <- colData[colData$strain == strain_name, ]
  strain_counts  <- counts_mat[, strain_samples$sample_id, drop = FALSE]

  # 过滤低表达
  smallest_group <- min(table(strain_samples$time))
  keep <- rowSums(strain_counts >= 10) >= smallest_group
  strain_counts <- strain_counts[keep, , drop = FALSE]
  cat(sprintf("  过滤后: %d genes\n", sum(keep)))

  # DESeq2
  strain_samples$time <- relevel(strain_samples$time, ref = "mock")
  dds <- DESeqDataSetFromMatrix(
    countData = strain_counts,
    colData   = strain_samples,
    design    = ~ time
  )
  dds <- DESeq(dds)

  # VST 归一化
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)

  # DEG 结果
  time_levels <- setdiff(levels(strain_samples$time), "mock")
  deg_list <- list()
  for (tp in time_levels) {
    res <- results(dds, contrast = c("time", tp, "mock"), alpha = PADJ_THRESHOLD)
    res_df <- as.data.frame(res) %>%
      rownames_to_column("Geneid") %>%
      mutate(comparison = sprintf("%s_%s_vs_mock", strain_name, tp),
             strain = strain_name, time = tp)
    deg_list[[tp]] <- res_df
  }
  all_degs <- bind_rows(deg_list)

  list(dds = dds, vsd = vsd, vst_mat = vst_mat,
       degs = all_degs, samples = strain_samples)
}

balb_res <- run_deseq_and_vst("BALB", counts_all, colData)
c57_res  <- run_deseq_and_vst("C57",  counts_all, colData)

# ---- 5. 选取 Top DEGs ----
cat("\n========================================\n")
cat("  选取 Top DEGs\n")
cat("========================================\n")

select_top_genes <- function(degs, top_n) {
  # 筛选显著基因
  sig <- degs %>%
    filter(!is.na(padj), padj < PADJ_THRESHOLD, abs(log2FoldChange) > LOG2FC_THRESHOLD) %>%
    arrange(comparison, padj) %>%
    group_by(comparison) %>%
    slice_head(n = top_n) %>%
    ungroup()
  # 返回去重基因名 (保持按比较组顺序)
  unique_genes <- unique(sig$Geneid)
  cat(sprintf("  各比较组top%d, 去重后: %d genes\n", top_n, length(unique_genes)))
  return(list(genes = unique_genes, full = sig))
}

cat("\nBALB/c:\n")
balb_top <- select_top_genes(balb_res$degs, TOP_N_PER_COMPARISON)
cat("C57BL/6:\n")
c57_top  <- select_top_genes(c57_res$degs,  TOP_N_PER_COMPARISON)

# ---- 6. 热图绘制函数 ----
make_heatmap <- function(strain_name, vst_mat, top_genes, sample_info, deg_info) {
  cat(sprintf("\n--- %s 热图 ---\n", strain_name))

  # 确保基因存在
  genes_present <- intersect(top_genes, rownames(vst_mat))
  cat(sprintf("  可用基因: %d / %d\n", length(genes_present), length(top_genes)))

  # 提取表达矩阵
  expr <- vst_mat[genes_present, , drop = FALSE]

  # 按时间点排列样本列
  sample_info$time_num <- as.numeric(sample_info$time)
  sample_order <- order(sample_info$time_num, sample_info$rep)
  sample_info  <- sample_info[sample_order, ]
  expr <- expr[, sample_info$sample_id, drop = FALSE]

  # Row Z-score
  expr_z <- t(scale(t(expr)))
  expr_z[is.na(expr_z)] <- 0
  expr_z[!is.finite(expr_z)] <- 0

  # 裁剪到 [-3, 3] 避免极端值拉偏色阶
  expr_z[expr_z > 3]  <- 3
  expr_z[expr_z < -3] <- -3

  # 列注释: 时间点
  ann_col <- data.frame(
    Time = sample_info$time,
    row.names = sample_info$sample_id,
    stringsAsFactors = FALSE
  )

  # 时间点颜色 (mock灰 → 3d浅 → 5d中 → 7d深 → 10d最深)
  time_colors <- c(
    "mock" = "#B0BEC5",
    "3d"   = "#81D4FA",
    "5d"   = "#29B6F6",
    "7d"   = "#0277BD",
    "10d"  = "#01579B"
  )
  ann_colors <- list(Time = time_colors)

  # 列间分隔线 (在每个时间点的最后一列后)
  time_counts <- table(sample_info$time)
  gaps_pos <- cumsum(time_counts)[-length(time_counts)]
  names(gaps_pos) <- NULL

  # 行注释: 基因来自哪个比较组 (只标注第一个命中的)
  gene_origin <- deg_info %>%
    filter(Geneid %in% genes_present) %>%
    arrange(padj) %>%
    distinct(Geneid, .keep_all = TRUE) %>%
    select(Geneid, comparison)

  # 保持 expr 中的行顺序
  gene_order <- rownames(expr_z)
  gene_origin <- gene_origin[match(gene_order, gene_origin$Geneid), ]

  # 比较组颜色
  comps <- unique(gene_origin$comparison)
  comp_colors <- setNames(
    colorRampPalette(brewer.pal(min(8, length(comps)), "Set2"))(length(comps)),
    comps
  )

  ann_row <- data.frame(
    DEG_Source = factor(gene_origin$comparison, levels = comps),
    row.names = gene_order
  )
  ann_colors$DEG_Source <- comp_colors

  # 标题
  main_title <- sprintf("Top DEGs Heatmap — %s Mice\n(Interferon Response, |log2FC|>%d, padj<%.2f)",
                        strain_name, LOG2FC_THRESHOLD, PADJ_THRESHOLD)

  # 色阶
  heat_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

  # 绘制
  ph <- pheatmap(
    expr_z,
    color             = heat_colors,
    annotation_col    = ann_col,
    annotation_row    = ann_row,
    annotation_colors = ann_colors,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 6.5,
    fontsize_col      = 9,
    fontsize          = 10,
    border_color      = NA,
    main              = main_title,
    gaps_col          = gaps_pos,
    gaps_row          = NULL,
    display_numbers   = FALSE,
    cellwidth         = 14,
    cellheight        = 9,
    treeheight_row    = 0,
    treeheight_col    = 0,
    angle_col         = 45,
    legend_breaks     = seq(-3, 3, 1),
    legend_labels     = c("-3","-2","-1","0","1","2","3")
  )

  return(list(plot = ph, expr_z = expr_z, expr_raw = expr,
              ann_col = ann_col, ann_row = ann_row))
}

# ---- 7. 生成热图 ----
cat("\n========================================\n")
cat("  绘制热图\n")
cat("========================================\n")

# BALB
balb_ht <- make_heatmap("BALB/c", balb_res$vst_mat,
                         balb_top$genes, balb_res$samples, balb_top$full)

# C57
c57_ht  <- make_heatmap("C57BL/6", c57_res$vst_mat,
                         c57_top$genes, c57_res$samples, c57_top$full)

# ---- 8. 保存矢量 PDF ----
cat("\n--- 保存 PDF ---\n")

pdf_file <- file.path(out_dir, "TopDEGs_Heatmap_C57_BALB.pdf")

# 根据基因数量动态调整高度
balb_h <- max(8, nrow(balb_ht$expr_z) * 0.22 + 3)
c57_h  <- max(8, nrow(c57_ht$expr_z)  * 0.22 + 3)
total_h <- balb_h + c57_h + 2

pdf(pdf_file, width = 14, height = total_h, useDingbats = FALSE)

# 两图上下排列
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(balb_h, c57_h), "inches"))))

# BALB 热图
pushViewport(viewport(layout.pos.row = 1))
grid.text(sprintf("BALB/c — %d Top DEGs", nrow(balb_ht$expr_z)),
          y = 0.98, gp = gpar(fontsize = 14, fontface = "bold"))
grid.draw(balb_ht$plot$gtable)
popViewport()

# C57 热图
pushViewport(viewport(layout.pos.row = 2))
grid.text(sprintf("C57BL/6 — %d Top DEGs", nrow(c57_ht$expr_z)),
          y = 0.98, gp = gpar(fontsize = 14, fontface = "bold"))
grid.draw(c57_ht$plot$gtable)
popViewport()

dev.off()
cat(sprintf("✓ PDF 已保存: %s\n", pdf_file))

# ---- 9. 单独保存每个品系的高分辨率热图 ----
for (strain in c("BALB", "C57")) {
  ht <- if (strain == "BALB") balb_ht else c57_ht
  strain_label <- if (strain == "BALB") "BALBc" else "C57BL6"
  single_pdf <- file.path(out_dir, sprintf("TopDEGs_Heatmap_%s.pdf", strain_label))

  h <- max(7, nrow(ht$expr_z) * 0.24 + 3)
  w <- max(10, ncol(ht$expr_z) * 0.6 + 4)

  pdf(single_pdf, width = w, height = h, useDingbats = FALSE)
  grid.newpage()
  grid.draw(ht$plot$gtable)
  dev.off()
  cat(sprintf("✓ 单图 PDF: %s\n", single_pdf))
}

# ---- 10. 导出 Excel ----
cat("\n--- 导出 Excel ---\n")

excel_file <- file.path(out_dir, "TopDEGs_Heatmap_Data.xlsx")
wb <- createWorkbook()

write_heatmap_sheet <- function(wb, sheet_prefix, ht, deg_full, vst_mat_all) {
  # Z-score 矩阵
  z_df <- as.data.frame(ht$expr_z) %>%
    rownames_to_column("Geneid") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  addWorksheet(wb, paste0(sheet_prefix, "_Zscore"))
  writeData(wb, paste0(sheet_prefix, "_Zscore"), z_df)
  freezePane(wb, paste0(sheet_prefix, "_Zscore"), firstRow = TRUE)

  # 原始 VST 值
  raw_df <- as.data.frame(ht$expr_raw) %>%
    rownames_to_column("Geneid") %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  addWorksheet(wb, paste0(sheet_prefix, "_VST"))
  writeData(wb, paste0(sheet_prefix, "_VST"), raw_df)
  freezePane(wb, paste0(sheet_prefix, "_VST"), firstRow = TRUE)

  # DEG 统计信息
  deg_info <- deg_full %>%
    filter(Geneid %in% rownames(ht$expr_z)) %>%
    select(Geneid, comparison, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(Geneid, padj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  addWorksheet(wb, paste0(sheet_prefix, "_DEGinfo"))
  writeData(wb, paste0(sheet_prefix, "_DEGinfo"), deg_info)
  freezePane(wb, paste0(sheet_prefix, "_DEGinfo"), firstRow = TRUE)
}

write_heatmap_sheet(wb, "BALB", balb_ht, balb_top$full)
write_heatmap_sheet(wb, "C57",  c57_ht,  c57_top$full)

# 汇总表
summary_sheet <- data.frame(
  Strain = c("BALB/c", "C57BL/6"),
  Total_Top_DEGs = c(length(balb_top$genes), length(c57_top$genes)),
  Genes_in_Heatmap = c(nrow(balb_ht$expr_z), nrow(c57_ht$expr_z)),
  Samples = c(ncol(balb_ht$expr_z), ncol(c57_ht$expr_z))
)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_sheet)

saveWorkbook(wb, excel_file, overwrite = TRUE)
cat(sprintf("✓ Excel 已保存: %s\n", excel_file))

# ---- 11. 完成 ----
cat("\n========================================\n")
cat("  热图分析完成！\n")
cat("========================================\n")
cat("\n输出文件:\n")
cat(sprintf("  1. %s\n", pdf_file))
cat(sprintf("  2. %s\n", file.path(out_dir, "TopDEGs_Heatmap_BALBc.pdf")))
cat(sprintf("  3. %s\n", file.path(out_dir, "TopDEGs_Heatmap_C57BL6.pdf")))
cat(sprintf("  4. %s\n", excel_file))
cat(sprintf("\n参数: top%d genes/comparison, |log2FC|>%d, padj<%.2f\n",
            TOP_N_PER_COMPARISON, LOG2FC_THRESHOLD, PADJ_THRESHOLD))
cat("========================================\n")
