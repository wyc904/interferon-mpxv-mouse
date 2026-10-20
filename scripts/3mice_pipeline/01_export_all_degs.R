# Run from repository root
source("scripts/3mice_pipeline/utils/config.R", encoding = "UTF-8")

# ============================================================================
# 导出所有显著差异基因 (用于挑选基因绘制热图)
# 输出: 完整DEG统计 + VST表达值 + 各比较组Sheet
# ============================================================================

library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)

options(scipen = 999)

PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1
out_dir <- OUT_DIR

# ---- 1. 读取数据 ----
cat("读取数据...\n")
read_fc <- function(fp) {
  raw <- read.table(fp, header = TRUE, sep = "\t", skip = 1,
                    comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  meta  <- raw[, 1:6]; colnames(meta) <- c("Geneid","Chr","Start","End","Strand","Length")
  counts <- raw[, -(1:6), drop = FALSE]
  colnames(counts) <- sapply(colnames(counts), function(nm) {
    sub("/Aligned\\.sortedByCoord\\.out\\.bam$", "", nm)
  })
  list(meta = meta, counts = counts)
}

data_dir <- DATA_DIR_3MICE
balb <- read_fc(file.path(data_dir, "BALB_F_counts_gene(1).txt"))
c57  <- read_fc(file.path(data_dir, "C57_F_counts_gene(1).txt"))

counts_all <- merge(
  data.frame(Geneid = balb$meta$Geneid, balb$counts, check.names = FALSE),
  data.frame(Geneid = c57$meta$Geneid,  c57$counts,  check.names = FALSE),
  by = "Geneid", all = FALSE
)
rownames(counts_all) <- counts_all$Geneid
counts_all$Geneid <- NULL

# ---- 2. 样本注释 ----
parse_sample <- function(colname) {
  parts <- strsplit(colname, "_")[[1]]
  n <- length(parts)
  rep_str  <- parts[n-1]
  time_str <- parts[n-2]
  prefix   <- paste(parts[1:(n-3)], collapse = "_")
  if (prefix == "") prefix <- parts[1]
  strain <- ifelse(grepl("BALB", prefix, ignore.case = TRUE), "BALB", "C57")
  tm <- c("mock"="mock","3d"="3d","5d"="5d","7d"="7d","10d"="10d")
  time_std <- tm[tolower(time_str)]
  if (is.na(time_std)) time_std <- tolower(time_str)
  c(strain = strain, time = time_std, rep = rep_str)
}

colData <- data.frame(sample_id = colnames(counts_all), stringsAsFactors = FALSE)
info <- do.call(rbind, lapply(colnames(counts_all), parse_sample))
colnames(info) <- c("strain","time","rep")
colData <- cbind(colData, as.data.frame(info, stringsAsFactors = FALSE))
colData$strain <- as.factor(colData$strain)
colData$time   <- factor(colData$time, levels = c("mock","3d","5d","7d","10d"))
colData$group  <- with(colData, factor(paste(strain, time, sep = "_")))

cat(sprintf("样本: %d genes × %d samples\n", nrow(counts_all), ncol(counts_all)))

# ---- 3. DESeq2 按品系 ----
run_degs <- function(strain_name, counts_mat, colData) {
  cat(sprintf("\n--- %s DESeq2 ---\n", strain_name))
  ss <- colData[colData$strain == strain_name, ]
  sc <- counts_mat[, ss$sample_id, drop = FALSE]

  smallest_group <- min(table(ss$time))
  keep <- rowSums(sc >= 10) >= smallest_group
  sc <- sc[keep, , drop = FALSE]

  ss$time <- relevel(ss$time, ref = "mock")
  dds <- DESeqDataSetFromMatrix(countData = sc, colData = ss, design = ~ time)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)

  time_levels <- setdiff(levels(ss$time), "mock")
  deg_list <- list()
  for (tp in time_levels) {
    res <- results(dds, contrast = c("time", tp, "mock"), alpha = PADJ_THRESHOLD)
    res_df <- as.data.frame(res) %>%
      rownames_to_column("Geneid") %>%
      mutate(
        comparison = sprintf("%s_%s_vs_mock", strain_name, tp),
        strain = strain_name,
        time   = tp
      )
    deg_list[[tp]] <- res_df
  }
  all_degs <- bind_rows(deg_list)
  list(degs = all_degs, vst_mat = vst_mat, dds = dds)
}

balb_res <- run_degs("BALB", counts_all, colData)
c57_res  <- run_degs("C57",  counts_all, colData)

all_degs <- bind_rows(balb_res$degs, c57_res$degs)

# VST 矩阵按基因合并 (两个品系独立过滤,基因集不同)
vst_balb_df <- as.data.frame(balb_res$vst_mat) %>% rownames_to_column("Geneid")
vst_c57_df  <- as.data.frame(c57_res$vst_mat)  %>% rownames_to_column("Geneid")
vst_full_df <- merge(vst_balb_df, vst_c57_df, by = "Geneid", all = TRUE)
rownames(vst_full_df) <- vst_full_df$Geneid
vst_full_df$Geneid <- NULL
# 用0填充缺失值(某品系被过滤但另一品系保留的基因)
vst_full_df[is.na(vst_full_df)] <- 0

cat(sprintf("\n总DEG条目: %d\n", nrow(all_degs)))

# ---- 4. 筛选显著DEG ----
sig_all <- all_degs %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD, abs(log2FoldChange) > LOG2FC_THRESHOLD) %>%
  mutate(
    direction = ifelse(log2FoldChange > 0, "Up", "Down"),
    negLog10padj = -log10(padj)
  ) %>%
  arrange(comparison, padj)

cat(sprintf("显著DEG条目: %d\n", nrow(sig_all)))
cat(sprintf("去重基因数:  %d\n", length(unique(sig_all$Geneid))))

# 各比较组统计
cat("\n各比较组显著DEG数:\n")
sig_all %>%
  group_by(comparison) %>%
  summarise(Total = n(), Up = sum(direction == "Up"), Down = sum(direction == "Down"),
            .groups = "drop") %>%
  print(n = 20)

# ---- 5. 导出 Excel ----
excel_file <- file.path(out_dir, "All_DEGs_for_Heatmap_Selection.xlsx")
wb <- createWorkbook()

# 5a. 目录/说明 sheet
desc <- data.frame(
  Item = c("分析日期", "比较组数", "显著DEG条目数", "去重基因总数",
           "显著性阈值", "|log2FC|阈值", "归一化方法",
           "Sheet说明",
           "All_Sig_DEGs", "DEGs_by_Comparison", "VST_Expression_AllGenes",
           "VST_Expression_DEGs", "DEG_Summary"),
  Value = c(format(Sys.time(), "%Y-%m-%d %H:%M"),
            length(unique(sig_all$comparison)),
            nrow(sig_all),
            length(unique(sig_all$Geneid)),
            sprintf("padj < %g", PADJ_THRESHOLD),
            sprintf("|log2FC| > %g", LOG2FC_THRESHOLD),
            "DESeq2 VST (variance stabilizing transformation)",
            "",
            "所有显著DEG的完整统计信息(含log2FC, padj, 表达方向等)",
            "按比较组拆分,方便按时间点/品系挑选基因",
            "全部基因的VST归一化表达值(33670 genes × 25 samples)",
            "仅显著DEG的VST表达值(可直接用于热图)",
            "各比较组DEG数量汇总")
)
addWorksheet(wb, "README")
writeData(wb, "README", desc)
setColWidths(wb, "README", cols = 1:2, widths = c(30, 65))

# 5b. 所有显著DEG (一条龙)
sig_export <- sig_all %>%
  select(Geneid, comparison, strain, time, direction,
         baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, negLog10padj) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))
addWorksheet(wb, "All_Sig_DEGs")
writeData(wb, "All_Sig_DEGs", sig_export)
freezePane(wb, "All_Sig_DEGs", firstRow = TRUE)
setColWidths(wb, "All_Sig_DEGs", cols = 1:13, widths = "auto")

# 5c. 按比较组拆分
comps <- unique(sig_all$comparison)
for (comp in comps) {
  comp_data <- sig_export %>%
    filter(comparison == comp) %>%
    arrange(padj)
  sheet_name <- substr(gsub("_vs_mock", "", comp), 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, comp_data)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:13, widths = "auto")
}

# 5d. 全部基因 VST 表达值
vst_all_df <- as.data.frame(vst_full_df) %>%
  rownames_to_column("Geneid") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
addWorksheet(wb, "VST_Expression_AllGenes")
writeData(wb, "VST_Expression_AllGenes", vst_all_df)
freezePane(wb, "VST_Expression_AllGenes", firstRow = TRUE)

# 5e. 仅显著DEG的VST表达值
deg_genes <- unique(sig_all$Geneid)
vst_deg_df <- vst_all_df %>% filter(Geneid %in% deg_genes)
addWorksheet(wb, "VST_Expression_DEGs")
writeData(wb, "VST_Expression_DEGs", vst_deg_df)
freezePane(wb, "VST_Expression_DEGs", firstRow = TRUE)

# 5f. 汇总
summary_df <- sig_all %>%
  group_by(strain, time) %>%
  summarise(
    Total_DEGs = n(),
    Up   = sum(direction == "Up"),
    Down = sum(direction == "Down"),
    .groups = "drop"
  ) %>%
  arrange(strain, factor(time, levels = c("3d","5d","7d","10d")))
addWorksheet(wb, "DEG_Summary")
writeData(wb, "DEG_Summary", summary_df)
setColWidths(wb, "DEG_Summary", cols = 1:5, widths = "auto")

saveWorkbook(wb, excel_file, overwrite = TRUE)
cat(sprintf("\n✓ Excel 已保存: %s\n", excel_file))
cat(sprintf("  包含 %d 个 sheet: README + All_Sig_DEGs + %d 个比较组 + VST表达矩阵 + 汇总\n",
            3 + length(comps), length(comps)))
cat("========================================\n")
