# Run from repository root
source("scripts/3mice_pipeline/utils/config.R", encoding = "UTF-8")

# ============================================================================
# GO 富集分析 + 气泡图: BALB/c & C57BL/6 干扰素小鼠
# 上调 DEGs → 上调 GO plot (红色系)
# 下调 DEGs → 下调 GO plot (蓝色系)
# 参考代码: D:\Desktop\干扰素小鼠\2-挑图\3-go\
# ============================================================================

library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(enrichplot)

options(scipen = 999)

out_dir <- file.path(OUT_DIR, "GO")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. 读取 DEG 数据 ----
cat("========================================\n")
cat("  读取 DEG 数据\n")
cat("========================================\n")

deg_file <- file.path(OUT_DIR, "All_DEGs_for_Heatmap_Selection.xlsx")
deg_data <- read.xlsx(deg_file, sheet = "All_Sig_DEGs")
cat(sprintf("DEG 条目: %d\n", nrow(deg_data)))
cat(sprintf("列名: %s\n", paste(colnames(deg_data), collapse = ", ")))

# 标准化列名
deg_data <- deg_data %>%
  mutate(
    strain = as.character(strain),
    time   = as.character(time),
    direction = as.character(direction),
    comparison = as.character(comparison)
  )

# 查看数据分布
cat("\n各比较组 DEG 数:\n")
deg_data %>%
  group_by(comparison, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  print(n = 20)

# ---- 2. 准备 GO 富集输入 (按品系+时间点+方向) ----
cat("\n========================================\n")
cat("  GO BP 富集分析 (compareCluster)\n")
cat("========================================\n")

# 构建分组列
deg_data <- deg_data %>%
  mutate(
    group = paste(strain, time, sep = "_"),
    dir_group = paste(strain, time, direction, sep = "_")
  )

# ---- 3. 分别对上调和下调做 GO 富集 ----
run_go_enrich <- function(deg_df, direction_filter, label) {
  cat(sprintf("\n--- %s GO 富集 ---\n", label))

  sub <- deg_df %>% filter(direction == direction_filter)
  cat(sprintf("  基因条目: %d\n", nrow(sub)))

  # 按比较组构建基因列表
  gene_list <- split(sub$Geneid, sub$group)
  # 移除空组
  gene_list <- gene_list[sapply(gene_list, length) > 0]
  cat(sprintf("  比较组数: %d\n", length(gene_list)))
  for (nm in names(gene_list)) {
    cat(sprintf("    %s: %d genes\n", nm, length(gene_list[[nm]])))
  }

  # GO BP 富集
  comp <- compareCluster(
    geneClusters = gene_list,
    fun          = "enrichGO",
    OrgDb        = org.Mm.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )

  if (is.null(comp) || nrow(as.data.frame(comp)) == 0) {
    cat("  ⚠ 无显著富集结果!\n")
    return(NULL)
  }

  # 转为数据框
  res_df <- as.data.frame(comp) %>%
    rename(Cluster = Cluster, ID = ID, Description = Description,
           GeneRatio = GeneRatio, BgRatio = BgRatio,
           pvalue = pvalue, p.adjust = p.adjust, qvalue = qvalue,
           geneID = geneID, Count = Count)

  # 解析 Cluster 为 strain + time
  res_df <- res_df %>%
    mutate(
      strain = str_extract(Cluster, "^(BALB|C57)"),
      time   = str_remove(Cluster, "^(BALB|C57)_"),
      mice   = strain  # 兼容参考代码的列名
    )

  # 计算 RichFactor 和 FoldEnrichment
  res_df <- res_df %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      BgRatio_num   = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      RichFactor    = GeneRatio_num,
      FoldEnrichment = GeneRatio_num / BgRatio_num,
      LogP = -log10(p.adjust)
    )

  cat(sprintf("  富集 GO terms: %d\n", nrow(res_df)))
  return(res_df)
}

go_up   <- run_go_enrich(deg_data, "Up",   "上调 (UP)")
go_down <- run_go_enrich(deg_data, "Down", "下调 (DOWN)")

# ---- 4. 气泡图绘制函数 ----
make_go_bubble <- function(go_data, direction_label, color_scheme, n_total = 15) {
  if (is.null(go_data) || nrow(go_data) == 0) {
    cat(sprintf("  ⚠ %s: 无数据可绘图\n", direction_label))
    return(NULL)
  }

  # 每个比较组取 Top N, 确保各品系×各时间点都有代表
  n_per_group <- ceiling(n_total / 4)  # 8组, 每组取 n_per_group 条

  go_sig <- go_data %>% filter(p.adjust < 0.05)

  top_terms <- go_sig %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_per_group) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()

  cat(sprintf("  %s: 每组Top%d → 去重后 %d 个 GO terms\n",
              direction_label, n_per_group, length(top_terms)))
  # 打印每个比较组贡献的term数
  coverage <- go_sig %>%
    filter(Description %in% top_terms) %>%
    group_by(Cluster) %>%
    summarise(n_terms = n(), .groups = "drop")
  print(as.data.frame(coverage))

  # 只保留这些 terms
  plot_data <- go_data %>%
    filter(Description %in% top_terms, p.adjust < 0.05)

  # 设置因子水平
  plot_data$time <- factor(plot_data$time, levels = c("3d", "5d", "7d", "10d"))
  plot_data$strain <- factor(plot_data$strain, levels = c("BALB", "C57"))
  plot_data$mice <- plot_data$strain

  # Y轴排序: 按出现的比较组排列
  plot_data <- plot_data %>%
    mutate(Description = fct_reorder(Description, as.numeric(time), .fun = min, .desc = TRUE))

  # 标题
  title_main <- sprintf("GO BP Enrichment — %s-Regulated DEGs\n(BALB/c & C57BL/6, Top%d per group → %d unique terms)",
                        direction_label, n_per_group, length(unique(plot_data$Description)))

  # 气泡图
  p <- ggplot(plot_data, aes(x = time, y = Description)) +
    geom_point(aes(size = Count, color = LogP), alpha = 0.85) +
    facet_grid(. ~ strain, scales = "free", space = "free") +
    scale_color_gradientn(
      colors = color_scheme,
      name   = expression(-log[10](P.adj)),
      guide  = guide_colorbar(
        title.position = "top",
        barwidth  = unit(5, "cm"),
        barheight = unit(0.5, "cm"),
        direction = "horizontal",
        order     = 1
      )
    ) +
    scale_size_continuous(
      range = c(2.5, 7.5),
      name  = "Gene Count",
      guide = guide_legend(
        title.position = "top",
        direction = "horizontal",
        order     = 2
      )
    ) +
    labs(
      title = title_main,
      x     = "Time Point",
      y     = "GO Terms"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y        = element_text(size = 9),
      strip.text         = element_text(face = "bold", size = 12),
      strip.background   = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position    = "bottom",
      legend.box         = "horizontal",
      legend.box.just    = "center",
      legend.title       = element_text(size = 10, face = "bold"),
      legend.text        = element_text(size = 9),
      legend.key.width   = unit(1, "cm"),
      legend.spacing.x   = unit(0.5, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  # 动态高度
  n_terms <- length(unique(plot_data$Description))
  dynamic_h <- max(7, n_terms * 0.3 + 2)

  list(plot = p, data = plot_data, height = dynamic_h, n_terms = n_terms)
}

# 颜色方案
red_colors  <- c("#FFF0F0", "#FFBBBB", "#FF7777", "#EE3333", "#BB0000", "#880000", "#550000")
blue_colors <- c("#F0F4FF", "#C5D5F0", "#8BAFDE", "#4B80CA", "#2255AA", "#113388", "#081D5C")

cat("\n========================================\n")
cat("  绘制气泡图\n")
cat("========================================\n")

N_TOP_TOTAL <- 15

up_plot   <- make_go_bubble(go_up,   "Up",  red_colors,  n_total = N_TOP_TOTAL)
down_plot <- make_go_bubble(go_down, "Down", blue_colors, n_total = N_TOP_TOTAL)

# ---- 5. 保存矢量 PDF ----
cat("\n--- 保存 PDF ---\n")

# 上调 GO 图
if (!is.null(up_plot)) {
  pdf_up <- file.path(out_dir, "GO_Bubble_Up_Regulated_v2.pdf")
  ggsave(pdf_up, up_plot$plot, width = 10, height = up_plot$height,
         device = "pdf", useDingbats = FALSE)
  cat(sprintf("✓ %s (%.1f × %.1f in)\n", pdf_up, 10, up_plot$height))
}

# 下调 GO 图
if (!is.null(down_plot)) {
  pdf_down <- file.path(out_dir, "GO_Bubble_Down_Regulated_v2.pdf")
  ggsave(pdf_down, down_plot$plot, width = 10, height = down_plot$height,
         device = "pdf", useDingbats = FALSE)
  cat(sprintf("✓ %s (%.1f × %.1f in)\n", pdf_down, 10, down_plot$height))
}

# ---- 6. 导出 Excel ----
cat("\n--- 导出 Excel ---\n")

# 6a. 绘图依据 Excel (仅图中用到的 GO terms)
plot_excel <- file.path(out_dir, "GO_Bubble_Plot_Data_v2.xlsx")
wb_plot <- createWorkbook()

if (!is.null(up_plot)) {
  addWorksheet(wb_plot, "Up_Regulated_PlotData")
  up_export <- up_plot$data %>%
    select(strain, time, ID, Description, GeneRatio, BgRatio,
           RichFactor, FoldEnrichment, pvalue, p.adjust, Count, LogP, geneID) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
    arrange(strain, time, p.adjust)
  writeData(wb_plot, "Up_Regulated_PlotData", up_export)
  freezePane(wb_plot, "Up_Regulated_PlotData", firstRow = TRUE)
  setColWidths(wb_plot, "Up_Regulated_PlotData", cols = 1:13, widths = "auto")
}

if (!is.null(down_plot)) {
  addWorksheet(wb_plot, "Down_Regulated_PlotData")
  down_export <- down_plot$data %>%
    select(strain, time, ID, Description, GeneRatio, BgRatio,
           RichFactor, FoldEnrichment, pvalue, p.adjust, Count, LogP, geneID) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
    arrange(strain, time, p.adjust)
  writeData(wb_plot, "Down_Regulated_PlotData", down_export)
  freezePane(wb_plot, "Down_Regulated_PlotData", firstRow = TRUE)
  setColWidths(wb_plot, "Down_Regulated_PlotData", cols = 1:13, widths = "auto")
}

# 添加说明
desc_df <- data.frame(
  Info = c("绘图数据说明",
           paste0("上调 GO: Top ", N_TOP_TOTAL, " terms (全局最显著)"),
           paste0("下调 GO: Top ", N_TOP_TOTAL, " terms (全局最显著)"),
           "来源: BALB/c & C57BL/6 干扰素小鼠 DEGs",
           "筛选标准: p.adjust < 0.05, |log2FC| > 1",
           "",
           "列说明:",
           "strain = 小鼠品系 (BALB/C57)",
           "time = 时间点 (3d/5d/7d/10d)",
           "Description = GO BP 通路名称",
           "GeneRatio = 该通路中差异基因/总差异基因",
           "RichFactor = 富集因子",
           "FoldEnrichment = 富集倍数",
           "p.adjust = BH校正后的p值",
           "Count = 该通路中的差异基因数",
           "LogP = -log10(p.adjust) 用于气泡图颜色映射")
)
addWorksheet(wb_plot, "README")
writeData(wb_plot, "README", desc_df)
setColWidths(wb_plot, "README", cols = 1:2, widths = c(20, 60))

saveWorkbook(wb_plot, plot_excel, overwrite = TRUE)
cat(sprintf("✓ 绘图数据: %s\n", plot_excel))

# 6b. 完整 GO 分析 Excel (供后续挑选通路)
go_excel <- file.path(out_dir, "GO_Enrichment_BALB_C57_Full_v2.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "README")
readme <- data.frame(
  Info = c("完整 GO 富集分析结果",
           "BALB/c & C57BL/6 干扰素小鼠",
           "",
           "包含所有显著 GO terms (p.adjust < 0.2)",
           "用于后续挑选通路进行定制热图/气泡图绘制",
           "",
           "Sheet 说明:",
           "  GO_Up_All   = 所有上调 DEGs 的 GO BP 富集结果",
           "  GO_Down_All = 所有下调 DEGs 的 GO BP 富集结果",
           "  UP_BALB_3d 等 = 按比较组拆分的详细结果",
           "  Summary     = 汇总统计")
)
writeData(wb, "README", readme)
setColWidths(wb, "README", cols = 1:2, widths = c(25, 65))

# UP GO results
if (!is.null(go_up)) {
  addWorksheet(wb, "GO_Up_All")
  writeData(wb, "GO_Up_All", go_up %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
  freezePane(wb, "GO_Up_All", firstRow = TRUE)

  for (grp in unique(go_up$Cluster)) {
    grp_data <- go_up %>% filter(Cluster == grp) %>% arrange(p.adjust)
    if (nrow(grp_data) > 0) {
      sn <- substr(paste0("UP_", grp), 1, 31)
      addWorksheet(wb, sn)
      writeData(wb, sn, grp_data %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
    }
  }
}

# DOWN GO results
if (!is.null(go_down)) {
  addWorksheet(wb, "GO_Down_All")
  writeData(wb, "GO_Down_All", go_down %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
  freezePane(wb, "GO_Down_All", firstRow = TRUE)

  for (grp in unique(go_down$Cluster)) {
    grp_data <- go_down %>% filter(Cluster == grp) %>% arrange(p.adjust)
    if (nrow(grp_data) > 0) {
      sn <- substr(paste0("DOWN_", grp), 1, 31)
      addWorksheet(wb, sn)
      writeData(wb, sn, grp_data %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
    }
  }
}

# 汇总
summary_df <- bind_rows(
  if (!is.null(go_up)) go_up %>% mutate(Direction = "Up") else NULL,
  if (!is.null(go_down)) go_down %>% mutate(Direction = "Down") else NULL
) %>%
  group_by(Direction, strain, time) %>%
  summarise(GO_Terms = n(), .groups = "drop")

if (nrow(summary_df) > 0) {
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_df)
}

saveWorkbook(wb, go_excel, overwrite = TRUE)
cat(sprintf("✓ 完整GO结果: %s\n", go_excel))

# ---- 7. 完成 ----
cat("\n========================================\n")
cat("  GO 富集分析 + 气泡图完成!\n")
cat("========================================\n")
cat("\n输出文件:\n")
if (!is.null(up_plot)) cat(sprintf("  1. GO_Bubble_Up_Regulated_v2.pdf (矢量, %d terms)\n", up_plot$n_terms))
if (!is.null(down_plot)) cat(sprintf("  2. GO_Bubble_Down_Regulated_v2.pdf (矢量, %d terms)\n", down_plot$n_terms))
cat(sprintf("  3. GO_Bubble_Plot_Data_v2.xlsx — 绘图依据数据\n"))
cat(sprintf("  4. GO_Enrichment_BALB_C57_Full_v2.xlsx — 完整GO结果(供挑选通路)\n"))
cat("========================================\n")
