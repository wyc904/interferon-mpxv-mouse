# ============================================================================
# common.R — 干扰素小鼠 3-1009 批次分析：公共配置与工具函数
# 用法：在仓库根目录运行 stage 脚本；stage 脚本会 source 本文件
# ============================================================================
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
})

## ---- 路径 ----------------------------------------------------------------
ROOT <- normalizePath(Sys.getenv("MPXV_PROJECT_ROOT", "."), winslash = "/", mustWork = FALSE)
DATA_DIR <- normalizePath(Sys.getenv("MPXV_DATA_DIR", file.path(ROOT, "data", "raw")),
                          winslash = "/", mustWork = FALSE)
META_DIR <- file.path(ROOT, "metadata")
RES_DIR  <- file.path(ROOT, "results", "1009")
CACHE    <- file.path(RES_DIR, "_cache")
LOG_DIR  <- file.path(RES_DIR, "logs")
SUB_DIR  <- file.path(RES_DIR, "00-汇总")

## 计数矩阵：默认在 data/raw/ 下查找；也可用 MPXV_COUNTS_FILE 指定
COUNTS_FILE <- Sys.getenv("MPXV_COUNTS_FILE", "")
if (!nzchar(COUNTS_FILE)) {
  cand <- list.files(DATA_DIR, pattern = "^MPXV_RB_counts_gene.*\\.txt$", full.names = TRUE)
  COUNTS_FILE <- if (length(cand)) cand[1] else NA_character_
}
for (d in c(RES_DIR, SUB_DIR, CACHE, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- 分析参数 --------------------------------------------------------------
PADJ_CUT           <- 0.05   # padj 显著性阈值
LFC_CUT            <- 1.0    # |log2FC| 阈值
MIN_COUNT          <- 10     # 低表达过滤: 某基因 count >= 10 的样本数须 >= 最小分组样本数
TOP_GO_PER_CLUSTER <- 15     # GO 气泡图每个比较组取 Top N
TOP_HEAT_PER_COMP  <- 50     # 热图每个比较组取 Top N
VOLCANO_LABEL_N    <- 10     # 火山图每面板标注 Top N 基因

## ---- 分组定义 --------------------------------------------------------------
GROUP_LEVELS <- c("G1", "G2", "G3", "G4")
GROUP_DESC   <- c(G1 = "uninfected", G2 = "infected", G3 = "treatment 1", G4 = "treatment 2")
GROUP_COLORS <- c(G1 = "#708090", G2 = "#A52A2A", G3 = "#1B9E77", G4 = "#7570B3")
UP_COL   <- "#C75D69"    # 上调 浅红
DOWN_COL <- "#4682B4"    # 下调 钢蓝
NS_COL   <- "#D3D3D3"    # 不显著 浅灰

## ---- 数据集（品系 × 组织）--------------------------------------------------
DATASETS <- list(
  list(id = "BALBc_Lung",   strain = "BALB/c", tissue = "Lung",   dir = "BALBc-肺",  label = "BALB/c — Lung"),
  list(id = "BALBc_Spleen", strain = "BALB/c", tissue = "Spleen", dir = "BALBc-脾",  label = "BALB/c — Spleen"),
  list(id = "AG129_Lung",   strain = "AG129",  tissue = "Lung",   dir = "AG129-肺",  label = "AG129 — Lung"),
  list(id = "AG129_Spleen", strain = "AG129",  tissue = "Spleen", dir = "AG129-脾",  label = "AG129 — Spleen")
)
names(DATASETS) <- sapply(DATASETS, `[[`, "id")

## ---- 6 个比较 --------------------------------------------------------------
COMPARISONS <- list(
  list(id = "G2vsG1", trt = "G2", ref = "G1"),
  list(id = "G3vsG1", trt = "G3", ref = "G1"),
  list(id = "G3vsG2", trt = "G3", ref = "G2"),
  list(id = "G4vsG1", trt = "G4", ref = "G1"),
  list(id = "G4vsG2", trt = "G4", ref = "G2"),
  list(id = "G4vsG3", trt = "G4", ref = "G3")
)
names(COMPARISONS) <- sapply(COMPARISONS, `[[`, "id")
COMP_ORDER  <- names(COMPARISONS)
comp_label  <- function(cid) gsub("vs", " vs ", cid)
COMP_LABELS <- setNames(sapply(COMP_ORDER, comp_label), COMP_ORDER)

## ---- 样本表 ----------------------------------------------------------------
read_sample_sheet <- function() {
  ss <- read.csv(file.path(META_DIR, "sample_sheet.csv"), fileEncoding = "UTF-8-BOM", check.names = FALSE)
  ss$col         <- ss$fc_col
  ss$label_short <- paste0(ss$group, "-", ss$rep)
  ss$strain_c    <- ss$strain
  ss
}

## ---- 计数矩阵 --------------------------------------------------------------
read_counts_matrix <- function() {
  if (is.na(COUNTS_FILE) || !file.exists(COUNTS_FILE)) {
    stop("未找到计数矩阵：请把 MPXV_RB_counts_gene.txt 放入 data/raw/，或设置 MPXV_COUNTS_FILE")
  }
  fp <- COUNTS_FILE
  raw <- read.table(fp, header = TRUE, sep = "\t", skip = 1,
                    comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE,
                    row.names = 1)
  gene_meta <- raw[, 1:5, drop = FALSE]
  colnames(gene_meta) <- c("Chr", "Start", "End", "Strand", "Length")
  cnt <- as.matrix(raw[, -(1:5), drop = FALSE])
  mode(cnt) <- "integer"
  attr(cnt, "gene_meta") <- gene_meta
  cnt
}

## ---- 结果分类 Up/Down/NS ---------------------------------------------------
classify_change <- function(df) {
  df$change <- with(df, ifelse(!is.na(padj) & padj < PADJ_CUT & log2FoldChange >  LFC_CUT, "Up",
                        ifelse(!is.na(padj) & padj < PADJ_CUT & log2FoldChange < -LFC_CUT, "Down", "NS")))
  df$change <- factor(df$change, levels = c("Down", "NS", "Up"))
  df
}

## ---- 保存图（PDF 矢量 + PNG 300dpi）----------------------------------------
save_fig <- function(p, base, w, h) {
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = "pdf", useDingbats = FALSE)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  invisible(base)
}

## ---- 公共主题 --------------------------------------------------------------
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
          legend.position = "bottom")
}

## ---- 日志 ------------------------------------------------------------------
open_log <- function(name) {
  con <- file(file.path(LOG_DIR, paste0(name, ".log")), open = "wt", encoding = "UTF-8")
  assign("..log_con..", con, envir = .GlobalEnv)
  sink(con, split = TRUE)
  sink(con, type = "message")
  invisible(con)
  cat("== log start:", format(Sys.time()), "==\n")
  cat("R:", R.version.string, "| DESeq2:", as.character(packageVersion("DESeq2")), "\n")
}
close_log <- function() {
  cat("== log end:", format(Sys.time()), "==\n")
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  if (exists("..log_con..", envir = .GlobalEnv)) {
    close(get("..log_con..", envir = .GlobalEnv))
    rm("..log_con..", envir = .GlobalEnv)
  }
}

## ---- 缓存路径 --------------------------------------------------------------
cache_path  <- function(id) file.path(CACHE, paste0(id, ".rds"))
ds_res_dir  <- function(ds) file.path(RES_DIR, ds$dir)
