## =============================
## 0) Packages
## =============================
library(tidyverse)
library(openxlsx)
library(ComplexUpset)
library(readr)     # 用 parse_number 更稳
library(svglite)

## =============================
## 1) Read data
## =============================
df_raw <- read.xlsx("绘制.xlsx", sheet = 1)

# 统一列类型：把需要计算的列转成 numeric
num_cols <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

df <- df_raw %>%
  mutate(
    time = as.character(time),
    mice = as.character(mice),
    gene = as.character(gene)
  ) %>%
  mutate(
    across(all_of(num_cols), ~ parse_number(as.character(.), locale = locale(decimal_mark = ".")))
  )

# 快速检查：哪些行的 log2FoldChange 仍然是 NA（说明原始单元格不是数字）
cat("log2FoldChange 转换后 NA 数:", sum(is.na(df$log2FoldChange)), "\n")
cat("padj 转换后 NA 数:", sum(is.na(df$padj)), "\n")

## =============================
## 2) Significant DEGs (Up+Down)
## =============================
padj_threshold   <- 0.05
log2fc_threshold <- 1

sig_df <- df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  filter(padj < padj_threshold, abs(log2FoldChange) > log2fc_threshold) %>%
  distinct(time, mice, gene, .keep_all = TRUE)

cat("显著差异基因总数:", nrow(sig_df), "\n")

## =============================
## 3) Define sets
## =============================
time_points_want   <- c("3d", "5d", "7d", "10d")
mouse_strains_want <- c("c57", "8bc")

time_points  <- time_points_want[time_points_want %in% unique(sig_df$time)]
mouse_strains <- mouse_strains_want[mouse_strains_want %in% unique(sig_df$mice)]

if (length(time_points) == 0 || length(mouse_strains) == 0) {
  stop("数据里找不到指定 time/mice；请检查 time/mice 字段取值是否为 3d/5d/7d/10d 与 c57/8bc。")
}

sets_order <- as.vector(outer(mouse_strains, time_points, paste, sep = "_"))
cat("绘图集合顺序:", paste(sets_order, collapse = ", "), "\n")

## =============================
## 4) Build membership matrix (Gene × sets)
## =============================
all_genes <- unique(sig_df$gene)
mat <- tibble(Gene = all_genes)

for (s in sets_order) {
  parts <- str_split(s, "_", simplify = TRUE)
  m0 <- parts[1]; t0 <- parts[2]
  genes_s <- sig_df %>% filter(mice == m0, time == t0) %>% pull(gene) %>% unique()
  mat[[s]] <- mat$Gene %in% genes_s  # TRUE/FALSE
}

pretty_set_names <- function(x) {
  x %>%
    str_replace("^c57_", "C57-") %>%
    str_replace("^8bc_", "8BC-") %>%
    str_replace("_", "-")
}

## =============================
## 5) Plot (Cell/Nature-like)
## =============================
p_pub <- ComplexUpset::upset(
  data = mat,
  intersect = sets_order,
  name = "DEGs",
  min_size = 5,
  n_intersections = 30,
  sort_sets = FALSE,
  sort_intersections = "descending",
  sort_intersections_by = c("degree", "cardinality"),
  width_ratio  = 0.23,
  height_ratio = 0.60,
  matrix = intersection_matrix(
    geom = geom_point(size = 2.0),
    segment = geom_segment(linewidth = 0.35),
    outline_color = list(active = "black", inactive = "grey80")
  ),
  set_sizes = upset_set_size(geom = geom_bar(width = 0.65)),
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts = TRUE,
      text = list(size = 2.6, vjust = -0.25)
    ) + ylab("Intersection size (genes)")
  ),
  labeller = pretty_set_names
) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) +
  ggtitle("Differentially expressed genes across strains and timepoints")

ggsave("UpSet_Total_DEGs_All_Timepoints_pub.pdf", p_pub,
       width = 180, height = 130, units = "mm", device = grDevices::cairo_pdf)

ggsave("UpSet_Total_DEGs_All_Timepoints_pub.svg", p_pub,
       width = 180, height = 130, units = "mm", device = svglite::svglite)

cat("Done: PDF + SVG 已输出。\n")
