# ============================================================================
# packages.R — 安装本项目所需的 R 包
# 测试环境：R 4.4.0 (2024-04-24) + Bioconductor 3.19
# 用法：在仓库根目录运行 Rscript environment/packages.R
# ============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "DESeq2",          # 1.44.0
  "clusterProfiler", # 4.12.6
  "org.Mm.eg.db",    # 3.19.1
  "enrichplot",      # 1.24.4
  "ggVennDiagram"    # 1.5.7
)

cran_pkgs <- c(
  "tidyverse",       # 元包；02_upset_analysis.R 使用
  "ggplot2",         # 4.0.2
  "dplyr",           # 1.2.0
  "tidyr",           # 1.3.2
  "tibble",          # 3.3.1
  "openxlsx",        # 4.2.8.1
  "writexl",         # 1.5.4
  "pheatmap",        # 1.0.13
  "UpSetR",          # 1.4.0
  "ragg",            # 1.5.2
  "viridis",         # 0.6.5
  "RColorBrewer",    # 1.1.3
  "stringr",         # 1.6.0
  "forcats",         # 1.0.1
  "patchwork",       # 1.3.2
  "cowplot",         # 1.2.0
  "ggrepel"          # 用于火山图/韦恩图标注
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(setdiff(cran_pkgs, rownames(installed.packages())),
                 repos = "https://cloud.r-project.org")

cat("\n安装完成。可用 environment/sessionInfo.txt 对照版本。\n")
