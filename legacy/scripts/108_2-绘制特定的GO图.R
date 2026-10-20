#=============================
# 共分析4个特定GO通路 & 可视化
#=============================

library(tidyverse)  # 已经在前面library(tidyverse)就不用重复了，这里为了自洽

# 1. 设置GO结果所在文件夹
go_dir <- "D:/Desktop/LRT_TimeCourse_Analysis"

# 2. 定义感兴趣的GO term（必须与Description列中的文字完全一致）
terms_of_interest <- c(
  "cellular response to interferon-beta",
  "response to interferon-beta",
  "regulation of innate immune response",
  "double-strand break repair"
)

# 3. 读取文件夹下所有GO结果csv文件
go_files <- list.files(go_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(go_files) == 0) {
  stop("在目录中没有找到任何csv文件，请检查路径或文件名。")
}

# 4. 从每个文件中筛选出这4个GO term，并合并
go_selected <- purrr::map_dfr(go_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  if (!"Description" %in% colnames(df)) {
    warning("文件没有 Description 列，已跳过: ", f)
    return(NULL)
  }
  
  df_sub <- df %>%
    dplyr::filter(Description %in% terms_of_interest) %>%
    dplyr::mutate(
      Source = tools::file_path_sans_ext(basename(f))  # 用文件名标记来源，如 Cluster1_xxx
    )
  
  df_sub
})

if (nrow(go_selected) == 0) {
  stop("在任何一个GO结果文件中都没有找到指定的4个GO通路，请检查Description的拼写是否一致。")
}

# 5. 整理因子顺序：Y轴按你指定的顺序，X轴按文件名顺序
go_selected$Description <- factor(go_selected$Description,
                                  levels = terms_of_interest)

go_selected$Source <- factor(go_selected$Source,
                             levels = sort(unique(go_selected$Source)))

# 6. 如果没有Count列，但有geneID列，则根据geneID拆分计算Count
if (!"Count" %in% colnames(go_selected) && "geneID" %in% colnames(go_selected)) {
  go_selected <- go_selected %>%
    dplyr::mutate(
      Count = sapply(strsplit(geneID, "/"), length)
    )
}

# 7. 计算 -log10(p.adjust)，若没有p.adjust但有pvalue，就暂时用pvalue代替
if (!"p.adjust" %in% colnames(go_selected) && "pvalue" %in% colnames(go_selected)) {
  message("未找到 p.adjust 列，将临时使用 pvalue 代替 p.adjust。")
  go_selected$p.adjust <- go_selected$pvalue
}

if (!"p.adjust" %in% colnames(go_selected)) {
  stop("既没有 p.adjust 也没有 pvalue，无法做显著性可视化。")
}

go_selected <- go_selected %>%
  dplyr::mutate(
    negLog10Padj = -log10(p.adjust)
  )

#=============================
# 图1：气泡图（不同分析中4个GO通路的显著性对比）
#=============================
p_bubble <- ggplot(go_selected,
                   aes(x = Source, y = Description)) +
  geom_point(aes(size = Count, color = negLog10Padj)) +
  scale_color_gradient(
    low  = "lightblue",
    high = "red",
    name = "-log10(p.adjust)"
  ) +
  scale_size_continuous(name = "基因数(Count)") +
  theme_bw() +
  labs(
    title = "4个关键免疫相关GO通路在不同分析中的显著性对比",
    x = "分析组别 / 文件（Source）",
    y = "GO Biological Process"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存气泡图
ggsave(file.path(go_dir, "Key4_GO_terms_BubblePlot.pdf"),
       p_bubble, width = 9, height = 5)

#=============================
# 图2：热图风格（-log10(p.adjust)）
#=============================
p_heat <- ggplot(go_selected,
                 aes(x = Source, y = Description, fill = negLog10Padj)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low  = "white",
    high = "red",
    name = "-log10(p.adjust)"
  ) +
  # 在每个格子里标注p.adjust，大概三位小数
  geom_text(aes(label = signif(p.adjust, 3)), size = 3) +
  theme_bw() +
  labs(
    title = "4个关键免疫相关GO通路的显著性热图",
    x = "分析组别 / 文件（Source）",
    y = "GO Biological Process"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(go_dir, "Key4_GO_terms_Heatmap.pdf"),
       p_heat, width = 9, height = 5)

message("✅ 4个GO通路的共分析已完成：\n",
        "  - 气泡图: ", file.path(go_dir, "Key4_GO_terms_BubblePlot.pdf"), "\n",
        "  - 热图:   ", file.path(go_dir, "Key4_GO_terms_Heatmap.pdf"))
