# ===============================
# 0. 加载必要包
# ===============================

library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(grid)

# ===============================
# 1. 读取GO数据
# ===============================

go_data <- read.csv(
  "D:/Desktop/compareCluster_Down_GO_BP.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ===============================
# 2. 数值列转换
# ===============================

numeric_columns <- c(
  "RichFactor",
  "p.adjust",
  "Count",
  "FoldEnrichment"
)

for(col in numeric_columns){
  
  if(col %in% colnames(go_data)){
    
    go_data[[col]] <- as.numeric(
      as.character(go_data[[col]])
    )
  }
}

# ===============================
# 3. 数据过滤
# ===============================

filtered_go_data <- go_data %>%
  
  filter(
    if_all(
      all_of(numeric_columns),
      ~ !is.na(.)
    )
  ) %>%
  
  filter(p.adjust < 0.05)

# 计算-log10(p.adjust)

filtered_go_data$LogP <- -log10(
  filtered_go_data$p.adjust
)

# ===============================
# 4. 设置比较组顺序
# ===============================

expected_comparisons <- c(
  "DMSO_vs_mock",
  "Tecovirimat_vs_mock",
  "Tecovirimat_vs_DMSO"
)

available_comparisons <- unique(
  filtered_go_data$Cluster
)

available_comparisons <- available_comparisons[
  !is.na(available_comparisons)
]

actual_comparisons <- intersect(
  expected_comparisons,
  available_comparisons
)

if(length(actual_comparisons) == 0){
  
  actual_comparisons <- available_comparisons
  
  cat("警告：未找到预期比较组\n")
}

filtered_go_data$Cluster <- factor(
  filtered_go_data$Cluster,
  levels = actual_comparisons
)

filtered_go_data <- filtered_go_data %>%
  filter(!is.na(Cluster))

# ===============================
# 5. 提取Top15
# ===============================

top_n_per_cluster <- 15

top_terms_list <- filtered_go_data %>%
  
  group_by(Cluster) %>%
  
  arrange(p.adjust) %>%
  
  slice_head(n = top_n_per_cluster) %>%
  
  ungroup() %>%
  
  pull(Description) %>%
  
  unique()

plot_data <- filtered_go_data %>%
  
  filter(
    Description %in% top_terms_list
  )

cat(
  "共筛选GO通路数量:",
  length(top_terms_list),
  "\n"
)

# ===============================
# 6. 自定义GO term顺序
# ===============================

custom_order <- c(
  
  "wound healing",
  
  "response to peptide",
  
  "regulation of T cell activation",
  
  "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction",
  
  "positive regulation of antigen processing and presentation",
  
  "peptidyl-serine modification",
  
  "negative regulation of supramolecular fiber organization",
  
  "negative regulation of organelle organization",
  
  "memory",
  
  "exploration behavior",
  
  "exocytosis",
  
  "cellular response to peptide hormone stimulus",
  
  "cellular response to peptide",
  
  "cell maturation",
  
  "blood coagulation",
  
  "antigen processing and presentation",
  
  "actin filament organization",
  
  "lymphocyte mediated immunity",
  
  "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
  
  "nucleoside monophosphate metabolic process",
  
  "deoxyribose phosphate metabolic process",
  
  "deoxyribonucleotide metabolic process",
  
  "deoxyribonucleoside monophosphate metabolic process",
  
  "adult locomotory behavior",
  
  "2'-deoxyribonucleotide metabolic process"
  
)

# 仅保留存在于数据中的term

custom_order <- custom_order[
  custom_order %in% plot_data$Description
]

# 设置Y轴顺序

plot_data$Description <- factor(
  plot_data$Description,
  levels = rev(custom_order)
)

# ===============================
# 7. 颜色设置
# ===============================

blue_colors <- c(
  "#EAF3FA",
  "#BCD7EA",
  "#86B4D4",
  "#315F8A",
  "#1D3F5F"
)

# ===============================
# 8. X轴标签
# ===============================

label_mapping <- c(
  "DMSO_vs_mock" = "DMSO_vs_mock",
  "Tecovirimat_vs_mock" = "Tecovirimat_vs_mock",
  "Tecovirimat_vs_DMSO" = "Tecovirimat_vs_DMSO"
)

actual_labels <- label_mapping[
  levels(plot_data$Cluster)
]

actual_labels[is.na(actual_labels)] <-
  names(actual_labels[is.na(actual_labels)])

# ===============================
# 9. 绘图
# ===============================

p_basic <- ggplot(
  plot_data,
  aes(
    x = Cluster,
    y = Description
  )
) +
  
  geom_point(
    aes(
      size = Count,
      color = LogP
    ),
    alpha = 0.85
  ) +
  
  scale_color_gradientn(
    
    colors = blue_colors,
    
    name = "-Log10(P.adj)",
    
    guide = guide_colorbar(
      
      title.position = "top",
      
      title.hjust = 0.5,
      
      barwidth = unit(5, "cm"),
      
      barheight = unit(0.5, "cm"),
      
      direction = "horizontal",
      
      order = 1
    )
  ) +
  
  scale_size_continuous(
    
    range = c(3, 8),
    
    name = "Gene Count",
    
    guide = guide_legend(
      
      title.position = "top",
      
      title.hjust = 0.5,
      
      direction = "horizontal",
      
      nrow = 1,
      
      label.position = "bottom",
      
      order = 2,
      
      override.aes = list(alpha = 1)
    )
  ) +
  
  scale_x_discrete(
    labels = actual_labels
  ) +
  
  labs(
    
    title = paste(
      "Top",
      top_n_per_cluster,
      "Enriched GO Terms per Group"
    ),
    
    x = NULL,
    
    y = NULL
  ) +
  
  theme_minimal() +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    # 竖直参考线
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    
    # 水平参考线
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    
    legend.box = "horizontal",
    
    legend.box.just = "center",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.width = unit(1, "cm"),
    
    legend.margin = margin(
      t = 5,
      r = 0,
      b = 0,
      l = 0
    ),
    
    legend.spacing.x = unit(0.5, "cm"),
    
    legend.box.margin = margin(
      t = 10,
      r = 0,
      b = 5,
      l = 0
    ),
    
    plot.margin = unit(
      c(1,1,1.5,1),
      "cm"
    )
  )

# ===============================
# 10. 显示图形
# ===============================

print(p_basic)

# ===============================
# 11. 保存PDF
# ===============================

output_dir <- "D:/Desktop/1"

if(!dir.exists(output_dir)){
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

output_pdf <- file.path(
  output_dir,
  "GO_bubble_plot_down_blue.pdf"
)

ggsave(
  
  filename = output_pdf,
  
  plot = p_basic,
  
  width = 16,
  
  height = 10,
  
  units = "in",
  
  device = "pdf",
  
  useDingbats = FALSE
)
