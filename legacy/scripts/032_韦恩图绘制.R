suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
})

# ==============================================================================
# 7. 用 ggVennDiagram 绘制：Tecovirimat/DMSO 与 mock 的 Up/Down 韦恩图
# ==============================================================================

# 7.1 提取两个对比的 Up/Down 基因集合（你已有 get_deg_lists）
deg_tec_vs_mock  <- get_deg_lists(dds, "Tecovirimat", "mock", FC_cutoff = FC_cutoff, P_cutoff = P_cutoff)
deg_dmso_vs_mock <- get_deg_lists(dds, "DMSO",        "mock", FC_cutoff = FC_cutoff, P_cutoff = P_cutoff)

up_sets_2 <- list(
  Tecovirimat_vs_mock = deg_tec_vs_mock$up,
  DMSO_vs_mock        = deg_dmso_vs_mock$up
)

down_sets_2 <- list(
  Tecovirimat_vs_mock = deg_tec_vs_mock$down,
  DMSO_vs_mock        = deg_dmso_vs_mock$down
)

cat("\n>>> Tecovirimat vs mock: Up=", length(deg_tec_vs_mock$up),
    " Down=", length(deg_tec_vs_mock$down), "\n")
cat(">>> DMSO vs mock:       Up=", length(deg_dmso_vs_mock$up),
    " Down=", length(deg_dmso_vs_mock$down), "\n")

# 7.2（可选）导出交集/特异基因列表
up_intersect   <- intersect(up_sets_2$Tecovirimat_vs_mock, up_sets_2$DMSO_vs_mock)
up_tec_only    <- setdiff(up_sets_2$Tecovirimat_vs_mock, up_sets_2$DMSO_vs_mock)
up_dmso_only   <- setdiff(up_sets_2$DMSO_vs_mock, up_sets_2$Tecovirimat_vs_mock)

down_intersect <- intersect(down_sets_2$Tecovirimat_vs_mock, down_sets_2$DMSO_vs_mock)
down_tec_only  <- setdiff(down_sets_2$Tecovirimat_vs_mock, down_sets_2$DMSO_vs_mock)
down_dmso_only <- setdiff(down_sets_2$DMSO_vs_mock, down_sets_2$Tecovirimat_vs_mock)

write.table(up_intersect,   "Venn_Up_intersect_Tecovirimat_DMSO_vs_mock.txt",   quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(up_tec_only,    "Venn_Up_Tecovirimat_only_vs_mock.txt",             quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(up_dmso_only,   "Venn_Up_DMSO_only_vs_mock.txt",                    quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(down_intersect, "Venn_Down_intersect_Tecovirimat_DMSO_vs_mock.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(down_tec_only,  "Venn_Down_Tecovirimat_only_vs_mock.txt",           quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(down_dmso_only, "Venn_Down_DMSO_only_vs_mock.txt",                  quote=FALSE, row.names=FALSE, col.names=FALSE)

cat(">>> Up 交集:", length(up_intersect), " | Tec-only:", length(up_tec_only), " | DMSO-only:", length(up_dmso_only), "\n")
cat(">>> Down 交集:", length(down_intersect), " | Tec-only:", length(down_tec_only), " | DMSO-only:", length(down_dmso_only), "\n")

# 7.3 颜色（沿用你示例的三色风格，这里两组就取蓝+橙；Tec 用棕色也可以）
group_colors_two <- c("#8B4513", "#708090")  # Tec 深棕，DMSO 灰（你也可换成蓝橙）

# 7.4 上调 Venn（红色系渐变填充）
p_up <- ggVennDiagram(
  up_sets_2,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#FDEDEC", high = "#C0392B") +  # 红色系（Up）
  scale_color_manual(values = group_colors_two) +
  labs(title = "", subtitle = "") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, "mm")
  ) +
  coord_flip()

# 7.5 下调 Venn（蓝色系渐变填充）
p_down <- ggVennDiagram(
  down_sets_2,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#EBF5FB", high = "#21618C") +  # 蓝色系（Down）
  scale_color_manual(values = group_colors_two) +
  labs(title = "", subtitle = "") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, "mm")
  ) +
  coord_flip()

print(p_up)
print(p_down)

# 7.6 保存（可选）
ggsave("Venn_Up_Tecovirimat_DMSO_vs_mock_ggVennDiagram.pdf",   p_up,   width = 7, height = 6)
ggsave("Venn_Down_Tecovirimat_DMSO_vs_mock_ggVennDiagram.pdf", p_down, width = 7, height = 6)

cat("\n>>> 已输出：Venn_Up_*_ggVennDiagram.pdf 与 Venn_Down_*_ggVennDiagram.pdf\n")
