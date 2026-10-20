library(ggplot2)

# =========================
# 1) 数据
# =========================
df <- data.frame(
  Comparison = c("Tecovirimat_vs_mock","Tecovirimat_vs_mock",
                 "DMSO_vs_mock","DMSO_vs_mock",
                 "Tecovirimat_vs_DMSO","Tecovirimat_vs_DMSO"),
  Direction  = c("down","up","down","up","down","up"),
  Count      = c(1334,447,495,465,1307,640)
)

# 控制显示顺序：每组里先 down 后 up（更像示例图）
df$Direction <- factor(df$Direction, levels = c("down", "up"))
df$Comparison <- factor(df$Comparison,
                        levels = c("Tecovirimat_vs_mock","DMSO_vs_mock","Tecovirimat_vs_DMSO"))

# =========================
# 2) 颜色（Down 红、Up 蓝）
# =========================
fill_cols <- c("down" = "#F26B6B",  # 红
               "up"   = "#4E79A7")  # 蓝

# =========================
# 3) 作图
# =========================
p <- ggplot(df, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = fill_cols,
    breaks = c("down","up"),
    labels = c("Down-regulated","Up-regulated"),
    name = "Regulation Direction"
  ) +
  labs(
    x = "Comparison",
    y = "Number of Differentially Expressed Genes"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  # 给柱顶文字留空间
  expand_limits(y = max(df$Count) * 1.15)

print(p)

# 可选：保存
# ggsave("DEG_barplot_up_down.pdf", p, width = 7.5, height = 5, useDingbats = FALSE)
# ggsave("DEG_barplot_up_down.png", p, width = 7.5, height = 5, dpi = 300)
