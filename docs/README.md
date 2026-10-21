# 文档索引

建议按以下顺序阅读：

| 顺序 | 文档 | 内容 | 适合谁 |
|---:|---|---|---|
| 1 | [`../README.md`](../README.md) | 项目简介、研究设计、快速开始、目录结构、图件索引 | 所有人 |
| 2 | [`analysis_notes.md`](analysis_notes.md) | 详细分析说明：样本、参数、每步方法、结果速览、注意事项 | 想复现或审稿的人 |
| 3 | [`figure_map.md`](figure_map.md) | 每张图 → 脚本 → 输入/输出 → 预览 | 想找某张图的代码 |
| 4 | [`data_dictionary.md`](data_dictionary.md) | 样本表字段、计数矩阵格式、输出表字段、分析参数 | 要接入数据或改代码的人 |
| 5 | [`../scripts/README.md`](../scripts/README.md) | 代码分类、运行顺序、路径变量 | 要运行代码的人 |
| 6 | [`../legacy/README.md`](../legacy/README.md) | 历史脚本归档说明 + `MANIFEST.tsv` 用法 + 阶段分类 | 要追溯旧版图件的人 |
| 7 | [`../legacy/FAMILIES.md`](../legacy/FAMILIES.md) | 133 个历史脚本的逐条分类索引（编号/文件名/主分类/原路径） | 要快速定位旧脚本的人 |
| 8 | [`../CITATION.cff`](../CITATION.cff) | 机器可读的论文引用信息（作者、期刊、DOI） | 要引用本项目的人 |
| 9 | [`../LICENSE`](../LICENSE) | MIT 许可证 | 要复用代码的人 |

## 图件预览

`docs/figures/` 下是低分辨率预览（不是投稿原图）：

- `upset_DEGs.pdf` — UpSet 图
- `heatmap_topDEGs_BALBc.pdf`、`heatmap_topDEGs_C57BL6.pdf`、`heatmap_topDEGs_C57_BALB.pdf` — Top DEG 热图
- `GO_bubble_up.pdf`、`GO_bubble_down.pdf` — GO 气泡图
- `volcano_panel_c.png` — 火山图网格
- `composite_BALBc_Lung.png`、`composite_BALBc_Spleen.png`、`composite_AG129_Lung.png`、`composite_AG129_Spleen.png` — 四套数据集合成图

## 数据

本仓库**不包含原始数据**。需要的数据文件、放置位置和命名规则见 [`../data/README.md`](../data/README.md)。
