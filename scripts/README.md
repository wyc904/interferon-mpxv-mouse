# scripts/ 代码分类与运行顺序

所有脚本都从**仓库根目录**运行。数据路径由 `utils/common.R` 或 `utils/config.R` 统一管理；默认从 `data/raw/` 读取，也可用环境变量覆盖。

---

## 1. `1009_pipeline/` — 主流程（4 套数据集 × 6 个比较）

适用数据：66 个样本（BALB/c 38 + AG129 28，肺 + 脾），按「品系 × 组织」拆成 4 套独立 DESeq2 分析。

| 顺序 | 脚本 | 作用 | 依赖 | 主要输出 |
|---:|---|---|---|---|
| 1 | `01_deseq2.R` | DESeq2 差异分析、VST 归一化、6 个比较 | 原始计数矩阵 + `metadata/sample_sheet.csv` | `results/1009/<dataset>/00-表格/DEG_全部比较_*.xlsx`、`_cache/<id>.rds` |
| 2 | `02_basic_figures.R` | PCA、DEG 数量柱状图、韦恩图、火山图 | 步骤 1 的缓存 | `01-PCA/`、`02-DEG数量/`、`03-韦恩图/`、`04-火山图/`、`_cache/plots_<id>.rds` |
| 3 | `03_go_enrichment.R` | GO 富集（旧版 Top15），并为合成图生成缓存 | 步骤 1 的缓存 | `05-GO富集/`、`_cache/plots_go_all.rds` |
| 4 | `03b_go_top30.R` | GO 富集**正式版**：每组 Top5、每图 ≤30 条 | 步骤 1 的缓存 | `05-GO富集/`、`00-汇总/GO富集_Top30_汇总.xlsx` |
| 5 | `04_heatmap.R` | 热图旧版：每组 Top50 | 步骤 1 的缓存 | `06-热图/` |
| 6 | `04b_heatmap_top50.R` | 热图**正式版**：每组 Top8、每图 ≤48 个基因 | 步骤 1 的缓存 | `06-热图/`、`00-表格/热图数据_*.xlsx` |
| 7 | `05_composite.R` | 模板式 a–f 合成图 | 步骤 2 的 `plots_<id>.rds` + 步骤 3 的 `plots_go_all.rds` | `07-合成图/Composite_*.pdf/png` |

> **顺序不能乱**：`05_composite.R` 依赖 `02` 和 `03` 的缓存；`03b` 会覆盖 `03` 的正式 GO 图，所以 `03b` 放在 `03` 之后；`04b` 放在 `04` 之后。
>
> 辅助文件：`utils/common.R`（路径、阈值、分组、数据集、配色、日志、保存函数）。

---

## 2. `3mice_pipeline/` — 三品系图件（BALB/c vs C57BL/6）

适用数据：`data/raw/3mice/` 下的 `BALB_F_counts_gene(1).txt` 和 `C57_F_counts_gene(1).txt`。

| 顺序 | 脚本 | 作用 | 输出 |
|---:|---|---|---|
| 1 | `01_export_all_degs.R` | 导出所有显著 DEG + VST 表达，供后续挑选 | `results/3mice/All_DEGs_for_Heatmap_Selection.xlsx` |
| 2 | `02_upset_analysis.R` | UpSet 图 + 交集基因列表 | `Upset_DEGs_Interferon_Mice.pdf`、`DEGs_Upset_Results.xlsx` |
| 3 | `03_heatmap_top_genes.R` | Top DEG 热图（BALB/c、C57BL/6、合并） | `TopDEGs_Heatmap_*.pdf`、`TopDEGs_Heatmap_Data.xlsx` |
| 4 | `04_go_enrich_bubble.R` | GO 富集 + 上/下调气泡图 | `GO/GO_Bubble_*_v2.pdf`、`GO_Bubble_Plot_Data_v2.xlsx`、`GO_Enrichment_BALB_C57_Full_v2.xlsx` |

`01` 是 `02`–`04` 的共同输入，必须先跑。辅助文件：`utils/config.R`。

---

## 3. `utilities/` — Python 辅助工具

这些脚本不直接出图，用于准备和核对元数据：

| 脚本 | 用途 |
|---|---|
| `build_sample_sheet.py` | 从组织样本对照表 + 计数矩阵表头构建 `sample_sheet.csv`，并做 66/66 一致性核对 |
| `inspect_source_meta.py` | 检查源数据表头、分组、缺失样本 |
| `extract_template_docx.py` | 从分析模板 docx 中提取文字/图片，供参考 |
| `_heat_summary.py` | 汇总热图选基因信息（内部辅助） |

---

## 4. 路径与数据约定

| 变量 | 默认值 | 说明 |
|---|---|---|
| `MPXV_PROJECT_ROOT` | 当前目录 | 仓库根目录；从其他目录运行时设置 |
| `MPXV_DATA_DIR` | `data/raw/` | 主流程计数矩阵目录 |
| `MPXV_COUNTS_FILE` | 自动查找 `MPXV_RB_counts_gene*.txt` | 直接指定计数矩阵文件 |
| `MPXV_DATA_DIR_3MICE` | `data/raw/3mice/` | 三品系计数矩阵目录 |
| `MPXV_OUT_3MICE` | `results/3mice/` | 三品系输出目录 |

**原始数据、缓存和结果都不会提交到 Git**：`data/raw/`、`data/processed/`、`results/` 均被 `.gitignore` 排除。

---

## 5. 与 `legacy/` 的区别

- `scripts/` = 经过路径重构、参数集中、文档化的**可复现流程**。
- `legacy/` = 首版投稿/二次调图/探索阶段的**历史脚本**，含旧的绝对路径和手工步骤，仅用于追溯，不能直接运行。
- 两者不重复：已收录到 `scripts/` 的源版本没有再放进 `legacy/`。

---

## 6. 快速语法检查

```bash
Rscript -e 'for (f in list.files("scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE)) { parse(f); cat("OK", f, "\n") }'
```
