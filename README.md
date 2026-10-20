# MPXV Clade Ib 小鼠感染模型 RNA-seq 分析代码

**Pathogenic Signatures and Antiviral Evaluation of MPXV Clade Ib in Mouse Models**

本仓库包含该研究转录组部分的**可复现分析代码**：DESeq2 差异分析、PCA、DEG 数量柱状图、韦恩图、UpSet 图、火山图、GO 富集气泡图、热图以及投稿用合成图。

> **本仓库不包含原始数据。** 计数矩阵、RDS 缓存、Excel 结果表、论文稿件、参考文献 PDF 均未上传。仓库只保留代码、样本表、分析说明和图件预览。

**研究背景**：猴痘病毒（MPXV）Clade Ib 是近年传播的毒株分支。本项目用 **BALB/c（免疫正常）** 与 **AG129（I 型干扰素受体缺陷）** 小鼠比较 MPXV 感染后的肺、脾转录组反应，并评估两个治疗组（G3/G4）相对感染组的回调效应。分析重点是干扰素/抗病毒通路、DEG 数量与方向、GO 富集和关键基因热图。

---

## 1. 研究设计

- **样本**：66 个 = BALB/c 38 个 + AG129 28 个；每只小鼠取**肺（Lung）**和**脾（Spleen）**。
- **分组**：
  - `G1`：未感染（uninfected）
  - `G2`：感染（infected）
  - `G3` / `G4`：两个治疗组（treatment 1 / treatment 2；具体药物名待实验方确认）
- **数据来源**：featureCounts v2.0.3 基因水平计数（`-p -B -C -t exon -g gene_name`；参考基因组 GRCm39）。
- **分析策略**：按「品系 × 组织」拆成 4 套数据集，每套独立做 DESeq2，再统一出图。

| 数据集 | G1 | G2 | G3 | G4 | 合计 |
|---|---:|---:|---:|---:|---:|
| BALB/c — 肺 | 5 | 5 | 5 | 5 | 20 |
| BALB/c — 脾 | 5 | 5 | 4 | 4 | 18 |
| AG129 — 肺 | 3 | 4 | 4 | 4 | 15 |
| AG129 — 脾 | 3 | 4 | 2 | 4 | 13 |

缺失样本：BALB/c G3 脾-1、G4 脾-5；AG129 G1 肺-1/脾-2、G3 脾-1/脾-3。

**6 个比较组**：`G2vsG1`、`G3vsG1`、`G3vsG2`、`G4vsG1`、`G4vsG2`、`G4vsG3`。

---

## 2. 仓库结构

```text
interferon-mpxv-mouse/
├─ README.md
├─ .gitignore
├─ environment/
│  ├─ packages.R              # 安装 R 包
│  └─ sessionInfo.txt         # 测试环境版本
├─ metadata/
│  ├─ sample_sheet.csv        # 66 个样本的分组/组织/重复信息
│  └─ sample_sheet_notes.txt  # 样本表构建与核对报告
├─ scripts/
│  ├─ README.md              # 代码分类与运行顺序（先读这个）
│  ├─ 1009_pipeline/          # 主流程：4 套数据集 × 6 个比较
│  │  ├─ 01_deseq2.R
│  │  ├─ 02_basic_figures.R
│  │  ├─ 03_go_enrichment.R   # GO Top15（旧版，供合成图取缓存）
│  │  ├─ 03b_go_top30.R       # GO Top5（正式版）
│  │  ├─ 04_heatmap.R         # 热图 Top50（旧版）
│  │  ├─ 04b_heatmap_top50.R  # 热图 Top8（正式版）
│  │  ├─ 05_composite.R
│  │  └─ utils/common.R       # 公共参数、路径、函数
│  ├─ 3mice_pipeline/         # BALB/c vs C57BL/6 图件流程
│  │  ├─ 01_export_all_degs.R
│  │  ├─ 02_upset_analysis.R
│  │  ├─ 03_heatmap_top_genes.R
│  │  ├─ 04_go_enrich_bubble.R
│  │  └─ utils/config.R
│  └─ utilities/              # 样本表构建、源文件检查等 Python 小工具
├─ docs/
│  ├─ README.md              # 文档索引与阅读顺序
│  ├─ analysis_notes.md       # 详细分析说明（推荐先读）
│  ├─ figure_map.md           # 图件 → 脚本 → 输入/输出对照
│  ├─ data_dictionary.md      # 样本表、计数矩阵、输出表字段说明
│  └─ figures/                # 少量低分辨率图件预览（非原始数据）
├─ data/
│  ├─ raw/                    # 放原始计数矩阵（.gitignore 已排除）
│  └─ processed/              # 放手工整理的表（.gitignore 已排除）
├─ results/                   # 运行输出（.gitignore 已排除）
└─ legacy/                    # 133 个历史脚本（按内容去重）+ MANIFEST.tsv + FAMILIES.md
```

### 2.1 代码分类

| 类别 | 目录 | 用途 | 入口 |
|---|---|---|---|
| 主分析流程 | `scripts/1009_pipeline/` | 4 套数据集 × 6 个比较的 DESeq2、PCA、DEG、韦恩、火山、GO、热图、合成图 | `01_deseq2.R` → `05_composite.R` |
| 三品系图件 | `scripts/3mice_pipeline/` | BALB/c vs C57BL/6 的 UpSet、Top DEG 热图、GO 气泡图 | `01_export_all_degs.R` → `04_go_enrich_bubble.R` |
| 工具脚本 | `scripts/utilities/` | 样本表构建、源文件检查、模板提取 | 按需运行 |
| 元数据 | `metadata/` | 66 个样本的分组/组织/重复；构建与核对报告 | 不需要运行 |
| 文档 | `docs/` | 分析说明、图件对照、数据字典、图件预览 | 先读 `docs/README.md` |
| 运行环境 | `environment/` | `packages.R` + `sessionInfo.txt` | `Rscript environment/packages.R` |
| 历史归档 | `legacy/` | 首版投稿/二次调图的 133 个去重脚本 + `MANIFEST.tsv` | 仅追溯，不直接运行 |

> 每个脚本的用途、依赖和运行顺序见 [`scripts/README.md`](scripts/README.md)。

---

## 3. 快速开始

### 3.1 环境

- R 4.4.0（测试环境；Bioconductor 3.19）
- 依赖包及版本见 `environment/sessionInfo.txt`
- 安装：

```bash
Rscript environment/packages.R
```

### 3.2 放入数据

主流程默认从 `data/raw/` 读取计数矩阵，文件名匹配 `MPXV_RB_counts_gene*.txt`：

```text
data/raw/MPXV_RB_counts_gene.txt
```

样本表已包含在 `metadata/sample_sheet.csv`，无需修改。若文件在别处，可设置环境变量：

```bash
# Windows PowerShell
$env:MPXV_COUNTS_FILE = "G:/path/to/MPXV_RB_counts_gene.txt"
$env:MPXV_DATA_DIR   = "G:/path/to/data"
# macOS / Linux
export MPXV_COUNTS_FILE="/path/to/MPXV_RB_counts_gene.txt"
export MPXV_DATA_DIR="/path/to/data"
```

### 3.3 运行主流程（4 套数据集）

**从仓库根目录运行**，脚本用相对路径查找 `common.R`：

```bash
Rscript scripts/1009_pipeline/01_deseq2.R
Rscript scripts/1009_pipeline/02_basic_figures.R
Rscript scripts/1009_pipeline/03_go_enrichment.R   # 旧版 Top15；生成合成图所需的 plots_go_all.rds
Rscript scripts/1009_pipeline/03b_go_top30.R       # 正式版 GO：每组 Top5，每图 ≤30 条
Rscript scripts/1009_pipeline/04b_heatmap_top50.R  # 正式版热图：每组 Top8，每图 ≤48 个基因
Rscript scripts/1009_pipeline/05_composite.R
```

> `03_go_enrichment.R` 与 `03b_go_top30.R` 都会写 `05-GO富集/`，正式版以 **03b** 的输出为准；`03` 还需要为 `05_composite.R` 提供 `_cache/plots_go_all.rds`，所以不要跳过。
>
> `04_heatmap.R` 是每组 Top50 的旧版；正式版用 **04b**。两者输出同名，后运行的会覆盖前者。

输出在 `results/1009/`（结构与原项目 `4-结果/` 一致）。

### 3.4 运行三品系图件流程（BALB/c vs C57BL/6）

数据放到 `data/raw/3mice/`，文件名需与脚本一致：

```text
data/raw/3mice/BALB_F_counts_gene(1).txt
data/raw/3mice/C57_F_counts_gene(1).txt
```

```bash
Rscript scripts/3mice_pipeline/01_export_all_degs.R
Rscript scripts/3mice_pipeline/02_upset_analysis.R
Rscript scripts/3mice_pipeline/03_heatmap_top_genes.R
Rscript scripts/3mice_pipeline/04_go_enrich_bubble.R
```

输出在 `results/3mice/`。`01_export_all_degs.R` 产出的 `All_DEGs_for_Heatmap_Selection.xlsx` 是后续热图和 GO 气泡图的输入。

---

## 4. 主流程说明

| 步骤 | 方法 | 脚本 | 主要输出 |
|---|---|---|---|
| 差异分析 | DESeq2，design `~ group`（G1 为参照）；低表达过滤：count ≥ 10 的样本数 ≥ 最小分组样本数；阈值 `padj < 0.05` 且 `\|log2FC\| > 1` | `01_deseq2.R` | `00-表格/DEG_全部比较_*.xlsx`、`00-汇总/汇总_DEG数量_四套.xlsx` |
| PCA / DEG 数量 / 韦恩 / 火山 | VST 取方差最大 500 基因；6 个比较；火山图标 Top15、y 轴 99.5% 分位截断 | `02_basic_figures.R` | `01-PCA/`、`02-DEG数量/`、`03-韦恩图/`、`04-火山图/` |
| GO 富集（旧版 Top15） | `clusterProfiler::compareCluster + enrichGO`（BP，`org.Mm.eg.db`，BH，`p < 0.05`，`q < 0.2`） | `03_go_enrichment.R` | `05-GO富集/`、`_cache/plots_go_all.rds` |
| GO 富集（正式版 Top5） | 同上；每组取 Top5，每图合计 ≤30 条 | `03b_go_top30.R` | `05-GO富集/GO_Bubble_*.pdf/png`、`00-汇总/GO富集_Top30_汇总.xlsx` |
| 热图（旧版 Top50） | VST → 行 Z-score（±3 截断）→ pheatmap | `04_heatmap.R` | `06-热图/` |
| 热图（正式版 Top8） | 每组 Top8 显著 DEG 并集，每图 ≤48 个基因 | `04b_heatmap_top50.R` | `06-热图/`、`00-表格/热图数据_*.xlsx` |
| 合成图 | 模板式 a–f 组合（PCA / DEG 数量 / 韦恩×2 / 火山网格 / GO×2） | `05_composite.R` | `07-合成图/Composite_*.pdf/png` |

详细说明、参数和结果速览见 [`docs/analysis_notes.md`](docs/analysis_notes.md)。

---

## 5. 图件索引

| 图件族 | 脚本 | 输出 / 预览 |
|---|---|---|
| UpSet：BALB/c + C57BL/6 各时间点 DEG 交集 | `scripts/3mice_pipeline/02_upset_analysis.R` | `results/3mice/Upset_DEGs_Interferon_Mice.pdf`；预览 `docs/figures/upset_DEGs.pdf` |
| Top DEG 热图（BALB/c、C57BL/6、合并） | `scripts/3mice_pipeline/03_heatmap_top_genes.R` | `results/3mice/TopDEGs_Heatmap_*.pdf`；预览 `docs/figures/heatmap_topDEGs_*.pdf` |
| GO 气泡图（BALB/c + C57BL/6，上/下调） | `scripts/3mice_pipeline/04_go_enrich_bubble.R` | `results/3mice/GO/GO_Bubble_*_v2.pdf`；预览 `docs/figures/GO_bubble_*.pdf` |
| 火山图网格（panel c） | `legacy/` 中 `8bc-2次调图.R` 等首版脚本 | 预览 `docs/figures/volcano_panel_c.png` |
| 4 套数据集整套合成图 | `scripts/1009_pipeline/05_composite.R` | `results/1009/*/07-合成图/Composite_*.pdf/png`；预览 `docs/figures/composite_*.png` |

> 论文正文/补充材料的最终图号请以投稿版本为准；本表按「图件族」组织。更细的对照见 [`docs/figure_map.md`](docs/figure_map.md)。

---

## 6. 数据与参数

- **计数矩阵格式**：featureCounts 输出，第 1 列基因名，随后 6 列注释（Geneid、Chr、Start、End、Strand、Length），再接 66 个样本列；`common.R::read_counts_matrix()` 会自动识别前 5 列注释 + 样本列。
- **样本表**：`metadata/sample_sheet.csv`，字段说明见 [`docs/data_dictionary.md`](docs/data_dictionary.md)。
- **显著阈值**：`padj < 0.05`、`|log2FC| > 1`，可在 `scripts/1009_pipeline/utils/common.R` 修改 `PADJ_CUT` / `LFC_CUT`。
- **随机性**：DESeq2 结果确定；PCA/热图基于 VST，无额外随机种子。
- **运行环境**：R 4.4.0、DESeq2 1.44.0、clusterProfiler 4.12.6（完整版本见 `environment/sessionInfo.txt`）。

---

## 7. 注意事项

1. **AG129 G3 脾只有 2 个样本**，统计效力有限；BALB/c 脾也各缺 1 个重复。解读该组结果时请谨慎。
2. **G1–G4 是代号**；G3/G4 的具体药物名尚未写入脚本。拿到药名后，改 `common.R` 里的 `GROUP_DESC` 即可统一替换图例。
3. **GO 正式版为每组 Top5（每图 ≤30 条）**，热图正式版为每组 Top8（每图 ≤48 个基因）；`07-合成图` 里的 GO 面板沿用旧版 Top15，正式 GO 图请看 `05-GO富集/`。
4. **火山图**：`padj` 为 NA 的基因不显示；y 轴按 99.5% 分位截断，个别超高显著性点会被压到顶端。
5. **路径**：主流程从仓库根目录运行；如从其他目录运行，设置 `MPXV_PROJECT_ROOT` 指向仓库根。
6. **原始数据未上传**：请勿把 `data/raw/`、`results/` 下的文件提交到公开仓库；`.gitignore` 已做排除，但提交前仍建议 `git status` 确认。
7. **历史脚本**：`legacy/` 里的 133 个脚本按内容去重，保留了首版投稿和二次调图的中间版本；它们含旧的绝对路径，不能直接运行，仅用于追溯。对照表见 `legacy/MANIFEST.tsv`。

---

## 8. 数据与代码可用性

- **原始数据**：本仓库不包含。原始计数矩阵和测序数据请按论文中的数据可用性声明获取。
- **代码许可**：暂未指定开源许可证；在论文发表前请联系作者确认使用范围。
- **引用**：如果本仓库代码对你有帮助，请引用对应论文：
  *Pathogenic Signatures and Antiviral Evaluation of MPXV Clade Ib in Mouse Models*（作者、期刊信息以正式发表版本为准）。

---

## 9. 联系

问题、复现困难或数据请求，请通过 GitHub Issues 联系仓库维护者。
