# 图件 → 脚本 → 输入/输出对照

本表先给出已发表论文的图号对应，再按「图件族」列出仓库内的完整流程。

---

## 0. 已发表论文图号 → 脚本

论文：Han H, Wang X, Wu Y, et al. *Pathogenic signatures and therapeutic evaluation of emergent MPXV Clade Ib in low-susceptibility and immunocompromised mouse models.* **BMC Microbiology.** 2026;26:723. doi:[10.1186/s12866-026-05273-4](https://doi.org/10.1186/s12866-026-05273-4)

| 论文图 | 内容 | 对应脚本 | 说明 |
|---|---|---|---|
| Fig. 1 | 四种小鼠的易感性、体重/生存、病毒载量、IHC/H&E | 无 | 湿实验，不涉及转录组代码 |
| Fig. 2a | 采样示意图 | 无 | 示意图 |
| Fig. 2b | UpSet：C57BL/6 vs BALB/c，3/5/7/10 dpi | `scripts/3mice_pipeline/02_upset_analysis.R` | 输出 `Upset_DEGs_Interferon_Mice.pdf` |
| Fig. 2c | 火山图 2×4（两品系 × 四时间点） | `legacy/scripts/015_8bc-2次调图.R` | 交互式脚本，需预先载入 `count_matrix` 和 `sample_info`；只生成单张火山图，2×4 拼图未在仓库中自动化 |
| Fig. 2d | Top DEG 热图 | `scripts/3mice_pipeline/03_heatmap_top_genes.R` | 输出 `TopDEGs_Heatmap_C57_BALB.pdf` |
| Fig. 2e | GO 气泡图（上/下调） | `scripts/3mice_pipeline/04_go_enrich_bubble.R` | 输出 `GO_Bubble_Up/Down_Regulated_v2.pdf`；具体通路条数随版本/参数变化 |
| Fig. 3 | MPXV Ib vs IIb 在 STAT1−/− 小鼠中的致病性 | 无 | 湿实验 |
| Fig. 4 | Tecovirimat 药效（C57BL/6、AGB6、STAT1−/−） | 无 | 湿实验；3-1009 的 G3/G4 是另一批 BALB/c + AG129 数据，不对应 Fig. 4 |

> **3-1009 批次的位置**：该批次（BALB/c + AG129，肺 + 脾，G1–G4）是独立的 RNA-seq 实验，**未包含在 BMC Microbiology 2026 论文的 Fig. 1–4 中**；图号待定（如用于后续论文，请在此处补充）。下面的第 1 节是它的完整流程。

---

## 1. 主流程：4 套数据集（BALB/c / AG129 × 肺 / 脾）

运行目录：仓库根目录。数据集 ID：`BALBc_Lung`、`BALBc_Spleen`、`AG129_Lung`、`AG129_Spleen`。

| 图件 | 脚本 | 输入 | 输出（`results/1009/`） |
|---|---|---|---|
| PCA | `scripts/1009_pipeline/02_basic_figures.R` | `_cache/<id>.rds`（`01_deseq2.R` 生成） | `<dataset>/01-PCA/PCA_<id>.pdf/.png` |
| DEG 数量柱状图 | 同上 | 同上 | `<dataset>/02-DEG数量/DEGcount_<id>.pdf/.png` |
| 韦恩图（每治疗组 All/Up/Down） | 同上 | 同上 | `<dataset>/03-韦恩图/Venn_{G3,G4}_{All,Up,Down}_<id>.pdf/.png` |
| 火山图（单图 + 2×3 网格） | 同上 | 同上 | `<dataset>/04-火山图/Volcano_<comp>_<id>.pdf`、`Volcano_Grid_<id>.pdf/.png` |
| GO 气泡图（旧版 Top15） | `scripts/1009_pipeline/03_go_enrichment.R` | `_cache/<id>.rds` | `<dataset>/05-GO富集/GO_Bubble_{Up,Down}_<id>.pdf/.png`、`_cache/plots_go_all.rds` |
| GO 气泡图（正式版 Top5） | `scripts/1009_pipeline/03b_go_top30.R` | `_cache/<id>.rds` | `<dataset>/05-GO富集/GO_Bubble_{Up,Down}_<id>.pdf/.png`、`00-汇总/GO富集_Top30_汇总.xlsx` |
| 热图（旧版 Top50） | `scripts/1009_pipeline/04_heatmap.R` | `_cache/<id>.rds` | `<dataset>/06-热图/TopDEGs_Heatmap_<id>.pdf/.png` |
| 热图（正式版 Top8） | `scripts/1009_pipeline/04b_heatmap_top50.R` | `_cache/<id>.rds` | `<dataset>/06-热图/TopDEGs_Heatmap_<id>.pdf/.png` |
| 合成图（a–f 模板） | `scripts/1009_pipeline/05_composite.R` | `_cache/plots_<id>.rds` + `_cache/plots_go_all.rds` | `<dataset>/07-合成图/Composite_<id>.pdf/.png` |
| DEG 明细表 | `scripts/1009_pipeline/01_deseq2.R` | `data/raw/MPXV_RB_counts_gene*.txt` + `metadata/sample_sheet.csv` | `<dataset>/00-表格/DEG_全部比较_<id>.xlsx`、`00-汇总/汇总_DEG数量_四套.xlsx` |
| 绘图依据表 | `scripts/1009_pipeline/02_basic_figures.R` | `_cache/<id>.rds` | `<dataset>/00-表格/绘图依据_PCA-DEG数量-韦恩_<id>.xlsx` |
| GO 明细表 | `scripts/1009_pipeline/03*.R` | `_cache/<id>.rds` | `<dataset>/00-表格/GO富集_<id>.xlsx` |
| 热图数据表 | `scripts/1009_pipeline/04*.R` | `_cache/<id>.rds` | `<dataset>/00-表格/热图数据_<id>.xlsx` |

> **缓存依赖**：`02_basic_figures.R` 会写 `_cache/plots_<id>.rds`；`03_go_enrichment.R` 会写 `_cache/plots_go_all.rds`。`05_composite.R` 依赖这两类缓存，因此不能跳过 `02` 和 `03`（即使用 03b 作为正式 GO 版）。

---

## 2. 三品系图件流程（BALB/c vs C57BL/6）

输入：`data/raw/3mice/BALB_F_counts_gene(1).txt`、`C57_F_counts_gene(1).txt`。输出：`results/3mice/`。

| 步骤 | 脚本 | 输出 |
|---|---|---|
| 导出全部显著 DEG + VST 表达 | `scripts/3mice_pipeline/01_export_all_degs.R` | `All_DEGs_for_Heatmap_Selection.xlsx` |
| UpSet 图 | `scripts/3mice_pipeline/02_upset_analysis.R` | `Upset_DEGs_Interferon_Mice.pdf`、`DEGs_Upset_Results.xlsx` |
| Top DEG 热图 | `scripts/3mice_pipeline/03_heatmap_top_genes.R` | `TopDEGs_Heatmap_BALBc.pdf`、`TopDEGs_Heatmap_C57BL6.pdf`、`TopDEGs_Heatmap_C57_BALB.pdf`、`TopDEGs_Heatmap_Data.xlsx` |
| GO 气泡图 | `scripts/3mice_pipeline/04_go_enrich_bubble.R` | `GO/GO_Bubble_Up_Regulated_v2.pdf`、`GO_Bubble_Down_Regulated_v2.pdf`、`GO_Bubble_Plot_Data_v2.xlsx`、`GO_Enrichment_BALB_C57_Full_v2.xlsx` |

`01_export_all_degs.R` 的 `All_DEGs_for_Heatmap_Selection.xlsx` 是 `03` 和 `04` 的共同输入，必须先跑。

---

## 3. 历史脚本（`legacy/`）

`legacy/scripts/` 下是首版投稿与二次调图的代码，按文件内容去重后保留 133 个；完整原始路径见 `legacy/MANIFEST.tsv`，按阶段/关键词分类见 `legacy/README.md` 和 `legacy/FAMILIES.md`。这些脚本包含旧的绝对路径（如 `D:/Desktop/...`），不能直接运行，仅用于追溯每张图的生成过程。

主要图件族与脚本关键词：

| 图件族 | 关键词 / 代表脚本 |
|---|---|
| 8BC 株系分析、火山图 | `8bc.R`、`8bc-2次调图.R`、`8bc-counts` 相关脚本 |
| AGB6 用药分析 + GO + 韦恩 | `AGB6-用药.R`、`AGB6-用药-go上调/下调.R`、`绘制AGB6--韦恩.R` |
| STAT1 用药分析 + GO + 韦恩 | `STAT1.R`、`STAT1-用药.R`、`stat1-用药-go/成功绘制-*.R`、`绘制stat1--韦恩.R` |
| C57 热图 / GO | `c57.R`、`C57-热图绘制.R`、`C57-跑通代码go分析*.R` |
| C57-8BC 通路 GO（二次调图） | `成功绘制-c57-8bc通路-终-上调.R`、`成功绘制-c57-8bc通路-下调-终.R` |
| 三小鼠 7 天热图 / 柱状图 | `1214-3小鼠-7天特制基因热图.R`、`5-3小鼠-7天特制基因热图.R` |
| 三小鼠 UpSet | `终-upset调图.R`、`美化的所有时间点的upset图.R` |
| 三小鼠 GO | `1215-绘制3小鼠的go通路-上调/下调.R` |

---

## 4. 预览图

`docs/figures/` 下是低分辨率预览，仅用于 README 展示，不是投稿原图：

- `upset_DEGs.pdf`
- `heatmap_topDEGs_BALBc.pdf` / `_C57BL6.pdf` / `_C57_BALB.pdf`
- `GO_bubble_up.pdf` / `GO_bubble_down.pdf`
- `volcano_panel_c.png`
- `composite_BALBc_Lung.png` / `_BALBc_Spleen.png` / `_AG129_Lung.png` / `_AG129_Spleen.png`
