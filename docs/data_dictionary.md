# 数据字典

## 1. `metadata/sample_sheet.csv`

66 行 × 8 列，一行一个样本（每只小鼠的肺或脾）。文件为 UTF-8 with BOM；第一列 `col_index` 是计数矩阵中的样本列序号（1–66）。

| 列名 | 含义 | 示例 |
|---|---|---|
| `col_index` | 计数矩阵样本列序号，1–66 | `1` |
| `fc_col` | 计数矩阵中的原始列名 | `1-lung-1/Aligned.sortedByCoord.out.bam` |
| `run_id` | 实验方样本编号 | `1-肺-1` |
| `seq_name` | 测序/分组编号 | `BALB/-G1-Lung-1` |
| `strain` | 品系：`BALB/c` 或 `AG129` | `BALB/c` |
| `group` | 分组：`G1`–`G4` | `G1` |
| `tissue` | 组织：`Lung` 或 `Spleen` | `Lung` |
| `rep` | 组内重复编号 | `1` |

分组定义：

| 组 | 含义 |
|---|---|
| `G1` | 未感染（uninfected） |
| `G2` | 感染（infected） |
| `G3` | 治疗组 1（treatment 1） |
| `G4` | 治疗组 2（treatment 2） |

> G3/G4 的具体药物名尚未写入脚本；更新 `scripts/1009_pipeline/utils/common.R` 中的 `GROUP_DESC` 即可统一替换图例。

---

## 2. 原始计数矩阵

主流程默认读取 `data/raw/MPXV_RB_counts_gene.txt`（或同目录下匹配 `MPXV_RB_counts_gene*.txt` 的文件）。

- 来源：featureCounts v2.0.3
- 参数：`-p -B -C -t exon -g gene_name`
- 参考基因组：Cell Ranger `refdata-gex-GRCm39-2024-A`（小鼠 GRCm39）
- 规模：33,670 个基因 × 66 个样本
- 格式：

```text
Geneid  Chr  Start  End  Strand  Length  <66 个样本列>
Gene1   chr1 ...    ...  +       ...     0 12 35 ...
...
```

第 1 行是 featureCounts 的命令行注释，`common.R::read_counts_matrix()` 会 `skip = 1` 跳过；前 5 列注释写入 `gene_meta`，其余列作为 count 矩阵。

三品系流程使用另外两个计数矩阵，放在 `data/raw/3mice/`：

```text
BALB_F_counts_gene(1).txt
C57_F_counts_gene(1).txt
```

---

## 3. 主流程输出表

### 3.1 `DEG_全部比较_<dataset>.xlsx`

| Sheet | 内容 |
|---|---|
| `README` | 数据集、样本数、分组、过滤参数、阈值、生成时间 |
| `All_Sig_DEGs` | 所有比较的显著 DEG（`padj < 0.05` 且 `\|log2FC\| > 1`）合并表 |
| `G2vsG1` … `G4vsG3` | 每个比较组的完整 DESeq2 结果 |
| `DEG_Summary` | 每个比较的 Up / Down / Total_Sig / NotSig / NA_padj 计数 |
| `VST_Expression_AllGenes` | 全部基因的 VST 归一化表达矩阵 |
| `VST_Expression_DEGs` | 仅显著 DEG 的 VST 表达矩阵 |

DEG 表主要字段：

| 字段 | 含义 |
|---|---|
| `Geneid` | 基因名 |
| `comparison` | 比较组 ID（如 `G2vsG1`） |
| `comparison_label` | 可读标签（如 `G2 vs G1`） |
| `change` | `Up` / `Down` / `NS` |
| `baseMean` | DESeq2 标准化平均表达 |
| `log2FoldChange` | log2 倍数变化（G1 为参照） |
| `lfcSE` | log2FC 标准误 |
| `stat` | Wald 统计量 |
| `pvalue` | 原始 p 值 |
| `padj` | BH 校正后 p 值 |

### 3.2 `绘图依据_PCA-DEG数量-韦恩_<dataset>.xlsx`

保存 PCA 坐标、DEG 数量、韦恩图各集合的基因列表，便于核对图中数字。

### 3.3 `GO富集_<dataset>.xlsx`

保存每个比较组的 GO 富集结果；`Plot_Up` / `Plot_Down` 等 sheet 是气泡图的绘图数据（通路、Count、p.adjust、FoldEnrichment、RichFactor 等）。

### 3.4 `热图数据_<dataset>.xlsx`

保存每个比较组选中的 Top 基因、VST 表达值、行 Z-score 和样本注释，便于复现热图。

### 3.5 `00-汇总/`

| 文件 | 内容 |
|---|---|
| `汇总_DEG数量_四套.xlsx` | 4 套数据集 × 6 个比较的 Up/Down/Total 汇总 + 分组结构 |
| `DEG数量_四套汇总.pdf/.png` | 四套数据集合并柱状图 |
| `GO富集_Top30_汇总.xlsx` | GO 气泡图（正式版 Top5）的绘图数据 |

---

## 4. 三品系流程输出

| 文件 | 内容 |
|---|---|
| `All_DEGs_for_Heatmap_Selection.xlsx` | 所有显著 DEG + VST 表达，供热图和 GO 挑选 |
| `DEGs_Upset_Results.xlsx` | UpSet 各集合的基因列表与交集计数 |
| `TopDEGs_Heatmap_Data.xlsx` | Top DEG 热图的表达矩阵与注释 |
| `GO_Bubble_Plot_Data_v2.xlsx` | GO 气泡图绘图数据 |
| `GO_Enrichment_BALB_C57_Full_v2.xlsx` | BALB/c 与 C57BL/6 的完整 GO 富集结果 |

---

## 5. 分析参数

| 参数 | 默认值 | 位置 |
|---|---|---|
| 显著阈值 `PADJ_CUT` | 0.05 | `scripts/1009_pipeline/utils/common.R` |
| 倍数变化阈值 `LFC_CUT` | 1.0（即 `\|log2FC\| > 1`） | 同上 |
| 低表达过滤 `MIN_COUNT` | 10 | 同上 |
| GO 每图 Top N | 5（正式版） | `scripts/1009_pipeline/03b_go_top30.R` |
| 热图每组 Top N | 8（正式版） | `scripts/1009_pipeline/04b_heatmap_top50.R` |
| 火山图标注数 | 15（02 脚本内） | `scripts/1009_pipeline/02_basic_figures.R` |
| 火山图 y 轴截断 | 99.5% 分位 | 同上 |

> 修改参数后建议清空 `results/1009/_cache/` 再重跑，避免旧缓存混入。
