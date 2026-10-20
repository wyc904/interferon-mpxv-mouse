# data/ 目录说明

**本目录不包含任何原始数据。** `.gitignore` 已排除 `data/raw/*` 与 `data/processed/*`，只保留 `.gitkeep`。

## 主流程需要的文件

```text
data/raw/
└─ MPXV_RB_counts_gene.txt     # featureCounts 基因计数矩阵（33,670 基因 × 66 样本）
```

- 文件名匹配 `MPXV_RB_counts_gene*.txt` 即可；也可以用环境变量 `MPXV_COUNTS_FILE` 指向任意位置。
- 格式说明见 [`../docs/data_dictionary.md`](../docs/data_dictionary.md) 第 2 节。
- 计数矩阵来自论文对应的测序数据；请按论文的数据可用性声明获取，或联系作者。
- 样本表 `metadata/sample_sheet.csv` 已包含在仓库中，无需另外准备。

## 三品系流程需要的文件

```text
data/raw/3mice/
├─ BALB_F_counts_gene(1).txt
└─ C57_F_counts_gene(1).txt
```

文件名必须与脚本中的读取路径一致（脚本按 `(1)` 命名）。如果文件名不同，可在 `scripts/3mice_pipeline/utils/config.R` 中调整 `DATA_DIR_3MICE`，或直接修改脚本里的文件名。

## data/processed/

预留给手工整理的表（例如早期 `1214-3小鼠-7天特制基因热图.R` 用到的 `热图绘制.xlsx`）。主流程不依赖此目录；该目录同样被 `.gitignore` 排除。

## 检查是否误提交数据

```bash
git status --short
git ls-files data/
```

正常情况下，`git ls-files data/` 只应看到 `data/README.md`、`data/raw/.gitkeep`、`data/processed/.gitkeep`。
