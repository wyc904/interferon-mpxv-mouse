# -*- coding: utf-8 -*-
"""从「肺与脾组织样本编号对照表.xlsx」生成 tidy sample_sheet.csv（按矩阵列顺序），
并与 MPXV_RB_counts_gene(1).txt 表头逐列核对、统计每个样本的文库大小/检出基因数。

输出（均在 2-script/）：
  sample_sheet.csv       col_index, fc_col, run_id, seq_name, strain, group, tissue, rep
  sample_sheet_notes.txt 核对与 QC 报告
"""
import csv
import re
import statistics
from pathlib import Path

import openpyxl

SRC = Path(r"D:\Desktop\临时保存\5-干扰素小鼠\3-1009\1-源")
OUTDIR = SRC.parent / "2-script"
CSV_OUT = OUTDIR / "sample_sheet.csv"
NOTES = OUTDIR / "sample_sheet_notes.txt"

notes = []


def n(s=""):
    notes.append(str(s))


# ---------------- 1. 解析对照表 ----------------
wb = openpyxl.load_workbook(SRC / "肺与脾组织样本编号对照表.xlsx", data_only=True)
ws = wb["Sheet1"]

pairs = []
for row in ws.iter_rows(values_only=True):
    row = list(row)
    while len(row) < 6:
        row.append(None)
    a, b, e, f = row[0], row[1], row[4], row[5]
    if a is not None and str(a).strip() and b is not None and str(b).strip():
        pairs.append((str(a).strip(), str(b).strip()))
    if e is not None and str(e).strip() and f is not None and str(f).strip():
        pairs.append((str(e).strip(), str(f).strip()))

RUN_RE = re.compile(r"^([0-9]+|[A-Za-z])-(肺|脾)-([0-9]+)$")
SEQ_RE = re.compile(r"^(.*?)-G([0-9])-(Lung|Spleen)-([0-9]+)$")
TISSUE = {"肺": "Lung", "脾": "Spleen"}
GROUP_OF = {"A": "G1", "B": "G2", "C": "G3", "D": "G4"}

problems = []
rows = []
for run_id, seq_name in pairs:
    m = RUN_RE.match(run_id)
    if not m:
        problems.append(f"编号解析失败: {run_id!r}")
        continue
    mid, tis_cn, rep = m.group(1), m.group(2), int(m.group(3))
    if mid.isdigit():
        strain, group = "BALB/c", f"G{int(mid)}"
    else:
        strain, group = "AG129", GROUP_OF.get(mid.upper(), "?")
    tissue = TISSUE[tis_cn]
    fc_col = f"{mid}-{tissue.lower()}-{rep}/Aligned.sortedByCoord.out.bam"

    sm = SEQ_RE.match(seq_name)
    if not sm:
        problems.append(f"测序名解析失败: {seq_name!r} (编号 {run_id})")
    else:
        pfx, g, t, r = sm.groups()
        want = {"BALB/", "BALB/-", "BALB/c"} if strain == "BALB/c" else {"AG129"}
        if pfx not in want or f"G{g}" != group or t != tissue or int(r) != rep:
            problems.append(f"编号-测序名不一致: {run_id} vs {seq_name} (pfx={pfx!r})")

    rows.append(dict(run_id=run_id, seq_name=seq_name, strain=strain, group=group,
                     tissue=tissue, rep=rep, fc_col=fc_col, mouse=mid))

n("样本表构建与核对报告 (build_sample_sheet.py)")
n()
n("[1] 对照表解析")
n(f"- 解析出样本数: {len(rows)}  (BALB/c {sum(r['strain']=='BALB/c' for r in rows)} + "
  f"AG129 {sum(r['strain']=='AG129' for r in rows)})")
n(f"- 解析问题: {problems if problems else '无'}")

# ---------------- 2. 与矩阵表头核对 ----------------
txt = SRC / "MPXV_RB_counts_gene(1).txt"
with open(txt, encoding="utf-8", errors="replace") as fh:
    header = None
    for ln in fh:
        if ln.startswith("#"):
            continue
        header = ln.rstrip("\n")
        break
cols = header.split("\t")
assert cols[0] == "Geneid", cols[:6]
meta_cols, sample_cols = cols[:6], cols[6:]
by_fc = {}
for r in rows:
    by_fc.setdefault(r["fc_col"], []).append(r)

n()
n("[2] 与计数矩阵表头核对")
n(f"- 矩阵列: {cols[:6]} + 样本列 {len(sample_cols)} 个")
miss = [c for c in sample_cols if c not in by_fc]
extra = [c for c in by_fc if c not in set(sample_cols)]
dups = [c for c, v in by_fc.items() if len(v) > 1]
n(f"- 矩阵有而对照表无: {miss if miss else '无'}")
n(f"- 对照表有而矩阵无: {extra if extra else '无'}")
n(f"- 重复 fc_col: {dups if dups else '无'}")
n(f"- 结果: {'完全一致' if not (miss or extra or dups) else '!! 存在不一致'}")

# ---------------- 3. 读计数，算文库大小 ----------------
sums, nz = [0] * len(sample_cols), [0] * len(sample_cols)
n_genes = 0
with open(txt, encoding="utf-8", errors="replace") as fh:
    seen_header = False
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        if not seen_header:
            seen_header = True
            continue
        parts = ln.rstrip("\n").split("\t")
        if len(parts) != 6 + len(sample_cols):
            problems.append(f"数据行列数异常({len(parts)}): {parts[0]!r}")
            continue
        n_genes += 1
        for i, v in enumerate(parts[6:]):
            x = int(v)
            sums[i] += x
            if x:
                nz[i] += 1

# ---------------- 4. 写 CSV（按矩阵列顺序） ----------------
order = {c: i for i, c in enumerate(sample_cols)}
final = sorted(rows, key=lambda r: order[r["fc_col"]])
with open(CSV_OUT, "w", newline="", encoding="utf-8-sig") as fh:
    wr = csv.writer(fh)
    wr.writerow(["col_index", "fc_col", "run_id", "seq_name", "strain", "group", "tissue", "rep"])
    for i, r in enumerate(sample_cols, 1):
        r0 = by_fc[r][0]
        wr.writerow([i, r, r0["run_id"], r0["seq_name"], r0["strain"], r0["group"], r0["tissue"], r0["rep"]])

# ---------------- 5. 报告：结构 + QC ----------------
n()
n("[3] 分组结构（以矩阵为准）")
n(f"- 基因行数 = {n_genes}, 唯一基因名 = {n_genes}")
idx = {c: i for i, c in enumerate(sample_cols)}
grid = {}
for r in rows:
    grid.setdefault((r["strain"], r["group"], r["tissue"]), []).append(r["rep"])
for strain in ("BALB/c", "AG129"):
    for g in ("G1", "G2", "G3", "G4"):
        for tis in ("Lung", "Spleen"):
            reps = sorted(grid.get((strain, g, tis), []))
            full = list(range(1, 6)) if strain == "BALB/c" else list(range(1, 5))
            missing = [f"{tis}-{x}" for x in full if x not in reps]
            n(f"  {strain:6s} {g} {tis:6s} n={len(reps)} reps={reps} 缺={missing if missing else '无'}")

n()
n("[4] 文库大小 (featureCounts 总计数) 与检出基因数")
stat = []
for i, c in enumerate(sample_cols):
    r0 = by_fc[c][0]
    stat.append((sums[i], nz[i], f"{r0['strain']}/{r0['group']}/{r0['tissue']}/{r0['rep']}", c))
stat.sort()
n(f"- 全部 {len(stat)} 个样本: min={stat[0][0]:,}  median={int(statistics.median(s[0] for s in stat)):,}  max={stat[-1][0]:,}")
n("- 最小 5 个:")
for s in stat[:5]:
    n(f"    {s[2]:32s} sum={s[0]:>12,}  检出基因={s[1]:,}")
n("- 最大 5 个:")
for s in stat[-5:]:
    n(f"    {s[2]:32s} sum={s[0]:>12,}  检出基因={s[1]:,}")
n("- 各分组中位文库大小:")
for strain in ("BALB/c", "AG129"):
    for g in ("G1", "G2", "G3", "G4"):
        for tis in ("Lung", "Spleen"):
            vals = [sums[idx[r["fc_col"]]] for r in rows
                    if r["strain"] == strain and r["group"] == g and r["tissue"] == tis]
            if vals:
                n(f"    {strain:6s} {g} {tis:6s} n={len(vals)}  median={int(statistics.median(vals)):>12,}"
                  f"  [{min(vals):,} ~ {max(vals):,}]")

n()
n("[5] 备注 / 待确认")
n("- seq_name 前缀字面为 'BALB/-' (疑为 BALB/c 的录入变体), 已在 strain 列规范为 BALB/c")
n("- docx 仅说明: 数字开头=BALB/c, 字母开头=AG129; G1-G4 代表什么(时间点/剂量/处理)文档未说明")
n("- counts 来源: featureCounts v2.0.3, -p -B -C -t exon -g gene_name, GTF=Cell Ranger refdata-gex-GRCm39-2024-A")
if problems:
    n()
    n("!! 解析/核对问题汇总:")
    for p in problems:
        n("  - " + p)

NOTES.write_text("\n".join(notes), encoding="utf-8")
print("WROTE:", CSV_OUT)
print("WROTE:", NOTES)
