# -*- coding: utf-8 -*-
"""
Step 4a: 把 AiCE 提名的变体写成 Boltz 输入 YAML（准备结构验证）
用法:
  python step4_make_yamls.py --top 20                 # top 20 单突变
  python step4_make_yamls.py --combos results/aice_multi_combos.csv --n 5
"""
import argparse
import os
import re

import pandas as pd

from puf12_config import (PUF12_SEQ, RNA_CONTEXT_17, parse_repeats,
                          target_rna, write_boltz_yaml)


def apply_mut(seq, mutstr):
    """mutstr 形如 'E42K+A100V'，1-based。"""
    s = list(seq)
    for m in mutstr.split("+"):
        mo = re.fullmatch(r"([A-Z])(\d+)([A-Z])", m)
        assert mo, f"突变格式错误: {m}"
        wt, pos, mt = mo.group(1), int(mo.group(2)), mo.group(3)
        assert s[pos - 1] == wt, f"{m}: 第{pos}位是{s[pos-1]}不是{wt}"
        s[pos - 1] = mt
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranked", default="results/aice_single_ranked.csv")
    ap.add_argument("--combos", default=None)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--outdir", default="workdir/validate_yamls")
    ap.add_argument("--rna", default=None,
                    help="验证用 RNA; 默认用 17nt 上下文序列(含 C388), "
                         "与 AF3 WT baseline 可比")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 默认 17nt (ACAUGGAGGACGUGCGC): 与 wildtype_seed1 的 AF3 任务同序列,
    # ipTM 可直接和 WT 0.91 对比
    rna = args.rna or RNA_CONTEXT_17
    jobs = []

    # 野生型对照
    write_boltz_yaml(PUF12_SEQ, rna, os.path.join(args.outdir, "WT.yaml"))
    jobs.append(("WT", PUF12_SEQ))

    if args.combos and os.path.exists(args.combos):
        cdf = pd.read_csv(args.combos).head(args.n)
        for i, row in cdf.iterrows():
            s = apply_mut(PUF12_SEQ, row["combo"])
            nm = f"combo{i+1}"
            write_boltz_yaml(s, rna, os.path.join(args.outdir, f"{nm}.yaml"))
            jobs.append((nm, s))
    else:
        df = pd.read_csv(args.ranked)
        df = df[df["aice_votes"] > 0].head(args.top)
        for i, row in df.iterrows():
            mut = f"{row['wt']}{row['pos_1based']}{row['ligandmpnn_mut']}"
            s = apply_mut(PUF12_SEQ, mut)
            write_boltz_yaml(s, rna, os.path.join(args.outdir, f"{mut}.yaml"))
            jobs.append((mut, s))

    print(f"已生成 {len(jobs)} 个 YAML -> {args.outdir}/")
    print("下一步: bash step4_validate.sh")


if __name__ == "__main__":
    main()
