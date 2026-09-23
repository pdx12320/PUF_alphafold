# -*- coding: utf-8 -*-
"""
Step 5: 汇总 Boltz 验证结果 -> 最终排名表
- 解析每个变体的 confidence json: ipTM / pLDDT / pTM
- 解析复合物 cif: 蛋白-RNA 接触数 (CA 距 RNA 重原子 < 8A 的残基数)
- 与 WT 对比, 输出 final_ranked.csv
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd


def complex_metrics(cif):
    """从 mmCIF 提取 RNA 接触残基数。"""
    import warnings
    warnings.filterwarnings("ignore")
    from Bio.PDB import MMCIFParser
    s = MMCIFParser(QUIET=True).get_structure("m", cif)[0]
    chains = {c.id: c for c in s}
    pa = "A" if "A" in chains else list(chains)[0]
    ra = "R" if "R" in chains else list(chains)[-1]
    rna_coords = np.array([a.coord for a in chains[ra].get_atoms()])
    n_contact = 0
    for r in chains[pa]:
        if r.id[0] != " " or "CA" not in r:
            continue
        d = np.sqrt(((rna_coords - r["CA"].coord) ** 2).sum(1)).min()
        if d < 8.0:
            n_contact += 1
    return n_contact


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz_dir", default="workdir/validate_boltz")
    ap.add_argument("--ranked", default="results/aice_single_ranked.csv")
    ap.add_argument("--out", default="results/final_ranked.csv")
    args = ap.parse_args()

    rows = []
    for conf in glob.glob(os.path.join(
            args.boltz_dir, "boltz_results_*", "predictions", "*",
            "confidence_*.json")):
        name = conf.split("/predictions/")[-1].split("/")[0]
        c = json.load(open(conf))
        cif = glob.glob(os.path.join(os.path.dirname(conf),
                                     f"{name}_model_0.cif"))
        nc = complex_metrics(cif[0]) if cif else None
        rows.append({
            "variant": name,
            "iptm": round(c.get("iptm", 0), 4),
            "ptm": round(c.get("ptm", 0), 4),
            "plddt": round(c.get("complex_plddt", 0), 4),
            "rna_contacts_lt8A": nc,
        })
    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("未找到 boltz 结果, 检查 --boltz_dir")

    wt = df[df["variant"] == "WT"]
    wt_iptm = float(wt["iptm"].iloc[0]) if len(wt) else np.nan
    wt_nc = float(wt["rna_contacts_lt8A"].iloc[0]) if len(wt) else np.nan
    df["delta_iptm_vs_WT"] = (df["iptm"] - wt_iptm).round(4)
    df["delta_contacts_vs_WT"] = df["rna_contacts_lt8A"] - wt_nc

    # 合并 AiCE 筛选信息 (单突变)
    if os.path.exists(args.ranked):
        r = pd.read_csv(args.ranked)
        r["variant"] = r.apply(
            lambda x: f"{x['wt']}{x['pos_1based']}{x['ligandmpnn_mut']}", axis=1)
        df = df.merge(r[["variant", "pos_1based", "rank_score", "consensus",
                         "ligandmpnn_f_max", "proteinmpnn_f_max",
                         "at_rna_interface"]],
                      on="variant", how="left")

    # 综合分: AiCE 出现率 + Boltz 验证
    df["rank_score"] = df.get("rank_score", pd.Series(0, index=df.index))
    df["final_score"] = (df["rank_score"].fillna(0)
                         + 2.0 * df["delta_iptm_vs_WT"].fillna(0)
                         + 0.02 * df["delta_contacts_vs_WT"].fillna(0))
    df = df.sort_values("final_score", ascending=False)
    df.to_csv(args.out, index=False)
    print(df.to_string(index=False))
    print(f"\n-> {args.out}")
    print("优先送实验: final_score 高且 consensus=True 且 code_safe 的变体")


if __name__ == "__main__":
    main()
