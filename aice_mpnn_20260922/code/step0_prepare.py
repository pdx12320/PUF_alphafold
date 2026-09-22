# -*- coding: utf-8 -*-
"""
Step 0: 解析 PUF12 重复单元 -> 反推 RNA 靶点 -> 生成 Boltz 输入 YAML
用法: python step0_prepare.py [--outdir workdir]
"""
import argparse
import csv
import os

from puf12_config import (APOE4_MRNA, PUF12_SEQ, parse_repeats, target_rna,
                          write_boltz_yaml, save_config)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="workdir")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cfg = save_config(args.outdir)
    repeats = parse_repeats()
    rna = target_rna(repeats)

    print(f"PUF12 长度: {len(PUF12_SEQ)} aa, 解析出 {len(repeats)} 个 repeat")
    print(f"{'repeat':>6} {'pos(1-based)':>12} {'triplet':>8} {'base':>5}")
    for r in repeats:
        print(f"{r['repeat']:>6} {r['pos12']+1:>12} {r['triplet']:>8} {r['base']:>5}")
    print(f"\nRNA 靶点 (5'->3'): {rna}")
    idx = APOE4_MRNA.find(rna)
    if idx >= 0:
        print(f"位于 APOE4 mRNA 第 {idx+1}-{idx+len(rna)} 位; "
              f"C388 位点(1-based 388)在靶点下游第 {388-(idx+len(rna))} nt")
    else:
        print("警告: 靶序列未在 APOE4 mRNA 中找到, 请检查密码表/序列!")

    # repeat 解析结果落盘, 供人工核对
    with open(os.path.join(args.outdir, "repeat_parse.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["repeat", "pos12_1based", "pos13_1based", "pos16_1based",
                    "triplet", "base"])
        for r in repeats:
            w.writerow([r["repeat"], r["pos12"] + 1, r["pos13"] + 1,
                        r["pos16"] + 1, r["triplet"], r["base"]])

    # Boltz 输入: 全长 PUF12 + 12nt 靶 RNA (各取前后 2nt 保护碱基可选)
    write_boltz_yaml(PUF12_SEQ, rna,
                     os.path.join(args.outdir, "puf12_complex.yaml"))
    # 另出一份 17nt 上下文版 (含 C388), 用于变体验证, 与 AF3 WT 任务可比
    from puf12_config import RNA_CONTEXT_17
    write_boltz_yaml(PUF12_SEQ, RNA_CONTEXT_17,
                     os.path.join(args.outdir, "puf12_complex_17nt.yaml"))
    print(f"\n已生成: {args.outdir}/puf12_complex.yaml  (Boltz 输入, 12nt)")
    print(f"已生成: {args.outdir}/puf12_complex_17nt.yaml  (验证用, 17nt 含 C388)")
    print(f"已生成: {args.outdir}/design_config.json, repeat_parse.csv")
    print("下一步: bash step1_fold_complex.sh  (或放入你已有的复合物 PDB)")


if __name__ == "__main__":
    main()
