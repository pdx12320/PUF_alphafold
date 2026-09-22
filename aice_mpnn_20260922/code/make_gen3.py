# -*- coding: utf-8 -*-
"""
第三代 PUF12 构建清单生成器
- 从同目录 puf12_config.py 读 WT 序列
- 每个突变先对 WT 做断言校验，再落地 fasta + csv
输出: ../results/gen3_constructs.fasta, ../results/gen3_constructs.csv
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puf12_config import PUF12_SEQ  # noqa: E402

OUT_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))

# ---- 突变包定义 ------------------------------------------------------------
GVE  = ["S288G", "Y289V"]                       # 湿实验赢家 R8 三联体 (第一批 #31)
NTQ  = ["Y325T"]                                # R9: NYQ->NTQ (备选高保留骨架)
SYVIRR = ["E256R"]                              # R7 局部 252-257 SYVIER->SYVIRR
CORE = ["T350V", "D374E", "M458L", "L274I", "T278I", "V366I"]  # AiCE 共识核心(避开 R8 的 V294I)
GEO  = ["N12S", "S30A", "L41R"]                 # C388 几何区 (N端帽子)
IFACE = ["A438G", "H392D"]                      # 界面亲和力

CONSTRUCTS = [
    ("G0_P8GVE_ctrl",   GVE,                                  "对照: 湿实验赢家骨架 (第一批 #31)"),
    ("G1_GVE+R45T",     GVE + ["R45T"],                       "几何: C388 朝向单点"),
    ("G2_GVE+core",     GVE + CORE,                           "亲和力/稳定性核心包"),
    ("G3_GVE+R45T+core", GVE + ["R45T"] + CORE,               "几何+核心"),
    ("G4_GVE+iface",    GVE + IFACE,                          "界面亲和力"),
    ("G5_GVE+geo",      GVE + GEO,                            "C388 几何包 (N12S/S30A/L41R)"),
    ("G6_GVE+all",      GVE + ["R45T"] + CORE + GEO + IFACE,  "全叠加 14 突变 (风险最高, 观察用)"),
    ("G7_NTQ+core",     NTQ + ["R45T"] + CORE,                "备选骨架 P9-NTQ (C388 保留 0.94)"),
    ("G8_NTQ+iface",    NTQ + ["R45T"] + IFACE,               "备选骨架 + 界面"),
    ("G9_GVE+NTQ+core", GVE + NTQ + ["R45T"] + CORE,          "双选择性骨架叠加"),
    ("G10_SYVIRR+core", SYVIRR + CORE,                        "高 C388 骨架 (湿实验保留率 1.10) + 核心包"),
]


def apply_muts(seq, muts):
    seq = list(seq)
    for m in muts:
        wt, pos, new = m[0], int(m[1:-1]), m[-1]
        assert seq[pos - 1] == wt, f"{m}: WT 第 {pos} 位是 {seq[pos-1]}, 不是 {wt}"
        assert wt != new, f"{m}: 无效突变"
        seq[pos - 1] = new
    return "".join(seq)


def main():
    rows, fasta = [], []
    for name, muts, note in CONSTRUCTS:
        s = apply_muts(PUF12_SEQ, muts)
        assert len(s) == 493
        rows.append([name, "+".join(muts), len(muts), note, len(s)])
        fasta.append(f">{name} n_mut={len(muts)}\n{s}\n")
        print(f"[ok] {name:20s} {len(muts):2d} muts: {'+'.join(muts)}")

    with open(os.path.join(OUT_DIR, "gen3_constructs.fasta"), "w") as f:
        f.write("".join(fasta))
    with open(os.path.join(OUT_DIR, "gen3_constructs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["构建", "突变", "n_mut", "设计意图", "序列长度"])
        w.writerows(rows)
    print(f"\n已生成 {len(rows)} 个构建 -> gen3_constructs.fasta / gen3_constructs.csv")


if __name__ == "__main__":
    main()
