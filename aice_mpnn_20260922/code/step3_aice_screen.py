# -*- coding: utf-8 -*-
"""
Step 3: AiCE 筛选 (核心脚本, 复刻 Cell 2025 AiCE 流程 + PUF 特异性保护)

AiCE_single (论文 STAR Methods / 官方 repo 推荐阈值):
  AiCE1 = {i | f_max(i) >= beta=0.8}                     全局高出现率
  AiCE2 = {i | i 在柔性区(DSSP 非螺旋非折叠) 且 f_max(i) >= gamma=0.5}
  AiCE_single = AiCE1 ∪ AiCE2

AiCE_multi (EC + LD):
  EC: Rivoire 加权协方差, 取全局 EC 前 10% 的位点对
  LD: 采样序列的连锁不平衡 (官方流程用 plink 对伪 DNA 算 r^2;
      此处用氨基酸层面的 Cramér's V 近似, 阈值 0.5)

PUF 特有约束:
  - 12 个 repeat 的识别位点 (每个 repeat 第 12/13/16 位) 若突变,
    新三联体必须仍在同义密码表内 (保持识别原碱基), 否则剔除。

用法: python step3_aice_screen.py [--beta 0.8] [--gamma 0.5] [--outdir results]
"""
import argparse
import glob
import json
import os
import re
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from puf12_config import (PUF12_SEQ, parse_repeats, specificity_positions,
                          SYNONYMOUS_CODES)

AA = "ACDEFGHIKLMNPQRSTVWY"
AAI = {a: k for k, a in enumerate(AA)}
# UniProt 天然氨基酸背景频率 (q_a, EC 权重用)
Q_BG = dict(zip("ARNDCQEGHILKMFPSTWYV",
                [0.0825, 0.0553, 0.0406, 0.0545, 0.0137, 0.0393, 0.0675,
                 0.0707, 0.0227, 0.0593, 0.0966, 0.0584, 0.0242, 0.0386,
                 0.0470, 0.0656, 0.0534, 0.0108, 0.0292, 0.0687]))
Q = np.array([Q_BG[a] for a in AA])


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------
def read_sampled_fastas(pattern):
    """读取 MPNN 输出 fasta, 跳过 native 记录。
    兼容两种 header: ProteinMPNN 原版 'sample=1' / LigandMPNN 'id=1';
    native 记录 header 含 'num_res=' 且无 id/sample 编号。"""
    seqs = []
    for fn in glob.glob(pattern):
        name, seq = None, []

        def flush():
            if name is None or not seq:
                return
            s = "".join(seq)
            if len(s) != len(PUF12_SEQ):
                return
            if ("sample=" in name) or re.search(r"\bid=\d+", name):
                seqs.append(s)

        for line in open(fn):
            line = line.strip()
            if line.startswith(">"):
                flush()
                name, seq = line[1:], []
            else:
                seq.append(line)
        flush()
    if not seqs:
        raise SystemExit(f"未读到采样序列: {pattern}")
    arr = np.array([[AAI.get(c, 20) for c in s] for s in seqs])
    print(f"  读取 {len(seqs)} 条序列 x {arr.shape[1]} 位点 <- {pattern}")
    return arr


# ---------------------------------------------------------------------------
# DSSP 柔性区 (优先 mkdssp, 退化到 pydssp)
# ---------------------------------------------------------------------------
def flexible_mask(pdb_path, chain="A", L=None):
    flex = None
    try:
        from Bio.PDB import PDBParser, DSSP
        s = PDBParser(QUIET=True).get_structure("s", pdb_path)[0]
        dssp = DSSP(s, pdb_path, dssp="mkdssp")
        m = []
        for r in s[chain]:
            if r.id[0] != " ":
                continue
            ss = dssp.get((chain, r.id), (None, "C", None))[1]
            m.append(ss not in ("H", "E", "G", "I", "B"))  # 非规则二级结构=柔性
        flex = np.array(m)
        print("  DSSP: mkdssp")
    except Exception:
        pass
    if flex is None:
        try:
            import pydssp
            from Bio.PDB import PDBParser
            from Bio.PDB.Polypeptide import is_aa
            s = PDBParser(QUIET=True).get_structure("s", pdb_path)[0]
            res = [r for r in s[chain] if is_aa(r, standard=True)]
            bb = np.array([[r[a].coord for a in ("N", "CA", "C", "O")]
                           for r in res])          # (L, 4, 3)
            if bb.ndim == 3 and bb.shape[1] == 4:
                # pydssp.assign 接受 (N, L, 4, 3) 或 (L, 4, 3), 返回每残基 H/E/C
                ss = pydssp.assign(bb, out_type="c3")
                ss = np.atleast_1d(ss)
                flex = np.array([c == "C" for c in ss])
                print("  DSSP: pydssp (C=柔性)")
        except Exception as e:
            print(f"  DSSP 不可用({e}), 柔性区全部记 False -> 只用全局阈值 beta")
    if flex is not None and L is not None and len(flex) != L:
        print(f"  警告: DSSP 长度 {len(flex)} != 序列长度 {L}, 忽略柔性信息")
        flex = None
    return flex


# ---------------------------------------------------------------------------
# 界面位点: 蛋白 CA 距 RNA 任意重原子 < cutoff
# ---------------------------------------------------------------------------
def interface_mask(pdb_path, prot_chain="A", rna_chain="R", cutoff=10.0):
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("s", pdb_path)[0]
    rna_atoms = [a.coord for a in s[rna_chain].get_atoms()]
    if not rna_atoms:
        print("  警告: PDB 中无 RNA 链 R, 界面信息缺失")
        return None
    rna_atoms = np.array(rna_atoms)
    mask = []
    for r in s[prot_chain]:
        if r.id[0] != " ":
            continue
        ca = r["CA"].coord if "CA" in r else list(r.get_atoms())[0].coord
        d = np.sqrt(((rna_atoms - ca) ** 2).sum(1)).min()
        mask.append(d < cutoff)
    m = np.array(mask)
    print(f"  界面位点数 (CA 距 RNA<{cutoff}A): {m.sum()}")
    return m


# ---------------------------------------------------------------------------
# appearance rate / f_max
# ---------------------------------------------------------------------------
def appearance(arr):
    """arr: (M, L) int -> per-position freq 矩阵 (L,20), f_wt, f_max, x_mut"""
    M, L = arr.shape
    freq = np.zeros((L, 20))
    for k in range(20):
        freq[:, k] = (arr == k).mean(0)
    wt_idx = np.array([AAI[c] for c in PUF12_SEQ])
    f_wt = freq[np.arange(L), wt_idx]
    freq_wo_wt = freq.copy()
    freq_wo_wt[np.arange(L), wt_idx] = -1
    x_mut = freq_wo_wt.argmax(1)
    f_max = freq_wo_wt[np.arange(L), x_mut]
    # 只保留 mut 出现率 > wt 的 (AiCE 定义)
    valid = f_max > f_wt
    return freq, f_wt, f_max, x_mut, valid


# ---------------------------------------------------------------------------
# EC score (Rivoire 加权协方差) 与 LD (Cramér's V)
# ---------------------------------------------------------------------------
def ec_ld(arr, positions):
    """在候选位点子集上计算 EC 与 LD。"""
    P = len(positions)
    sub = arr[:, positions]                      # (M, P)
    M = sub.shape[0]
    f1 = np.zeros((P, 20))
    for k in range(20):
        f1[:, k] = (sub == k).mean(0)
    f1 = np.clip(f1, 1e-9, 1 - 1e-9)
    phi = np.log(f1 * (1 - Q)[None, :] / ((1 - f1) * Q[None, :]))
    EC = np.zeros((P, P))
    V = np.zeros((P, P))
    for i in range(P):
        xi = sub[:, i]
        for j in range(i + 1, P):
            xj = sub[:, j]
            fab = np.zeros((20, 20))
            for a in range(20):
                ma = xi == a
                if ma.any():
                    fab[a] = [(xj[ma] == b).mean() * ma.mean() for b in range(20)]
            Cab = fab - np.outer(f1[i], f1[j])
            EC[i, j] = EC[j, i] = np.abs(phi[i][:, None] * phi[j][None, :] * Cab).sum()
            # Cramér's V
            chi2 = (fab - np.outer(f1[i], f1[j])) ** 2 / \
                   np.clip(np.outer(f1[i], f1[j]), 1e-12, None)
            n_eff = M * chi2.sum()
            V[i, j] = V[j, i] = np.sqrt(max(n_eff / M, 0) / 19.0)
    return EC, V


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, default=0.8, help="全局阈值 (官方推荐 0.8)")
    ap.add_argument("--gamma", type=float, default=0.5, help="柔性区阈值 (官方推荐 0.5)")
    ap.add_argument("--pdb", default="workdir/complex.pdb")
    ap.add_argument("--lmpnn", default="workdir/sample_ligand_mpnn/seqs/*.fa")
    ap.add_argument("--pmpnn", default="workdir/sample_protein_mpnn/seqs/*.fa")
    ap.add_argument("--config", default="workdir/design_config.json")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--top_combo", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cfg = json.load(open(args.config))
    spec_pos = set(cfg["specificity_positions_0based"])
    repeats = cfg["repeats"]
    zone_lo, zone_hi = cfg.get("c388_zone_1based", [1, 45])
    L = len(PUF12_SEQ)

    print("[1/5] 读取采样序列")
    arrs = {}
    arrs["ligandmpnn"] = read_sampled_fastas(args.lmpnn)
    try:
        arrs["proteinmpnn"] = read_sampled_fastas(args.pmpnn)
    except SystemExit:
        print("  (无 ProteinMPNN 结果, 只用 LigandMPNN)")

    print("[2/5] 结构约束: DSSP 柔性区 + RNA 界面")
    flex = flexible_mask(args.pdb, "A", L)
    iface = interface_mask(args.pdb, "A", "R")

    print("[3/5] appearance rate / f_max / AiCE_single")
    stats = {}
    for name, arr in arrs.items():
        freq, f_wt, f_max, x_mut, valid = appearance(arr)
        stats[name] = dict(freq=freq, f_wt=f_wt, f_max=f_max,
                           x_mut=x_mut, valid=valid)

    # PUF 特异性保护: 识别位点的突变三联体必须是同义密码
    def code_safe(pos, aa_idx):
        if pos not in spec_pos:
            return True
        # 找到所属 repeat, 用新氨基酸替换对应位后检查三联体
        for r in repeats:
            for key, offset in (("pos12", 0), ("pos13", 1), ("pos16", 4)):
                if r[key] == pos:
                    t = list(PUF12_SEQ[r["pos12"]] + PUF12_SEQ[r["pos13"]] +
                             PUF12_SEQ[r["pos16"]])
                    t[{0: 0, 1: 1, 4: 2}[offset]] = AA[aa_idx]
                    return "".join(t) in SYNONYMOUS_CODES.get(r["base"], set())
        return True

    rows = []
    for i in range(L):
        row = {"pos_1based": i + 1, "wt": PUF12_SEQ[i],
               "in_flexible": bool(flex[i]) if flex is not None else None,
               "at_rna_interface": bool(iface[i]) if iface is not None else None,
               # C388 毗邻区(N端帽子区): 此处的突变会改变脱氨酶融合端相对
               # C388 的位置几何 -> 影响 C388 选择性, 与亲和力位点区分开
               "near_c388_zone": zone_lo <= i + 1 <= zone_hi,
               "is_specificity_pos": i in spec_pos}
        aice_votes = 0
        for name in arrs:
            st = stats[name]
            mut = AA[st["x_mut"][i]]
            fm = st["f_max"][i]
            row[f"{name}_mut"] = mut
            row[f"{name}_f_max"] = round(float(fm), 4)
            row[f"{name}_f_wt"] = round(float(st["f_wt"][i]), 4)
            in_aice = st["valid"][i] and (
                fm >= args.beta or
                (flex is not None and flex[i] and fm >= args.gamma))
            safe = code_safe(i, st["x_mut"][i])
            row[f"{name}_aice_single"] = bool(in_aice and safe)
            row[f"{name}_code_safe"] = bool(safe)
            aice_votes += int(in_aice and safe)
        row["consensus"] = aice_votes == len(arrs)
        row["aice_votes"] = aice_votes
        rows.append(row)

    df = pd.DataFrame(rows)
    df["rank_score"] = (
        df[[f"{n}_f_max" for n in arrs]].max(axis=1)
        + 0.2 * df["consensus"].astype(float)
        + 0.1 * df["at_rna_interface"].fillna(False).astype(float)
    )
    df = df.sort_values("rank_score", ascending=False)
    df.to_csv(os.path.join(args.outdir, "aice_single_ranked.csv"), index=False)
    nom = df[df["aice_votes"] > 0]
    print(f"  AiCE_single 提名 {len(nom)} 个位点 "
          f"(其中两模型共识 {df['consensus'].sum()} 个)")

    print("[4/5] EC + LD (AiCE_multi)")
    cand = nom["pos_1based"].values - 1
    combos = []
    if len(cand) >= 2:
        arr = arrs["ligandmpnn"]
        EC, V = ec_ld(arr, cand)
        ec_vals = EC[np.triu_indices(len(cand), 1)]
        ec_thr = np.quantile(ec_vals, 0.9) if len(ec_vals) else np.inf
        pairs = []
        for a in range(len(cand)):
            for b in range(a + 1, len(cand)):
                if EC[a, b] >= ec_thr or V[a, b] >= 0.5:
                    pairs.append((cand[a], cand[b], EC[a, b], V[a, b]))
        # 贪心组合: 按 EC+LD 归一化分数排序, 互不冲突地组合 2-4 个位点
        pairs.sort(key=lambda x: -(x[2] + x[3]))
        used_sets = []
        for p in pairs:
            merged = False
            for s in used_sets:
                if p[0] in s[0] or p[1] in s[0]:
                    if len(s[0]) < 4:
                        s[0].update((p[0], p[1]))
                        s[1] += p[2] + p[3]
                    merged = True
                    break
            if not merged:
                used_sets.append([{p[0], p[1]}, p[2] + p[3]])
        for s, sc in used_sets[:args.top_combo]:
            muts = []
            for pos in sorted(s):
                muts.append(f"{PUF12_SEQ[pos]}{pos+1}"
                            f"{AA[stats['ligandmpnn']['x_mut'][pos]]}")
            combos.append({"combo": "+".join(muts), "n_mut": len(muts),
                           "ec_ld_score": round(sc, 3)})
        print(f"  高 EC/LD 位点对 {len(pairs)} 个 -> 组合 {len(combos)} 组")
    pd.DataFrame(combos).to_csv(
        os.path.join(args.outdir, "aice_multi_combos.csv"), index=False)

    print("[5/5] 输出")
    print(f"  {args.outdir}/aice_single_ranked.csv   (单突变提名, 已排序)")
    print(f"  {args.outdir}/aice_multi_combos.csv    (组合突变提名)")
    print("下一步: 取 top N 变体 -> bash step4_validate.sh")


if __name__ == "__main__":
    main()
