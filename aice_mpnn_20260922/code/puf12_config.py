# -*- coding: utf-8 -*-
"""
PUF12 / APOE4 项目的共享配置与序列解析模块
序列来源: APOE4 mRNA+PUF12 Protein.docx
所有 step 脚本从这里 import，保证位点编号一致。
"""
import json
import os

# ---------------------------------------------------------------------------
# 序列
# ---------------------------------------------------------------------------
PUF12_SEQ = (
    "GRSRLLEDFRNNRYPNLQLREIAGHIMEFSQDQHGSYFIELKLERATPAERQLVFNEILQAAYQLMVDVF"
    "GSYVIRKFFEFGSLEQKLALAERIRGHVLSLALQMYGCRVIQKALEFIPSDQQNEMVRELDGQVFALSTHPY"
    "GSYVIERILEHCLPDQTLPILEELHQHTEQLVQDQYGSYVIEHVLEHGRPEDKSKIVAEIRGNVLVLSQHKF"
    "ACRVVQKCVTHASRTERAVLIDELDGQVFALSTHPYGSYVIERILEHCLPDQTLPILEELHQHTEQLVQDQY"
    "GSYVIEHVLEHGRPEDKSKIVAEIRGNVLVLSQHKFANYVVQKCVTHASRTERAVLIDKIRPHTEQLVQDQY"
    "GCRVVQHVLEHGRPEDKSKIVAEIRGNVLVLSQHKFASYVIRKCVTHASRTERAVLIDEVCTMNDGPHSALY"
    "TMMKDQYACRVVQKMIDVAEPGQRKIVMHKIRPHIATLRKYTYGKHILAKLEKYYMKNGVDLG"
)

APOE4_MRNA = (
    "AUGAAGGUUCUGUGGGCUGCGUUGCUGGUCACAUUCCUGGCAGGAUGCCAGGCCAAGGUGGAGCAAGCGG"
    "UGGAGACAGAGCCGGAGCCCGAGCUGCGCCAGCAGACCGAGUGGCAGAGCGGCCAGCGCUGGGAACUGGCAC"
    "UGGGUCGCUUUUGGGAUUACCUGCGCUGGGUGCAGACACUGUCUGAGCAGGUGCAGGAGGAGCUGCUCAGCU"
    "CCCAGGUCACCCAGGAACUGAGGGCGCUGAUGGACGAGACCAUGAAGGAGUUGAAGGCCUACAAAUCGGAAC"
    "UGGAGGAACAACUGACCCCGGUGGCGGAGGAGACGCGGGCACGGCUGUCCAAGGAGCUGCAGGCGGCGCAGG"
    "CCCGGCUGGGCGCGGACAUGGAGGACGUGCGCGGCCGCCUGGUGCAGUACCGCGGCGAGGUGCAGGCCAUGC"
    "UCGGCCAGAGCACCGAGGAGCUGCGGGUGCGCCUCGCCUCCCACCUGCGCAAGCUGCGUAAGCGGCUCCUCC"
    "GCGAUGCCGAUGACCUGCAGAAGCGCCUGGCAGUGUACCAGGCCGGGGCCCGCGAGGGCGCCGAGCGCGGCC"
    "UCAGCGCCAUCCGCGAGCGCCUGGGGCCCCUGGUGGAACAGGGCCGCGUGCGGGCCGCCACUGUGGGCUCCC"
    "UGGCCGGCCAGCCGCUACAGGAGCGGGCCCAGGCCUGGGGCGAGCGGCUGCGCGCGCGGAUGGAGGAGAUGG"
    "GCAGCCGGACCCGCGACCGCCUGGACGAGGUGAAGGAGCAGGUGGCGGAGGUGCGCGCCAAGCUGGAGGAGC"
    "AGGCCCAGCAGAUACGCCUGCAGGCCGAGGCCUUCCAGGCCCGCCUCAAGAGCUGGUUCGAGCCCCUGGUGG"
    "AAGACAUGCAGCGCCAGUGGGCCGGGCUGGUGGAGAAGGUGCAGGCUGCCGUGGGCACCAGCGCCGCCCCUG"
    "UGCCCAGCGACAAUCACUGA"
)

# ---------------------------------------------------------------------------
# PUF 识别密码表 (triplet -> base)
# 经典密码 (Wang et al. / Filipovska et al.): 每个 repeat 的第 12/13/16 位
# 可按需要扩展（如 TDA 识别的宽松密码、双碱基识别 repeat 等）
# ---------------------------------------------------------------------------
PUF_CODE = {
    "SYE": "G", "SNE": "G", "CKE": "G",
    "SYR": "C", "SNR": "C",
    "CRQ": "A", "CKQ": "A", "SQQ": "A",
    "NYQ": "U", "NKQ": "U", "NYC": "U",
}

# 带上下文的 17nt RNA (mRNA 374-390, 含 C388=第15位), 与 AF3 wildtype_seed1
# 任务保持一致, 用于变体验证 (Boltz/AF3) 使指标可与 WT baseline (ipTM 0.91) 对比
RNA_CONTEXT_17 = "ACAUGGAGGACGUGCGC"
C388_NUCLEOTIDE_INDEX = 15      # 在 17nt RNA 中的位置 (1-based)

# C388 毗邻区: 结构上 C388 贴近 PUF 的 N 端帽子区 (6A 内残基 F9/R10/N12/P15/
# H34/Y37/F38/L41), 脱氨酶融合端的位置几何由这一带决定
C388_ZONE_1BASED = (1, 45)

# 同义密码: 识别同一碱基所允许的 triplet 集合（筛选时用于保护特异性）
SYNONYMOUS_CODES = {}
for _trip, _base in PUF_CODE.items():
    SYNONYMOUS_CODES.setdefault(_base, set()).add(_trip)


def parse_repeats(seq=PUF12_SEQ, max_gap=48, min_gap=28):
    """扫描序列中所有 PUF 识别三联体 (位置 i, i+1, i+4)，按间距筛出连续的 repeat 链。
    返回 list[dict]: repeat 序号、triplet、三联体残基 0-based 索引、识别碱基。
    """
    cands = []
    for i in range(len(seq) - 4):
        trip = seq[i] + seq[i + 1] + seq[i + 4]
        if trip in PUF_CODE:
            cands.append((i, trip))
    sel = []
    for pos, trip in cands:
        if not sel or min_gap <= pos - sel[-1][0] <= max_gap:
            sel.append((pos, trip))
    repeats = []
    for k, (pos, trip) in enumerate(sel, 1):
        repeats.append({
            "repeat": k,
            "triplet": trip,
            "base": PUF_CODE[trip],
            "pos12": pos,       # 0-based, repeat 内第 12 位
            "pos13": pos + 1,   # 第 13 位
            "pos16": pos + 4,   # 第 16 位
        })
    return repeats


def target_rna(repeats):
    """PUF repeat 从 N 到 C 依次结合 RNA 的 3' 到 5'，故反转得到 5'->3' 靶序列。"""
    return "".join(r["base"] for r in repeats)[::-1]


def specificity_positions(repeats):
    """所有识别位点的 0-based 索引（采样后筛选时必须保护的位置）。"""
    out = []
    for r in repeats:
        out += [r["pos12"], r["pos13"], r["pos16"]]
    return out


def write_boltz_yaml(protein_seq, rna_seq, path, name="puf12_complex"):
    """Boltz-1/2 输入: 蛋白链 A + RNA 链 R。"""
    text = (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {protein_seq}\n"
        "  - rna:\n"
        "      id: R\n"
        f"      sequence: {rna_seq}\n"
    )
    with open(path, "w") as f:
        f.write(text)
    return path


def save_config(outdir):
    repeats = parse_repeats()
    rna = target_rna(repeats)
    idx = APOE4_MRNA.find(rna)
    cfg = {
        "puf12_length": len(PUF12_SEQ),
        "n_repeats": len(repeats),
        "repeats": repeats,
        "rna_target_5to3": rna,
        "rna_context_17nt": RNA_CONTEXT_17,
        "c388_nucleotide_index": C388_NUCLEOTIDE_INDEX,
        "c388_zone_1based": C388_ZONE_1BASED,
        "mrna_1based_start": idx + 1 if idx >= 0 else None,
        "specificity_positions_0based": specificity_positions(repeats),
    }
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "design_config.json"), "w") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    return cfg


if __name__ == "__main__":
    cfg = save_config(".")
    print(json.dumps(cfg, indent=2, ensure_ascii=False))
