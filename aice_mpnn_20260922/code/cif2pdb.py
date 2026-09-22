# -*- coding: utf-8 -*-
"""
mmCIF -> PDB 转换器 (Boltz / AlphaFold3 输出通用)
保留蛋白链和 RNA 链, 蛋白链命名为 A, RNA 链命名为 R (LigandMPNN 要求)。
用法: python cif2pdb.py input.cif output.pdb
"""
import sys
import warnings

warnings.filterwarnings("ignore")

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa


def main():
    cif, out = sys.argv[1], sys.argv[2]
    s = MMCIFParser(QUIET=True).get_structure("cx", cif)[0]

    # 识别链类型: 含标准氨基酸残基 -> 蛋白; 含 RNA 残基 -> RNA
    RNA_NAMES = {"A", "U", "G", "C", "RA", "RU", "RG", "RC",
                 " DA", " DT", " DG", " DC"}
    prot, rna = [], []
    for ch in s:
        names = {r.resname.strip() for r in ch}
        if any(is_aa(r) for r in ch):
            prot.append(ch)
        elif names & {n.strip() for n in RNA_NAMES}:
            rna.append(ch)

    if not prot:
        raise SystemExit("未找到蛋白链!")

    # 重命名: 第一条蛋白链 -> A, 第一条 RNA 链 -> R
    prot[0].id = "A"
    if rna:
        rna[0].id = "R"
    else:
        print("警告: 未找到 RNA 链, LigandMPNN 将没有配体上下文!")

    class Keep(Select):
        def accept_chain(self, chain):
            return chain.id in ("A", "R")

        def accept_residue(self, res):
            return res.id[0] == " "  # 去掉水/离子等 HETATM

    io = PDBIO()
    io.set_structure(s)
    io.save(out, select=Keep())

    from Bio.PDB import PDBParser
    st = PDBParser(QUIET=True).get_structure("c", out)
    for ch in st.get_chains():
        n = len([r for r in ch if r.id[0] == " "])
        print(f"chain {ch.id}: {n} residues")
    print("written:", out)


if __name__ == "__main__":
    main()
