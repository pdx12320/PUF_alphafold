#!/usr/bin/env bash
# ============================================================================
# Step 2: AiCE 核心步骤 1 —— 逆折叠模型大批量采样
#
# AiCE 论文设置 (Cell 2025, STAR Methods):
#   - 每个模型生成 10,000 条逆折叠序列
#   - 采样温度 (sequence diversity) = 0.5
#   - 模型: ProteinMPNN / LigandMPNN / (ESM-IF1 可选)
#
# 本脚本在同一结构上分别跑两个模型:
#   1) ligand_mpnn  —— RNA 链 R 作为配体上下文 (界面位点更可靠)
#   2) protein_mpnn —— 只看蛋白骨架 (核心/稳定性位点更可靠)
#
# 用法: bash step2_sample_mpnn.sh [n_seqs=10000] [temp=0.5]
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"

N=${1:-10000}
TEMP=${2:-0.5}
BATCH=50
NB=$(( N / BATCH ))
PDB=${PDB:-workdir/complex.pdb}
LMPNN_DIR=${LMPNN_DIR:-LigandMPNN}

[ -f "$PDB" ] || { echo "缺少 $PDB, 先运行 bash step1_fold_complex.sh"; exit 1; }
[ -d "$LMPNN_DIR" ] || { echo "缺少 LigandMPNN 目录, 先运行 bash setup_env.sh"; exit 1; }

for MODEL in ligand_mpnn protein_mpnn; do
  OUT="workdir/sample_${MODEL}"
  mkdir -p "$OUT"
  if [ "$MODEL" = ligand_mpnn ]; then
    CKPT="$LMPNN_DIR/model_params/ligandmpnn_v_32_010_25.pt"
    EXTRA="--checkpoint_ligand_mpnn $CKPT"
  else
    CKPT="$LMPNN_DIR/model_params/proteinmpnn_v_48_020.pt"
    EXTRA="--checkpoint_protein_mpnn $CKPT"
  fi

  echo "==> [$MODEL] 采样 ${N} 条序列, T=${TEMP}, 只设计蛋白链 A"
  echo "    (RNA 链 R 不是标准氨基酸, 会自动进入配体上下文, 无需额外参数)"
  python "$LMPNN_DIR/run.py" \
    --model_type "$MODEL" \
    $EXTRA \
    --seed 111 \
    --pdb_path "$PDB" \
    --chains_to_design "A" \
    --out_folder "$OUT" \
    --number_of_batches "$NB" \
    --batch_size "$BATCH" \
    --temperature "$TEMP" \
    --save_stats 1
done

echo
echo "采样完成:"
echo "  LigandMPNN : workdir/sample_ligand_mpnn/seqs/"
echo "  ProteinMPNN: workdir/sample_protein_mpnn/seqs/"
echo "下一步: python step3_aice_screen.py"
