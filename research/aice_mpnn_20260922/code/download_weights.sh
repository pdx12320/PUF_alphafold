#!/usr/bin/env bash
# ============================================================================
# LigandMPNN 权重下载（替代官方 get_model_params.sh: 带进度/断点续传/重试）
#
# 用法:
#   bash download_weights.sh          # 只下载本管线必需的 2 个权重 (~10 MB)
#   bash download_weights.sh --all    # 下载全部权重
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"

BASE="http://files.ipd.uw.edu/pub/ligandmpnn"
DIR="LigandMPNN/model_params"
mkdir -p "$DIR"

# 本管线必需 (step2_sample_mpnn.sh 用到)
FILES="proteinmpnn_v_48_020.pt ligandmpnn_v_32_010_25.pt"

# --all 时追加的其余模型
if [ "${1:-}" = "--all" ]; then
  FILES="$FILES proteinmpnn_v_48_002.pt proteinmpnn_v_48_010.pt \
ligandmpnn_v_32_005_25.pt ligandmpnn_v_32_020_25.pt \
ligandmpnn_sc_v_32_002_16.pt \
per_residue_label_membrane_mpnn_v_48_020.pt \
global_label_membrane_mpnn_v_48_020.pt solublempnn_v_48_020.pt"
fi

for f in $FILES; do
  if [ -s "$DIR/$f" ]; then
    echo "[已有] $f, 跳过"
    continue
  fi
  echo "[下载] $f"
  if command -v aria2c &>/dev/null; then
    aria2c -x 8 -c -d "$DIR" -o "$f" "$BASE/$f" || rm -f "$DIR/$f"
  else
    wget -c --show-progress --tries=10 --timeout=60 -P "$DIR" "$BASE/$f" \
      || { echo "[失败] $f (检查网络/代理)"; rm -f "$DIR/$f"; }
  fi
done

echo
ls -lh "$DIR"
echo "完成。继续: bash step0_prepare.py (或你中断前的下一步)"
