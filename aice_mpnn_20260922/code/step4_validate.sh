#!/usr/bin/env bash
# ============================================================================
# Step 4b: Boltz-1 批量重折叠变体 + 收集 ipTM 打分
#
# 判断标准 (建议):
#   - ipTM >= 0.6 且不低于 WT 太多: 复合物结合构象保持
#   - 蛋白 pLDDT 保持: 折叠稳定
#   - 最终以 ipTM 变化量 (Δ vs WT) 排序
#
# 用法: bash step4_validate.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"

YDIR=${YDIR:-workdir/validate_yamls}
ODIR=${ODIR:-workdir/validate_boltz}
mkdir -p "$ODIR"

shopt -s nullglob
YAMLS=("$YDIR"/*.yaml)
[ ${#YAMLS[@]} -gt 0 ] || { echo "缺少 YAML, 先运行 python step4_make_yamls.py"; exit 1; }

for y in "${YAMLS[@]}"; do
  name=$(basename "$y" .yaml)
  if ls "$ODIR/boltz_results_${name}/predictions/${name}/"*model_0.cif >/dev/null 2>&1; then
    echo "[跳过] $name 已完成"
    continue
  fi
  echo "[Boltz-1] $name"
  boltz predict "$y" --out_dir "$ODIR" \
    --recycling_steps 3 --diffusion_samples 1 \
    --accelerator gpu --devices 1 --use_msa_server \
    || boltz predict "$y" --out_dir "$ODIR" \
         --recycling_steps 3 --diffusion_samples 1 \
         --accelerator gpu --devices 1
done

python step5_collect_results.py --boltz_dir "$ODIR" \
  --ranked results/aice_single_ranked.csv \
  --out results/final_ranked.csv

echo
echo "最终排名: results/final_ranked.csv"
