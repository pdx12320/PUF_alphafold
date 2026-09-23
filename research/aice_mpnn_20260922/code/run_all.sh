#!/usr/bin/env bash
# ============================================================================
# PUF12 × AiCE 一键运行
#
#   bash run_all.sh smoke    # 冒烟测试: 每模型 200 条, ~10 分钟
#   bash run_all.sh          # 正式运行: 每模型 10000 条 (AiCE 论文设置), ~1-2 小时
#
# 内置: LigandMPNN 代码 + 两个模型权重 + AF3 复合物结构 (workdir/complex.pdb)
# 全程不需要外网 (首次建 conda 环境装 python 包除外)
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
MODE=${1:-full}
ENV=puf_aice
PIP_INDEX=${PIP_INDEX:-https://pypi.tuna.tsinghua.edu.cn/simple}

echo "==> [1/5] 检查 conda 环境 ($ENV)"
command -v conda >/dev/null || { echo "未找到 conda"; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^${ENV} "; then
  echo "    创建环境..."
  conda create -n $ENV python=3.10 -y -q
fi
conda activate $ENV

# 缺啥装啥 (可重复执行, 已有则秒过)
if ! python -c "import numpy, scipy, pandas, Bio, prody, ml_collections, tree, torch" 2>/dev/null; then
  echo "    安装 python 依赖 (镜像: $PIP_INDEX)..."
  pip install -q -i "$PIP_INDEX" numpy scipy pandas biopython prody pydssp pyyaml ml_collections dm-tree filelock fsspec networkx sympy jinja2
  python -c "import torch" 2>/dev/null || \
    pip install -q torch --index-url https://download.pytorch.org/whl/cu121
fi

python - <<'EOF'
import torch
print(f"    torch {torch.__version__} | CUDA 可用: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print("    GPU:", torch.cuda.get_device_name(0))
EOF

echo "==> [2/5] 检查权重与结构"
for f in LigandMPNN/model_params/proteinmpnn_v_48_020.pt \
         LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt \
         workdir/complex.pdb; do
  [ -s "$f" ] || { echo "缺少 $f —— 压缩包不完整, 请重新解压"; exit 1; }
done
echo "    OK"

echo "==> [3/5] step0: 解析 PUF12 repeat / 生成设计配置"
python step0_prepare.py --outdir workdir

N=10000
[ "$MODE" = "smoke" ] && N=200
echo "==> [4/5] step2: MPNN 采样 (每模型 $N 条, T=0.5)"
bash step2_sample_mpnn.sh $N 0.5

echo "==> [5/5] step3: AiCE 筛选 (单突变 + 组合突变)"
python step3_aice_screen.py --pdb workdir/complex.pdb

cat <<'DONE'

============================================================
 全部完成! 查看结果:
   results/aice_single_ranked.csv   单突变提名 (按 rank_score 排序)
   results/aice_multi_combos.csv    组合突变提名

 优先送实验: consensus=True 且 code_safe 的位点;
   near_c388_zone=True 的位点影响 C388 几何选择性

 可选结构验证 (需外网装 boltz):
   python step4_make_yamls.py --top 20 && bash step4_validate.sh
============================================================
DONE
