#!/usr/bin/env bash
# ============================================================================
# PUF12-AiCE 环境安装（4090 服务器, CUDA 12.x）
# 用法: bash setup_env.sh
# ============================================================================
set -euo pipefail

ENV_NAME=${ENV_NAME:-puf_aice}

# ---- conda 环境 -------------------------------------------------------------
if ! command -v conda &>/dev/null; then
  echo "未找到 conda，请先安装 Miniconda/Mambaforge"; exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -n "$ENV_NAME" python=3.10 -y
conda activate "$ENV_NAME"

# ---- PyTorch (CUDA 12.1, 适配 4090) -----------------------------------------
pip install torch --index-url https://download.pytorch.org/whl/cu121

# ---- 科学计算 + 结构工具 ------------------------------------------------------
pip install numpy scipy pandas biopython prody tqdm pydssp pyyaml ml_collections dm-tree filelock fsspec networkx sympy jinja2
# 官方 AiCE 用 mkdssp 4.4.7；pydssp 是纯 Python 近似。若想与论文完全一致:
#   conda install -c salilab dssp -y   # 提供 mkdssp 二进制

# ---- LigandMPNN (内含 ProteinMPNN 权重与 run.py) -----------------------------
if [ ! -d LigandMPNN ]; then
  git clone https://github.com/dauparas/LigandMPNN.git
fi
# 官方 get_model_params.sh 用 wget -q 静默下载, 国内直连极易"假卡住";
# 改用带进度/断点续传的下载脚本 (默认只下本管线必需的 2 个权重)
bash download_weights.sh

# ---- 官方 AiCE 仓库（可选, 用于 plink LD / pySCA EC 的完全复刻流程）-----------
if [ ! -d AiCE_official ]; then
  git clone https://github.com/ScorpioLea/AiCE.git AiCE_official || true
fi

# ---- Boltz-1（复合物建模与验证, 支持 蛋白+RNA）--------------------------------
pip install boltz -U || pip install boltz==1.*

echo
echo "=================================================="
echo " 环境安装完成: conda activate $ENV_NAME"
echo " 下一步: python step0_prepare.py"
echo "=================================================="
