# AiCE × PUF12：MPNN 反向折叠筛选与湿实验对照（2026-09-22）

模仿 **AiCE**（AI-informed Codon/sequence Engineering，高通量蛋白进化筛选框架）的思路，把 PUF12–APOE4 mRNA 复合物喂给反向折叠模型，对全蛋白 493 个位点做单突变扫描，再用既有两批湿实验数据（84 个变体记录）做回顾性对照，并给出第三代构建清单。

## 思路

1. **结构输入**：AF3 预测的 PUF12(493 aa) + 17nt APOE4 mRNA 上下文（`ACAUGGAGGACGUGCGC`，含目标位点 C388）复合物 → `inputs/complex.pdb`（链 A = 蛋白，链 R = RNA）。
2. **反向折叠采样**：只设计蛋白链，RNA 作为上下文。两个模型各采样 10,000 条序列（T=0.5）：
   - **ProteinMPNN**：把 RNA 链当固定配体骨架；
   - **LigandMPNN**：显式把 RNA 原子作为配体上下文（`num_ligand_res=89`，确认 RNA 参与特征化）。
3. **AiCE 式打分**：对每个位点，比较模型给 WT 氨基酸的概率 `f_wt` 与最优替代氨基酸的频率 `f_max`；`f_max >> f_wt` 视为"模型认为 WT 不是最优"的提名位点。识别三联体位点（每个 PUF repeat 的第 12/13/16 位）单独标记，不与非识别位点混排。
4. **湿实验对照**：把两批 84 条变体记录（C295/C388/C871 编辑率）与模型频率逐点算 Spearman 相关。
5. **第三代设计**：湿实验赢家骨架 × AiCE 非识别位点提名，叠加成 11 个构建（G0–G10）。

## 结果

### 1. AiCE 提名（`results/aice_single_ranked.csv`，493 位点全扫描）

- 共提名 **81 个单突变**（`results/recommended_mutations.csv`），其中 **17 个双模型共识**；
- **识别位点 0 提名**——模型认为 WT 的三联体（SYE/CRQ/NYQ…）已经最优，与 PUF 密码表一致；
- 共识核心：T350V、D374E、V294I、M458L、L274I、T278I、V366I、R93K、R45T 等；
- C388 几何区（N 端帽子，1–45 位）：R45T、S30A、L41R、N12S；界面：A438G、H392D；
- 组合突变提名见 `results/aice_multi_combos.csv`（top：T247K+Q283K+Q352E…）。

### 2. 与湿实验的相关性（`results/wetlab_vs_prediction.csv`，84 条记录）

| 对照 | Spearman ρ | p | 结论 |
|---|---|---|---|
| ProteinMPNN 联合频率 vs C388 保留率 | ≈0.23–0.27 | 弱 | 方向对但强度低 |
| LigandMPNN 联合频率 vs C388 保留率 | 0.013 | 0.911 | **完全无相关** |
| LigandMPNN 联合频率 vs C871/C388 | 0.101 | 0.391 | 无相关 |
| LigandMPNN 边缘频率 log vs C871/C388 | **0.289** | **0.012** | 弱显著：模型越不喜欢 → 旁观者编辑越低 |

关键反例：湿实验最优变体 **P8-R6-GVE**（C871/C388=0.153，两批最优，C388 保留 89%）在两个模型里的联合采样频率都是 **0**；VFQ/SYVIRR/ETD/LRD/VWH 同样为 0，NTQ 仅 0.0001。

**结论（分工定论）**：
- **识别三联体（repeat 第 12/13/16 位）的筛选归湿实验 + PUF 密码表**，MPNN 类模型在此不可信——它们把实验赢家全部判了死刑；
- **非识别位点（稳定性/亲和力/C388 几何）归 AiCE/MPNN 提名**，与识别位点正交叠加；
- 加 RNA 上下文（LigandMPNN）并没有救回识别位点的预测力，与 AiCE 论文自身对 LigandMPNN 偏弱的观察一致；
- 弱显著的 ρ=0.289 方向支持"亲和力气窗"假说：旁观者编辑高往往因为结合太弱，而不是识别错误。

### 3. 第三代构建清单（`results/gen3_constructs.fasta/.csv`）

湿实验赢家骨架 × AiCE 非识别位点提名，11 个构建，每个突变都经过 WT 序列断言校验：

| 构建 | 突变数 | 设计意图 |
|---|---|---|
| G0_P8GVE_ctrl | 2 | 对照：湿实验赢家骨架 |
| G1_GVE+R45T | 3 | C388 几何单点 |
| G2_GVE+core | 8 | 共识核心包（刻意避开 R8 的 V294I） |
| G3_GVE+R45T+core | 9 | 几何+核心 |
| G4_GVE+iface | 4 | 界面亲和力（A438G+H392D） |
| G5_GVE+geo | 5 | C388 几何包（N12S/S30A/L41R） |
| G6_GVE+all | 14 | 全叠加，风险最高，仅观察用 |
| G7_NTQ+core | 8 | 备选骨架 P9-NTQ（C388 保留 0.94） |
| G8_NTQ+iface | 4 | 备选骨架+界面 |
| G9_GVE+NTQ+core | 10 | 双选择性骨架叠加 |
| G10_SYVIRR+core | 7 | 高 C388 骨架（湿实验保留率 1.10）+核心包 |

注：R7 局部序列 252–257 野生型为 SYVIER（257 位本就是 R），故 SYVIRR 只需 E256R 单突变。

## 哪些 P / 残基重要

### Repeat 层面（湿实验 84 条记录按被突变的 repeat 分组，`results/repeat_importance.csv`）

| 分组 | Repeat | 依据 |
|---|---|---|
| **选择性改造主战场** | **P8、P9、P7** | 两批全部最优变体都出在这里：P8-GVE（C871/C388=0.153，全场最优）、P9-NTQ/NPS/NSG/NFP/NWP（保留率 0.88+）、P7-SYVIRR（保留率 0.98）；其中 P8 只有 GVE 一个变体，远未饱和 |
| **C295 旁观者开关** | **P3** | P3 被突变时 C295 中位数直接掉到 0；同时出了 VFQ（C871/C388=0.123），但保留率中位 0.67，属高风险高收益 |
| **敏感别动** | **P2、P5** | 突变后 C388 保留率中位仅 0.69/0.73，且几乎不出选择性改善；P5 的 14 个变体破识别最多 |
| **耐受但无收益** | P1、P4、P6、P10、P11 | 保留率 0.80–0.93，动了不坏死但也不改善 C871/C388（中位 0.69–0.87） |
| **空白** | P12 | 两批均无数据，且 R12 三联体在 439 位、远离常规 repeat 间距，待补 |

注意样本不均衡：P5 有 14 个变体、P8 只有 1 个，分组中位数的置信度差异很大。

### 残基层面（AiCE 双模型）

**该动的**（17 个双模型共识，按功能分）：
- 稳定性/亲和力核心：T350V、D374E、V294I、M458L、L274I、T278I、L166I、V366I、R93K、E351L、R462L
- C388 几何区（1–45 位，脱氨酶融合端朝向）：**R45T** 是唯一位于几何区的共识，另 S30A/L41R/N12S 为单模型提名
- 界面：A438G、H392D（H212D 系列为单模型提名）

**别动的**（模型 f_wt≈1.0，即"WT 已最优"）：各 repeat 交界处的甘氨酸 G35/G71/G107/G179/G251/G287、界面 V290/V362、P348、L101 等——这批是 PUF 骨架的保守铰链，与"敏感别动"的 P2/P5 结论互洽。

**模型信号要反着看的**：36 个识别位点（repeat 第 12/13/16 位）里 27 个 f_wt<0.5，模型普遍"想改"三联体，但湿实验证明 WT 三联体最优、实验赢家模型频率全为 0——识别位点完全交给湿实验+密码表。

## 后续计划

1. **第三代湿实验**：按 G0–G10 构建、测 C295/C388/C871 三点编辑率，检验"识别归实验、非识别归 AiCE"的叠加策略是否可加性成立；优先补 P8（除 GVE 外的疏水/极性三联体）和 P12 空白；
2. **第三轮迭代**：第三代数据回来后重新做湿实验 vs 预测对照，更新核心包/几何包的取舍；
3. **模型侧**：若继续用 MPNN 类模型，只用于非识别位点；识别位点的计算探索改用几何/能量类方法或直接饱和突变。

## 复现

```bash
cd code
bash setup_env.sh                 # conda 环境: torch/prody/pandas/biopython 等
bash download_weights.sh          # ProteinMPNN + LigandMPNN 权重
# LigandMPNN 源码需单独 clone: https://github.com/dauparas/LigandMPNN
# （已知兼容性问题见下方"限制"）
cp ../inputs/complex.pdb workdir/ # 或 bash step0 后自行折叠
bash run_all.sh                   # step0->step3 全流程
python make_gen3.py               # 重新生成第三代构建清单
```

## 限制与诚实声明

- 相关性分析是**回顾性**的：84 条湿实验记录并非为验证 MPNN 而设计，位点分布不均；
- ρ=0.289 仅边缘显著，不能当作预测工具使用，只作为"气窗"假说的方向性证据；
- 服务器无外网环境下 LigandMPNN 依赖需要逐个手装（`ml_collections`、`dm-tree`、`ProDy`、`pydssp`），且其捆绑的 openfold 代码与 numpy≥1.24 不兼容（`np.int` 已移除），需补丁或降版本；
- G6（14 突变）为高风险观察构建，预期可能不折叠/不表达；
- 本目录所有"提名"均为计算候选，**没有独立新实验确证**之前不应视为有效设计。

## 文件清单

| 路径 | 内容 |
|---|---|
| `code/` | 自写流水线 step0–step5、run_all.sh、环境脚本、make_gen3.py |
| `inputs/complex.pdb` | AF3 PUF12+17nt RNA 复合物（链 A 蛋白 / 链 R RNA） |
| `inputs/design_config.json`、`repeat_parse.csv` | repeat 解析与设计配置 |
| `results/aice_single_ranked.csv` | 493 位点全扫描双模型打分表 |
| `results/recommended_mutations.csv` | 81 个提名突变（含优先级分类） |
| `results/aice_multi_combos.csv` | 组合突变提名 |
| `results/wetlab_vs_prediction.csv` | 84 条湿实验记录 × 双模型频率对照总表 |
| `results/repeat_importance.csv` | 按 repeat 分组的湿实验敏感性统计 |
| `results/mpnn_samples_{proteinmpnn,ligandmpnn}.fa.gz` | 各 10,000 条采样序列（gzip） |
| `results/gen3_constructs.{fasta,csv}` | 第三代 11 个构建 |
