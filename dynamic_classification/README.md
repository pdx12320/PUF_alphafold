# 三个位点动态阈值分类（2026-09-08）

取消 C388 >50% 入组限制，使用本批全部 30 个有实验数据的突变体。WT 仅作为参考，不进入训练。C388、C871、C295 独立建模；编辑率达到位点分界定义为高活性。

## 主要结果

以下为训练折内动态选择模型分数阈值后的平衡准确率（BA）。编辑率分界和模型先经本批数据探索性筛选，因此仍有选择偏差，尚无独立测试验证。

| 位点 | 候选编辑率分界 | 模型 | 留一 BA | 按突变位置组合留出 BA |
|---|---:|---|---:|---:|
| C388 | 30% | 接触面 CP + 逻辑回归 | 76.0% | 78.0% |
| C871 | 45% | ContactSeek CCR + 随机森林 | 78.5% | 79.4% |
| C295 | 12.83675% | ContactSeek CCR + 随机森林 | 82.5% | 52.5% |

C388、C871 可作下一批验证候选；C295 跨突变位置泛化不足。C388 高活性用于保留目标编辑，C871/C295 低活性可用于筛选降低相应编辑的构建。每个位点模型都使用全部 30 个突变体。

编辑率的等价分界区间分别为 C388 (25.33%,32.484%]、C871 (41.619%,45.5995%]、C295 (12.468%,13.2055%]。区间内所有分界产生同样的本批标签，不能解读成精确生物学界限。

## 文件导航

- [完整 HTML 报告](Dynamic_threshold_classification_report.html)：下载后用浏览器打开。
- [混淆矩阵](results/confusion_matrices.png)；[PDF](results/confusion_matrices.pdf)。
- [全部探索性阈值扫描](results/exploratory_threshold_scan.csv)：编辑率分界×6 个固定模型；其最佳分数阈值依据全体 OOF 标签选择，属于探索性上限。
- [训练内分数阈值验证](results/training_only_score_threshold_validation.csv)。
- [逐构建验证预测](results/validated_predictions.csv)、[候选标签](results/candidate_labels.csv)。
- [模型特征权重](results/candidate_feature_weights_exploratory.csv)：全数据拟合，仅供探索。
- `results/candidate_models.joblib`：探索性模型、WT CP、接触筛选与 CCR 定义。
- `inputs/`：已聚合的实验、CP 与结构特征输入，支持直接复现。
- `sources/`、`source_provenance.json`、`source_checksums.json`：参考源码及版本来源。
- `validation_checks.json`：180 条外层预测唯一性、分类判定和指标复核。

## 方法

借鉴 [PUF_alphafold 联合阈值方法](../c388_threshold/README.md) 和 ContactSeek Top3/CCR。扫描编辑率 5%–95%（步长 5 个百分点）及相邻观测均值中点，合并等价标签，要求每类至少 5 个构建。比较总 CP、接触面 CP、S4/S12/S24、19 项 CP+结构、ContactSeek Top3、CCR。

LR 使用 C=0.1 和 balanced 权重；RF 使用 100 棵树、深度 2、叶节点至少 3、balanced 权重。每个位点选探索性 BA 最佳组合，再用训练集内层 3 折预测选分数阈值，外层分别做留一构建和突变位置组合留出。接触残基、CCR 和标准化在相应训练折中拟合。未联合嵌套选择编辑率分界和模型，不能消除事后候选筛选偏差。

93 个构建×位点组合（含 WT），97 个 seed 任务；每 seed 的 5 个 CP 一致，未计成独立实验样本。蛋白裁去前 7 aa，使用 512×15 蛋白–RNA CP。WT 每个位点仅一个 seed；AYE_C871 的 seed 敏感性需复核。原表部分零值代表缺少 C→T 记录，按原表保留，未做排除这些记录的敏感性验证。

## 复现

Python 3.12。从仓库根目录运行：

```bash
python -m pip install -r dynamic_classification/requirements.txt
cd dynamic_classification
python classify.py
python report.py
```

已包含建模所需输入，无需重新运行 AlphaFold。`extract.py` 用于从原始压缩包重提取，需提供同级 `upload/` 和前次 `work_cp_modeling/` 汇总文件；这些原始输入不随本目录重复上传。`report.py` 会重建 HTML、图和上级目录的结果 ZIP，并生成简版 README；本归档 README 提供更完整的方法说明。
