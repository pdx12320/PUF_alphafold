# 重现说明

所有命令在仓库根目录、Python 3.12环境运行；依赖版本见requirements.txt。脚本是按分析阶段保存的研究代码；没有隐式自动提交或发布动作。

## 已缓存数据重现（推荐）

`python tools/reproduce_current.py`读取以下已提供文件：

- architecture_validation/main_features.csv、all_features.csv
- architecture_validation/all_arrays.npz、all_residue_mapping.csv
- combined12/manifest.json

这些数据足以重现精确置换、局部候选发现、三种验证、重要性、图表、报告与最简模型。固定随机种子2026。运行将覆盖相应结果文件，建议先在独立checkout运行。

主要统计定义：24个construct的C(24,4)=10,626种标签排列；保持architecture内成功数量的360种条件排列；BH提供家族与全局校正。核心minimum pLDDT为“残基先平均，再取核心最小值”；高置信CP阈值作用于聚合矩阵，原CP summary阈值则逐model计算。定义详见特征字典。

## 其他分析阶段

| 命令 | 用途 |
|---|---|
| `python combined12/train.py` | 从两批缓存特征重训24种PUF12固定模型 |
| `python rf_tuning/tune.py` | 嵌套LOCO的RF参数搜索 |
| `python rf_tuning/report.py` | 调参结果报告/对照表 |
| `python architecture_validation/stats.py` | 精确置换、组内控制、CCR/contact分析 |
| `python architecture_validation/models.py` | LOCO/整组留出及敏感性分析 |

代码保留与当时结果的对应关系。旧模型或旧CCR结果的解释以最新报告为准。

## 原始AlphaFold输入

原始ZIP不入Git，以免重复提交大型MSA/template。将输入放到原目录：

- inputs/folds_2026_09_06_03_44.zip
- inputs/folds_2026_09_06_03_46.zip
- upload/folds_2026_09_07_02_48.zip
- upload/folds_2026_09_07_02_49.zip
- upload/91-96.zip

`docs/raw_input_inventory.json`标记当前文件大小、SHA256与ZIP有效性。原始03_46.zip当前本地副本曾被截断；最新分析复用了此前验证完成的14构建缓存。若重跑原始提取，应提供完整原ZIP，不能使用清单中zip_valid=false的截断副本。

`extended/extract.py`用于初始批；`new_batch/analyze.py`与`new_batch/map_repeats.py`现已限定第二批日期文件名，避免意外读取TRM ZIP；`trm_validation/analyze.py`读取91–96；`architecture_validation/extract.py`读取新批和TRM，随后prepare.py与原14种缓存合并。

`architecture_validation/prepare.py`中与旧CP特征的数值比对可验证衔接。不同repeat位置采用P1–P12编号，source采用R1–R8；不要将不同来源repeat的同一位置直接解释成相同残基。

## 文件核验与安全读取

`python tools/validate_snapshot.py`校验已归档文件的SHA256、主集标签及各验证预测覆盖范围。重现后浮点/序列化字节可能随环境变化，原SHA256用于已发布快照，不要求所有重跑二进制逐字节一致。

joblib模型只能读取可信来源文件。模型加载审计是可读性验证；泛化指标来自heldout_predictions.csv，不来自重新拟合的训练评分。
