# C388 local / nonlocal optimization：阶段性结果补充

2026-09-19 收录。来源：`c388_local_nonlocal_optimization_partial.zip`。本次补充不重新训练模型，不重筛 work 样本，不覆盖原仓库的 architecture、C295 或 WT-relative 分析。

**归档范围：本目录提供可浏览、可核验的结果快照，原始94文件压缩包未整体嵌入仓库。** 原包含全部脚本、输入缓存、14,784条预测、672行性能指标以及9组PNG/PDF；完整重跑需该原包。下方 CSV 为明确标识的原表摘录，数值原样复制。

[原包报告](SOURCE_README.md) · [比较指标摘录](results/comparison_metrics.csv) · [关键方案逐样本预测](results/key_predictions.csv) · [本次独立核验](INTEGRATION_AUDIT.json)

## 1. 固定研究问题

沿用原包的22个 work mutants：16个 preserved、6个 strongly reduced。WT仅作参考，C388编辑率为72.561%。正类定义：`ΔC388 = mutant − WT ≥ −15个百分点`。本终点允许最多15个百分点下降，不能称作编辑增强。

与原仓库不同阶段的样本范围、标签和验证方式分开报告，所有旧模型与结论保持原状。

## 2. 核心结果

| 特征与模型 | score协议 | 分组留出BA | 留一BA | 分组留出AUC |
|---|---|---:|---:|---:|
| Interface_CP / RF_depth2 | 固定0.5 | 0.937500 | 0.854167 | 0.947917 |
| Total_CP / LR_C0.1 | 固定0.5 | 0.687500 | 0.718750 | 0.947917 |
| Total_CP / LR_C1 | 固定0.5 | 0.750000 | 0.750000 | 0.927083 |
| Total_CP / RF_depth2 | 固定0.5 | 0.708333 | 0.708333 | 0.796875 |
| CP_summary / LR_C1 | 固定0.5 | 0.875000 | 0.791667 | 0.947917 |
| C388_pm0 / LR_C0.1 | 固定0.5 | 0.687500 | 0.687500 | 见CSV |
| Interface_CP / RF_depth2 | 训练内选score阈值 | 0.583333 | 0.854167 | 0.947917 |
| Interface_CP / nested_model_selection | 训练内选模型与score阈值 | 0.500000 | 0.770833 | 0.875000 |
| nested_best_available | 训练内选特征族、模型与score阈值 | 0.500000 | 0.656250 | 0.687500 |

以完整CSV为准。固定0.5的Interface_CP/RF_depth2：分组留出TP=14、TN=6、FP=0、FN=2，accuracy=0.909091、MCC=0.810093；留一TP=14、TN=5、FP=1、FN=2。

**解释边界：** 0.9375来自多候选探索性比较，仍有模型选择偏差。Total_CP/LR_C0.1的分组留出AUC也为0.947917，因此该比较的BA差异不能概括为所有指标上的优势。模型和阈值选择会明显改变分类表现；这里尚无独立新突变体测试。

## 3. Interface_CP的准确含义

保留原字段名，以便与历史结果对应。其唯一输入为：

```text
sqrt(sum((mutant protein–RNA CP − site-WT protein–RNA CP)^2) / (519 × 15))
```

这是全protein–RNA接触矩阵相对WT的RMS变化，分母沿用历史定义。缓存蛋白512个残基对应全长519的8–519位。**该指标不构成相邻repeat protein–protein界面或非局部蛋白通讯的直接证据。** `Total_CP`同样来自protein–RNA CP总和变化。

当前支持的表述：在这22个已入组突变体中，PR矩阵整体扰动是C388保留分类的一个探索性候选特征。局部RNA窗口尚未显示稳定的额外分类收益。CP关联、特征重要性和分类分数均不能独立确立因果机制。

## 4. 验证方式与剩余缺口

分组沿用P4、P5、P6、P7、P8、P4+P7。组合组与单位置组仍可能共享组成位置；验证覆盖位置组合留出，尚未严格隔离全部组成位置。多数测试组为单类，组内BA/AUC记缺失。P5有9个样本、2个假阴性，组内BA=0.75。

尚未完成mutation-centered >4/8/12/16/24、严格global >4/8/12/16/24，以及mutation→C388相关区域的protein–protein coupling。原包缺少WT及这22个mutants对应的完整PP CP矩阵、token/sequence映射和seed对应。旧S4/S12/S24使用`|i−j|≥k`、CP>0.5的接触数/蛋白长度，不能替代这些待计算特征。

## 5. 文件与核验

`results/comparison_metrics.csv`：从原672行性能表摘录12种配置、两种验证，共24行。保留模型与阈值协议。

`results/key_predictions.csv`：从原14,784行表摘录3种关键配置、两种验证，共132行；省略训练成员等长字段。config_id映射见`verify_results.py`。原包全表的训练/测试成员已在本次集成核验中检查。

`results/best_feature_set.csv`、`best_available_per_group_metrics.csv`、`frozen_cohort.csv`、`baseline_scalar_features.csv`、`baseline_feature_definitions.json`、`input_audit.json`、`validation_checks.json`直接保留原文件。`provenance.json`记录原分析来源；`INTEGRATION_AUDIT.json`另行记录本次集成，二者语境分开。

本次核对原包93个manifest条目全部通过；全部672组指标由原预测复算，最大绝对误差1.11e-16。原包状态仍为`full_task_complete=false`。

仅核验本快照的关键结果，无第三方依赖：

```bash
python c388_local_nonlocal_optimization/verify_results.py
```

恢复本地完整原包（须提供用户已有ZIP；脚本不会联网下载或自动提交）：

```bash
python c388_local_nonlocal_optimization/import_source_archive.py \
  /path/to/c388_local_nonlocal_optimization_partial.zip \
  --out ./c388_source_snapshot
cd c388_source_snapshot/c388_local_nonlocal_optimization
python -m pip install -r requirements.txt
python validate_available.py
```

完整重跑命令见原包README。完整原包是进一步扩展PP矩阵分析的输入基础；本次提交不包含新增PP计算或最终部署权重。
