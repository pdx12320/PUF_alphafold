# C388 local / nonlocal optimization — 阶段性结果，未全部完成

固定使用附件 work_mutant_selection.csv 中标记的22个work mutants，未重新筛选。WT只作参考（72.561%）。主终点为 ΔC388≥−15 pp：16 preserved / 6 strongly reduced。

已完成：原动态流程完整重跑；固定终点的8个baseline族、5个RNA窗口、WT Top5/10/20/30×4窗口、3个已有global标量及5个局部+已有global CP3组合。每族比较 balanced LR C=0.1/1 和 RF 200树、深度2/3、叶≥2。仅最后三项自适应策略在训练折内选择feature family。新输出与原分析完全分离。

未完成：mutation-centered >4/8/12/16/24、严格global >4/8/12/16/24、mutation→C388 protein–protein coupling，以及这些特征的组合、解释与真正的coupling heatmap。唯一剩余输入缺口是 **WT及这22个mutants的C388完整protein–protein CP矩阵（含token/sequence编号和seed对应），或AF3原始ZIP**。不能由protein–RNA矩阵反推出这些值。

## Baseline复现

重新执行了每一个outer-training fold内的动态标签/特征/模型/score阈值选择。44条外层预测全部匹配：True，最大分数误差 1.11e-16。全训练探索性赢家仍为Total_CP、−15 pp。

| validation | balanced_accuracy | MCC | AUC |
|---|---|---|---|
| LOOCV | 0.5729 | 0.1550 | 0.7708 |
| mutation_group_out | 0.8583 | 0.7258 | 0.6833 |

这些动态BA并非固定−15 pp的Total_CP成绩。固定终点的Total_CP如下：

| model | score_policy | balanced_accuracy_group_out | balanced_accuracy_loocv |
|---|---|---|---|
| LR_C0.1 | fixed_0.5 | 0.6875 | 0.7188 |
| LR_C0.1 | training_only_threshold | 0.5000 | 0.7188 |
| nested_model_selection | training_only_threshold | 0.5000 | 0.7188 |

## 暂时的已测结果

按group-out BA、group-out MCC、LOOCV BA排序。下表使用常规固定score=0.5；为简洁每个family仅展示外层排行榜中最佳参数，故整个表是探索性比较，不能据此把winner当作独立验证的已选模型。所有参数的完整结果保存在CSV。不能仅凭此表宣称泛化改善。

| feature_set | model | balanced_accuracy_group_out | MCC_group_out | balanced_accuracy_loocv |
|---|---|---|---|---|
| Interface_CP | RF_depth2 | 0.9375 | 0.8101 | 0.8542 |
| CP_summary | LR_C1 | 0.8750 | 0.6708 | 0.7917 |
| ContactSeek_Top3 | LR_C1 | 0.7500 | 0.4629 | 0.7500 |
| Total_CP | LR_C1 | 0.7500 | 0.4629 | 0.7500 |
| ContactSeek_CCR | LR_C1 | 0.7500 | 0.4629 | 0.6979 |
| WT_Top5_pm0 | LR_C0.1 | 0.7188 | 0.4183 | 0.7500 |
| WT_Top5_pm1 | LR_C1 | 0.7188 | 0.4183 | 0.6042 |
| protein_RNA_CP | LR_C1 | 0.6875 | 0.3750 | 0.8125 |
| C388_pm0 | LR_C0.1 | 0.6875 | 0.3750 | 0.6875 |
| C388_pm1 | LR_C0.1 | 0.6875 | 0.3750 | 0.6875 |

训练内选参数和score阈值的结果独立列出。`nested_best_local`每折选择5种窗口，`nested_best_WT_contact`选择WT集合/窗口，`nested_best_available`选择所有已测family；外层标签始终固定−15 pp。这套更复杂的选择策略在当前样本上不稳定：

| feature_set | balanced_accuracy_group_out | MCC_group_out | balanced_accuracy_loocv |
|---|---|---|---|
| Interface_CP | 0.5000 | 0.0000 | 0.7708 |
| Total_CP | 0.5000 | 0.0000 | 0.7188 |
| nested_best_available | 0.5000 | 0.0000 | 0.6562 |
| nested_best_local | 0.4062 | -0.2433 | 0.6979 |
| nested_best_WT_contact | 0.4062 | -0.2433 | 0.3021 |

## 原请求的13个问题

1. C388 local CP 是否优于 Total_CP？默认score=0.5下，最佳已测局部窗口 C388_pm0/LR_C0.1 的 group-out BA=0.6875，Total_CP/LR C=0.1为0.6875，差值 +0.00 pp。这是全数据排行榜的探索性候选比较；不能视为独立模型选择评估。训练内窗口/模型/score阈值选择策略另列，未建立可靠改善。
2. target only、±1（LR C=0.1/1）、±2（LR C=1）在指定BA/MCC/LOOCV排序指标上并列；不能宣称唯一最佳窗口。若按窗口更窄优先，target only是更简洁候选。窗口均在 RNA 末端截断。
3. Top5/10/20 的最佳已测组合：WT_Top5_pm0，group-out BA=0.7188。这也包含窗口比较，详见全部表。
4. mutation-centered >4/8/12/16/24：未计算，缺少完整 protein–protein CP。
5. global 与 mutation-centered 比较：尚不能回答。已有 S4/S12/S24 是 |i-j|>=k、CP>0.5 的接触数/蛋白长度，不等同于本轮要求的严格 >k 或按合格pair数归一化密度。
6. 局部+已有global CP3 的最佳已测组合是 C388_pm1_plus_existing_CP3，group-out BA=0.6875；不能替代局部+mutation-centered nonlocal。
7. mutation→C388 contact CP 的额外信息：未计算，缺少 protein–protein CP。
8. 默认0.5下按指定排序的暂时最佳已测候选 Interface_CP/RF_depth2：LOOCV BA=0.8542。此处最佳来自多候选排行榜，不是独立验证的最终赢家。
9. 同一暂时最佳已测族的 group-out BA=0.9375。
10. 对固定Total_CP/LR C=0.1、同样score=0.5：group-out BA 差 +25.00 pp，LOOCV 差 +13.54 pp。不同阈值协议不混在同一差值中。
11. 暂时最佳族相对同协议 Total_CP 的group-out 错→对：P5-R6-SWD, P5-R6-TWD, P6-R7-QFW, P6-R7-QRW, P6-R7-TRQ, P6-R7-WFD, P7-R5-SYVIRR, P8-R6-GVE。对→错：无。
12. 最难预测的group：见下方各组错误数。单类小组的 BA/AUC 未定义，不能用填零的BA排序。
13. 最终最小 feature set：尚不能推荐。mutation-centered、严格global和communication特征未取得；候选赢家及训练内选择策略只是现有数据的阶段性结果。

| mutation_group | n | accuracy | balanced_accuracy | FP | FN |
|---|---|---|---|---|---|
| P4 | 2 | 1.0000 | nan | 0 | 0 |
| P4+P7 | 4 | 1.0000 | 1.0000 | 0 | 0 |
| P5 | 9 | 0.7778 | 0.7500 | 0 | 2 |
| P6 | 5 | 1.0000 | nan | 0 | 0 |
| P7 | 1 | 1.0000 | nan | 0 | 0 |
| P8 | 1 | 1.0000 | nan | 0 | 0 |

## 定义与边界

- 用户确认C388为 `ACAUGGAGGACGUGC` 第15位。target/±1/±2/±3/±4实际使用15、14–15、13–15、12–15、11–15，不存在下游nt，不补0。
- 缓存蛋白512残基对应全长519的8–519位；序列差异重新核验，mutation_mapping_audit.csv同时提供两种编号。独立样本是construct，seed先等权平均CP；非线性local汇总施加于平均矩阵。
- WT target最大CP=0.09。>0.1/0.2/0.3/0.5集合均为空，输出缺失值并跳过模型。target仅16个非零接触，因此Top20/30包含零CP残基；并列按序列位置排序，不能都称为强接触。
- local模型6维：Δsum、Δmax、Δtop5 mean、Δtop10 mean、mean|ΔCP|、RMS ΔCP。top-k指各construct最高k个CP值的均值，再减WT，不是最高k个Δ。
- WT集合模型4维：Δsum、RMS ΔCP、下降比例、上升比例。集合一律按WT target定义，固定用于不同窗口和所有mutants。
- baseline CP_summary保留原14维，CCR保留原9维以公平复现；新增特征族≤6维。原Total_CP是全protein–RNA CP总和，不是全蛋白CP总和。
- preprocessing和ContactSeek选择仅在对应训练折。内层为原流程3折GroupKFold；score阈值按BA、MCC、接近0.5优先选。新模型参数按inner BA/MCC后固定配置顺序，adaptive family并列按名称，不使用外层成绩打破平局。
- group-out保持原P4/P5/P6/P7/P8/P4+P7位置组合分组；P4+P7与单P4/P7仍可能共享位置。这不是严格所有组成位置隔离的验证。
- 固定score=0.5和训练内score阈值均完整输出、分开报告。分数未经校准。固定候选的排行榜及family赢家是在多个模型上探索比较，不声称显著改善或可部署。图01/02/03/04每族显示探索性最佳默认0.5候选；图05最后两柱采用训练内选择/阈值，其余为固定候选，属于协议对照。
- feature_importance.csv是暂时最佳族的训练折系数或RF impurity importance；不同importance类型不得直接混合排序。local_interpretation.csv是事后描述关联，不参与特征筛选。
- 04图展示已有>=k global摘要，08图展示RNA-local ΔCP；文件标题明确注明，二者不冒充缺少矩阵的请求结果。

## 复现

在安装requirements.txt所列版本后运行：

```text
python run_baselines.py --stage reproduce
python run_baselines.py --stage fixed --include-local
python report_available.py
```

baseline_reference保留附件的所需输入、ContactSeek源码及许可、原classify.py和原C388结果。run_baselines.py只写本目录results，不执行原脚本的输出入口。复现允许覆盖本新目录的结果，原仓库与附件不变。

固定阈值的逐构建预测、9组可生成的PNG/PDF、原baseline复现误差、定义和暂时最佳族均已保存；缺少的nonlocal结果不生成空白CSV冒充完成。补齐PP矩阵后仍需完成其余比较、secondary dynamic新特征分析、最小signature选择及最终结论。
