# PUF AlphaFold CP analysis

利用 AlphaFold 接触概率、结构置信度和几何特征，研究 PUF repeat 排列的构建成败与 TRM 变化后的位点编辑表现。骨架功能、C388活性及C295/C871旁观者编辑分别评估。

## 2026-09-22 更新

**[最新结果总览](results_20260922/README.md)** · **[第一模块 Wiki 中文](docs/wiki/PUF_Module1_ZH.md)** · **[English Wiki module](docs/wiki/PUF_Module1_EN.md)**

| 本次发布内容 | 结果与范围 | 文件 |
|---|---|---|
| 新补随机森林：初始骨架 | 初始14个留一验证14/14正确；仅12个PUF12为12/12 | [RF结果](results_20260922/scaffold_RF/results/metrics.csv) |
| 新补随机森林：扩展24骨架 | CP＋结构RF留一22/24正确；排列留出AUC0.9375、BA0.75 | [第一模块](docs/wiki/PUF_Module1_ZH.md) |
| v4动态阈值、全/非局部/界面CP | 全81个PR可用TRM留一：C295 57/81、C388 36/81、C871 49/81正确 | [全部81构建](results_20260922/v4/all81_LOCO.csv) |
| 指定五构建分别留出 | 保留同P其他构建；三个C共13/15命中 | [五构建结果](results_20260922/five_construct_holdout/predictions.csv) |
| 十构建同时留出 | 五候选＋五中等表现比较构建；预测Top5含3候选，综合Spearman 0.103 | [完整排名](results_20260922/ten_construct_holdout/ranking.csv) |
| 固定旧配方的两批TRM对比 | 同一十构建测试集，加入第二批没有使每个终点都改善 | [独立分析结果](results_20260922/wiki_two_batch/TRM_heldout_metrics.csv) |

这些都是既有实验记录上的回顾性验证。早期RF的100%限于初始数据；五构建13/15不能代替全体或十构建同时留出的成绩。不同标签、队列和验证协议的最高值不作为同一个“最佳模型”。v4非局部分支的高分另配[其自身动态阈值的多数类基线](results_20260922/v4/nonlocal_own_threshold_baseline.csv)。

本次补跑的代码和完整小型特征输入已发布；原始大体积AF/CP ZIP未直接入Git，见[来源清单与恢复方式](results_20260922/README.md#what-is-actually-included)。未删除下方历史分析。

## 2026-09-20 历史分析入口

[干实验 DBTL：设计、构建、测试与学习](docs/DRY_LAB_DBTL.md) · [历史复现说明](docs/REPRODUCIBILITY.md)

| 任务 | 历史数据与结果 | 入口 |
|---|---|---|
| 骨架跨排列验证 | 24个PUF12，4成功/20失败；三密度RF跨排列AUC0.86875，CP＋结构RF0.9375 | [报告](architecture_validation/PUF_architecture_validation_report.md) |
| C388联合动态阈值 | 31构建；40%标签、CP＋结构LR；组留出BA80.7%，LOCO67.9% | [报告](c388_threshold/README.md) |
| C871/C295动态标签分类 | 固定22突变体；位置组合选出候选BA87.5%/73.2%，含探索性选择偏差 | [报告](trm22_offtarget/REPORT.md) |
| C388活性保留 | 22突变体、ΔC388≥−15pp；固定候选组留出BA93.75%，完整自适应选择BA50% | [阶段性报告](c388_local_nonlocal_optimization/README.md) |

历史和新结果的样本、任务与验证方式见各自报告。当前仍无独立新实验确证的通用联合筛选器。骨架模型共享前部R123设计背景；功能标签也不等同于编辑特异性。

## 复现

新增RF分类器从小型特征表复现：

```bash
python results_20260922/scaffold_RF/retrain.py
```

历史检查命令保留：

```bash
python -m pip install -r requirements.txt
python tools/validate_snapshot.py
python c388_local_nonlocal_optimization/verify_results.py
```

历史三密度骨架评分：

```bash
python architecture_validation/predict.py --features new_construct_features.csv --out new_construct_scores.csv
```

输入需要按同一蛋白范围和model→seed→construct规则得到的 `pp_nonlocal4_high_per_res`、`pp_nonlocal12_high_per_res`、`pp_nonlocal24_high_per_res`。分数未校准。

[历史全部骨架指标](architecture_validation/validation_metrics.csv) · [历史折外预测](architecture_validation/heldout_predictions.csv) · [构建标签](architecture_validation/construct_registry.csv) · [特征定义](architecture_validation/CP_summary_dictionary.csv)

从缓存重现旧统计和图表的 `python tools/reproduce_current.py` 会覆盖旧输出，应在独立checkout中运行。初始阶段曾按2026-09-20整理要求删除的目录仍可从固定历史commit追溯；[原始输入清单](docs/raw_input_inventory.json)与新的[上传ZIP哈希](results_20260922/source_archives.json)分别保留其来源。
