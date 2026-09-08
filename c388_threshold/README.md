# C388 实验标签阈值 × 模型分数阈值联合扫描

31 个独立构建，以每个构建的两次 C388 编辑率均值定义标签；模型/seed 已先聚合。分数 >= 分类阈值时预测 work。改变 C388 阈值会改变预测任务，不能理解为同一任务性能的直接提升。

## 范围与方法

- C388 5%–95%（5 个百分点步长）及全部相邻观测均值的中点；等价标签分组只训练一次，优先使用整 5% 分界值。保留两类各至少 2 个构建的 28 种标签划分。小于 2 个的极端类别未建模。
- CP3 / CP+Structure9 两组低维特征，各使用固定参数 Logistic Regression 和 shallow RF。
- LOCO 与 mutation-group-out。后者按 TRM 改变的位置组合留出，不等于 architecture-out。
- 分数阈值遍历留出分数相邻中点，另含 0、0.5 和全阴性边界。分别优化 accuracy、balanced accuracy；joint.py 中平局先比较另一指标，再取接近 0.5 者。
- 所有特征预处理只在训练折拟合；但是联合扫描最佳结果依据全部留出标签选阈值/模型，属于探索性上限，仍有选择偏差。
- 补充验证：在固定 C388 30%、35%、40%、50%、55% 和固定 LR 模型下，分类阈值仅通过外层训练集的内层 3 折 OOF 分数选择，目标 balanced accuracy；再预测外层测试构建。此步骤防止分类阈值直接使用测试标签，不能消除事后挑选 C388 阈值和模型的偏差。

## 结果

两类均至少 5 个构建时：

| 选择目标 | 验证 | C388 阈值 | 模型 | 分数阈值 | Accuracy | Balanced accuracy |
|---|---|---:|---|---:|---:|---:|
| 最高 accuracy | LOCO | 30% | CP+Structure9 LR | 0.376386 | 90.32% | 70.00% |
| 最高 balanced accuracy | LOCO | 35% | CP3 LR | 0.453027 | 83.87% | 79.46% |
| 最高 balanced accuracy | mutation-group-out | 40% | CP+Structure9 LR | 0.458533 | 80.65% | 82.88% |

若重点关注分组留出表现，下一轮可预先固定测试的候选为 C388 >= 40%，CP+Structure9 LR 分数 >= 0.46；该四舍五入阈值在本次分组留出预测上的分类与精确扫描阈值相同。40% 将样本划分为 23 work / 8 non-work；当前样本中任意 (37.332%, 44.306%] 的 C388 截断都会产生相同标签，数据不支持把 40% 解释为精确生物学边界。

该分组留出组合的探索性混淆矩阵：TP18、TN7、FP1、FN5，precision 94.74%、recall 78.26%。使用训练折内选分类阈值后，accuracy 77.42%、balanced accuracy 80.71%，TP17、TN7、FP1、FN6；训练所得阈值中位数 0.480867、范围 0.429228–0.520244。普通 LOCO 下同模型、同 C388 40%、训练内选分类阈值的 accuracy 64.52%、balanced accuracy 67.93%，提示阈值稳定性不足。

保留原 C388 50% 定义时，CP3 LR 训练折内动态阈值的 LOCO accuracy 77.42%、balanced accuracy 77.53%；mutation-group-out 为 70.97%、72.98%。

如果允许极小类别，扫描最高 accuracy 为 93.55%（多组并列），例如 C388 5%、CP+Structure9 LR、LOCO 分数阈值 0.235101。此时 29 work / 2 non-work，全部猜 work 也有 93.55%。这两个零值来自原表缺少变体记录，不等于已确认编辑率精确为零；因此不建议以 5% 作为有效最佳阈值。极端类别还可产生很高 balanced accuracy，其估计同样不稳定。

## 结论边界

没有已经独立验证的单一“最优阈值”。若用于决定实验是否足够有效，C388 标准应依据实验目标预先确定；不能只因为更容易预测就改变成功定义。40% / 0.46 是本次分组留出扫描的探索性候选，不是已校准的成功概率或经过新批次验证的部署规则。

## 文件

- joint_threshold_results.csv：全部标签阈值 × 模型 × 验证方式 × 优化目标的最佳结果，含等价 C388 阈值区间。
- joint_predictions.csv：所有外层留出构建的连续分数及真实标签。
- joint_inner_threshold_validation.csv / joint_inner_threshold_predictions.csv：训练内选分类阈值的补充验证及逐构建结果。
- frozen_joint_results.csv：历史冻结模型的探索性联合阈值扫描；不参与新模型最优推荐。
- apparent_best_thresholds.csv / dynamic_threshold_validation.csv：固定 C388 50% 的前一轮结果。
- joint.py、validate_joint.py、frozen_joint.py、thresholds.py：计算代码。
- training_data.csv、manifest.json、predictions.csv：复现输入；代码默认从 c388_analysis 对应路径读取。

LR: C=0.1, class_weight=balanced, solver=liblinear；RF: 100 trees, max_depth=2, min_samples_leaf=3, max_features=1.0, class_weight=balanced。随机种子 2026。没有新增 RF 参数搜索。
