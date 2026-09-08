# C388≥50%：PUF work/non-work 重新定义与模型验证

## 结论
旧三密度RF无法直接迁移：AUC=0.331，9个non-work中误报8个。重新训练出现中等排序信号，但模型选择不稳定：固定CP＋结构Logistic的LOCO/突变位置组留出AUC为0.778/0.798；把模型选择纳入嵌套验证后为0.576/0.702。当前结果不足以支持可靠的二分类预测。

若需要保守、可解释的研究基线，可保留三密度Logistic（LOCO AUC=0.722，balanced accuracy=0.753，TP16/FN6/TN7/FP2）。9维CP＋结构Logistic作为排序对照，不能把本数据上选择出的最高AUC当独立测试性能。

## 数据、匹配与标签
- Excel数值为0–1比例；0.5=50%。同名两条记录平均后，≥0.5为work，低于0.5为non-work。31种有标签蛋白，22 work/9 non-work，两条重复均未跨阈值。
- 仅使用C388配对结构提取输入。全部RNA为ACAUGGAGGACGUGC（15 nt）。C295/C871结构不作为额外独立样本。
- 33个C388折叠job：32种蛋白，其中P5-R6-SYQ有2个seed，其余1个seed；每seed5个model。主集31种蛋白对应160个model；先model均值、再seed等权均值，每蛋白最终一行。
- P4-R5-VTF有结构、无实验标签，排除训练与评估，单列未标注预测。97.zip补齐P6-R7-MHQ C388。
- 蛋白519 aa，统一去除额外N端SGSETPG七个残基，按512 aa范围计算旧模型特征。这沿用上轮TRM核对方案，并非针对当前标签调节。完整519 aa评分作为敏感性对照，旧RF AUC=0.462，仍无有效区分。
- 核心12×36残基映射为上传序列41–472位；逐构建真实突变、序列哈希与seed见sequence_mapping.csv。保留源标签中P7-P5-YYE拼写，不静默重命名。
- 原模型使用repeat arrangement构建work标签；本轮使用C388高编辑标签。虽然都简称work，两种终点不同。本轮低编辑不自动等同于蛋白未成功构建或折叠失败。

## 验证设计
冻结测试使用原24种PUF12训练的已保存CP_density3与S12_only模型，未用新标签重新拟合或反转分数。主模型未包含此前TRM敏感性样本，但WT核心已在原架构训练集中出现；另报告排除WT结果，不能称整批为完全独立结构验证。旧RNA背景为13 nt，本轮为15 nt，存在RNA上下文迁移。

重新训练采用固定参数LR（C=0.1、balanced）和浅层RF（100树、depth2、叶节点至少3、balanced）。比较CP3、CP14、结构6、CP＋结构9、RNA4、CP＋RNA7、CP＋界面6。阈值固定0.5；填补和标准化仅在训练折内拟合。

LOCO一次留出一蛋白。突变位置组留出将P4、P4+P7、P5、P6、P7、P8、WT分别整组留出；它检验新突变位置组迁移，不能称为不同repeat architecture验证。所有蛋白共享同一来源排列，无法做跨architecture泛化评估。

嵌套验证在外层训练集内，以三折内层验证从5个特征集合×LR/RF共10候选中按AUC选择，再预测外层留出样本。位置组留出对应内层也按位置分组，避免内外层分组混用。外层固定模型比较为探索性；嵌套结果用于审查模型选择偏差。所有逐构建预测均单列，禁止用重拟合训练分数替代CV。

## 主要性能
| family | validation | model | AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| frozen | transfer | CP_density3_RF | 0.331 | 0.396 | 0.652 | 0.682 | 1 | 8 | 7 | 15 |
| retrained | LOCO | CP3_LR | 0.722 | 0.753 | 0.889 | 0.727 | 7 | 2 | 6 | 16 |
| retrained | LOCO | CP3_RF | 0.631 | 0.697 | 0.842 | 0.727 | 6 | 3 | 6 | 16 |
| retrained | LOCO | CP_Structure9_LR | 0.778 | 0.606 | 0.800 | 0.545 | 6 | 3 | 10 | 12 |
| retrained | LOCO | CP_Structure9_RF | 0.657 | 0.687 | 0.818 | 0.818 | 5 | 4 | 4 | 18 |
| retrained | LOCO | RNA4_LR | 0.566 | 0.619 | 0.789 | 0.682 | 5 | 4 | 7 | 15 |
| retrained | mutation_group_out | CP3_LR | 0.697 | 0.674 | 0.833 | 0.682 | 6 | 3 | 7 | 15 |
| retrained | mutation_group_out | CP3_RF | 0.674 | 0.641 | 0.800 | 0.727 | 5 | 4 | 6 | 16 |
| retrained | mutation_group_out | CP_Structure9_LR | 0.798 | 0.707 | 0.875 | 0.636 | 7 | 2 | 8 | 14 |
| retrained | mutation_group_out | CP_Structure9_RF | 0.662 | 0.609 | 0.773 | 0.773 | 4 | 5 | 5 | 17 |
| retrained | mutation_group_out | RNA4_LR | 0.556 | 0.573 | 0.765 | 0.591 | 5 | 4 | 9 | 13 |
| nested_selection | LOCO | inner_select_10_candidates | 0.576 | 0.606 | 0.800 | 0.545 | 6 | 3 | 10 | 12 |
| nested_selection | mutation_group_out | inner_select_10_candidates | 0.702 | 0.684 | 0.867 | 0.591 | 7 | 2 | 9 | 13 |

AP、梯形PR-AUC、Brier及全部模型见metrics.csv。22/31为阳性，因此全部判work也有71.0% accuracy，但balanced accuracy只有0.5；仅看普通准确率会误导。

## 结构解释与不确定性
S12在本轮work/non-work均值约1.31965/1.32379，方向与旧架构数据不同，相关性较弱（r=-0.186，原始置换p约0.318）。旧规则“非局部CP越高越work”没有迁移到C388阈值终点。不能据此直接解释为降低接触能提高编辑。

本轮最强单变量候选为normalized Rg：work/non-work均值4.31282/4.34900，Hedges g=-1.075，r=-0.460，Monte Carlo置换p=0.0077，16候选BH q=0.123。没有候选达到本轮BH q<0.05。RNA接触总量有一定原始组差异，但RNA4模型LOCO AUC约0.566，单独预测较弱。不能因终点二分类有差异而声称能预测连续编辑率。

只有9个non-work，多数蛋白只有1个seed；SYQ的两个seed不足以证明全体CP稳定。5000次按标签分层重采样的AUC区间见conditional_score_bootstrap.csv，仅描述固定留出分数的条件不确定性，未重训，不能作为完整泛化置信区间。

## 敏感性结果
两个零值构建P4-R5-ETD、P4-R5-YTH+P7-R5-YTH的备注为variant表未检出C→T记录。主分析遵循表中0值；缺少覆盖度完整证据，不能把“未检出”当作严格测得0。剔除后29种、22 work/7 non-work：CP3 LR LOCO AUC=0.656，CP＋结构LR=0.708，提示部分区分来自低端样本。

排除WT后30种：CP3 LR LOCO AUC=0.720，CP＋结构LR=0.772。阈值由用户预先规定50%，未搜索更有利阈值。重复范围不跨50%，无标签冲突。

## 输出与使用
per_construct_comparison.csv：31种蛋白的实验均值/范围、真实标签、旧RF与新模型的逐构建留出分数。

unlabelled_predictions.csv：P4-R5-VTF使用全部31种重拟合模型的探索分数，真实标签未知，不能用其计算性能。

推荐暂作研究基线的CP3_LR.joblib与辅助CP_Structure9_LR.joblib已保存；系数、训练均值与标准差在logistic_scoring_coefficients.csv。分数均未校准；现在不建议自动给新设计作可靠work/non-work判定。更有价值的下一步是冻结模型，在新增TRM构建上独立测C388，尤其补充50%附近和不同突变位置的阴性。
