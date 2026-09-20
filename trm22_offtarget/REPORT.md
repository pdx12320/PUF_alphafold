# 22个固定入组TRM突变体：C871与C295动态标签及分类方法比较

## 数据与目标

名单逐项沿用C388阶段的 frozen_cohort.csv，22个突变体全部进入两个位点分析。未按本次C388、C871或C295编辑率重新入组；WT仅作对应位点参考。实验率为两条实验记录的均值，Δ=突变体−位点WT，单位为百分点。C871 WT=67.5015%，C295 WT=32.958%。

正类定义为 Δ < t，表示下降超过指定幅度。编辑标签在−5至−60个百分点、5个百分点步长的固定网格扫描；全队列两类各至少3个才计算，正式候选摘要要求各至少5个。这个网格未穷举相邻观测中点。原C295 −20标签仅3个阳性，留作小类敏感性结果。所有结论均限定于这22个已入组突变体。

## 方法

比较8类输入：蛋白–RNA总CP差值、全PR矩阵RMS差值、14项CP汇总、19项仓库CP/结构特征、突变局部ΔCP、序列变化、序列+局部ΔCP、ContactSeek CCR。局部序列距离窗为≤4、≤18、19–54、≥55残基及全蛋白，保留各RNA碱基的差值；这些窗口没有指定编辑中心。序列特征包含P位置、氨基酸变化与TRM12/13/16位变化。没有加入新的完整PP矩阵或nonlocal coupling。

Interface_RMS沿用C388的519×15分母，缓存实际512×15，先等权平均seed再相对同位点WT求RMS。它描述protein–RNA矩阵整体扰动。与旧报告每seed统计后汇总的字段可能在多seed构建上不同，输入审计明确保留seed数。

9套分类配置：L2逻辑回归C=0.01/0.1/1；RF深度1/2；ExtraTrees深度2；自动收缩LDA；GaussianNB；RBF-SVM C=1。树模型100棵、叶节点至少3；RF/ET/LR/SVM使用balanced类权重，LDA/NB使用均衡先验。缺失值填补和标准化均在训练折拟合，CCR仅由相应训练折重新发现。SVM输出为decision_function的sigmoid变换，仅用作排序分数，未经概率校准。

外层验证包含留一构建、位置组合留出、严格位置留出。组合组是P4/P5/P6/P7/P8/P4+P7；严格位置每次排除含指定P的所有构建，双突变可能被测试两次，按预测次数倒数加权。严格位置仍允许双突变的另一个P出现在训练中。内层LOO流程使用最多3折分层，少数类只有1例时使用内层LOO；组合验证使用3折GroupKFold；严格位置使用训练内严格位置留出。内层单类训练折用常数预测并记录，不删除困难折。

模型分数阈值依据内层折外预测最大化BA，候选含0、0.5、相邻分数中点和全负边界；平局优先接近0.5。固定0.5作为对照。

每个固定编辑标签/特征/模型的外层成绩单独保存。另在每个外层训练集中同时选择编辑标签、特征、模型及分数阈值；标签要求该训练集两类各至少3例，平局优先下降幅度较小的标签。最终外层测试标签根据当折训练选定的编辑界限生成。

## 标签划分

| site | threshold_pp | absolute_editing_pct | positive | negative | eligible_headline |
| --- | --- | --- | --- | --- | --- |
| C871 | -5 | 62.5015 | 19 | 3 | False |
| C871 | -10 | 57.5015 | 15 | 7 | True |
| C871 | -15 | 52.5015 | 13 | 9 | True |
| C871 | -20 | 47.5015 | 10 | 12 | True |
| C871 | -25 | 42.5015 | 9 | 13 | True |
| C871 | -30 | 37.5015 | 7 | 15 | True |
| C871 | -35 | 32.5015 | 4 | 18 | False |
| C295 | -5 | 27.9580 | 18 | 4 | False |
| C295 | -10 | 22.9580 | 14 | 8 | True |
| C295 | -15 | 17.9580 | 8 | 14 | True |
| C295 | -20 | 12.9580 | 3 | 19 | False |

## 探索性最佳固定方案

每种验证分别扫描得到的最高值存在多候选选择偏差；不同标签定义了不同预测任务。不能把改变标签后的BA升高单独归因于模型改善。

| site | validation | label_threshold_pp | family | model | BA | AUC | precision | recall | TP | TN | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | LOOCV | -15 | Repository | RF_depth2 | 0.8846 | 0.8376 | 1.0000 | 0.7692 | 10.0000 | 9.0000 | 0.0000 | 3.0000 |
| C871 | Position_group_out | -20 | Repository | RF_depth1 | 0.8750 | 0.7917 | 0.7692 | 1.0000 | 10.0000 | 9.0000 | 3.0000 | 0.0000 |
| C871 | Strict_position_out | -30 | CP_summary | RF_depth1 | 0.8571 | 0.8714 | 1.0000 | 0.7143 | 5.0000 | 15.0000 | 0.0000 | 2.0000 |
| C295 | LOOCV | -15 | Sequence_Local | RF_depth1 | 0.7768 | 0.7500 | 0.8333 | 0.6250 | 5.0000 | 13.0000 | 1.0000 | 3.0000 |
| C295 | Position_group_out | -10 | Local_deltaCP | LR_0.01 | 0.7321 | 0.5536 | 0.8333 | 0.7143 | 10.0000 | 6.0000 | 2.0000 | 4.0000 |
| C295 | Strict_position_out | -10 | ContactSeek_CCR | LR_1 | 0.8571 | 0.7299 | 0.8710 | 0.9643 | 13.5000 | 6.0000 | 2.0000 | 0.5000 |

## 按位置组合成绩选出的同一个候选，在三种验证下的表现

| site | validation | policy | label_threshold_pp | family | model | BA | AUC | precision | recall |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | LOOCV | fixed_0.5 | -20 | Repository | RF_depth1 | 0.7833 | 0.8167 | 0.6923 | 0.9000 |
| C871 | LOOCV | training_only_threshold | -20 | Repository | RF_depth1 | 0.6917 | 0.8167 | 0.6154 | 0.8000 |
| C871 | Position_group_out | fixed_0.5 | -20 | Repository | RF_depth1 | 0.7750 | 0.7917 | 0.7273 | 0.8000 |
| C871 | Position_group_out | training_only_threshold | -20 | Repository | RF_depth1 | 0.8750 | 0.7917 | 0.7692 | 1.0000 |
| C871 | Strict_position_out | fixed_0.5 | -20 | Repository | RF_depth1 | 0.8083 | 0.7708 | 0.7037 | 0.9500 |
| C871 | Strict_position_out | training_only_threshold | -20 | Repository | RF_depth1 | 0.8083 | 0.7708 | 0.7037 | 0.9500 |
| C295 | LOOCV | fixed_0.5 | -10 | Local_deltaCP | LR_0.01 | 0.5982 | 0.5357 | 0.7273 | 0.5714 |
| C295 | LOOCV | training_only_threshold | -10 | Local_deltaCP | LR_0.01 | 0.5625 | 0.5357 | 0.7000 | 0.5000 |
| C295 | Position_group_out | fixed_0.5 | -10 | Local_deltaCP | LR_0.01 | 0.5536 | 0.5536 | 0.7143 | 0.3571 |
| C295 | Position_group_out | training_only_threshold | -10 | Local_deltaCP | LR_0.01 | 0.7321 | 0.5536 | 0.8333 | 0.7143 |
| C295 | Strict_position_out | fixed_0.5 | -10 | Local_deltaCP | LR_0.01 | 0.5536 | 0.5982 | 0.7143 | 0.3571 |
| C295 | Strict_position_out | training_only_threshold | -10 | Local_deltaCP | LR_0.01 | 0.4375 | 0.5982 | 0.5385 | 0.2500 |

各外层折训练内选择的分数阈值：

| site | validation | score_threshold_median | score_threshold_min | score_threshold_max |
| --- | --- | --- | --- | --- |
| C871 | LOOCV | 0.3906 | 0.1881 | 0.6186 |
| C871 | Position_group_out | 0.5000 | 0.4260 | 0.5641 |
| C871 | Strict_position_out | 0.5000 | 0.4838 | 0.5573 |
| C295 | LOOCV | 0.4889 | 0.3893 | 0.5912 |
| C295 | Position_group_out | 0.5639 | 0.3943 | 0.6500 |
| C295 | Strict_position_out | 0.4666 | 0.2699 | 0.5608 |

## 完整训练内动态选择流程

以下评估整个自适应流程。不同外层折可能选择不同编辑标签，因此其BA、AUC不能解释为某一个统一下降标准的固定终点表现；原始逐折阈值与预测均单独提供。

| site | validation | BA | accuracy | AUC | precision | recall | TP | TN | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C295 | LOOCV | 0.5729 | 0.6818 | 0.8021 | 0.4000 | 0.3333 | 2.0000 | 13.0000 | 3.0000 | 4.0000 |
| C295 | Position_group_out | 0.5104 | 0.5909 | 0.4271 | 0.7333 | 0.6875 | 11.0000 | 2.0000 | 4.0000 | 5.0000 |
| C295 | Strict_position_out | 0.4643 | 0.4773 | 0.5631 | 0.6522 | 0.5000 | 7.5000 | 3.0000 | 4.0000 | 7.5000 |
| C871 | LOOCV | 0.5128 | 0.5455 | 0.4444 | 0.4286 | 0.3333 | 3.0000 | 9.0000 | 4.0000 | 6.0000 |
| C871 | Position_group_out | 0.3482 | 0.4091 | 0.4375 | 0.5333 | 0.5714 | 8.0000 | 1.0000 | 7.0000 | 6.0000 |
| C871 | Strict_position_out | 0.6417 | 0.6364 | 0.6385 | 0.7000 | 0.5833 | 7.0000 | 7.0000 | 3.0000 | 5.0000 |

## 不同训练方法的探索性比较

每行在位置组合留出中选择该方法最好的标签/输入，仅作为探索性比较。

| site | model | label_threshold_pp | family | BA | AUC | precision | recall |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | RF_depth1 | -20 | Repository | 0.8750 | 0.7917 | 0.7692 | 1.0000 |
| C871 | RF_depth2 | -20 | Repository | 0.8750 | 0.7833 | 0.7692 | 1.0000 |
| C871 | ExtraTrees_depth2 | -20 | ContactSeek_CCR | 0.8500 | 0.8000 | 1.0000 | 0.7000 |
| C871 | LDA_shrinkage | -30 | CP_summary | 0.7524 | 0.5714 | 0.8000 | 0.5714 |
| C871 | LR_1 | -30 | Local_deltaCP | 0.7333 | 0.6190 | 0.4667 | 1.0000 |
| C295 | LR_0.01 | -10 | Local_deltaCP | 0.7321 | 0.5536 | 0.8333 | 0.7143 |
| C295 | LDA_shrinkage | -15 | Sequence | 0.7232 | 0.7143 | 0.5385 | 0.8750 |
| C295 | RF_depth2 | -10 | Interface_RMS | 0.7232 | 0.6875 | 0.8889 | 0.5714 |
| C295 | RF_depth1 | -10 | Interface_RMS | 0.7232 | 0.6027 | 0.8889 | 0.5714 |
| C871 | LR_0.01 | -25 | Total_CP | 0.7137 | 0.6752 | 0.5714 | 0.8889 |
| C871 | LR_0.1 | -25 | Total_CP | 0.6966 | 0.6667 | 0.5833 | 0.7778 |
| C295 | ExtraTrees_depth2 | -10 | Total_CP | 0.6964 | 0.3393 | 0.8182 | 0.6429 |
| C295 | SVM_RBF | -15 | Total_CP | 0.6875 | 0.5000 | 0.5000 | 0.8750 |
| C871 | GaussianNB | -10 | Sequence | 0.6857 | 0.6429 | 0.8000 | 0.8000 |
| C295 | GaussianNB | -10 | Sequence | 0.6786 | 0.6786 | 0.7500 | 0.8571 |
| C295 | LR_1 | -10 | ContactSeek_CCR | 0.6518 | 0.6607 | 0.8571 | 0.4286 |
| C295 | LR_0.1 | -15 | Sequence | 0.6518 | 0.6071 | 0.4667 | 0.8750 |
| C871 | SVM_RBF | -25 | Total_CP | 0.6368 | 0.6923 | 0.5000 | 0.8889 |

## 解释边界与复现

所有候选使用同一批既有实验，尚无独立新突变体测试。WT与多数突变体seed有限，WT及实验重复的不确定性没有传播；没有全流程置换检验。22个名单自身来自C388入组，结论不能外推到全部TRM突变体。旧30个突变体的指标不能与本轮直接当作同任务的提升。没有用全数据重拟合成绩报告泛化，也未发布部署模型。

运行 `python run.py` 重算模型，默认4个CPU进程，可用TRM_JOBS调整。运行 `python report.py` 重算指标、验证拆分并生成本报告。inputs包含复现缓存；sources保留ContactSeek源码及许可证；results/provenance.json记录版本和输入哈希。保存的逐样本分数可复核所有指标。run.log记录训练过程。

训练方法遵循scikit-learn技能中的训练折内预处理与嵌套验证流程。软件工作流参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. https://doi.org/10.48550/arXiv.2609.00065 。该文献用于记录分析工具来源，不作为本项目生物学结论的证据。


## 严格位置验证选出的候选，在全部验证下的表现

同样属于事后探索性选择，需与前文位置组合选出的候选一起阅读。

| site | label_threshold_pp | family | model | validation | BA | AUC | precision | recall | TP | TN | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | -30 | CP_summary | RF_depth1 | LOOCV | 0.8286 | 0.9333 | 0.6667 | 0.8571 | 6.0000 | 12.0000 | 3.0000 | 1.0000 |
| C871 | -30 | CP_summary | RF_depth1 | Position_group_out | 0.8571 | 0.8190 | 1.0000 | 0.7143 | 5.0000 | 15.0000 | 0.0000 | 2.0000 |
| C871 | -30 | CP_summary | RF_depth1 | Strict_position_out | 0.8571 | 0.8714 | 1.0000 | 0.7143 | 5.0000 | 15.0000 | 0.0000 | 2.0000 |
| C295 | -10 | ContactSeek_CCR | LR_1 | LOOCV | 0.6429 | 0.6518 | 0.7333 | 0.7857 | 11.0000 | 4.0000 | 4.0000 | 3.0000 |
| C295 | -10 | ContactSeek_CCR | LR_1 | Position_group_out | 0.6518 | 0.6607 | 0.8571 | 0.4286 | 6.0000 | 7.0000 | 1.0000 | 8.0000 |
| C295 | -10 | ContactSeek_CCR | LR_1 | Strict_position_out | 0.8571 | 0.7299 | 0.8710 | 0.9643 | 13.5000 | 6.0000 | 2.0000 | 0.5000 |