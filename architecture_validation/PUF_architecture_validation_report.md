# PUF12跨repeat排列：结构机制与泛化复核

主结论：当前模型同时利用了接触完整性相关信号和repeat排列背景。高non-local CP density有跨两个已标注混合组的同方向证据，但尚不足以认定为所有PUF12 architecture通用的成功机制。原CCR01的增强/协同变化未在扩充数据中复现。建议继续保留3个高CP密度特征的浅层RF作保守基线，分数用于优先级排序，不能作为可靠成功概率。

## 1. 数据与分组

主集为两批24种PUF12构建：4 success、20 failure，96个seed、480个预测model。第一批原14构建中2个PUF11另做敏感性分析；第二批ID9/10/11为15-repeat而排除。仅改TRM的WT_PUF12-9与P8-R6-GVE另外加入26样本敏感性分析，同一蛋白的三个RNA先分别聚合，再等权汇总；它们不增加独立骨架类别。

每个model先提特征，model→seed→construct等权聚合。原批已经核验的缓存用于复现；新批重新提取。主集14项CP summary与此前合并表比对通过（validation.json）。模型中每种构建一行，seed/model没有作为独立训练样本。

按去除明确MNDGPHS插入段后的36-aa repeat core映射位置，利用忽略TRM位点的source模块匹配确定source/order。P1–P12指排列位置；R1–R8指来源repeat。source不同的同一P位置不能直接视作同一氨基酸。各构建原始编号、氨基酸和source均已输出。第二批ID16按实际序列来源顺序归类，与名字暗示的顺序存在差异。

| architecture | source_order | n_construct | n_success | n_failure |
| --- | --- | --- | --- | --- |
| G01 | 1,2,3,4,5,6,6,5,4,3,2,8 | 2 | 0 | 2 |
| G02 | 1,2,3,4,5,6,7,3,4,5,6,8 | 2 | 0 | 2 |
| G03 | 1,2,3,5,6,7,5,6,7,6,7,8 | 10 | 3 | 7 |
| G04 | 1,2,3,5,6,7,6,7,5,6,7,8 | 3 | 1 | 2 |
| G05 | 1,2,3,5,6,7,8,5,6,7,8,8 | 2 | 0 | 2 |
| G06 | 1,2,3,5,6,7,8,6,7,8,4,4 | 4 | 0 | 4 |
| G08 | 1,2,3,6,5,6,6,5,1,3,2,8 | 1 | 0 | 1 |

G03对应R123/R567/R567/R678类，G04对应R123/R567/R6R7/R5678类，G05对应R123/R5678/R56788类，G06对应含末端R4R4类。G07为PUF11组，仅敏感性分析使用。其余按表中实际顺序定义。G03的loop变体全部同组。

完整逐construct清单：construct_registry.csv。所有主集构建仍共享前部R123，因此本次评价范围为R123背景下的后部repeat排列迁移；缺少替换前部R123的真实数据，不能估计该情形性能。

## 2. 特征定义、范围与重要性

设L为分析蛋白残基数，N为RNA长度，CP为接触概率。原CP summary共14项，6项蛋白内部特征与8项蛋白–RNA特征。高接触阈值均为CP>0.5；原CP summary先逐model计算再平均，阈值不施加在先平均的CP上。

| feature | definition | minimum | maximum | LOCO_RF_MDI | architecture_out_RF_MDI |
| --- | --- | --- | --- | --- | --- |
| pr_cp_sum | sum of all protein-RNA CP_ij | 42.56 | 54.39 | 0.04391 | 0.04711 |
| pr_cp_per_nt | protein-RNA CP sum / RNA nucleotide count | 3.274 | 4.184 | 0.03716 | 0.04014 |
| pr_top20_mean | mean of 20 largest protein-RNA CP values | 0.9279 | 0.9959 | 0.02139 | 0.02203 |
| pr_rna_max_mean | mean across RNA nucleotides of maximum CP over protein residues | 0.7021 | 0.8444 | 0.05289 | 0.05424 |
| pr_rna_max_min | minimum across RNA nucleotides of maximum CP over protein residues | 0.0825 | 0.21 | 0.008 | 0.01086 |
| pr_rna_coverage_05 | fraction of RNA nucleotides with maximum protein CP>0.5 | 0.6923 | 0.8077 | 0.0333 | 0.03727 |
| pr_protein_coverage_05 | fraction of protein residues with maximum RNA CP>0.5 | 0.05408 | 0.07749 | 0.03005 | 0.02374 |
| pr_high_contacts_per_nt | count protein-RNA CP>0.5 / RNA nucleotide count | 2.654 | 4.058 | 0.0192 | 0.0185 |
| pp_nonlocal4_per_res | sum CP_ij / L, protein i<j and j-i>=4 | 2.345 | 2.44 | 0.05102 | 0.05648 |
| pp_nonlocal4_high_per_res | count(CP_ij>0.5) / L, protein i<j and j-i>=4 | 2.213 | 2.387 | 0.1967 | 0.1943 |
| pp_nonlocal12_per_res | sum CP_ij / L, protein i<j and j-i>=12 | 1.286 | 1.36 | 0.06751 | 0.0749 |
| pp_nonlocal12_high_per_res | count(CP_ij>0.5) / L, protein i<j and j-i>=12 | 1.189 | 1.333 | 0.1704 | 0.1636 |
| pp_nonlocal24_per_res | sum CP_ij / L, protein i<j and j-i>=24 | 1.104 | 1.162 | 0.05158 | 0.05092 |
| pp_nonlocal24_high_per_res | count(CP_ij>0.5) / L, protein i<j and j-i>=24 | 1.027 | 1.128 | 0.2169 | 0.206 |

完整结构/界面字典：feature_dictionary_ranges_importance.csv。范围为24种构建的观测范围，不能作为新设计的硬性生物学边界。

新增并统计：核心平均/最小pLDDT、全核心pair平均PAE、contact-weighted PAE、各repeat-pair PAE、全长Rg及normalized Rg、核心anisotropy、compactness、置信度过滤接触密度、每个repeat内部密度、全部66个repeat-pair CP、相邻与非相邻聚合指标、局部残基/接触/CCR。

- 为统一所有批次，minimum pLDDT定义为“每个残基先跨model/seed平均，再取核心残基最小值”；未称作全长逐model最小值平均。
- global_core_PAE使用所有核心i<j的对称PAE平均；模型同时保留原全长、gap≥12的contact-weighted PAE。
- 高置信接触密度定义为：核心平均CP>0.5、对称平均PAE<5 Å、两端平均pLDDT≥80且native gap≥12的pair数/核心长度。这是聚合后阈值指标，与原summary的逐model阈值定义分开。
- repeat internal density：native gap≥4的独特pair CP之和/36。repeat interface CP：两个core模块间native gap≥4的pair CP平均。非相邻定义为P位置相差至少2。
- compactness=L/(4πRg³/3)，仅为几何代理；与normalized Rg单调冗余，因此未同时塞入低维模型。统计表仍保留两者。
- 原K204/A235和CCR01采用参考成功构建编号；新构建对应残基可能改变，不能直接据此指定突变。

RF图为训练fold的平均杂质重要性，相关特征会分摊重要性。feature_importance_summary.csv另外提供留出样本特征置为训练中位数后的log-loss变化；正值表示替换后变差，负值表示当前特征使用未带来该诊断下的收益。该替换诊断也可能受到特征相关性影响，不能当作因果效应。

## 3. 精确置换与architecture混杂

所有单变量检验以construct为单位。对24种构建枚举C(24,4)=10,626种成功标签分配，以绝对组均值差为统计量，得到双侧精确p；输出Pearson/point-biserial、Spearman、Hedges g。BH同时提供按特征家族校正q_family和全部1,541项校正q_global。没有任何无条件关联达到本轮家族BH q<0.05；这不等于证明所有效应不存在。

对“feature ~ success + architecture”，用architecture固定效应残差化计算success回归系数和partial r，并枚举保持各architecture成功数不变的360种标签分配。只有G03和G04含两类样本；其余组无法估计组内success−failure差异。architecture_R2衡量特征变异被组别解释的比例，不是因果归因。

| feature | success_mean | failure_mean | delta | hedges_g | pearson_r | spearman_r | p_exact | q_family | architecture_R2 | adjusted_success_beta | p_within_architecture_exact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pp_nonlocal12_high_per_res | 1.323 | 1.275 | 0.04799 | 1.286 | 0.4604 | 0.5815 | 0.0191 | 0.08915 | 0.7446 | 0.01928 | 0.06389 |
| P5_P6_CP | 0.03895 | 0.03758 | 0.001363 | 0.7273 | 0.2814 | 0.323 | 0.1619 | 0.9071 | 0.9872 | 0.000194 | 0.2528 |
| K204_slot_density | 3.189 | 2.282 | 0.907 | 0.8433 | 0.3219 | 0.5653 | 0.1023 | 0.1534 | 0.9908 | 0.1769 | 0.04167 |
| A235_slot_density | 3.344 | 2.71 | 0.6333 | 0.862 | 0.3282 | 0.4522 | 0.08743 | 0.1526 | 0.9902 | 0.08949 | 0.1222 |
| CCR01_fixed | 3.067 | 3.107 | -0.03946 | -0.08938 | -0.03601 | -0.09691 | 0.9145 | 0.9145 | 0.9975 | -0.03066 | 0.04167 |

组内方向检查：

| feature | architecture | n_success | n_failure | delta |
| --- | --- | --- | --- | --- |
| pp_nonlocal12_high_per_res | G03 | 3 | 7 | 0.01178 |
| pp_nonlocal12_high_per_res | G04 | 1 | 2 | 0.04288 |
| P5_P6_CP | G03 | 3 | 7 | 0.0001096 |
| P5_P6_CP | G04 | 1 | 2 | 0.0004599 |
| K204_slot_density | G03 | 3 | 7 | 0.09762 |
| K204_slot_density | G04 | 1 | 2 | 0.4267 |
| A235_slot_density | G03 | 3 | 7 | 0.02888 |
| A235_slot_density | G04 | 1 | 2 | 0.2804 |
| CCR01_fixed | G03 | 3 | 7 | -0.03153 |
| CCR01_fixed | G04 | 1 | 2 | -0.02792 |

S12高接触密度仍是优先候选之一，但没有证据认定它唯一最稳定：S4、S12、S24相关且统计表现接近；RF三密度模型的architecture折平均重要性分别约0.359、0.358、0.283。S12成功均值1.3230、失败1.2750，g=1.286，Pearson r=0.460，原始p=0.0191，CP家族BH q=0.0892；组内控制architecture后p=0.0639。两个混合组方向均正，有有限跨组一致性，但只有两个背景提供这类证据。

P5–P6界面CP约98.7%的变异可由architecture解释，组内p=0.253；优先标记为architecture-associated feature。K204/A235位置密度同样约99%的变异随architecture变化；K204组内p=0.0417，但多重校正未通过。它们可作为实验候选，尚不能称为通用必要位点。

## 4. repeat界面、接触对与CCR

每种construct的P1–P2至P11–P12与全部非相邻pair已输出，重复对CP变化另见repeat_pair_delta_reference.csv。参考为原512-aa成功构建；subtract-reference只作可视化，固定参考相减不改变单变量相关系数。

接触pair采用标签无关筛选：aligned gap≥12、跨构建最大CP≥0.1、方差非零，共928个pair。完整结果见univariate_contact.csv；前100个映射到参考和每个构建原始编号，见top_contact_pairs_mapped.csv及top_contact_native_mapping.csv。最靠前的参考421–458与204–239接触增强，原始p约0.0010/0.0012，但接触家族BH q均约0.568。不能称为显著功能接触。全体repeat界面/残基检验也未通过家族BH校正。

CCR以非局部CP残基密度为profile，在每个P模块内寻找连续3–8个残基，要求内部两两Pearson及Spearman均≥0.6。新发现5个候选区域；候选生成基于无标签协同变化，成功关联另行检验。它们代表本次操作定义下的密度协同区域，不等同于已验证RNA结合位点。模型每轮只在训练construct内部重新发现并选1个CCR，避免将全数据筛出的区域直接放进验证。

| feature | reference_residues | repeat_slot | min_internal_Pearson | min_internal_Spearman | pearson_r | p_exact | q_family | p_within_architecture_exact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CCR01_fixed | 178,179,180 | 5 | -0.8922 | -0.1001 | -0.03601 | 0.9145 | 0.9145 | 0.04167 |
| CCR_new_01 | 162,163,164 | 4 | 0.7405 | 0.7548 | -0.05276 | 0.659 | 0.659 | 0.002778 |
| CCR_new_02 | 195,196,197,198 | 5 | 0.7622 | 0.6139 | -0.1915 | 0.5305 | 0.659 | 0.2222 |
| CCR_new_03 | 411,412,413 | 11 | 0.6435 | 0.7722 | 0.3762 | 0.07086 | 0.1856 | 0.9889 |
| CCR_new_04 | 414,415,416 | 11 | 0.7828 | 0.6144 | 0.3657 | 0.07425 | 0.1856 | 0.7 |
| CCR_new_05 | 447,448,449 | 12 | 0.6513 | 0.709 | 0.2521 | 0.2372 | 0.3954 | 0.075 |

原CCR01（参考178/179/180，P5起始）与成功关联r=-0.036，p=0.915；内部最小Pearson降到-0.892，已不满足本次协同区域要求。它在两个混合组内的success−failure差异均为负。原正向CCR结论未复现，应从通用机制表述中撤回。新CCR也未建立经过BH校正的成功关联。

## 5. 模型和验证设计

固定小样本模型：LR C=0.1、balanced、liblinear；RF 300棵树、max_depth=2、min_samples_leaf=2、max_features=sqrt、balanced；XGBoost 100棵depth=1树、eta=.05、lambda=10、alpha=1、训练集内类别权重。随机种子2026；填补及标准化仅拟合训练集；阈值固定0.5。

低维模型集合：

- CP_summary14：完整原14特征对照。
- CP_density3：S4/S12/S24三个高CP密度，作为主要简约基线。
- structural_only：核心平均/最小pLDDT、global core PAE、contact-weighted PAE、normalized Rg、anisotropy，共6维。
- CP_structural：三密度＋6结构，共9维。
- CP_interface：三密度＋相邻interface CP均值/最小值＋非相邻CP均值，共6维。
- full：三密度＋6结构＋3界面聚合＋训练内选1残基、1CCR、1contact，共15维。没有flatten全CP直接训练；完整候选池仅用于发现和统计。
- S12_only：1维附加基线。

以上采用预先固定的低维代表指标作模型比较；大量逐repeat、PAE和HC候选保留在单变量审计中，未全部堆入分类器。“full”表示覆盖各特征类别，并非输入全部1,541项。候选类型/模型比较本身仍有探索性选择偏差。

A：24折LOCO。B：按source顺序忽略loop插入位置定义的7组architecture整组留出。C：source-repeat-order整组留出。本数据定义下B与C完全相同，性能表保留两种名称，不能将它们算作两份独立验证证据。

各折均有训练阳性，但留出G03时训练中仅剩1个成功construct，模型不稳定风险很高。只有G03/G04测试组可单独计算AUC；纯失败组AUC、balanced accuracy记NA，报告实际误报数。总表AUC来自拼接所有留出预测，包含组间区分贡献；不是每个新组的预期AUC。

RF主要结果：

| feature_set | model | validation | AUC | AP | PR_AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CP_density3 | RF | LOCO | 0.9125 | 0.5845 | 0.4798 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| CP_interface | RF | LOCO | 0.925 | 0.622 | 0.5193 | 0.8 | 0.5 | 0.75 | 17 | 3 | 1 | 3 |
| CP_structural | RF | LOCO | 0.9125 | 0.7875 | 0.7688 | 0.85 | 0.75 | 0.75 | 19 | 1 | 1 | 3 |
| CP_summary14 | RF | LOCO | 0.85 | 0.4444 | 0.3608 | 0.7 | 0.5 | 0.5 | 18 | 2 | 2 | 2 |
| S12_only | RF | LOCO | 0.875 | 0.4833 | 0.3994 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| full | RF | LOCO | 0.85 | 0.6104 | 0.5701 | 0.725 | 0.6667 | 0.5 | 19 | 1 | 2 | 2 |
| structural_only | RF | LOCO | 0.75 | 0.7917 | 0.7871 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |
| CP_density3 | RF | architecture_out | 0.8688 | 0.475 | 0.6542 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_interface | RF | architecture_out | 0.9062 | 0.6625 | 0.7479 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_structural | RF | architecture_out | 0.9375 | 0.8611 | 0.8524 | 0.75 | 1 | 0.5 | 20 | 0 | 2 | 2 |
| CP_summary14 | RF | architecture_out | 0.7625 | 0.4196 | 0.3182 | 0.65 | 0.3333 | 0.5 | 16 | 4 | 2 | 2 |
| S12_only | RF | architecture_out | 0.8 | 0.3409 | 0.583 | 0.725 | 0.3333 | 0.75 | 14 | 6 | 1 | 3 |
| full | RF | architecture_out | 0.85 | 0.7019 | 0.6822 | 0.625 | 1 | 0.25 | 20 | 0 | 3 | 1 |
| structural_only | RF | architecture_out | 0.7625 | 0.7935 | 0.7888 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |

所有LR/RF/XGB结果：

| feature_set | model | validation | AUC | AP | PR_AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CP_density3 | LR | LOCO | 0.9125 | 0.5845 | 0.4798 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_density3 | RF | LOCO | 0.9125 | 0.5845 | 0.4798 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| CP_density3 | XGB | LOCO | 0.9 | 0.6012 | 0.4131 | 0.675 | 0.4 | 0.5 | 17 | 3 | 2 | 2 |
| CP_interface | LR | LOCO | 0.925 | 0.6917 | 0.65 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_interface | RF | LOCO | 0.925 | 0.622 | 0.5193 | 0.8 | 0.5 | 0.75 | 17 | 3 | 1 | 3 |
| CP_interface | XGB | LOCO | 0.9 | 0.6012 | 0.4131 | 0.675 | 0.4 | 0.5 | 17 | 3 | 2 | 2 |
| CP_structural | LR | LOCO | 0.8375 | 0.8088 | 0.8028 | 0.625 | 0.2308 | 0.75 | 10 | 10 | 1 | 3 |
| CP_structural | RF | LOCO | 0.9125 | 0.7875 | 0.7688 | 0.85 | 0.75 | 0.75 | 19 | 1 | 1 | 3 |
| CP_structural | XGB | LOCO | 0.65 | 0.3738 | 0.2761 | 0.55 | 0.25 | 0.25 | 17 | 3 | 3 | 1 |
| CP_summary14 | LR | LOCO | 0.825 | 0.4132 | 0.3309 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_summary14 | RF | LOCO | 0.85 | 0.4444 | 0.3608 | 0.7 | 0.5 | 0.5 | 18 | 2 | 2 | 2 |
| CP_summary14 | XGB | LOCO | 0.8875 | 0.5179 | 0.4131 | 0.675 | 0.4 | 0.5 | 17 | 3 | 2 | 2 |
| S12_only | LR | LOCO | 0.9 | 0.65 | 0.6077 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| S12_only | RF | LOCO | 0.875 | 0.4833 | 0.3994 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| S12_only | XGB | LOCO | 0.675 | 0.4417 | 0.2996 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| full | LR | LOCO | 0.725 | 0.3125 | 0.2394 | 0.675 | 0.2727 | 0.75 | 12 | 8 | 1 | 3 |
| full | RF | LOCO | 0.85 | 0.6104 | 0.5701 | 0.725 | 0.6667 | 0.5 | 19 | 1 | 2 | 2 |
| full | XGB | LOCO | 0.3375 | 0.1655 | 0.1261 | 0.375 | 0 | 0 | 15 | 5 | 4 | 0 |
| structural_only | LR | LOCO | 0.7875 | 0.6568 | 0.6201 | 0.6 | 0.2143 | 0.75 | 9 | 11 | 1 | 3 |
| structural_only | RF | LOCO | 0.75 | 0.7917 | 0.7871 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |
| structural_only | XGB | LOCO | 0.75 | 0.7917 | 0.7871 | 0.75 | 1 | 0.5 | 20 | 0 | 2 | 2 |
| CP_density3 | LR | architecture_out | 0.9 | 0.7042 | 0.6646 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_density3 | RF | architecture_out | 0.8688 | 0.475 | 0.6542 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_density3 | XGB | architecture_out | 0.8938 | 0.5179 | 0.5298 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| CP_interface | LR | architecture_out | 0.8875 | 0.7409 | 0.7205 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_interface | RF | architecture_out | 0.9062 | 0.6625 | 0.7479 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_interface | XGB | architecture_out | 0.9062 | 0.5679 | 0.5548 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| CP_structural | LR | architecture_out | 0.925 | 0.85 | 0.8417 | 0.75 | 0.2857 | 1 | 10 | 10 | 0 | 4 |
| CP_structural | RF | architecture_out | 0.9375 | 0.8611 | 0.8524 | 0.75 | 1 | 0.5 | 20 | 0 | 2 | 2 |
| CP_structural | XGB | architecture_out | 0.8812 | 0.5 | 0.5119 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| CP_summary14 | LR | architecture_out | 0.8 | 0.4623 | 0.3606 | 0.7 | 0.25 | 1 | 8 | 12 | 0 | 4 |
| CP_summary14 | RF | architecture_out | 0.7625 | 0.4196 | 0.3182 | 0.65 | 0.3333 | 0.5 | 16 | 4 | 2 | 2 |
| CP_summary14 | XGB | architecture_out | 0.8875 | 0.5179 | 0.4131 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| S12_only | LR | architecture_out | 0.8875 | 0.7409 | 0.7205 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| S12_only | RF | architecture_out | 0.8 | 0.3409 | 0.583 | 0.725 | 0.3333 | 0.75 | 14 | 6 | 1 | 3 |
| S12_only | XGB | architecture_out | 0.8125 | 0.5088 | 0.6562 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| full | LR | architecture_out | 0.875 | 0.5417 | 0.4375 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| full | RF | architecture_out | 0.85 | 0.7019 | 0.6822 | 0.625 | 1 | 0.25 | 20 | 0 | 3 | 1 |
| full | XGB | architecture_out | 0.5438 | 0.2917 | 0.3165 | 0.55 | 0.25 | 0.25 | 17 | 3 | 3 | 1 |
| structural_only | LR | architecture_out | 0.8875 | 0.8269 | 0.8197 | 0.725 | 0.2667 | 1 | 9 | 11 | 0 | 4 |
| structural_only | RF | architecture_out | 0.7625 | 0.7935 | 0.7888 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |
| structural_only | XGB | architecture_out | 0.3812 | 0.375 | 0.3393 | 0.625 | 1 | 0.25 | 20 | 0 | 3 | 1 |

AP为阶梯加权average precision；PR_AUC为precision–recall曲线梯形面积。并列分数较多时二者可能有明显差异，不能互换。TN/FP/FN/TP均按success为正类；分数尚未校准。

## 6. 对泛化的直接回答

1. 最重要信号属于高CP非局部接触密度家族，S4/S12/S24都应考虑；具体训练重要性取决于特征集合，不能将S12视作唯一确定机制。
2. S12在两个同时有成功和失败的architecture内方向一致，但组内精确检验和BH结果仍不足以证明普适。
3. PAE/pLDDT/geometry提供部分额外排序信息：RF三密度＋结构architecture-out AUC=0.9375，AP=0.8611；但0.5阈值下TP2/FN2/FP0/TN20，漏掉一半成功。6维结构RF可TP3/FP0，但AUC仅0.7625，仍存在极低分成功构建。不能仅凭一个指标判定稳定优越。
4. 局部特征没有一致改善泛化：full RF architecture-out AUC=0.8500、TP1/FN3，比三密度基线漏判更多；full XGB明显不稳定。原CCR01失效，新候选多受architecture影响。
5. 现有证据更符合“结构完整性相关模式与architecture背景混合”。无法从这些观察数据因果拆分两者。具体局部信号明显依赖排列；non-local density相对更有跨组一致性。
6. 可量化的迁移变化：三密度RF AUC 0.9125→0.86875，下降0.04375；balanced accuracy 0.825→0.750，下降0.075；precision 0.600→0.375；recall维持0.75，误报2→5。其他模型并非统一下降：9维CP＋结构RF AUC 0.9125→0.9375，增加0.025，但recall 0.75→0.50。这些是当前留出任务的实际差异，不能外推为任意新排列的固定下降幅度。
7. 实用最简基线仍推荐3密度浅层RF（3维、无RNA与局部筛选）。若只要最简连续排序，可用1维S12 Logistic：architecture-out AUC=0.8875，召回4/4，但误报9/20，适合宽松初筛。9维结构RF作为辅助排序对照保留研究价值。

留出G03时，三密度RF组内AUC=0.7857（TP3/FP5），低于拼接全部组的0.86875。留出G04时仅1个成功样本，组内AUC=1但该成功评分未到0.5，被漏判。后者清楚说明AUC=1不能解释为模型成熟。

从固定留出预测中逐一去掉一个成功样本（不重训），三密度RF architecture-out AUC可在0.850–0.925间变化，9维RF可在0.9167–1.000间变化。这只是少阳性敏感性诊断，不是置信区间。未宣称性能变化达到统计显著。

## 7. PUF11与TRM敏感性

| subset | feature_set | model | validation | AUC | AP | PR_AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| plus_PUF11 | CP_density3 | RF | LOCO | 0.9205 | 0.5845 | 0.4798 | 0.8295 | 0.6 | 0.75 | 20 | 2 | 1 | 3 |
| plus_PUF11 | S12_only | RF | LOCO | 0.8864 | 0.4833 | 0.3994 | 0.8295 | 0.6 | 0.75 | 20 | 2 | 1 | 3 |
| plus_PUF11 | full | RF | LOCO | 0.8977 | 0.6528 | 0.6108 | 0.7273 | 0.6667 | 0.5 | 21 | 1 | 2 | 2 |
| plus_PUF11 | CP_density3 | RF | architecture_out | 0.8807 | 0.475 | 0.6542 | 0.7614 | 0.375 | 0.75 | 17 | 5 | 1 | 3 |
| plus_PUF11 | S12_only | RF | architecture_out | 0.8182 | 0.3409 | 0.583 | 0.7386 | 0.3333 | 0.75 | 16 | 6 | 1 | 3 |
| plus_PUF11 | full | RF | architecture_out | 0.9318 | 0.85 | 0.8417 | 0.5 | 0 | 0 | 22 | 0 | 4 | 0 |
| plus_TRM | CP_density3 | RF | LOCO | 0.925 | 0.6829 | 0.6088 | 0.8667 | 0.7143 | 0.8333 | 18 | 2 | 1 | 5 |
| plus_TRM | S12_only | RF | LOCO | 0.8833 | 0.569 | 0.5058 | 0.8667 | 0.7143 | 0.8333 | 18 | 2 | 1 | 5 |
| plus_TRM | full | RF | LOCO | 0.925 | 0.8357 | 0.8238 | 0.8083 | 0.8 | 0.6667 | 19 | 1 | 2 | 4 |
| plus_TRM | CP_density3 | RF | architecture_out | 0.8875 | 0.6042 | 0.7566 | 0.7917 | 0.5 | 0.8333 | 15 | 5 | 1 | 5 |
| plus_TRM | S12_only | RF | architecture_out | 0.8167 | 0.4557 | 0.6792 | 0.7667 | 0.4545 | 0.8333 | 14 | 6 | 1 | 5 |
| plus_TRM | full | RF | architecture_out | 0.975 | 0.9444 | 0.941 | 0.75 | 1 | 0.5 | 20 | 0 | 3 | 3 |

加入PUF11后，三密度RF architecture-out AUC约0.8807，接近主集；full RF虽然AUC升到0.9318，却漏掉全部4个成功样本，进一步表明小样本评分阈值的不稳定。此结果未纳入主PUF12结论。

加入TRM两种蛋白后，26种记录有6个成功，但其中WT核心与原成功骨架重复，独立成功architecture仍只有2组。TRM已按前次核对去掉额外N端7 aa，并平均三个RNA背景；属于事后预处理敏感性检查。所有同source/order记录整组留出；分组权重虽不泄漏，重复骨架仍可能改变样本权重，不能据此声称新增独立阳性证据。

## 8. 评分与交付

推荐RF文件CP_density3_RF.joblib采用全部24种主集重拟合，包含填补/标准化流程。输入必须按同一蛋白范围、同一CP定义从model→seed→construct聚合；端部标签与长度变化可能显著改变分数。完整训练范围及评分脚本可检查明显的边界外输入；范围内不代表已验证可迁移。

附加单特征Logistic公式：

p = sigmoid[-0.0744141 + 0.4747608 × (S12 − 1.2830477)/0.03884619]

该p为未校准分数。RF默认以0.5分类；不能根据已知结果临时裁剪序列或调阈值，再声称是独立测试性能。

输出：feature_importance、success_failure_distributions、architecture_feature_distributions、repeat_pair_heatmap、adjacent_interface_heatmap、CCR_heatmap、validation_comparison、heldout_probabilities，均含PNG和矢量PDF。图中C01等编号映射见plot_ID_mapping.csv；每张图数据均在对应表/NPZ中。完整特征、p/q、分组、逐构建预测、局部筛选日志、敏感性结果和可复现脚本一并提供。

所有结论限于现有小样本和当前AF预测流程；构建成功标签不能直接替代编辑效率/脱靶标签。
