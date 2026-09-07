# 分析沿革与解释边界

1. **原14种构建**：3种成功，11种失败。其中2种为PUF11。成功标签为R123/R567/R567/R678、R123/R567/R6R7/R5678、R123/R567/R5loop67/loopR678。原CP summary和CCR发现结果保存于previous/results与extended。
2. **第二批**：实际15种设计，12种为PUF12，3种为15-repeat。用户补充实验标签：该批PUF12仅ID1成功。new_batch中的早期冻结评分先于这次标签更新，不能按其分数认定实验成功。
3. **合并24种PUF12**：成功4、失败20；每构建4 seed×5 model。combined12为原固定模型与24折LOCO；rf_tuning为内层3折调参、外层LOCO评估。高维/更深模型没有整体改善。
4. **91–96 TRM验证**：6个复合物实际为2种蛋白×3个RNA。WT核心与旧成功骨架一致，突变体相对本批WT为S304G/Y305V；两者均有额外N端7 aa。剔除额外片段的评分变化属于事后敏感性诊断；不同RNA和model不能视为独立成功构建。
5. **最新跨architecture验证**：24种PUF12为主集，PUF11与TRM分别做敏感性分析。实际source/order划分7组；成功仅出现在其中2组。相同source/order的loop变体整组留出。旧CCR01协同/正向信号未复现。

## 不应混淆的结果

- LOCO、整组留出、训练集重拟合分数是不同对象。
- 原14种、合并24种、加PUF11、加TRM的样本量与独立背景数不同。
- architecture-out与source-order-out当前完全相同。
- AP使用阶梯加权；PR_AUC使用梯形面积，并列分数可使二者不同。
- 全阳性TRM验证不能计算有效ROC-AUC或特异度。
- R123仍为共同前部背景；没有验证替换前部R123后的性能。
- 全数据探索出的残基与CCR不能直接加入LOCO。当前模型脚本在每轮训练集内重新筛选。
- 高CP/pLDDT与构建成功的关联不能证明因果机制，也不能直接预测RNA编辑效率。

## 模型文件状态

`docs/model_load_audit.json`记录当前保存模型的读取检查。历史版本保留供溯源；预测新设计优先使用architecture_validation中的模型。不能只依据旧模型文件名推断其训练集或特征定义。

历史extended/final_summary_RF.joblib存在截断。仓库使用已重建且与保存的第二批RF评分逐项一致的new_batch/restored_frozen_RF.joblib恢复该文件；旧字节另存legacy_corrupt_final_summary_RF.joblib.bin。详见model_repair.json。此修复未改变当前24构建模型或预测结果。
