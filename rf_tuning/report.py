import json,zipfile
from pathlib import Path
import pandas as pd
O=Path('rf_tuning');m=pd.read_csv(O/'nested_LOCO_metrics.csv');p=pd.read_csv(O/'nested_LOCO_predictions.csv');b=pd.read_csv('combined12/model_metrics.csv');b=b[(b.validation=='combined_LOCO')&b.model.str.endswith('_RF')].copy();b['method']='Original fixed RF';b['threshold']='fixed_0.5';m['method']='Nested tuned RF';cmp=pd.concat([b,m],ignore_index=True);cols=['method','model','threshold','AUC','AP','balanced_accuracy','accuracy','TP','FN','FP','TN'];cmp[cols].to_csv(O/'baseline_vs_tuned.csv',index=False)
oldp=pd.read_csv('combined12/heldout_predictions.csv');oldp=oldp[(oldp.validation=='combined_LOCO')&oldp.model.str.endswith('_RF')];wide=p.pivot(index=['construct','batch','design_id','success'],columns='scope',values='score').reset_index() if False else p[['construct','batch','design_id','success']].drop_duplicates().merge(p.pivot(index='construct',columns='scope',values='score'),on='construct')
for s in ['CP_density3','CP_summary14','CP_structure17']:
 wide=wide.merge(oldp[oldp.model==s+'_RF'][['construct','score']].rename(columns={'score':'baseline_'+s}),on='construct')
wide.to_csv(O/'per_construct_comparison.csv',index=False)
sel=p.groupby(['scope','candidate_id','feature_set']).size().reset_index(name='selected_outer_folds');sel.to_csv(O/'selection_frequency.csv',index=False)
def tab(df):
 def fmt(v):return f'{v:.4f}' if isinstance(v,float) else str(v)
 return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)])
manifest=json.load(open(O/'final_model_manifest.json'));joint=next(x for x in manifest if x['scope']=='joint')
text='''# 两批12-repeat PUF：Random Forest调参结果

## 主要结论

本轮搜索未整体超过原浅层RF。3个CP密度特征调参后AUC从0.9125变为0.9000，分类计数保持TP3/FN1/FP2/TN18；14个CP summary特征AUC从0.8500升到0.8750，成功检出从2/4升到3/4；17个特征AUC从0.9250降到0.8125。内层联合选择特征与参数的AUC为0.8625。训练内选阈值未改善任何流程的balanced accuracy。

本批1号仍正确识别，4号与7号仍误报。建议继续保留原CP_density3_RF（300棵树、深度2、叶节点最少2个样本、sqrt特征抽样、balanced类别权重）作为当前基线。调参模型保留作研究比较，不声称优于原模型。以上变化为描述性比较，未证明统计显著差异。

## 数据与验证

共24种构建，4种work、20种不work。上一批12种中3种成功；本批12种仅1号成功。每种构建先聚合多个seed/model，训练中每种构建权重相同。沿用前次特征提取定义与标签，不包含11或15-repeat构建。

- 外层：24折LOCO，每次整种构建留出。
- 内层：仅外层训练集，3折分层CV。每折至少包含一个阳性。48个候选组合（16组RF参数×3组特征），以平均内层ROC AUC选择，平分时以AP再以候选编号决胜。
- 参数范围：深度1/2/3/4/不限，min_samples_leaf 1/2/3，max_features sqrt/全部，class_weight balanced/无/balanced_subsample；采用预先指定的16种组合，并非完整笛卡尔积。完整组合见search_protocol.json。
- 内层搜索每森林100棵树；外层预测及最终模型使用300棵树；随机种子2026。搜索树数与重拟合树数的差异属于本次固定训练流程。
- 分别报告每个特征组内调参，以及在内层联合选择特征组和参数的joint流程。joint评估包含特征组选择。
- 固定阈值0.5与内层OOF选择阈值分别报告；后者在0.10–0.80、步长0.05网格上最大化balanced accuracy，平分选最接近0.5者。没有用外层留出标签选阈值。
- 本次未新增CCR；沿用3个接触密度、14个CP summary、17个CP+结构特征，便于单独评估RF参数调整。

## 原模型与调参后比较

'''+tab(cmp[cols])+'''

TP=成功且预测成功；FN=成功但漏判；FP=失败但误判为成功；TN=失败且正确排除。
全部预测失败也能得到83.3% accuracy，因此应同时看成功召回、误报与balanced accuracy。阈值变化不改变AUC/AP。

## 本批逐构建LOCO评分

'''+tab(wide[wide.batch=='new'].sort_values('design_id'))+'''

上表均为外层留出评分，可与实验标签比较；分数尚未校准为成功概率。

## 最终部署模型

所有24种构建内层搜索后重新拟合，供新设计使用。最终全数据选定参数不等同于每个LOCO折都选择相同参数；最终模型自身没有新的独立测试成绩。

联合选择最终配置：

```json
'''+json.dumps(joint,ensure_ascii=False,indent=2)+'''
```

每组特征的独立模型配置见final_model_manifest.json；各外层折选中参数见nested_LOCO_predictions.csv和selection_frequency.csv。完整内层搜索日志同时保留。

## 使用

```bash
python rf_tuning/predict.py --features new_construct_features.csv --scope joint --out predictions.csv
```

默认0.5阈值；使用训练内选择阈值可添加`--threshold inner`。输入必须是按原流程聚合到构建级的同定义特征；不能将单个seed/model当作独立构建。要求12-repeat且与训练域相近。旧特征提取脚本保留在此前PUF_12repeat_combined_training.zip中。

## 解释边界

只有4个独立成功构建，超参数、特征组和阈值选择容易不稳定。嵌套LOCO避免本次参数搜索直接利用外层测试标签，但两批数据已用于此前分析及特征设计，结果仍属探索性内部验证。任何新模型优势都应在后续新批次锁定参数后复核。

feature_importance文件来自全数据最终RF，属于训练模型的杂质重要性；相关CP特征会分摊重要性，不能直接证明因果结构机制。本次不根据全数据调参最高分报告泛化性能。
'''
(O/'RF_tuning_report.md').write_text(text)
with zipfile.ZipFile('PUF_RF_tuning.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in O.iterdir():
  if f.is_file() and f.suffix!='.log':z.write(f,str(f))
 for f in ['combined12/training_data.csv','combined12/model_metrics.csv','combined12/heldout_predictions.csv','extended/requirements.txt']:z.write(f,f)
with zipfile.ZipFile('PUF_RF_tuning.zip') as z:assert z.testzip() is None
print(cmp[cols].to_string(index=False));print('FINAL',json.dumps(joint));print(wide[wide.batch=='new'].sort_values('design_id').to_string(index=False))
