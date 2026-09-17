import pathlib,json,hashlib,sys,shutil,zipfile
import numpy as np,pandas as pd,joblib
from sklearn.metrics import roc_auc_score,confusion_matrix,average_precision_score
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=pathlib.Path('c388_analysis');a=pd.read_csv(O/'training_data.csv');m=pd.read_csv(O/'metrics.csv');p=pd.read_csv(O/'predictions.csv');f=pd.read_csv(O/'construct_features.csv');manifest=json.load(open(O/'manifest.json'));rng=np.random.default_rng(20260908)
# Main per-construct table contains only held-out new-model predictions.
b=a[['construct','mutation_group','mutations','C388_mean','C388_min','C388_max','work','zero_variant_record']].copy()
configs=[('frozen','transfer','CP_density3_RF','old_RF'),('retrained','LOCO','CP3_LR','new_CP3_LR_LOCO'),('retrained','LOCO','CP_Structure9_LR','new_CPstruct_LR_LOCO'),('retrained','mutation_group_out','CP_Structure9_LR','new_CPstruct_LR_groupout'),('nested_selection','LOCO','inner_select_10_candidates','nested_LOCO'),('nested_selection','mutation_group_out','inner_select_10_candidates','nested_groupout')]
for fam,val,model,col in configs:
 t=p[(p.family==fam)&(p.validation==val)&(p.model==model)&(p.scope=='matched512')].set_index('construct');b[col+'_score']=b.construct.map(t.score);b[col+'_prediction']=(b[col+'_score']>=.5).astype(int)
b.sort_values('C388_mean',ascending=False).to_csv(O/'per_construct_comparison.csv',index=False)
# Conditional bootstrap uncertainty of existing held-out scores; no retraining or independence claim.
boots=[]
for fam,val,model,col in configs:
 vals=b[col+'_score'].to_numpy();y=b.work.to_numpy();stats=[]
 for _ in range(5000):
  ix=np.r_[rng.choice(np.where(y==1)[0],sum(y==1),replace=True),rng.choice(np.where(y==0)[0],sum(y==0),replace=True)];stats.append(roc_auc_score(y[ix],vals[ix]))
 boots.append(dict(method=col,AUC=roc_auc_score(y,vals),conditional_bootstrap_low=np.quantile(stats,.025),conditional_bootstrap_high=np.quantile(stats,.975),n_resamples=5000))
pd.DataFrame(boots).to_csv(O/'conditional_score_bootstrap.csv',index=False)
# Direct-transfer sensitivity excludes WT already represented as a core in historical training.
froz=[]
for sub,mask in [('exclude_WT',b.construct!='PUF12-9'),('exclude_zero_records',~b.zero_variant_record)]:
 yy=b.loc[mask,'work'];vv=b.loc[mask,'old_RF_score'];tn,fp,fn,tp=confusion_matrix(yy,vv>=.5,labels=[0,1]).ravel();froz.append(dict(subset=sub,n=sum(mask),AUC=roc_auc_score(yy,vv),TN=tn,FP=fp,FN=fn,TP=tp))
pd.DataFrame(froz).to_csv(O/'frozen_sensitivity.csv',index=False)
# Unlabelled construct: full-data fit scores, clearly distinct from evaluation.
u=f[(f.scope=='matched512')&(~f.construct.isin(a.construct))].copy();ur=[]
for _,r in u.iterrows():
 out={'construct':r.construct,'true_label':'unknown'}
 for key,typ in [('CP3','LR'),('CP_Structure9','LR')]:out[key+'_'+typ+'_refit_score']=float(joblib.load(O/(key+'_'+typ+'.joblib')).predict_proba(r[manifest['feature_sets'][key]].to_numpy(float)[None,:])[0,1])
 ur.append(out)
pd.DataFrame(ur).to_csv(O/'unlabelled_predictions.csv',index=False)
# Fitted LR formula and effect direction; refitted coefficients are descriptive, not CV importance.
coef=[]
for key in ['CP3','CP_Structure9']:
 mod=joblib.load(O/(key+'_LR.joblib'))
 for nm,mu,sd,w in zip(manifest['feature_sets'][key],mod[1].mean_,mod[1].scale_,mod[-1].coef_[0]):coef.append(dict(model=key+'_LR',feature=nm,mean=mu,scale=sd,standardized_coefficient=w,intercept=mod[-1].intercept_[0]))
pd.DataFrame(coef).to_csv(O/'logistic_scoring_coefficients.csv',index=False)
# Descriptive mutation-position performance, using held-out predictions.
gs=[]
for fam,val,model,col in configs:
 for g,t in b.groupby('mutation_group'):
  yy=t.work.to_numpy();vv=t[col+'_score'].to_numpy();tn,fp,fn,tp=confusion_matrix(yy,vv>=.5,labels=[0,1]).ravel();gs.append(dict(method=col,group=g,n=len(t),n_work=int(sum(yy)),AUC=roc_auc_score(yy,vv) if len(set(yy))==2 else np.nan,TN=tn,FP=fp,FN=fn,TP=tp))
pd.DataFrame(gs).to_csv(O/'per_mutation_group_metrics.csv',index=False)
# Scoped scientific figure settings; native source data is packaged.
with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42}):
 fig,axs=plt.subplots(1,3,figsize=(14,4.6),layout='constrained');cols=['old_RF_score','new_CPstruct_LR_LOCO_score','nested_LOCO_score'];titles=['Frozen architecture RF','Refit CP + structure LR: LOCO','Nested model selection: LOCO']
 for ax,col,title in zip(axs,cols,titles):
  for yy,color,marker in [(0,'#D55E00','x'),(1,'#0072B2','o')]:
   t=b[b.work==yy];ax.scatter(t.C388_mean*100,t[col],c=color,marker=marker,s=48,label='Non-work (<50%)' if yy==0 else 'Work (>=50%)')
  ax.axvline(50,color='gray',ls=':');ax.axhline(.5,color='gray',ls=':');ax.set(xlim=(-3,103),ylim=(-.03,1.03),xlabel='Mean C388 editing (%)',ylabel='Held-out / frozen model score',title=title);ax.legend(fontsize=8,loc='lower right')
 fig.suptitle('31 constructs; 22 work / 9 non-work; one point per protein',fontsize=13);fig.savefig(O/'C388_prediction_comparison.png',dpi=250);fig.savefig(O/'C388_prediction_comparison.pdf');plt.close(fig)
 fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained');keys=['CP3_LR','CP3_RF','CP_Structure9_LR','CP_Structure9_RF','RNA4_LR','RNA4_RF','inner_select_10_candidates'];labels=['CP3 LR','CP3 RF','CP + structure LR','CP + structure RF','RNA LR','RNA RF','Nested selection']
 for j,val in enumerate(['LOCO','mutation_group_out']):
  t=m[(m.validation==val)&m.model.isin(keys)&m.family.isin(['retrained','nested_selection'])].set_index('model').loc[keys];ax.scatter(t.AUC,np.arange(len(keys))+(j-.5)*.16,c=['#0072B2','#D55E00'][j],marker=['o','s'][j],label=val,s=65)
 ax.axvline(.5,ls=':',c='gray');ax.set(yticks=range(len(keys)),yticklabels=labels,xlim=(0,1),xlabel='ROC-AUC',title='New-label validation: model selection remains unstable');ax.invert_yaxis();ax.legend(loc='lower left');fig.savefig(O/'C388_validation.png',dpi=250);fig.savefig(O/'C388_validation.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(12,4.5),layout='constrained')
 for ax,feat,title in zip(axs,['pp_nonlocal12_high_per_res','Rg_full_length_normalized','pr_cp_per_nt'],['High non-local CP density (S12)','Normalized radius of gyration','Protein-RNA CP per nucleotide']):
  for yy,color,marker in [(0,'#D55E00','x'),(1,'#0072B2','o')]:
   v=a.loc[a.work==yy,feat];xx=np.full(len(v),yy)+rng.uniform(-.09,.09,len(v));ax.scatter(xx,v,c=color,marker=marker,s=44);ax.plot([yy-.18,yy+.18],[v.mean()]*2,c='black',lw=2)
  ax.set(xticks=[0,1],xticklabels=['Non-work (9)','Work (22)'],title=title)
 fig.suptitle('Construct means; horizontal bars show class means');fig.savefig(O/'C388_feature_distributions.png',dpi=250);fig.savefig(O/'C388_feature_distributions.pdf');plt.close(fig)
# Compact markdown report
sel=m[((m.family=='frozen')&(m.model=='CP_density3_RF')&(m.scope=='matched512'))|((m.family=='retrained')&(m.model.isin(['CP3_LR','CP3_RF','CP_Structure9_LR','CP_Structure9_RF','RNA4_LR'])))|(m.family=='nested_selection')].copy()
def mdtable(t):
 return '| '+' | '.join(t.columns)+' |\n| '+' | '.join(['---']*len(t.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in t.to_numpy())
mt=sel[['family','validation','model','AUC','balanced_accuracy','precision','recall','TN','FP','FN','TP']].copy()
for c in ['AUC','balanced_accuracy','precision','recall']:mt[c]=mt[c].map(lambda v:f'{v:.3f}')
report='''# C388≥50%：PUF work/non-work 重新定义与模型验证

## 固定阈值模型清理（2026-09-17）

按现有折外结果清理可复用权重。采用保守标准：C388≥50%、模型分数阈值0.5时，LOCO与突变位置分组留出的AUC均低于0.60，删除对应权重。该标准用于本次文件整理，不代表统计显著性或部署门槛。

| 处理 | 模型 | LOCO AUC / BA | 位置分组留出 AUC / BA |
|---|---|---|---|
| 删除权重 | RNA4_LR | 0.566 / 0.619 | 0.556 / 0.573 |
| 删除权重 | RNA4_RF | 0.439 / 0.374 | 0.500 / 0.452 |
| 删除权重 | Structure6_RF | 0.561 / 0.485 | 0.535 / 0.540 |
| 保留研究基线 | CP3_LR | 0.722 / 0.753 | 0.697 / 0.674 |
| 保留排序对照 | CP_Structure9_LR | 0.778 / 0.606 | 0.798 / 0.707 |

其余5个C388权重保留为中等或不稳定信号的探索性候选，逐项指标见 [model_retention.json](model_retention.json)。原architecture模型在其自身任务上有预测信号，连同用于迁移复核的reference_models一起保留。C295固定相对WT标签中的局部ΔCP保留；其他C295方案仅有比较记录，未发现单独的最终模型权重。

完整指标、折外预测、失败记录、数据、动态阈值结果与复现代码继续保留。analyze.py仍评估全部原候选，最终权重导出遵循model_retention.json，避免改变历史比较与嵌套选择范围。删除前文件可从Git提交历史恢复。所有保留结果仍受小样本、模型筛选及缺少独立验证的限制。

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
'''+mdtable(mt)+'''

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
'''
(O/'C388_analysis_report.md').write_text(report)
# Freeze reproducibility dependencies rather than relying on previous scratch state.
ref=O/'reference_models';ref.mkdir(exist_ok=True)
old=pathlib.Path('architecture_validation')
for q in old.glob('*_only_*.joblib'):shutil.copyfile(q,ref/q.name)
for q in old.glob('CP_density3_*.joblib'):shutil.copyfile(q,ref/q.name)
for nm in ['main_features.csv']:
 if (old/nm).exists():shutil.copyfile(old/nm,ref/nm)
if pathlib.Path('combined12/manifest.json').exists():shutil.copyfile('combined12/manifest.json',ref/'old_feature_manifest.json')
versions={};
for nm in ['numpy','pandas','scipy','sklearn','matplotlib','joblib','openpyxl']:
 mod=__import__(nm);versions[nm]=mod.__version__
json.dump(versions,open(O/'environment_versions.json','w'),indent=2)
raw=[]
for q in pathlib.Path('upload').glob('*'):
 if q.suffix not in ['.zip','.xlsx']:continue
 h=hashlib.sha256();
 with open(q,'rb') as fp:
  for block in iter(lambda:fp.read(1024*1024),b''):h.update(block)
 raw.append(dict(file=q.name,bytes=q.stat().st_size,sha256=h.hexdigest()))
if raw:json.dump(raw,open(O/'raw_input_manifest.json','w'),ensure_ascii=False,indent=2)
print('REPORT DONE');print(pd.DataFrame(boots).to_string(index=False));print('Unlabelled',ur)
