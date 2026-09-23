"""Independent arithmetic/split audit and Chinese report from saved predictions."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix,roc_auc_score,average_precision_score
P=Path(__file__).resolve().parent;O=P/'results'
pr=pd.read_csv(O/'outer_predictions.csv');perf=pd.read_csv(O/'performance.csv')
source=pd.read_csv(P/'inputs/construct_site_summary.csv').set_index(['site','construct'])
cohort=pd.read_csv(P/'inputs/frozen_cohort.csv');members=set(cohort.construct_key)
assert len(members)==22
assert not pr.duplicated(['site','label_threshold_pp','validation','family','model','policy','fold','construct']).any()
assert set(pr.construct)==members
assert ((pr.score>=pr.threshold).astype(int)==pr.prediction).all()
rate=np.array([source.loc[(s,c),'delta_editing_pp'] for s,c in zip(pr.site,pr.construct)])
assert ((rate<pr.label_threshold_pp).astype(int)==pr.y).all()
for a in json.loads((O/'split_audit.json').read_text()):
    tr=set(a['train']);te=set(a['test']);assert tr.isdisjoint(te) and tr|te==members
    for b in a['inner']:assert set(b['train']).isdisjoint(b['test']) and set(b['train'])|set(b['test'])==tr
    if a['validation']=='Strict_position_out':
        import re
        pos=a['fold'][1:]
        assert all(pos not in re.findall(r'(?:^|_plus_)p(\d+)_',c) for c in tr)
        assert all(pos in re.findall(r'(?:^|_plus_)p(\d+)_',c) for c in te)
keys=['site','label_threshold_pp','validation','family','model','policy'];maxerr=0.
for key,g in pr.groupby(keys):
    assert set(g.construct)==members
    w=1/g.construct.map(g.construct.value_counts()).values
    tn,fp,fn,tp=confusion_matrix(g.y,g.prediction,labels=[0,1],sample_weight=w).ravel()
    check={'BA':.5*(tp/(tp+fn)+tn/(tn+fp)),'accuracy':(tp+tn)/22,'TP':tp,'TN':tn,'FP':fp,'FN':fn,'AUC':roc_auc_score(g.y,g.score,sample_weight=w),'AP':average_precision_score(g.y,g.score,sample_weight=w)}
    q=perf
    for k,v in zip(keys,key):q=q[q[k]==v]
    assert len(q)==1
    maxerr=max(maxerr,max(abs(q.iloc[0][k]-v) for k,v in check.items()))
assert maxerr<1e-10
labels=[]
for site in ['C871','C295']:
    d=source.loc[site].loc[sorted(members)];x=d.delta_editing_pp
    for t in sorted(perf[perf.site==site].label_threshold_pp.unique(),reverse=True):
        n=int((x<t).sum());labels.append(dict(site=site,threshold_pp=t,absolute_editing_pct=float(d.wt_editing_pct.iloc[0]+t),positive=n,negative=22-n,equivalent_lower_exclusive=float(x[x<t].max()),equivalent_upper_inclusive=float(x[x>=t].min()),eligible_headline=min(n,22-n)>=5))
lab=pd.DataFrame(labels);lab.to_csv(O/'label_scan.csv',index=False)
pf=perf.merge(lab[['site','threshold_pp','positive','negative','eligible_headline']],left_on=['site','label_threshold_pp'],right_on=['site','threshold_pp'])
dynamic=pf[(pf.policy=='training_only_threshold')&pf.eligible_headline&pf.model.ne('nested_model_selection')]
best=[];compare=[]
for site in ['C871','C295']:
    for validation in ['LOOCV','Position_group_out','Strict_position_out']:
        q=dynamic[(dynamic.site==site)&(dynamic.validation==validation)].sort_values(['BA','AUC','precision','family','model','label_threshold_pp'],ascending=[False,False,False,True,True,False],kind='stable')
        best.append(q.iloc[0].to_dict())
    chosen=next(b for b in best if b['site']==site and b['validation']=='Position_group_out')
    for validation in ['LOOCV','Position_group_out','Strict_position_out']:
        for policy in ['fixed_0.5','training_only_threshold']:
            q=pf[(pf.site==site)&(pf.label_threshold_pp==chosen['label_threshold_pp'])&(pf.family==chosen['family'])&(pf.model==chosen['model'])&(pf.validation==validation)&(pf.policy==policy)]
            compare.extend(q.to_dict('records'))
pd.DataFrame(best).to_csv(O/'exploratory_best_by_validation.csv',index=False)
cmp=pd.DataFrame(compare);cmp.to_csv(O/'group_selected_candidate_comparison.csv',index=False)
selected_pr=[];throws=[]
for b in best:
    if b['validation']!='Position_group_out':continue
    q=pr[(pr.site==b['site'])&(pr.label_threshold_pp==b['label_threshold_pp'])&(pr.family==b['family'])&(pr.model==b['model'])]
    selected_pr.append(q)
    for val,g in q[q.policy=='training_only_threshold'].groupby('validation'):
        t=g.drop_duplicates('fold').threshold
        throws.append(dict(site=b['site'],validation=val,label_threshold_pp=b['label_threshold_pp'],family=b['family'],model=b['model'],score_threshold_median=t.median(),score_threshold_min=t.min(),score_threshold_max=t.max()))
pd.concat(selected_pr).to_csv(O/'selected_candidate_predictions.csv',index=False);pd.DataFrame(throws).to_csv(O/'selected_score_thresholds.csv',index=False)
joint=pd.read_csv(O/'joint_label_performance.csv');js=pd.read_csv(O/'joint_label_selection.csv')
jp=pd.read_csv(O/'joint_label_predictions.csv')
assert ((jp.score>=jp.threshold).astype(int)==jp.prediction).all()
assert ((np.array([source.loc[(s,c),'delta_editing_pp'] for s,c in zip(jp.site,jp.construct)])<jp.label_threshold_pp).astype(int)==jp.y).all()
for (site,val),g in jp.groupby(['site','validation']):
    assert set(g.construct)==members
    w=1/g.construct.map(g.construct.value_counts()).values
    tn,fp,fn,tp=confusion_matrix(g.y,g.prediction,labels=[0,1],sample_weight=w).ravel()
    expected=.5*(tp/(tp+fn)+tn/(tn+fp))
    assert abs(joint[(joint.site==site)&(joint.validation==val)].BA.iloc[0]-expected)<1e-10
for site,t in [('C871',-25),('C295',-20)]:
    f=O/f'{site}_labels.csv';d=pd.read_csv(f)
    d=d.rename(columns={'label_large_decrease':'reference_original_label'})
    d['reference_original_threshold_pp']=t;d.to_csv(f,index=False)
js.groupby(['site','validation','label_threshold_pp']).size().rename('n_folds').reset_index().to_csv(O/'joint_label_frequency.csv',index=False)
modelbest=dynamic[dynamic.validation=='Position_group_out'].sort_values(['BA','AUC'],ascending=False).groupby(['site','model'],sort=False).head(1)
modelbest.to_csv(O/'best_per_model_exploratory.csv',index=False)
audit={'unique_constructs':22,'prediction_rows':len(pr),'metric_rows':len(perf),'all_splits_disjoint':True,'labels_recomputed_from_experiment':True,'all_score_decisions_verified':True,'max_metric_error':maxerr,'no_new_experiments':True}
(O/'validation_checks.json').write_text(json.dumps(audit,indent=2))
def table(df,cols):
    def fmt(v):
        if isinstance(v,(float,np.floating)):return f'{v:.4f}'
        return str(v)
    return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df[cols].itertuples(index=False,name=None))
text='''# 22个固定入组TRM突变体：C871与C295动态标签及分类方法比较

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

'''+table(lab,['site','threshold_pp','absolute_editing_pct','positive','negative','eligible_headline'])+'\n\n## 探索性最佳固定方案\n\n每种验证分别扫描得到的最高值存在多候选选择偏差；不同标签定义了不同预测任务。不能把改变标签后的BA升高单独归因于模型改善。\n\n'+table(pd.DataFrame(best),['site','validation','label_threshold_pp','family','model','BA','AUC','precision','recall','TP','TN','FP','FN'])+'\n\n## 按位置组合成绩选出的同一个候选，在三种验证下的表现\n\n'+table(cmp,['site','validation','policy','label_threshold_pp','family','model','BA','AUC','precision','recall'])+'\n\n各外层折训练内选择的分数阈值：\n\n'+table(pd.DataFrame(throws),['site','validation','score_threshold_median','score_threshold_min','score_threshold_max'])+'\n\n## 完整训练内动态选择流程\n\n以下评估整个自适应流程。不同外层折可能选择不同编辑标签，因此其BA、AUC不能解释为某一个统一下降标准的固定终点表现；原始逐折阈值与预测均单独提供。\n\n'+table(joint,['site','validation','BA','accuracy','AUC','precision','recall','TP','TN','FP','FN'])+'\n\n## 不同训练方法的探索性比较\n\n每行在位置组合留出中选择该方法最好的标签/输入，仅作为探索性比较。\n\n'+table(modelbest,['site','model','label_threshold_pp','family','BA','AUC','precision','recall'])+'''

## 解释边界与复现

所有候选使用同一批既有实验，尚无独立新突变体测试。WT与多数突变体seed有限，WT及实验重复的不确定性没有传播；没有全流程置换检验。22个名单自身来自C388入组，结论不能外推到全部TRM突变体。旧30个突变体的指标不能与本轮直接当作同任务的提升。没有用全数据重拟合成绩报告泛化，也未发布部署模型。

运行 `python run.py` 重算模型，默认4个CPU进程，可用TRM_JOBS调整。运行 `python report.py` 重算指标、验证拆分并生成本报告。inputs包含复现缓存；sources保留ContactSeek源码及许可证；results/provenance.json记录版本和输入哈希。保存的逐样本分数可复核所有指标。run.log记录训练过程。

训练方法遵循scikit-learn技能中的训练折内预处理与嵌套验证流程。软件工作流参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. https://doi.org/10.48550/arXiv.2609.00065 。该文献用于记录分析工具来源，不作为本项目生物学结论的证据。
'''
strict_compare=[]
for b in best:
    if b['validation']!='Strict_position_out':continue
    q=pf[(pf.site==b['site'])&(pf.label_threshold_pp==b['label_threshold_pp'])&(pf.family==b['family'])&(pf.model==b['model'])&(pf.policy=='training_only_threshold')]
    strict_compare.extend(q.to_dict('records'))
sc=pd.DataFrame(strict_compare);sc.to_csv(O/'strict_selected_candidate_comparison.csv',index=False)
text+='\n\n## 严格位置验证选出的候选，在全部验证下的表现\n\n同样属于事后探索性选择，需与前文位置组合选出的候选一起阅读。\n\n'+table(sc,['site','label_threshold_pp','family','model','validation','BA','AUC','precision','recall','TP','TN','FP','FN'])
(P/'REPORT.md').write_text(text)
print(json.dumps(audit));print(pd.DataFrame(best).to_string(index=False));print(joint.to_string(index=False))
