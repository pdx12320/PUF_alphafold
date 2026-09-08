from pathlib import Path
import json, html
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent
out=root/'results'
data=pd.read_csv(root/'inputs/construct_site_summary.csv')
rep=pd.read_csv(root/'inputs/experimental_replicates.csv')
val=pd.read_csv(out/'training_only_score_threshold_validation.csv')
scan=pd.read_csv(out/'exploratory_threshold_scan.csv')
labels=pd.read_csv(out/'candidate_labels.csv')
pred=pd.read_csv(out/'validated_predictions.csv')
old=pd.read_csv(root/'audit_reference/exploratory_oof_scores.csv')
new=pd.read_csv(out/'exploratory_oof_scores.csv')
checks=[]
baseline=[]
for site,ss in data.groupby('site'):
    wt=float(rep.loc[rep.construct=='wt_puf12_9',site].mean()*100)
    assert np.allclose(ss.wt_editing_pct,wt)
    assert np.allclose(ss.editing_pct-ss.wt_editing_pct,ss.delta_editing_pp)
    actual=rep.groupby('construct')[site].mean()*100
    assert np.allclose(ss.editing_pct,actual.loc[ss.construct])
    baseline.append(dict(site=site,wt_editing_pct=wt,n=int(len(ss)),n_delta_ge_0=int((ss.delta_editing_pp>=0).sum()),n_delta_lt_0=int((ss.delta_editing_pp<0).sum()),delta_min_pp=ss.delta_editing_pp.min(),delta_max_pp=ss.delta_editing_pp.max(),zero_threshold_eligible=False))
    # Compare independently rerun scores for each identical label partition.
    om={}
    for (model,t),g in old[old.site==site].groupby(['model','label_threshold_pct']):
        g=g.sort_values('construct');om[(model,tuple(g.actual_class))]=g.score.to_numpy()
    nm={}
    for (model,t),g in new[new.site==site].groupby(['model','label_threshold_delta_pp']):
        g=g.sort_values('construct');key=(model,tuple(g.actual_class));nm[key]=g.score.to_numpy()
    assert set(om)==set(nm)
    maxdiff=max(float(np.abs(om[k]-nm[k]).max()) for k in om)
    assert maxdiff<1e-10
    checks.append(dict(site=site,matching_model_label_combinations=len(om),max_oof_score_difference=maxdiff,all_partitions_match=True))
base=pd.DataFrame(baseline);base.to_csv(out/'wt_baselines_and_zero_threshold.csv',index=False)
summary=[]
for site in ['C388','C871','C295']:
    v=val[val.site==site];r=v.iloc[0]
    s=scan[(scan.site==site)&(scan.model==r.model)&np.isclose(scan.label_threshold_delta_pp,r.label_threshold_delta_pp)].iloc[0]
    wt=base.set_index('site').loc[site,'wt_editing_pct']
    summary.append(dict(site=site,WT_editing_pct=wt,threshold_delta_pp=r.label_threshold_delta_pp,equivalent_lower_exclusive_pp=s.equivalent_delta_lower_exclusive_pp,equivalent_upper_inclusive_pp=s.equivalent_delta_upper_inclusive_pp,corresponding_absolute_threshold_pct=wt+r.label_threshold_delta_pp,model=r.model,n_delta_ge_threshold=int(r.n_high),n_delta_lt_threshold=int(r.n_low),LOOCV_BA=v.set_index('validation').loc['LOOCV','BA'],position_group_BA=v.set_index('validation').loc['Position_group_out','BA']))
summary=pd.DataFrame(summary);summary.to_csv(out/'wt_relative_summary.csv',index=False)
for frame in [labels,pred]:
    frame['class_1_meaning']='delta >= threshold: smaller decrease or increase'
    frame['class_0_meaning']='delta < threshold: larger decrease'
pred=pred.merge(summary[['site','threshold_delta_pp','WT_editing_pct']],on='site',validate='many_to_one')
assert np.all(pred.actual_class==(pred.delta_editing_pp>=pred.threshold_delta_pp).astype(int))
pred.to_csv(out/'wt_relative_per_construct_predictions.csv',index=False)
labels.to_csv(out/'candidate_labels.csv',index=False)
audit=dict(label_definition='delta_editing_pp = mutant_editing_pct - site_specific_WT_editing_pct; class 1 iff delta >= threshold_delta_pp',units='percentage points, not relative percent or fold change',replicate_means_verified=True,paired_delta_verified=True,all_30_constructs_included=True,score_comparison=checks,limitations=['Label threshold and model selected globally; exploratory selection bias remains.','Score thresholds and learned contact/CCR features fitted within training folds.','No independent validation or full-pipeline permutation test performed.','WT means treated as fixed; experimental baseline uncertainty not propagated.','Position groups follow prior pipeline, with double mutants grouped by combined positions; this does not guarantee every constituent position is unseen.'])
(out/'validation_checks.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
def table(d): return d.to_html(index=False,border=0,float_format=lambda x:f'{x:.5g}')
display=summary.copy()
display['LOOCV_BA']=display.LOOCV_BA.map(lambda x:f'{x:.1%}')
display['position_group_BA']=display.position_group_BA.map(lambda x:f'{x:.1%}')
display=display.rename(columns={'site':'位点','WT_editing_pct':'WT编辑率(%)','threshold_delta_pp':'Δ分界(pp)','equivalent_lower_exclusive_pp':'等价下界(不含)','equivalent_upper_inclusive_pp':'等价上界(含)','corresponding_absolute_threshold_pct':'对应绝对编辑率(%)','model':'模型','n_delta_ge_threshold':'Δ≥分界数','n_delta_lt_threshold':'Δ<分界数','LOOCV_BA':'留一平衡准确率','position_group_BA':'位置分组平衡准确率'})
body='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>PUF 相对 WT 动态阈值分类重分析</title><style>body{font:16px/1.7 system-ui,sans-serif;max-width:1400px;margin:36px auto;padding:0 24px;color:#172b40}h1,h2{color:#126b70}table{border-collapse:collapse;font-size:13px;width:100%;margin:20px 0}td,th{padding:9px;border-bottom:1px solid #d6e2e8;text-align:left}th{background:#ecf5f5}.box{background:#eef6fa;padding:16px;border-left:4px solid #187b85} .scroll{overflow:auto}</style><h1>PUF：相对 WT 编辑变化的动态阈值分类</h1>'''
body+='''<p class="box"><b>已重新运行三个位点全部30个突变体的分类扫描与验证。</b>标签改为 Δ编辑效率＝突变体平均编辑率−对应位点WT平均编辑率，单位为百分点（pp）。正值表示提高，负值表示下降。未使用相对百分比变化或倍数变化。</p><h2>结果</h2><p>类别1：Δ≥分界；类别0：Δ&lt;分界。当前候选分界均为负数，类别1包含下降较少和提高的构建，类别0表示下降较多。C388通常希望保留活性；C871/C295通常希望降低编辑。单个位点分类分数不等于联合特异性优化成功率。</p><div class="scroll">'''+table(display)+'</div>'
body+='''<h2>为何准确率可以与上轮相同</h2><p>同一位点的WT均值为固定常数。绝对编辑率≥t 与 Δ≥t−WT 给出相同标签。前后两轮都扫描了所有满足每类至少5个样本的相邻观测值分界，故候选分组集合相同。本轮在Δ尺度重新选择分界，并重新拟合模型、计算折外预测；按标签分组逐一比较两轮预测以验证这一点。阈值数值变化不代表模型性能提升。</p>'''+table(pd.DataFrame(checks))
body+='''<h2>与WT直接比较：Δ=0</h2><p>三个终点达到或超过WT的构建分别只有少数几个，均未达到每类至少5个的预设条件。因此未将Δ=0作为合格动态候选，也未给出稳健的“提高/下降”分类结论。</p>'''+table(base)
body+='''<h2>验证范围</h2><ul><li>模型与编辑变化分界仍在本批数据中筛选，表中成绩属于探索性结果。</li><li>模型分数阈值在训练折的内部验证中选择；接触残基筛选、CCR构建在训练折内完成。</li><li>沿用上一轮位置分组规则；双突变采用组合位置组，可能与其他组共享单个位置，不能视为严格的所有组成位置均未见验证。</li><li>WT与突变体分别对两条实验记录取均值；已复核输入，但本轮未传播实验误差，也未增加独立seed。</li><li>C295的位置分组成绩接近随机水平；C388/C871保留为后续验证候选。</li></ul><h2>验证指标与混淆矩阵</h2><div class="scroll">'''+table(val)+'</div>'
body+='<h2>逐构建结果</h2><p>表内类别均使用相对WT定义。prediction为类别1的判定；score为类别1分数，尚未校准为实验成功概率。</p><div class="scroll">'+table(pred)+'</div></html>'
(root/'WT_relative_classification_report.html').write_text(body)
(root/'README.md').write_text('''# 相对WT动态分类重分析

在此目录运行 `python classify.py`，完成后运行 `python report.py`。
依赖：numpy pandas scipy scikit-learn joblib threadpoolctl；来源版本见 versions.json。

Δ编辑效率（百分点）=突变体平均编辑率−该位点WT平均编辑率。
类别1：Δ≥动态分界；类别0：Δ<动态分界。阈值扫描为−100至100 pp的5 pp步长加全部相邻观测值中点，等价标签去重，每类至少5个构建。
模型：Total CP / Interface CP / S4-S12-S24 / Repository特征逻辑回归，ContactSeek Top3 / CCR浅层随机森林。
inputs保留已提取CP、实验记录和结构特征；sources保留本轮调用的ContactSeek函数与许可证；来源记录见source_provenance.json。源结构压缩包未包含；可从已提取输入复现本轮全部分类。
audit_reference保留上一轮折外分数，用于逐分组核验平移不变性。

编辑分界与模型经全体数据探索筛选，最终成绩未消除该层选择偏差。模型分数阈值和CCR在训练折内确定。组合突变的位置分组并非严格的所有组成位置隔离。
主要文件：results/wt_relative_summary.csv，results/wt_relative_per_construct_predictions.csv，results/validation_checks.json，WT_relative_classification_report.html。
''')
print(summary.to_string(index=False))
print(json.dumps(audit,ensure_ascii=False,indent=2))
