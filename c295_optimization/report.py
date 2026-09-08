from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix,roc_auc_score,average_precision_score
P=Path(__file__).resolve().parent;O=P/'results'
pr=pd.concat([pd.read_csv(O/'outer_predictions.csv'),pd.read_csv(O/'original_CCR_predictions.csv')],ignore_index=True)
rows=[]
for (validation,family),g in pr.groupby(['validation','family']):
    weights=1/g.construct.map(g.construct.value_counts()).values
    tn,fp,fn,tp=confusion_matrix(g.y,g.prediction,labels=[0,1],sample_weight=weights).ravel()
    rows.append(dict(validation=validation,family=family,BA=.5*(tp/(tp+fn)+tn/(tn+fp)),accuracy=(tn+tp)/(tn+fp+fn+tp),precision=tp/(tp+fp) if tp+fp else 0,recall=tp/(tp+fn),TN=tn,FP=fp,FN=fn,TP=tp,AUC=roc_auc_score(g.y,g.score,sample_weight=weights),AP=average_precision_score(g.y,g.score,sample_weight=weights),n_constructs=g.construct.nunique(),n_predictions=len(g)))
met=pd.DataFrame(rows);met.to_csv(O/'all_model_performance.csv',index=False);pr.to_csv(O/'all_outer_predictions.csv',index=False)
audit=json.loads((O/'sequence_and_window_audit.json').read_text());pos={r['construct']:set(r['positions']) for r in audit}
for fold in json.loads((O/'split_audit.json').read_text()):
    assert not set(fold['train'])&set(fold['test'])
    if fold['validation']=='Strict_position_out':assert all(int(fold['fold']) not in pos[c] for c in fold['train'])
for (val,fam),g in pr.groupby(['validation','family']):assert g.construct.nunique()==30
sel=pd.read_csv(O/'fold_model_selection.csv');sel.groupby(['validation','reported_family','family','parameters']).size().reset_index(name='n_selected').to_csv(O/'selection_frequency.csv',index=False)
order=['original_CCR_RF','sequence','local_CP','sequence_plus_CP','nested_family_selection']
names=['Original CCR RF','Sequence','Local ΔCP','Sequence + ΔCP','Nested input selection']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(12,5.8));fig.subplots_adjust(left=.09,right=.985,bottom=.24,top=.77)
x=np.arange(5);w=.32
for off,val,color,hatch,label in [(-w/2,'LOOCV','#237D86',None,'Leave-one-construct-out'),(w/2,'Strict_position_out','#435C96','//','Strict position-out')]:
    vals=met[met.validation==val].set_index('family').loc[order,'BA'].values*100
    bars=ax.bar(x+off,vals,w,color=color,hatch=hatch,edgecolor='white',label=label)
    for b,value in zip(bars,vals):ax.text(b.get_x()+w/2,value+1.5,f'{value:.1f}',ha='center',fontsize=11,weight='bold')
ax.axhline(50,color='#727984',lw=1,ls='--');ax.set(ylim=(0,100),xticks=x,xticklabels=names,ylabel='Balanced accuracy (%)');ax.yaxis.grid(alpha=.13);ax.set_axisbelow(True)
fig.suptitle('C295 optimization: does the signal transfer to a new position?',x=.09,ha='left',y=.96,fontsize=19,weight='bold')
fig.text(.09,.88,'Positive = editing decreases by >20 percentage points vs WT | 10 positive / 20 negative mutants',fontsize=11)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,frameon=False)
fig.text(.09,.035,'Tuning and score thresholds selected inside training folds. Strict position-out excludes every construct containing the held-out P.',fontsize=9)
fig.text(.09,.005,'Double mutants contribute to two position tests, with total weight 1 per mutant. Point estimates; no independent validation.',fontsize=9)
fig.savefig(P/'C295_optimization_comparison.png',dpi=300,facecolor='white');fig.savefig(P/'C295_optimization_comparison.pdf',facecolor='white');plt.close(fig)
summary=met[met.family.isin(order)].pivot(index='family',columns='validation',values='BA').loc[order]
strict=met[met.validation=='Strict_position_out'].set_index('family').loc[order]
disp=pd.DataFrame({'模型/输入':names,'留一BA':summary.LOOCV.values,'严格位置留出BA':summary.Strict_position_out.values,'严格位置precision':strict.precision.values,'严格位置recall':strict.recall.values})
display=disp.copy()
for col in display.columns[1:]:display[col]=display[col].map(lambda x:f'{x:.1%}')
txt='''# C295 相对WT分类优化

正类固定为 ΔC295 < −20个百分点；负类为Δ≥−20。30个突变体中正类10个、负类20个。WT为32.958%，该分界等价于突变体编辑率<12.958%。本轮没有扫描新的实验标签分界。

## 结果

'''+'| '+' | '.join(display.columns)+' |\n|'+'---|'*len(display.columns)+'\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in display.itertuples(index=False,name=None))+'''

BA是正负两类召回率的平均值；precision/recall均针对“下降超过20个百分点”这一正类。严格位置验证中双突变分别进入两个位置的测试，按每个构建预测次数的倒数加权，避免双突变在总成绩中被加倍计数；因此部分混淆矩阵计数可以为小数。n_unique=30，不能把全部预测行当成独立样本。

## 结论与位置分解

局部ΔCP是本轮三个输入族中表现较好的候选，严格位置BA=67.5%，高于同验证下原CCR的43.8%。但其precision=43.5%，低于一半；recall=100%来自当前小样本测试，不能外推。

局部ΔCP在P4测试的9个构建和P5测试的13个构建上均全部判为正类，组内BA均50%；P6的6个构建全部为负类且全部判对；P7的7个构建BA=70%；P8仅1个负类且判错。P6和P8不具备两类样本，单位置BA不可估计。整体改善不能直接解释为各位置内部均具备良好区分能力。

输入族也通过内层选择时，严格位置BA降为56.25%，留一BA=42.5%。参数选择不稳定：局部CP的5个严格位置折中，3折选L2 C=0.01，1折选L2 C=0.1，1折选Elastic-net C=1、l1_ratio=0.1。当前未选定一个可部署的最终模型，也未将全数据训练分数用于报告泛化。

建议保留局部ΔCP作为待验证候选，下一轮优先确认结构seed稳定性并扩充各P位置内部的正负构建，暂缓扩大模型搜索。当前改进尚未获得独立验证。

## 输入改动

序列特征由实际519 aa蛋白与WT逐残基比较得到：P位置、名称可解析的来源R、突变数、单/双位置、氨基酸组成差值、TRM第12/13/16位的分位氨基酸变化、Kyte–Doolittle疏水性差值、简化KR/DE电荷编码和侧链重原子数差值。侧链重原子数仅是大小代理，未计算物理体积。个别名称来源R字段不规范时记录缺失计数，不猜测修正。

局部CP用实际突变残基的序列距离划分：±4、±18、距离19–54、距离≥55及全蛋白。窗口为预设序列范围，未将其声称为精确repeat结构边界。每区统计ΔCP总和、绝对和、RMS、增强/减弱接触和、WT CP≥0.1接触面的Δ总和、每核苷酸ΔCP及每残基Top-3差值。缓存512个蛋白残基对应原输入8–519，已校准编号。15 nt RNA的每个碱基独立保留；未假定其中某个C就是编辑中心。

## 模型与参数

每个输入族比较26套：L2逻辑回归C=0.01/0.1/1/10；Elastic-net同样C×l1_ratio=0.1/0.5/0.9；自动收缩LDA；RF深度1/2/3×叶节点最少3/5/8，100棵树。分类器采用balanced权重或均衡先验。去常数特征、标准化和所有拟合均仅用训练折。

留一验证内部使用3折分层；严格位置验证内部再次留出包含某P的全部构建。内部折外分数用于选择参数和分数阈值，外层留出只用于最终评价。Nested input selection还把三个输入族的选择放入内部验证。原CCR RF以同一正类定义、相同严格位置分组及训练内阈值复核；CCR特征也在训练折内构建。

## 解释边界

严格位置验证覆盖P4/P5/P6/P7/P8，部分位置测试样本很少或只有一个类别，单位置BA会缺失。该实验只排除正在测试的P；一个双突变的另一P可能出现在训练集中，未要求全部突变位置同时未知。不同外层折模型分数未统一校准。

20 pp标签来自前期本批数据探索；即使本轮固定，仍有历史选择偏差。新增输入方案也基于已有结果提出，所有成绩均待新构建验证。每个C295突变体只有1个seed，未完成跨seed稳健性检验；未增加湿实验，未执行全流程置换或置信区间估计。单独降低C295不等于保留C388或同时降低C871。

## 文件与复现

运行 `python run.py` 进行新特征提取和嵌套模型比较；运行 `python original_baseline.py` 复核原CCR；最后运行 `python report.py` 生成汇总与图。依赖见 requirements.txt。

results/all_model_performance.csv 为统一对照结果；all_outer_predictions.csv 包含所有折外预测；fold_model_selection.csv 和 selection_frequency.csv 记录参数选择；split_audit.json 给出每折训练/测试名单；sequence_and_window_audit.json 给出真实氨基酸突变和窗口范围。inputs包含可复现输入缓存和从原始job_request提取的序列，sources保留ContactSeek调用源码与许可证。
'''
(P/'README.md').write_text(txt)
body='<html lang="zh"><meta charset="utf-8"><title>C295优化结果</title><style>body{max-width:1200px;margin:40px auto;font:16px/1.7 system-ui;color:#243447}table{border-collapse:collapse;width:100%}td,th{padding:10px;border-bottom:1px solid #ddd}img{width:100%}pre{white-space:pre-wrap}</style><h1>C295相对WT分类优化</h1>'+display.to_html(index=False,border=0)+'<img src="C295_optimization_comparison.png"><pre>'+__import__('html').escape(txt)+'</pre></html>'
(P/'C295_optimization_report.html').write_text(body)
import sklearn,scipy,joblib,threadpoolctl
(P/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in [('numpy',np.__version__),('pandas',pd.__version__),('scipy',scipy.__version__),('scikit-learn',sklearn.__version__),('joblib',joblib.__version__),('threadpoolctl',threadpoolctl.__version__),('matplotlib',__import__('matplotlib').__version__)])+'\n')
(O/'validation_checks.json').write_text(json.dumps({'label_positive':'delta_C295_pp < -20','n_constructs':30,'n_positive':10,'outer_construct_disjoint':True,'heldout_position_absent_from_training':True,'all_families_cover_all_30_constructs':True,'strict_aggregate_weight':'1 / number of tests for each construct','uncertainty':'point estimates only; no confidence intervals','no_new_seed_data':True},indent=2))
print(display.to_string(index=False))
