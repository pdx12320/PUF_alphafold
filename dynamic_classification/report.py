from pathlib import Path
import pandas as pd,json,html,shutil,zipfile,platform
p=Path('.');o=p/'results';s=pd.read_csv(o/'exploratory_threshold_scan.csv');v=pd.read_csv(o/'training_only_score_threshold_validation.csv');labels=pd.read_csv(o/'candidate_labels.csv');pr=pd.read_csv(o/'validated_predictions.csv')
best=s.sort_values(['BA','AUC','model','label_threshold_pct'],ascending=[False,False,True,True]).groupby('site').head(1).sort_values('site');best.to_csv(o/'exploratory_candidates.csv',index=False)
# All values are held-out predictions; candidate selection remains exploratory.
def tab(x):return x.to_html(index=False,float_format=lambda z:f'{z:.4f}',border=0)
summary=best[['site','label_threshold_pct','equivalent_lower_exclusive','equivalent_upper_inclusive','model','n_high','n_low','BA','accuracy','AUC']]
vv=v[['site','validation','model','label_threshold_pct','BA','accuracy','AUC','sensitivity','specificity','TN','FP','FN','TP','median_score_threshold','min_score_threshold','max_score_threshold','majority_accuracy']]
import joblib,numpy as np
bundles=joblib.load(o/'candidate_models.joblib');importance=[]
for site,b in bundles.items():
 f=b['feature'];md=b['model'];cols=b['repository_columns']
 if f in ['total','interface','density','repo']:
  cols= ['pr_cp_sum'] if f=='total' else [c for c in cols if c.startswith('pr_')] if f=='interface' else [f'pp_nonlocal{k}_high_per_res' for k in [4,12,24]] if f=='density' else cols
  for col,w in zip(cols,md[-1].coef_[0]):importance.append(dict(site=site,feature=col,weight=float(w),meaning='standardized_logistic_coefficient_for_high_activity'))
 else:
  im=md.feature_importances_
  regions=b['regions'] if f=='ccr' else [dict(positions=[i]) for i in np.where(b['contact_mask'])[0]]
  for j,r in enumerate(regions):importance.append(dict(site=site,feature='residues_'+','.join(str(int(i)+8) for i in r['positions']),weight=float(im[j*12:(j+1)*12].sum()),meaning='RF_importance_exploratory'))
imp=pd.DataFrame(importance);imp.to_csv(o/'candidate_feature_weights_exploratory.csv',index=False)
blocks=['<h2>候选模型使用的 CP 特征</h2><p>全数据拟合的探索性权重。逻辑回归正系数指向高活性，负系数指向低活性；相关特征的系数不能单独解释为机制。随机森林重要度仅表示模型使用程度。残基编号已换回输入蛋白 519 aa 编号。</p>'+tab(imp.assign(abs_weight=imp.weight.abs()).sort_values('abs_weight',ascending=False).groupby('site').head(8).drop(columns='abs_weight'))]
for site in ['C388','C871','C295']:
 z=labels[labels.site==site].copy();lo=pr[(pr.site==site)&(pr.validation=='LOOCV')][['construct','score','score_threshold','prediction']];z=z.merge(lo,on='construct').sort_values('editing_pct',ascending=False);blocks.append(f'<h2>{site}：实验标签与留一预测</h2>'+tab(z))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
fig,axes=plt.subplots(1,3,figsize=(11,3.8),layout='constrained')
for ax,site in zip(axes,['C388','C871','C295']):
 row=v[(v.site==site)&(v.validation=='LOOCV')].iloc[0]
 cm=np.array([[row.TN,row.FP],[row.FN,row.TP]],dtype=int)
 ax.imshow(cm,cmap='Blues',vmin=0,vmax=30)
 for i in range(2):
  for j in range(2):ax.text(j,i,str(cm[i,j]),ha='center',va='center',fontsize=19,color='white' if cm[i,j]>18 else '#172b3a')
 ax.set_xticks([0,1],['Low','High']);ax.set_yticks([0,1],['Low','High']);ax.set_xlabel('Predicted class');ax.set_ylabel('Observed class');ax.set_title(f'{site} | cutoff {row.label_threshold_pct:g}%\nBalanced accuracy {row.BA:.1%}',fontsize=11)
fig.suptitle('Training-only score thresholds: leave-one-construct-out\nExploratory candidate selection; n = 30 mutants per site',fontsize=12)
fig.savefig(o/'confusion_matrices.png',dpi=220);fig.savefig(o/'confusion_matrices.pdf');plt.close(fig)
figure='<img style="max-width:100%" alt="Confusion matrices" src="data:image/png;base64,'+base64.b64encode((o/'confusion_matrices.png').read_bytes()).decode()+'">'
text='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>PUF 动态阈值分类分析</title><style>body{font-family:system-ui,sans-serif;max-width:1400px;margin:36px auto;padding:0 24px;color:#172b3a;line-height:1.65}table{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto}th,td{padding:7px 9px;border-bottom:1px solid #d8e0e6;text-align:right}th{background:#e9f1f7}h1,h2{color:#133d59}.note{background:#fff4dc;padding:16px;border-left:4px solid #ac7300}</style><h1>PUF：三个 C 位点的动态阈值分类</h1>
<p>取消 C388 &gt;50% 入组条件，使用全部 30 个有实验数据的突变体；WT 仅提供结构与编辑率参考，不进入训练。每个构建两条实验记录取均值。C388、C871、C295 分别训练，标签统一为编辑率 ≥ 位点阈值的“高活性”（1），低于阈值的“低活性”（0）。C388 高活性有利于保留目标编辑；C871/C295 低活性可作为降低相应编辑的候选。三个任务均使用全体突变体，未按 C388 再次筛选。</p>
<p class="note">这是阈值与模型的探索性筛选。分数阈值的补充验证只使用训练折内预测，但候选编辑率阈值和模型已经依据本批数据选出，因此仍存在选择偏差。表中的最佳结果不能视为独立测试性能，也不能证明新的分界是精确生物学边界。换分类任务后，指标不能与此前回归 R² 直接比较。</p>
<h2>候选编辑率分界</h2><p>BA 为平衡准确率：高活性召回率与低活性召回率的平均，恒猜一个类别时为 0.5。下表 BA/accuracy 的模型分数阈值依据全部留出标签扫描，属于探索性上限。等价区间为（lower, upper]，此区间内任何编辑率分界产生相同分组。</p>'''+tab(summary)+'''<h2>训练折内选择模型分数阈值</h2><p>仅对每个位点探索性 BA 最佳组合复核。LOOCV 每次留出一个突变体；Position_group_out 按改变的 P 位置组合留出，训练内也按位置分组。接触残基筛选、CCR 定义、标准化全部仅使用相应训练折。分数阈值由内层 3 折预测最大化 BA，平局优先接近 0.5；外层测试标签不参与该阈值选择。下列模型分数未做概率校准。</p>'''+tab(vv)+figure+'''<h2>方法与数据边界</h2><ul>
<li>参考 PUF_alphafold 的 c388_threshold/thresholds.py 和 README：编辑率 5%–95% 每 5 个百分点及所有相邻观测均值中点，等价分组只训练一次。为避免极端小类，本次要求两类各至少 5 个突变体，优先保留整 5% 阈值。</li>
<li>比较 6 个固定模型：蛋白–RNA 总 CP LR；8 项接触面 CP 汇总 LR；蛋白内部 S4/S12/S24 LR；19 项 CP+结构特征 LR；ContactSeek 分区 Top3 RF；ContactSeek CCR RF。总 CP 指蛋白–RNA 子矩阵的求和，不代表整个复合体结合概率。所有特征仅来自对应 C 位点结构，无实验编辑率作为输入。</li>
<li>LR：C=0.1、balanced、liblinear；RF：100 棵树、深度 2、叶节点最少 3、balanced、全部特征参与候选；随机种子 2026。参数沿用参考代码，未另外搜索。</li>
<li>按参考 PUF 特征裁去输入蛋白前 7 aa，使用 512 aa × 15 nt。ContactSeek 对 RNA 1–5、6–10、11–15 分区，取每残基 Top3 CP 及 max(mutant)−max(WT)。先逐 seed 做非线性提取，再等权平均；同任务 5 份 CP 完全相同，未充当独立样本。</li>
<li>接触残基阈值 CP=0.15、差值=0.1；CCR 相关阈值 0.6、带宽 7。分类的动态阈值与这些结构筛选阈值含义不同。</li>
<li>共 93 个构建×位点组合（含 WT），97 个 seed 任务；表外 P4-R5-VTF 无对应实验标签，未训练。新补 P6-R7-MHQ C388 已纳入。</li>
<li>原表部分零值来自未检出 C→T 记录，本次按原表保留，不能解释为已确认精确零活性。重复类型未知；当前分析未估计标签测量误差，也未做排除零记录的敏感性验证。</li>
<li>WT 各位点仅一个 seed；AYE_C871 跨 seed 差异较大。分类分界附近的构建、分数阈值跨折变化及跨位置结果需要新批次实验验证。</li>
<li>训练内分数阈值验证未联合嵌套选择编辑率标签阈值和模型；未进行独立外部验证。提供全量扫描表，避免仅展示高分组合。</li></ul>
<h2>来源</h2><p><a href="https://github.com/pdx12320/PUF_alphafold/blob/553299d10213ec31b07c78874dbadb58ef3550e2/c388_threshold/README.md">PUF_alphafold 联合阈值方法</a>；<a href="https://github.com/menghaowei/ContactSeek/blob/7cadb91f8d629e62dabf3bdf939ce9688521eea5/ContactSeek/ContactExtraction.py">ContactSeek 接触特征</a>。来源快照及本次代码随结果包保存。</p>'''+''.join(blocks)+'</html>'
(p/'Dynamic_threshold_classification_report.html').write_text(text)
(p/'README.md').write_text('运行环境：Python 3，numpy pandas scipy scikit-learn joblib threadpoolctl。\n在本目录运行 python classify.py，然后 python report.py。inputs 为已提取的可复现输入，extract.py 可从同级 upload 原压缩包重新提取（还依赖前次 work_cp_modeling 的原始汇总表）。\n结果为探索性。训练内阈值验证仍受事后选择编辑率标签分界和模型的偏差影响。不得作为独立验证的部署规则。\n标签：editing_pct >= threshold 为1（高活性），低于为0。\n')
import sklearn,numpy,scipy
(p/'versions.json').write_text(json.dumps(dict(python=platform.python_version(),sklearn=sklearn.__version__,numpy=numpy.__version__,scipy=scipy.__version__),indent=2))
with zipfile.ZipFile('../Dynamic_threshold_classification_results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in p.rglob('*'):
  if f.is_file() and '__pycache__' not in str(f) and f.name!='run.log':z.write(f,'dynamic_classification/'+str(f))
print(summary.to_string(index=False));print(vv.to_string(index=False))
