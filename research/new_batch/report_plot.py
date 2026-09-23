import os,json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O='new_batch';a=pd.read_csv(O+'/ranked_results.csv').sort_values('design_id');s=pd.read_csv(O+'/seed_scores.csv');old=pd.read_csv('data/architecture_inputs/batch1_features.csv');audit=json.load(open(O+'/architecture_audit.json'));mp=pd.read_csv(O+'/repeat_mapping.csv')
a['design_id']=a.design_id.astype(int);a['domain']=a.design_id.map(lambda i:'15-repeat extrapolation' if i in [9,10,11] else ('new module order' if i>=12 else 'related 12-repeat'))
with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42}):
 fig,axes=plt.subplots(1,3,figsize=(14,7.5),layout='constrained',sharey=True)
 colors={'related 12-repeat':'#0072B2','15-repeat extrapolation':'#884DA7','new module order':'#A65A00'}
 for j,(_,r) in enumerate(a.iterrows()):
  v=s[s.construct==r.construct];c=colors[r.domain];axes[0].scatter(v.XGB_score,np.full(len(v),j),color=c,s=16,alpha=.4);axes[0].scatter(r.XGB_score,j,color=c,marker='D',s=50)
  f='pp_nonlocal4_high_per_res';axes[1].errorbar(r[f],j,xerr=v[f].std(),color=c,fmt='o',capsize=3)
  for typ,mark,col in [('XGB','D','#0072B2'),('RF','s','#884DA7'),('LR','x','#A65A00')]:axes[2].scatter(r[typ+'_score'],j,color=col,marker=mark,s=35,label=typ if j==0 else None)
 axes[0].set(xlim=(0,1),xlabel='Frozen XGBoost score',title='A  Construct score and individual seeds');axes[0].axvline(.5,ls='--',color='grey',lw=1)
 f='pp_nonlocal4_high_per_res';axes[1].axvspan(old[f].min(),old[f].max(),color='grey',alpha=.15,label='Old construct-mean range');axes[1].axvline(2.379752851,ls='--',color='grey',lw=1);axes[1].set(xlabel='High-CP pairs / residue (gap ≥ 4)',title='B  Contact density: mean ± seed SD');axes[1].legend(loc='lower right',fontsize=8)
 axes[2].set(xlim=(0,1),xlabel='Frozen model score',title='C  Model agreement');axes[2].axvline(.5,ls='--',color='grey',lw=1);axes[2].legend(loc='lower right')
 axes[0].set_yticks(range(len(a)),['ID %02d | %d repeats'%(r.design_id,r.repeat_count) for _,r in a.iterrows()]);axes[0].invert_yaxis()
 for ax in axes:ax.grid(axis='x',alpha=.12);ax.axhline(6.5,color='#aaaaaa',lw=.7);ax.axhline(9.5,color='#aaaaaa',lw=.7)
 fig.suptitle('New PUF designs: structural scores without experimental outcome labels',fontsize=16);fig.savefig(O+'/new_design_comparison.png',dpi=220,facecolor='white');fig.savefig(O+'/new_design_comparison.pdf',facecolor='white');plt.close(fig)

def md(x):
 def fmt(v):return f'{v:.4f}' if isinstance(v,(float,np.floating)) else str(v)
 return '\n'.join(['| '+' | '.join(x.columns)+' |','| '+' | '.join(['---']*len(x.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in x.itertuples(index=False,name=None)])
local=pd.read_csv(O+'/local_contact_means.csv');local=local.merge(a[['construct','design_id','domain']],on='construct');a.to_csv(O+'/final_summary.csv',index=False)
report='''# 2026-09-07 PUF新设计：冻结模型评分与结构核对

## 结论

这批包含15种设计、60个seed任务、300个结构模型；每种设计4个seed，每seed5个模型，编号缺6。训练标签未用于新增设计，未将新设计默认标为失败。所有预测RNA均为AUGGAGGACGUGC（13nt）。没有实验成功/失败标签，因此无法计算本批AUC或准确率。

在与旧训练骨架较接近的12-repeat设计中，7、8号可作为优先关注候选。7号构建均值XGBoost评分0.756，但仅3/4个seed过0.5；8号均值评分0.619，4/4个seed过0.5。1、4号RF和LR较高，XGBoost偏低，应保留为模型分歧候选。9–11号三种模型均高，但实际为15-repeat，属于明确外推，不能把其数值高直接解释为更可靠的实验成功预测。13、16号多模型均偏低，且平均pLDDT约87、CP加权PAE约2.26，弱于本批大部分设计。

## 完整结果

评分未校准为成功概率，阈值0.5沿用冻结模型。构建分数是先平均模型/seed的特征、再输入模型；seed通过数是分别对每个seed均值评分。树模型非线性，两者不必给出相同多数结论。

'''+md(a[['design_id','repeat_count','pp_nonlocal12_high_per_res','XGB_score','RF_score','LR_score','XGB_seed_positive_count','iptm','plddt_protein_mean']])+'''

## 重点设计

- **7号**：R123 R567 R567 Rloop67loop8，S12=1.32370；RF=0.737、LR=0.785。seed3的S4降到2.37572，XGBoost为0.267，其余三个seed为0.756。
- **8号**：R123 R567 R567 R6loop7loop8，S12=1.31647；RF=0.766、LR=0.736。四个seed的XGBoost为0.619–0.733。其seed分类更一致，但不代表已证明优于7号。
- **1号**：去掉loop的12-repeat设计，S12=1.33267，RF=0.780、LR=0.773；构建均值XGBoost=0.381。S4均值2.37871略低于冻结阈值2.37975，虽然三个seed为高分，seed1的低S4拉低均值。应视为边界与模型分歧案例。
- **4号**：S12=1.32324，RF=0.772、LR=0.764，XGBoost=0.381且四个seed都未过0.5。高S12不足以决定这个冻结XGBoost的分类，它还高度依赖S4。
- **9–11号**：分别620、620、627 aa，各15个repeat。S12=1.36976、1.35363、1.34370；全部四个seed为高分。但6个蛋白内部汇总特征都超过旧构建均值范围，旧模型也未训练过15-repeat。接触数除以长度无法保证消除repeat数与几何结构变化的影响。
- **12–16号**：模块顺序发生变化，旧数据没有覆盖新的前6个模块排列。12号模型分歧较大；14、15号偏低；13、16号三模型均低。应按新骨架单独积累验证数据。

## 输入序列核查

全部15种设计均未与旧14种构建的蛋白+RNA序列完全一致，旧实验标签没有直接继承。名称接近也不能等同于同一个构建。

**10号**名称含no8，但实际TRM计数为15，末端545–580位的核心序列仍高度匹配R8模块（排除TRM后的匹配率约97%）。请核对no8的原本含义及输入序列。

**16号**名称R12365673156loop8暗示来源顺序1,2,3,6,5,6,7,3,1,5,6,8；实际非TRM核心匹配支持1,2,3,6,5,6,6,5,1,3,2,8。主要不一致位于：

| 物理repeat | 文件名暗示来源 | 序列匹配来源 | 核心原始位置 |
|---|---|---|---|
| P7 | R7 | R6 | 250–285 |
| P8 | R3 | R5 | 286–321 |
| P10 | R5 | R3 | 358–393 |
| P11 | R6 | R2 | 394–436（包含loop） |

模块来源根据旧参考模块、排除TRM位点后最近序列匹配推断，存在边界氨基酸差异；完整匹配率及核心序列见repeat_mapping.csv。实际分数始终基于ZIP中的真实序列与结构，未按文件名改序列。

## 局部接触检查

沿用旧参考的K204–E239、L196–A235接触及CCR01（178/179/180）、CCR02（126/127/128）位置。仅1–11号（本批实际含1–5、7–11）的前六模块顺序与旧训练一致，适合直接对照；12–16号对应位置可能换成其他来源模块，不宜把原残基解释照搬过去。

'''+md(local[local.design_id<12][['design_id','CP_204_239','CP_196_235','CCR01_density','CCR02_density']].sort_values('design_id'))+'''

这些候选位点在上一轮没有通过多重检验校正，且加入CCR未稳定提升预测，因此局部对照仅帮助观察结构差异，不作为额外成功判据。本批未进行成功/失败组显著性检验。

## 方法与限制

直接读取contact_probs及PAE，Cα坐标用于回转半径，CIF原子置信度先平均到残基后求蛋白均值。构建均值按每seed五模型平均、再四seed等权平均。原14个CP summary特征顺序、标准化、模型参数、评分阈值均保持不变。

主要结果使用上一轮冻结XGBoost；LR及四项连续结构LR同样使用冻结文件。旧RF序列化文件无法完整读取，使用原14种训练构建、原特征及原固定参数（300树、深度2、叶节点最少2、balanced、random_state=2026）重建，仅用于模型一致性检查，没有使用新设计标签或重新调参。

S_d定义为完整蛋白中i<j、j−i≥d、CP>0.5的接触对数除以蛋白长度。冻结XGBoost主要依赖S4/S12/S24，因此评分是有限个离散值；不能用分数相同的候选推断严格生物学优劣。模型输入超出训练单变量最小/最大范围只是一种粗略外推标记，也不能据此宣布完全处于训练分布内。

本批ipTM约0.915–0.926，整体接近。高ipTM支持预测模型对蛋白–RNA相对排布有较高置信度，无法直接证明克隆、表达或编辑功能成功。旧成功标签本身的实验环节仍需要明确。

## 文件

- final_summary.csv：完整构建评分、seed稳定性及训练范围检查。
- seed_scores.csv / model_features.csv：60个seed和300个模型的完整特征。
- new_design_comparison.png/.pdf：构建与seed评分、S4密度、模型分歧。大菱形为构建均值特征评分，小点为seed评分，误差线为seed间SD。
- repeat_mapping.csv / architecture_audit.json：真实模块排列与来源匹配。
- local_contact_means.csv / local_contact_seed_values.csv：局部CP对照。
- sequences_and_audit.json：每种设计的真实序列、loop位置与旧序列匹配。

构建名映射：

'''+md(a[['design_id','construct','domain']])
open(O+'/new_PUF_analysis_report.md','w').write(report)
json.dump({'constructs':15,'seeds':60,'models':300,'experimental_labels_available':False,'missing_design_ids':[6],'validation':'frozen prediction only; no new-cohort AUC','figure_uncertainty':'SD of seed mean features, technical prediction variability','RF_reconstructed_from_old_training_only':True,'ranking':'XGBoost score ties are not resolved by biological evidence','alt_text':'Designs 7 and 8 have relatively high scores among related 12-repeat constructs. Designs 1 and 4 show model disagreement. Designs 9–11 are high-scoring but outside repeat-count training range. Designs 12–16 rearrange modules.'},open(O+'/provenance.json','w'),indent=2)
