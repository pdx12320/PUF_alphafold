import pandas as pd,numpy as np,json,os,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O='combined12';a=pd.read_csv(O+'/training_data.csv');m=pd.read_csv(O+'/model_metrics.csv');p=pd.read_csv(O+'/heldout_predictions.csv');cv=m[m.validation=='combined_LOCO'];new=p[(p.validation=='combined_LOCO')&(p.batch=='new')];table=new[new.model=='CP_density3_RF'][['design_id','success','score','prediction']].rename(columns={'score':'density_RF_score','prediction':'density_RF_prediction'}).merge(new[new.model=='CP_structure17_RF'][['design_id','score','prediction']].rename(columns={'score':'structure_RF_score','prediction':'structure_RF_prediction'}),on='design_id').sort_values('design_id');table.to_csv(O+'/new_batch_LOCO_comparison.csv',index=False)
with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
 fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained');types=['LR','RF','XGB'];sets=['CP_density3','CP_summary14','CP_structure17'];v=np.array([[cv[cv.model==k+'_'+t].AUC.iloc[0] for t in types] for k in sets]);im=axes[0].imshow(v,vmin=0,vmax=1,cmap='cividis');axes[0].set_xticks(range(3),types);axes[0].set_yticks(range(3),['3 CP densities','14 CP summaries','+ 3 structural features']);axes[0].set_title('Combined LOCO AUC: 24 constructs')
 for i in range(3):
  for j in range(3):axes[0].text(j,i,f'{v[i,j]:.3f}',ha='center',va='center',fontsize=14)
 fig.colorbar(im,ax=axes[0],fraction=.045)
 for i,(_,r) in enumerate(table.iterrows()):
  axes[1].scatter(r.density_RF_score,i,color='#0072B2',marker='o',label='3-density RF' if i==0 else None);axes[1].scatter(r.structure_RF_score,i,color='#B45B00',marker='s',label='17-feature RF' if i==0 else None)
 axes[1].set_yticks(range(len(table)),['ID %02d | %s'%(r.design_id,'work' if r.success else 'nonwork') for _,r in table.iterrows()]);axes[1].invert_yaxis();axes[1].set(xlim=(0,1),xlabel='Held-out score',title='New batch: only ID 01 works');axes[1].axvline(.5,color='grey',ls='--');axes[1].legend(loc='lower right');fig.savefig(O+'/combined12_results.png',dpi=220,facecolor='white');fig.savefig(O+'/combined12_results.pdf',facecolor='white');plt.close(fig)
def md(x):
 def f(v):return f'{v:.4f}' if isinstance(v,(float,np.floating)) else str(v)
 return '\n'.join(['| '+' | '.join(x.columns)+' |','| '+' | '.join(['---']*len(x.columns))+' |']+['| '+' | '.join(f(v) for v in row)+' |' for row in x.itertuples(index=False,name=None)])
report='''# 两批12-repeat PUF合并训练

## 数据冻结

将用户“12nt”按12-repeat构建设计处理。预测输入RNA实际仍为13nt AUGGAGGACGUGC，未裁剪RNA。上一批排除2种PUF11，保留12种（3 work / 9 nonwork）；这批排除9、10、11号15-repeat，保留12种，只有1号work（1 work / 11 nonwork）。共24个独立构建、4 work / 20 nonwork，来自96个seed任务、480个模型。构建作为独立样本；各seed内平均模型，再对seed等权平均。6号未出现在输入中。

这批所有真实序列均与上一批构建不完全相同。仍存在高度相近的骨架和loop变体，LOCO结果可能高于真正全新骨架的泛化表现。10、16号的名称与序列问题仍按上一报告保留，未擅自修改；10号由于15-repeat被排除，16号按实际结构及用户nonwork标签纳入。

## 方法

三组固定特征：S4/S12/S24三种全局高CP密度；原14个CP summary；14个CP summary加蛋白平均pLDDT、CP加权对称PAE、长度归一化Cα回转半径，共17个。未加入此前未验证的CCR。各模型参数沿用上一轮，未通过本轮LOCO调整阈值或超参数。LR C=0.1且balanced；RF 300树、深度2、叶节点最少2且balanced；XGBoost 100轮、深度1、学习率0.05、L1=1、L2=10、列采样0.8、按训练类别比例加权。额外提供深度1树作为旧简单规则基线。

24折LOCO：每折完整留出1种构建，标准化与训练仅使用另外23种；同构建的seed不跨训练与测试。阈值统一0.5，未用留出标签选择阈值。最终模型另用全部24种重新拟合并保存，最终拟合模型分数不得作为LOCO性能。

## 合并LOCO结果

'''+md(cv[['model','AUC','AP','accuracy','balanced_accuracy','precision','recall','F1','TP','FP','FN','TN']])+'''

推荐并列保留两个模型：

- **CP_density3_RF**：AUC=0.9125，21/24正确，work识别3/4，误报2个nonwork；平衡准确率0.825。作为当前默认分类基线。
- **CP_structure17_RF**：AUC=0.9250、AP=0.7679，20/24正确，work识别2/4、误报2个nonwork。排序更好，但固定阈值下漏判更多work，不能笼统称为所有方面更优。

全预测nonwork的准确率即20/24=83.33%，因此准确率必须与召回率、平衡准确率同时阅读。多模型比较中的最优结果有选择偏差，目前只有4个work，尚未获得下一批独立验证。

## 这批逐构建留出结果

'''+md(table)+'''

两个RF均在不使用1号本身标签的留出折里识别出1号；均将4、7号误判为work。三密度RF将7号排在1号前，17特征RF将1号排在所有24个构建留出分数的第一位。8号在两模型中都低于0.5。

17特征RF的另两个漏判为上一批成功构建R123_R567_R5loop67_loopR678与R123_R567_R6R7_R5678。因此，不能仅看新1号排第一就宣称问题已经解决。

## 旧批训练、新批测试

另仅使用上一批12种构建训练，再对这批12种已标注构建测试。这与合并LOCO不同，不让这批标签进入模型拟合。

'''+md(m[m.validation=='previous12_to_new12'][['model','AUC','accuracy','balanced_accuracy','precision','recall','TP','FP','FN','TN']])+'''

旧批12-repeat summary XGBoost漏判1号、误判7和8号，TP=0、FP=2、FN=1、TN=9，AUC=0.7727。上一轮实际发出的旧14构建冻结XGBoost在这批12-repeat上的结论也是上述混淆矩阵。这说明此前高CP候选推荐没有对应到本批真实work标签。

部分旧批RF/LR在这批只有1个阳性的测试中AUC达到1，但固定0.5阈值仍有多次误报。AUC=1仅说明唯一阳性的分数排在前面，不能解释为完全分类正确。

## 使用已训练模型

默认使用三密度RF：

```bash
python combined12/predict.py --zip new_folds.zip --out predictions.csv
```

ZIP须包含本批同格式AlphaFold-server full_data JSON与job_request JSON，目录名为*_rna_*_seedN。脚本只纳入识别出12个TRM的构建，避免将11/15-repeat混入。TRM识别针对当前骨架，若换用其他识别三肽需提供明确repeat映射。

其他模型可从特征CSV评分：

```bash
python combined12/predict.py --features features.csv --model CP_structure17_RF --out predictions.csv
```

需要的特征列、标签规则见manifest.json。所有评分均未校准为实验成功概率。模型文件保存后均重新加载并核对预测一致性。LR_formulas.csv保留各线性模型的均值、标准差、系数和截距；feature_importance_summary.csv提供跨LOCO训练折的重要性，属于描述性解释。

## 构建标签清单

'''+md(a[['construct','batch','design_id','success']])+'''

success=1表示work，0表示nonwork。新标签仅适用于这批纳入的12-repeat构建；本轮未训练或判断被排除的15-repeat构建。
'''
open(O+'/combined12_report.md','w').write(report)
