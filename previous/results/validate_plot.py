import os,json,itertools
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import balanced_accuracy_score,roc_auc_score,confusion_matrix
O=os.path.dirname(__file__);a=pd.read_csv(O+'/construct_features.csv');s=pd.read_csv(O+'/seed_features.csv');p=pd.read_csv(O+'/leave_construct_out_predictions.csv');y=a.success.values
cols=[c for c in a if c.startswith(('pr_','pp_','rna','protein_bin'))];X=a[cols].values
# A depth-one tree repeats feature and threshold selection using training constructs only.
def cv(labels,save=False):
 out=[]; records=[]
 for i in range(len(y)):
  train=np.arange(len(y))!=i;m=DecisionTreeClassifier(max_depth=1,class_weight='balanced',random_state=0).fit(X[train],labels[train]);out.append(m.predict(X[i:i+1])[0]);records.append({'construct':a.construct.iloc[i],'selected_feature':cols[m.tree_.feature[0]],'training_threshold':m.tree_.threshold[0],'actual':int(labels[i]),'predicted':int(out[-1])})
 if save:pd.DataFrame(records).to_csv(O+'/nested_stump_predictions.csv',index=False)
 return balanced_accuracy_score(labels,out),confusion_matrix(labels,out).ravel().tolist()
score,cm=cv(y,True);null=[]
for comb in itertools.combinations(range(len(y)),3):
 z=np.zeros(len(y),int);z[list(comb)]=1;null.append(cv(z)[0])
m=json.load(open(O+'/metrics.json'));m['nested_feature_threshold_selection_stump']={'LOCO_balanced_accuracy':score,'TN_FP_FN_TP':cm,'exact_permutation_p_balanced_accuracy':float(np.mean(np.array(null)>=score-1e-10))};json.dump(m,open(O+'/metrics.json','w'),indent=2)
f='pp_nonlocal12_high_per_res';order=a.sort_values(f,ascending=False).construct.tolist();mapping={c:'C%02d'%(i+1) for i,c in enumerate(order)}
a['ID']=a.construct.map(mapping);a[['ID','construct','success','protein_length',f,'pr_cp_sum']].sort_values('ID').to_csv(O+'/construct_summary.csv',index=False)
with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42}):
 fig,axes=plt.subplots(1,3,figsize=(15,6),layout='constrained',gridspec_kw={'width_ratios':[1.2,1,1]})
 for i,c in enumerate(order):
  row=a[a.construct==c].iloc[0];vals=s[s.construct==c][f].values;col='#0072B2' if row.success else '#B45B00';mark='D' if row.success else 'o'
  axes[0].scatter(vals,np.full(len(vals),i),s=15,color=col,alpha=.45);axes[0].errorbar(row[f],i,xerr=np.std(vals,ddof=1),fmt=mark,color=col,capsize=3,markersize=6)
  axes[1].scatter(row.pr_cp_sum,i,color=col,marker=mark,s=40)
  prob=p.loc[p.construct==c,'CP_summary_prob'].iloc[0];axes[2].scatter(prob,i,color=col,marker=mark,s=40)
 for ax in axes:ax.set_yticks(range(14));ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
 axes[0].set_yticklabels([mapping[c]+('  success' if a.loc[a.construct==c,'success'].iloc[0] else '  failure') for c in order]);axes[1].set_yticklabels([]);axes[2].set_yticklabels([])
 axes[0].set(xlabel='High-CP nonlocal pairs / protein residue',title='A  Protein internal contacts')
 axes[1].set(xlabel='Sum of protein–RNA CP',title='B  Protein–RNA contacts')
 axes[2].set(xlabel='Held-out model score (uncalibrated)',title='C  Held-out logistic regression',xlim=(0,1));axes[2].axvline(.5,color='#777777',ls='--',lw=1)
 fig.suptitle('PUF construction outcomes: 3 successes and 11 failures',fontsize=16)
 fig.savefig(O+'/CP_comparison.png',dpi=220,facecolor='white');fig.savefig(O+'/CP_comparison.pdf',facecolor='white');plt.close(fig)
 fig,ax=plt.subplots(figsize=(10,6),layout='constrained');mat=np.stack([a.loc[a.construct==c,['rna%02d_cp_sum'%j for j in range(1,14)]].values[0] for c in order]);im=ax.imshow(mat,aspect='auto',cmap='cividis',vmin=0,vmax=mat.max(),interpolation='nearest');ax.set_yticks(range(14),[mapping[c]+(' success' if a.loc[a.construct==c,'success'].iloc[0] else ' failure') for c in order]);ax.set_xticks(range(13),['%d %s'%(i+1,b) for i,b in enumerate('AUGGAGGACGUGC')]);ax.set(xlabel='RNA position and nucleotide',title='Protein–RNA CP summed over protein residues');fig.colorbar(im,ax=ax,label='Sum of CP');fig.savefig(O+'/RNA_CP_heatmap.png',dpi=220,facecolor='white');fig.savefig(O+'/RNA_CP_heatmap.pdf',facecolor='white');plt.close(fig)
print(json.dumps(m,indent=2));print(a[['ID','construct','success',f]].sort_values('ID').to_string(index=False));print(p[['construct','success','CP_summary_prob']].to_string(index=False))
