import pathlib,json
import numpy as np,pandas as pd,joblib
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import roc_auc_score,average_precision_score
O=pathlib.Path('c388_analysis');a=pd.read_csv(O/'training_data.csv');y=a.work.to_numpy();rng=np.random.default_rng(20260908)
features=['pp_nonlocal4_high_per_res','pp_nonlocal12_high_per_res','pp_nonlocal24_high_per_res','plddt_core_mean','plddt_core_min','global_core_PAE','pae_contact_weighted','Rg_full_length_normalized','anisotropy','pr_cp_per_nt','pr_rna_max_mean','pr_rna_coverage_05','pr_protein_coverage_05','adj_CP_mean','adj_CP_min','nonadj_CP_mean'];X=a[features].to_numpy();delta=X[y==1].mean(0)-X[y==0].mean(0);perms=np.array([rng.permutation(y) for _ in range(19999)]);stats=perms@X/sum(y)-(1-perms)@X/sum(1-y);p=((np.abs(stats)>=np.abs(delta)).sum(0)+1)/20000;order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);rows=[]
for i,f in enumerate(features):
 x=X[:,i];n1=sum(y);n0=sum(1-y);sd=np.sqrt(((n1-1)*x[y==1].var(ddof=1)+(n0-1)*x[y==0].var(ddof=1))/(n1+n0-2));rows.append(dict(feature=f,work_mean=x[y==1].mean(),nonwork_mean=x[y==0].mean(),delta=delta[i],hedges_g=(1-3/(4*(n1+n0)-9))*delta[i]/sd,point_biserial_r=pearsonr(x,y).statistic,spearman_C388=spearmanr(x,a.C388_mean).statistic,p_permutation_MC=p[i],q_BH=q[i],n_permutations=19999))
pd.DataFrame(rows).to_csv(O/'univariate_statistics.csv',index=False)
# Within-mutation-position directions are descriptive; only mixed groups are evaluable.
rows=[]
for g,b in a.groupby('mutation_group'):
 for f in features:rows.append(dict(group=g,feature=f,n_work=int(b.work.sum()),n_nonwork=int((1-b.work).sum()),delta=b.loc[b.work==1,f].mean()-b.loc[b.work==0,f].mean()))
pd.DataFrame(rows).to_csv(O/'within_position_directions.csv',index=False)
print(pd.read_csv(O/'univariate_statistics.csv').sort_values('p_permutation_MC').to_string(index=False))
