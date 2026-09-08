import ast,pathlib,numpy as np,pandas as pd
# Reuse only definitions/imports without launching the retraining sweep.
s=pathlib.Path('c388_threshold/joint.py').read_text();exec(s.split('u=np.unique(rate);thresholds=')[0])
p=pd.read_csv('c388_analysis/predictions.csv');out=[]
for name,b in p[(p.family=='frozen')&(p.scope=='matched512')&(p.validation=='transfer')].groupby('model'):
 r=b.C388_mean.to_numpy();u=np.unique(r)
 for c in np.unique(np.r_[np.arange(.05,.951,.05),(u[:-1]+u[1:])/2]):
  y=(r>=c).astype(int)
  if min(y.sum(),len(y)-y.sum())<2:continue
  ss=sweep(y,b.score.to_numpy());ss['distance']=(ss.score_threshold-.5).abs()
  for obj in ['accuracy','balanced_accuracy']:
   z=ss.sort_values([obj,'distance'],ascending=[False,True]).iloc[0].to_dict();out.append(dict(model=name,C388_threshold=c,n_work=int(y.sum()),n_nonwork=int(len(y)-y.sum()),objective=obj,AUC=roc_auc_score(y,b.score),**z))
pd.DataFrame(out).to_csv(O/'frozen_joint_results.csv',index=False)
