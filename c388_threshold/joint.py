import pathlib,json,numpy as np,pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
O=pathlib.Path('c388_threshold');a=pd.read_csv('c388_analysis/training_data.csv');man=json.load(open('c388_analysis/manifest.json'));rate=a.C388_mean.to_numpy();n=len(a)
def model(t):
 m=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026) if t=='LR' else RandomForestClassifier(n_estimators=100,max_depth=2,min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
 return make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),m)
def sweep(y,p):
 u=np.unique(p);ts=np.unique(np.r_[0,(u[:-1]+u[1:])/2,1.00000001,.5]);z=p[:,None]>=ts;tp=(z&y[:,None].astype(bool)).sum(0);tn=(~z&~y[:,None].astype(bool)).sum(0);fn=y.sum()-tp;fp=len(y)-y.sum()-tn
 return pd.DataFrame(dict(score_threshold=ts,accuracy=(tp+tn)/len(y),balanced_accuracy=(tp/y.sum()+tn/(len(y)-y.sum()))/2,TP=tp,TN=tn,FP=fp,FN=fn))
u=np.unique(rate);thresholds=np.unique(np.r_[np.arange(.05,.951,.05),(u[:-1]+u[1:])/2]);seen=set();jobs=[]
for c in sorted(thresholds,key=lambda c:(abs(c*20-round(c*20))>1e-7,c)):
 y=(rate>=c).astype(int);sig=tuple(y)
 if min(y.sum(),n-y.sum())<2 or sig in seen:continue
 seen.add(sig);jobs.append((c,y))
rows=[];pred=[]
for c,y in sorted(jobs):
 for key,typ in [('CP3','LR'),('CP_Structure9','LR'),('CP3','RF'),('CP_Structure9','RF')]:
  X=a[man['feature_sets'][key]].to_numpy()
  for val in ['LOCO','mutation_group_out']:
   groups=a.construct.to_numpy() if val=='LOCO' else a.mutation_group.to_numpy();p=np.full(n,np.nan)
   for g in sorted(set(groups)):
    tr=groups!=g;te=~tr
    if len(np.unique(y[tr]))<2:continue
    m=model(typ).fit(X[tr],y[tr]);p[te]=m.predict_proba(X[te])[:,1]
   if np.isnan(p).any():continue
   s=sweep(y,p);s['distance']=(s.score_threshold-.5).abs()
   for objective in ['accuracy','balanced_accuracy']:
    sec='balanced_accuracy' if objective=='accuracy' else 'accuracy';b=s.sort_values([objective,sec,'distance'],ascending=[False,False,True]).iloc[0].to_dict();rows.append(dict(C388_threshold=c,equivalent_cutoff_lower_exclusive=rate[y==0].max(),equivalent_cutoff_upper_inclusive=rate[y==1].min(),n_work=int(y.sum()),n_nonwork=int(n-y.sum()),model=key+'_'+typ,validation=val,objective=objective,AUC=roc_auc_score(y,p),majority_accuracy=max(y.sum(),n-y.sum())/n,**b))
   for i in range(n):pred.append(dict(C388_threshold=c,model=key+'_'+typ,validation=val,construct=a.construct.iloc[i],rate=rate[i],y=y[i],score=p[i]))
 pd.DataFrame(rows).to_csv(O/'joint_threshold_results.csv',index=False);pd.DataFrame(pred).to_csv(O/'joint_predictions.csv',index=False);print('done',c,int(y.sum()),flush=True)
