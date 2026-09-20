import json,pathlib
import numpy as np,pandas as pd
from sklearn.metrics import confusion_matrix,roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold,GroupKFold
O=pathlib.Path('c388_threshold');a=pd.read_csv('c388_threshold/inputs/training_data.csv');man=json.load(open('c388_threshold/inputs/manifest.json'));y=a.work.to_numpy();names=a.construct.to_numpy();groups=a.mutation_group.to_numpy()
def metrics(y,p,t):
 tn,fp,fn,tp=confusion_matrix(y,p>=t,labels=[0,1]).ravel();return dict(threshold=float(t),balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2,accuracy=(tp+tn)/len(y),precision=tp/(tp+fp) if tp+fp else 0,recall=tp/(tp+fn),specificity=tn/(tn+fp),TN=int(tn),FP=int(fp),FN=int(fn),TP=int(tp))
def sweep(y,p):
 u=np.unique(p);ts=np.unique(np.r_[0,(u[:-1]+u[1:])/2,1.00000001,.5]);return [metrics(y,p,t) for t in ts]
def best(y,p):return sorted(sweep(y,p),key=lambda r:(-r['balanced_accuracy'],abs(r['threshold']-.5),r['threshold']))[0]
def model(typ):
 est=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026) if typ=='LR' else RandomForestClassifier(n_estimators=100,max_depth=2,min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
 return make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),est)
# Fixed model families, threshold selected from training-only out-of-fold scores.
records=[];results=[]
for key,typ in [('CP3','LR'),('CP_Structure9','LR'),('CP3','RF'),('CP_Structure9','RF')]:
 X=a[man['feature_sets'][key]].to_numpy()
 for val in ['LOCO','mutation_group_out']:
  splitgroups=names if val=='LOCO' else groups;pp=np.zeros(len(y));tt=np.zeros(len(y))
  for g in sorted(set(splitgroups)):
   tr=np.where(splitgroups!=g)[0];te=np.where(splitgroups==g)[0];inner=list(StratifiedKFold(3,shuffle=True,random_state=2026).split(tr,y[tr])) if val=='LOCO' else list(GroupKFold(3).split(tr,y[tr],groups[tr]));ip=np.zeros(len(tr))
   for itr,ite in inner:
    mm=model(typ).fit(X[tr[itr]],y[tr[itr]]);ip[ite]=mm.predict_proba(X[tr[ite]])[:,1]
   t=best(y[tr],ip)['threshold'];mm=model(typ).fit(X[tr],y[tr]);pp[te]=mm.predict_proba(X[te])[:,1];tt[te]=t
   for i in te:records.append(dict(model=key+'_'+typ,validation=val,construct=names[i],work=int(y[i]),score=float(pp[i]),training_selected_threshold=t,prediction=int(pp[i]>=t)))
  # Per-fold thresholds yield binary decisions, aggregate once across outer holdouts.
  d=metrics(y,(pp>=tt).astype(float),.5);d.pop('threshold');results.append(dict(model=key+'_'+typ,validation=val,min_fold_threshold=min(tt),median_fold_threshold=np.median(tt),max_fold_threshold=max(tt),**d));print('NESTED',results[-1],flush=True)
pd.DataFrame(records).to_csv(O/'dynamic_threshold_predictions.csv',index=False);pd.DataFrame(results).to_csv(O/'dynamic_threshold_validation.csv',index=False)
