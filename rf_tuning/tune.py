import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,time,itertools,joblib
import numpy as np,pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score,average_precision_score,confusion_matrix,balanced_accuracy_score
from joblib import Parallel,delayed
from pathlib import Path
O=Path('rf_tuning');a=pd.read_csv('combined12/training_data.csv');y=a.success.to_numpy();sets=json.load(open('combined12/manifest.json'))['feature_sets']
# Fixed, modest search space; feature-set selection occurs inside each outer training fold.
params=[dict(max_depth=d,min_samples_leaf=l,max_features=f,class_weight=w) for d,l,f,w in [
(2,2,'sqrt','balanced'),(1,1,'sqrt','balanced'),(2,1,'sqrt','balanced'),(3,1,'sqrt','balanced'),
(4,1,'sqrt','balanced'),(None,1,'sqrt','balanced'),(2,3,'sqrt','balanced'),(3,2,'sqrt','balanced'),
(2,2,1.0,'balanced'),(3,1,1.0,'balanced'),(None,2,1.0,'balanced'),(2,1,'sqrt',None),
(3,2,'sqrt',None),(None,1,'sqrt',None),(2,2,'sqrt','balanced_subsample'),(3,2,1.0,'balanced_subsample')]]
candidates=[dict(id=len(params)*j+i,feature_set=s,**p) for j,s in enumerate(sets) for i,p in enumerate(params)]
json.dump({'candidates':candidates,'inner':'3-fold stratified CV; mean fold ROC AUC; seed2026','inner_trees':100,'outer_and_final_trees':300,'threshold_selection':'maximize inner out-of-fold balanced accuracy over fixed grid 0.10..0.80 step0.05; ties closest to 0.5','selection':'mean inner AUC, then AP, then candidate ID; full-data search is deployment only'},open(O/'search_protocol.json','w'),indent=2)
def model(c,n):return RandomForestClassifier(n_estimators=n,random_state=2026,n_jobs=1,**{k:c[k] for k in params[0]})
def metric(yy,p,t=.5):
 t=np.broadcast_to(t,len(yy));tn,fp,fn,tp=confusion_matrix(yy,p>=t,labels=[0,1]).ravel()
 return dict(AUC=roc_auc_score(yy,p),AP=average_precision_score(yy,p),accuracy=(tn+tp)/len(yy),balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2,TP=int(tp),FN=int(fn),FP=int(fp),TN=int(tn),precision=tp/max(tp+fp,1),recall=tp/(tp+fn))
def search(idx):
 yy=y[idx];splits=list(StratifiedKFold(3,shuffle=True,random_state=2026).split(idx,yy));rows=[];probs={}
 for c in candidates:
  X=a.iloc[idx][sets[c['feature_set']]].to_numpy();p=np.zeros(len(idx));au=[];ap=[]
  assert np.isfinite(X).all()
  for tr,va in splits:
   m=model(c,100).fit(X[tr],yy[tr]);p[va]=m.predict_proba(X[va])[:,1];au.append(roc_auc_score(yy[va],p[va]));ap.append(average_precision_score(yy[va],p[va]))
  rows.append(dict(candidate_id=c['id'],feature_set=c['feature_set'],inner_AUC=np.mean(au),inner_AP=np.mean(ap)));probs[c['id']]=p
 rows=sorted(rows,key=lambda r:(-r['inner_AUC'],-r['inner_AP'],r['candidate_id']))
 return rows,probs

def fold(i):
 t0=time.time();idx=np.delete(np.arange(len(y)),i);rows,probs=search(idx);out=[]
 for scope in ['joint']+list(sets):
  chosen=next(r for r in rows if scope=='joint' or r['feature_set']==scope);c=candidates[chosen['candidate_id']];X=a[sets[c['feature_set']]].to_numpy();m=model(c,300).fit(X[idx],y[idx]);p=float(m.predict_proba(X[i:i+1])[0,1]);grid=np.round(np.arange(.1,.801,.05),2);th=sorted(grid,key=lambda z:(-balanced_accuracy_score(y[idx],probs[c['id']]>=z),abs(z-.5),z))[0]
  out.append(dict(scope=scope,heldout=i,construct=a.construct.iloc[i],batch=a.batch.iloc[i],design_id=a.design_id.iloc[i],success=int(y[i]),score=p,inner_threshold=float(th),**chosen))
 pd.DataFrame(rows).to_csv(O/f'inner_search_fold_{i:02d}.csv',index=False);pd.DataFrame(out).to_csv(O/f'outer_fold_{i:02d}.csv',index=False)
 print(f'Completed outer fold {i+1}/24 in {time.time()-t0:.1f}s',flush=True);return out
if __name__=='__main__':
 results=Parallel(n_jobs=4)(delayed(fold)(i) for i in range(24));pred=pd.DataFrame(sum(results,[]));pred.to_csv(O/'nested_LOCO_predictions.csv',index=False);metrics=[]
 for scope,g in pred.groupby('scope'):
  for threshold in ['fixed_0.5','inner_selected']:
   metrics.append(dict(model=scope,threshold=threshold,**metric(g.success.to_numpy(),g.score.to_numpy(),.5 if threshold=='fixed_0.5' else g.inner_threshold.to_numpy())))
 pd.DataFrame(metrics).to_csv(O/'nested_LOCO_metrics.csv',index=False)
 rows,probs=search(np.arange(len(y)));pd.DataFrame(rows).to_csv(O/'full_data_inner_search.csv',index=False);final=[]
 for scope in ['joint']+list(sets):
  chosen=next(r for r in rows if scope=='joint' or r['feature_set']==scope);c=candidates[chosen['candidate_id']];features=sets[c['feature_set']];X=a[features].to_numpy();m=model(c,300).fit(X,y);p=O/f'RF_tuned_{scope}.joblib';joblib.dump(m,p,compress=3);assert np.allclose(m.predict_proba(X),joblib.load(p).predict_proba(X));pd.DataFrame({'feature':features,'importance':m.feature_importances_}).sort_values('importance',ascending=False).to_csv(O/f'feature_importance_{scope}.csv',index=False)
  grid=np.round(np.arange(.1,.801,.05),2);th=sorted(grid,key=lambda z:(-balanced_accuracy_score(y,probs[c['id']]>=z),abs(z-.5),z))[0]
  final.append(dict(scope=scope,features=features,parameters=c,n_estimators=300,random_state=2026,inner_selected_threshold=float(th),default_threshold=.5))
 json.dump(final,open(O/'final_model_manifest.json','w'),indent=2)
 print(pd.DataFrame(metrics).to_string(index=False),flush=True);print('DONE',flush=True)
