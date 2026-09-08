import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,io,contextlib,warnings,json,re
import numpy as np,pandas as pd,joblib
from sklearn.model_selection import LeaveOneOut,LeaveOneGroupOut,StratifiedKFold,GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix,roc_auc_score
from threadpoolctl import threadpool_limits
threadpool_limits(1)
sys.path.insert(0,str(Path('sources/menghaowei/ContactSeek').resolve()))
from ContactSeek.ContactExtraction import query_cp_with_top_n
from ContactSeek.FindContactResidue import find_contact_residues
from ContactSeek.CCRegionFinding import find_consensus_contact_regions
out=Path('results');out.mkdir(exist_ok=True)
data=pd.read_csv('inputs/construct_site_summary.csv');cons=sorted(data.construct.unique());n=len(cons);assert n==30
repo=pd.read_csv('inputs/repository_features.csv',index_col=0);mats=np.load('inputs/pr_by_seed.npz')
groups=np.array(['+'.join(re.findall(r'(?:^|_plus_)p(\d+)_',c)) for c in cons])
specs=[('Total_CP_LR','total','LR'),('Interface_CP_LR','interface','LR'),('S4_S12_S24_LR','density','LR'),('Repository_LR','repo','LR'),('ContactSeek_Top3_RF','regional','RF'),('ContactSeek_CCR_RF','ccr','RF')]
def model(kind):
 return make_pipeline(StandardScaler(),LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026)) if kind=='LR' else RandomForestClassifier(n_estimators=100,max_depth=2,min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
def fit(kind,x,y,xt):
 if len(np.unique(y))<2:return np.full(len(xt),y.mean())
 return model(kind).fit(x,y).predict_proba(xt)[:,1]
def met(y,p,t):
 tn,fp,fn,tp=confusion_matrix(y,p>=t,labels=[0,1]).ravel();return dict(BA=(tp/(tp+fn)+tn/(tn+fp))/2,accuracy=(tn+tp)/len(y),sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),TN=tn,FP=fp,FN=fn,TP=tp)
def best(y,p):
 u=np.unique(p);ts=np.unique(np.r_[0,.5,1.00000001,(u[:-1]+u[1:])/2]);rs=[dict(score_threshold=t,**met(y,p,t)) for t in ts];return sorted(rs,key=lambda r:(-r['BA'],abs(r['score_threshold']-.5),r['score_threshold']))[0]
scans=[];preds=[];validated=[];vp=[];finals={};labels=[]
for site in ['C388','C871','C295']:
 ss=data[data.site==site].set_index('construct').loc[cons];rate=ss.delta_editing_pp.values;assert np.allclose(rate,ss.editing_pct.values-ss.wt_editing_pct.values);assert ss.wt_editing_pct.nunique()==1;wt=mats['wt_puf12_9_'+site.lower()].mean(0)
 raw=[];dif=[]
 for c in cons:
  qs=[query_cp_with_top_n(m,list(range(1,513)),[(0,5),(5,10),(10,15)],top_n=3,ref_array=wt) for m in mats[c+'_'+site.lower()]]
  raw.append(np.mean([q['top_n_input_prob'] for q in qs],0));dif.append(np.mean([q['diff_prob'] for q in qs],0))
 raw=np.array(raw);dif=np.array(dif);regional=np.concatenate([raw,dif],2)
 rx=repo.loc[[c+'_'+site.lower() for c in cons]].copy();rx.index=cons;rx-=repo.loc['wt_puf12_9_'+site.lower()].values
 static=dict(total=rx[['pr_cp_sum']].values,interface=rx[[c for c in rx if c.startswith('pr_')]].values,density=rx[[f'pp_nonlocal{k}_high_per_res' for k in [4,12,24]]].values,repo=rx.values)
 cache={};defs={}
 def features(tr):
  key=tuple(tr)
  if key in cache:return cache[key]
  with warnings.catch_warnings(),contextlib.redirect_stdout(io.StringIO()):
   warnings.simplefilter('ignore');keep=find_contact_residues(list(raw[tr]),list(dif[tr]),min_cp_threshold=.15,min_diff_threshold=.1,verbose=False)
   if keep.sum()>1:_,_,regions=find_consensus_contact_regions([dict(y_g3=np.zeros(len(tr)),cp_raw_cas_nuc=list(raw[tr]))],keep,correlation_threshold=.6,band_width=7,max_merge_iterations=10,verbose=False,protein_name='PUF')
   else:regions=[dict(positions=np.where(keep)[0])] if keep.any() else []
  ccr=np.concatenate([regional[:,r['positions']].mean(1) for r in regions],1) if regions else np.zeros((n,1))
  cache[key]=dict(static,regional=regional[:,keep].reshape(n,-1) if keep.any() else np.zeros((n,1)),ccr=ccr);defs[key]=(keep,regions);return cache[key]
 # Coarse grid plus adjacent observed midpoints; collapse equivalent labels. >=5 per class.
 u=np.unique(rate);candidates=list(np.arange(-100,101,5))+list((u[:-1]+u[1:])/2);thresholds={}
 for t in candidates:
  y=(rate>=t).astype(int)
  if min(y.sum(),n-y.sum())>=5:thresholds.setdefault(tuple(y),float(t))
 splits=list(LeaveOneOut().split(rate));site_scan=[]
 for yy,t in thresholds.items():
  y=np.array(yy);lo=rate[y==0].max();hi=rate[y==1].min()
  for name,f,kind in specs:
   p=np.zeros(n)
   for tr,te in splits:
    xx=features(tr)[f];p[te]=fit(kind,xx[tr],y[tr],xx[te])
   result=dict(site=site,label_threshold_delta_pp=t,equivalent_delta_lower_exclusive_pp=lo,equivalent_delta_upper_inclusive_pp=hi,model=name,n_high=int(y.sum()),n_low=int(n-y.sum()),AUC=roc_auc_score(y,p),BA_fixed_score05=met(y,p,.5)['BA'],**best(y,p))
   scans.append(result);site_scan.append(result)
   for i,c in enumerate(cons):preds.append(dict(site=site,label_threshold_delta_pp=t,model=name,construct=c,editing_pct=ss.editing_pct.iloc[i],delta_editing_pp=rate[i],actual_class=y[i],score=p[i]))
  print('SCAN',site,t,flush=True)
 pd.DataFrame(scans).to_csv(out/'exploratory_threshold_scan.csv',index=False)
 # Only the top exploratory label/model is followed up; selection bias remains explicitly reported.
 chosen=sorted(site_scan,key=lambda r:(-r['BA'],-r['AUC'],r['model'],r['label_threshold_delta_pp']))[0]
 t=chosen['label_threshold_delta_pp'];name=chosen['model'];_,f,kind=next(s for s in specs if s[0]==name);y=(rate>=t).astype(int)
 for val,splits in [('LOOCV',list(LeaveOneOut().split(rate))),('Position_group_out',list(LeaveOneGroupOut().split(rate,groups=groups)))]:
  p=np.zeros(n);ts=np.zeros(n)
  for tr,te in splits:
   if val=='LOOCV':inner=list(StratifiedKFold(min(3,np.bincount(y[tr],minlength=2).min()),shuffle=True,random_state=2026).split(tr,y[tr]))
   else:inner=list(GroupKFold(min(3,len(np.unique(groups[tr])))).split(tr,y[tr],groups[tr]))
   ip=np.zeros(len(tr))
   for a,b in inner:
    ta=tr[a];tb=tr[b];xx=features(ta)[f];ip[b]=fit(kind,xx[ta],y[ta],xx[tb])
   th=best(y[tr],ip)['score_threshold'];xx=features(tr)[f];p[te]=fit(kind,xx[tr],y[tr],xx[te]);ts[te]=th
  validated.append(dict(site=site,model=name,label_threshold_delta_pp=t,validation=val,n_high=int(y.sum()),n_low=int(n-y.sum()),AUC=roc_auc_score(y,p),median_score_threshold=np.median(ts),min_score_threshold=min(ts),max_score_threshold=max(ts),majority_accuracy=max(y.mean(),1-y.mean()),**met(y,(p>=ts).astype(float),.5)))
  for i,c in enumerate(cons):vp.append(dict(site=site,model=name,validation=val,construct=c,editing_pct=ss.editing_pct.iloc[i],delta_editing_pp=rate[i],actual_class=y[i],score=p[i],score_threshold=ts[i],prediction=int(p[i]>=ts[i]),correct=bool((p[i]>=ts[i])==y[i])))
  print('VALIDATED',validated[-1],flush=True)
 allidx=np.arange(n);xx=features(allidx)[f];md=model(kind).fit(xx,y)
 # Deployment score threshold from full-training OOF, not an independent estimate.
 cp=pd.DataFrame(preds);cp=cp[(cp.site==site)&(cp.model==name)&(cp.label_threshold_delta_pp==t)].set_index('construct').loc[cons]
 finals[site]=dict(model=md,feature=f,label_threshold_delta_pp=t,score_threshold=best(y,cp.score.values)['score_threshold'],contact_mask=defs[tuple(allidx)][0],regions=defs[tuple(allidx)][1],wt_pr=wt,wt_editing_pct=float(ss.wt_editing_pct.iloc[0]),label_definition="delta_editing_pp >= label_threshold_delta_pp",repository_columns=list(rx),constructs=cons,exploratory=True)
 for i,c in enumerate(cons):labels.append(dict(site=site,construct=c,editing_pct=ss.editing_pct.iloc[i],delta_editing_pp=rate[i],wt_editing_pct=ss.wt_editing_pct.iloc[0],threshold_delta_pp=t,class_high=int(y[i]),position_group=groups[i]))
 np.savez_compressed(out/f'{site}_feature_arrays.npz',raw_top3=raw,delta_max=dif,**static)
 pd.DataFrame(validated).to_csv(out/'training_only_score_threshold_validation.csv',index=False)
pd.DataFrame(preds).to_csv(out/'exploratory_oof_scores.csv',index=False);pd.DataFrame(vp).to_csv(out/'validated_predictions.csv',index=False);pd.DataFrame(labels).to_csv(out/'candidate_labels.csv',index=False);joblib.dump(finals,out/'candidate_models.joblib',compress=3)
print('ALL DONE',flush=True)
