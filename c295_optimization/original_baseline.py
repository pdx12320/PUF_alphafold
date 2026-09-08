import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,io,contextlib,warnings
import numpy as np,pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold,LeaveOneOut
from sklearn.metrics import confusion_matrix
P=Path(__file__).resolve().parent;O=P/'results'
sys.path.insert(0,str(P/'sources'))
from ContactSeek.ContactExtraction import query_cp_with_top_n
from ContactSeek.FindContactResidue import find_contact_residues
from ContactSeek.CCRegionFinding import find_consensus_contact_regions
d=pd.read_csv(O/'labels.csv');cons=d.construct.tolist();y=d.label_large_decrease.values
audit=json.loads((O/'sequence_and_window_audit.json').read_text());pos=[set(r['positions']) for r in audit];mats=np.load(P/'inputs/pr_by_seed.npz');wt=mats['wt_puf12_9_c295'].mean(0)
raw=[];dif=[]
for c in cons:
 q=[query_cp_with_top_n(m,list(range(1,513)),[(0,5),(5,10),(10,15)],top_n=3,ref_array=wt) for m in mats[c+'_c295']]
 raw.append(np.mean([a['top_n_input_prob'] for a in q],0));dif.append(np.mean([a['diff_prob'] for a in q],0))
raw=np.array(raw);dif=np.array(dif);regional=np.concatenate([raw,dif],2);cache={}
def features(tr):
 key=tuple(tr)
 if key in cache:return cache[key]
 with warnings.catch_warnings(),contextlib.redirect_stdout(io.StringIO()):
  warnings.simplefilter('ignore');keep=find_contact_residues(list(raw[tr]),list(dif[tr]),min_cp_threshold=.15,min_diff_threshold=.1,verbose=False)
  if keep.sum()>1:_,_,regions=find_consensus_contact_regions([dict(y_g3=np.zeros(len(tr)),cp_raw_cas_nuc=list(raw[tr]))],keep,correlation_threshold=.6,band_width=7,max_merge_iterations=10,verbose=False,protein_name='PUF')
  else:regions=[dict(positions=np.where(keep)[0])] if keep.any() else []
 x=np.concatenate([regional[:,r['positions']].mean(1) for r in regions],1) if regions else np.zeros((len(y),1));cache[key]=x;return x
def fit(tr,te):
 if len(set(y[tr]))<2:return np.full(len(te),y[tr].mean())
 x=features(tr);m=RandomForestClassifier(n_estimators=100,max_depth=2,min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
 return m.fit(x[tr],y[tr]).predict_proba(x[te])[:,1]
def purged(idx):
 return [(np.array([j for j in idx if p not in pos[j]]),np.array([j for j in idx if p in pos[j]]),str(p)) for p in sorted(set.union(*(pos[j] for j in idx)))]
def best(yt,p,ww):
 u=np.unique(p);opts=[]
 for t in np.unique(np.r_[0,.5,1.00000001,(u[:-1]+u[1:])/2]):
  tn,fp,fn,tp=confusion_matrix(yt,p>=t,labels=[0,1],sample_weight=ww).ravel();ba=.5*(tp/(tp+fn)+tn/(tn+fp));opts.append((-ba,abs(t-.5),t))
 return min(opts)[2]
rows=[]
for val,folds in [('LOOCV',[(a,b,str(b[0])) for a,b in LeaveOneOut().split(y)]),('Strict_position_out',purged(np.arange(len(y))))]:
 for tr,te,g in folds:
  if val=='LOOCV':inner=[(tr[a],tr[b]) for a,b in StratifiedKFold(3,shuffle=True,random_state=2026).split(tr,y[tr])]
  else:inner=[(a,b) for a,b,_ in purged(tr) if len(a) and len(b)]
  ids=np.concatenate([b for a,b in inner]);pp=np.concatenate([fit(a,b) for a,b in inner]);ww=1/np.bincount(ids,minlength=len(y))[ids];t=best(y[ids],pp,ww);p=fit(tr,te)
  for i,score in zip(te,p):rows.append(dict(validation=val,fold=g,family='original_CCR_RF',construct=cons[i],y=int(y[i]),score=score,threshold=t,prediction=int(score>=t),delta_editing_pp=d.delta_editing_pp.iloc[i]))
  print('DONE',val,g,flush=True)
pd.DataFrame(rows).to_csv(O/'original_CCR_predictions.csv',index=False)
