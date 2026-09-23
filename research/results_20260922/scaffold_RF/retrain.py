#!/usr/bin/env python3
"""Reproduce fixed-recipe scaffold RF experiments without raw AF files."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import csv,json,platform
import numpy as np,sklearn,joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix,roc_auc_score,average_precision_score,matthews_corrcoef
A=Path(__file__).resolve().parent
for d in ['data','results','models','audit']:(A/d).mkdir(parents=True,exist_ok=True)
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def write(p,rows):
 ks=list(dict.fromkeys(k for r in rows for k in r));p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,ks);w.writeheader();w.writerows(rows)
cp=['pp_nonlocal4_high_per_res','pp_nonlocal12_high_per_res','pp_nonlocal24_high_per_res'];st=['plddt_core_mean','plddt_core_min','global_core_PAE','pae_contact_weighted','Rg_full_length_normalized','anisotropy']
full=read(A/'data/scaffolds24.csv');reg=read(A/'data/architecture_groups.csv')
orders={x['construct']:x['architecture'] for x in reg}
initial=[dict(x) for x in full if x['batch']=='previous']
initial14=read(A/'data/initial14_CP3.csv')
test=json.loads((A/'audit/protocol.json').read_text())['fixed4']
rf={'n_estimators':300,'max_depth':2,'min_samples_leaf':2,'max_features':'sqrt','class_weight':'balanced','random_state':2026,'n_jobs':1}
protocol={'source_commit':'0d68259e15c3d8d26b308bd938ef8acde16083f5','historical_commit':'958de2cc3567671c8f9c452cf819aa2133b630d5','random_state':2026,'labels':'historical scaffold success, not TRM activity','CP3':cp,'CP_structure9':cp+st,'RF':rf,'cutoff':0.5,'fixed4':test,'selection':'Reuse frozen four from previous Wiki; no resampling','frozen_before_training':True,'no_new_wetlab':True,'RF_S12_sensitivity':{'features':[cp[1]],'max_depth':1,'min_samples_leaf':1,'max_features':1.0},'early14_source':'12 old PUF12 rows plus two PUF11 CP3 records from pinned historical construct_features.csv','policy':'Fixed RF recipes, no selection from test scores; S12 sensitivity inherits historical feature-selection exposure.'}
(A/'audit/executed_protocol.json').write_text(json.dumps(protocol,indent=2))
pred=[];summ=[];splits=[]
def metric(y,p):
 q=np.asarray(p)>=.5;tn,fp,fn,tp=confusion_matrix(y,q,labels=[0,1]).ravel();both=len(set(y))==2
 return {'n':len(y),'correct':int(tp+tn),'accuracy':float(np.mean(q==y)),'BA':float(.5*(tp/(tp+fn)+tn/(tn+fp))) if both else None,'AUC':float(roc_auc_score(y,p)) if both else None,'AP':float(average_precision_score(y,p)) if both else None,'MCC':float(matthews_corrcoef(y,q)),'TP':int(tp),'TN':int(tn),'FP':int(fp),'FN':int(fn)}
def run(dataset,dsname,fs,columns,params,mode):
 ids=[x['construct'] for x in dataset];X=np.array([[float(x[c]) for c in columns] for x in dataset]);y=np.array([int(x['success']) for x in dataset]);b=np.array([x['batch'] for x in dataset]);idx=np.arange(len(ids))
 if mode=='LOCO':folds=[(ids[j],idx[idx!=j],np.array([j])) for j in idx]
 elif mode=='batch1_to_batch2':folds=[('batch2',idx[b=='previous'],idx[b=='new'])]
 elif mode=='fixed4':te=np.array([ids.index(x) for x in test]);folds=[('frozen4',np.array([j for j in idx if j not in te]),te)]
 elif mode=='architecture_out':folds=[(g,np.array([j for j in idx if orders[ids[j]]!=g]),np.array([j for j in idx if orders[ids[j]]==g])) for g in sorted(set(orders.values()))]
 else:raise ValueError(mode)
 ys=[];ps=[]
 for fold,tr,te in folds:
  assert not set(tr)&set(te)
  model=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),RandomForestClassifier(**params));model.fit(X[tr],y[tr]);scores=model.predict_proba(X[te])[:,1]
  for j,p in zip(te,scores):pred.append({'dataset':dsname,'validation':mode,'model':fs,'fold':fold,'construct':ids[j],'true_work':int(y[j]),'score':float(p),'predicted_work':int(p>=.5),'correct':int((p>=.5)==y[j]),'n_train':len(tr)})
  splits.append({'dataset':dsname,'validation':mode,'model':fs,'fold':fold,'train':[ids[k] for k in tr],'test':[ids[k] for k in te]})
  if mode in ['fixed4','batch1_to_batch2']:joblib.dump(model,A/'models'/f'{dsname}_{fs}_{mode}.joblib',compress=3)
  ys.extend(y[te]);ps.extend(scores)
 row={'dataset':dsname,'validation':mode,'model':fs,**metric(np.array(ys),ps)};summ.append(row);print(row,flush=True)
run(initial14,'initial14','RF_CP3',cp,rf,'LOCO');run(initial,'initial12_PUF12','RF_CP3',cp,rf,'LOCO')
run(initial14,'initial14','RF_S12_sensitivity',[cp[1]],dict(rf,max_depth=1,min_samples_leaf=1,max_features=1.0),'LOCO')
for fs,cols in [('RF_CP3',cp),('RF_CP_structure9',cp+st)]:
 for mode in ['LOCO','batch1_to_batch2','fixed4','architecture_out']:run(full,'PUF12_24',fs,cols,rf,mode)
write(A/'results/metrics.csv',summ);write(A/'results/predictions.csv',pred);(A/'audit/splits.json').write_text(json.dumps(splits,indent=2));(A/'audit/environment.json').write_text(json.dumps({'python':platform.python_version(),'sklearn':sklearn.__version__,'numpy':np.__version__},indent=2))
print('DONE',A,flush=True)
