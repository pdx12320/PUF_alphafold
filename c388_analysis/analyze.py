import os,json,pathlib,sys
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold,GroupKFold
from sklearn.metrics import roc_auc_score,average_precision_score,precision_recall_curve,auc,confusion_matrix,brier_score_loss
O=pathlib.Path('c388_analysis');old=O/'reference_models';lab=pd.read_csv(O/'labels.csv');features=pd.read_csv(O/'construct_features.csv');meta=pd.read_csv(O/'sequence_mapping.csv').drop_duplicates('construct');a=features[features.scope=='matched512'].merge(lab,on='construct').merge(meta[['construct','mutation_group','mutations','sequence_sha256']],on='construct').sort_values('construct').reset_index(drop=True);a.to_csv(O/'training_data.csv',index=False);y=a.work.to_numpy();groups=a.mutation_group.to_numpy();names=a.construct.to_numpy()
cp3=[f'pp_nonlocal{k}_high_per_res' for k in [4,12,24]];cp14=json.load(open(old/'old_feature_manifest.json'))['feature_sets']['CP_summary14'];st=['plddt_core_mean','plddt_core_min','global_core_PAE','pae_contact_weighted','Rg_full_length_normalized','anisotropy'];rna=['pr_cp_per_nt','pr_rna_max_mean','pr_rna_coverage_05','pr_protein_coverage_05'];inter=['adj_CP_mean','adj_CP_min','nonadj_CP_mean'];sets={'CP3':cp3,'CP14':cp14,'Structure6':st,'CP_Structure9':cp3+st,'RNA4':rna,'CP_RNA7':cp3+rna,'CP_Interface6':cp3+inter}
def model(typ):
 e=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026) if typ=='LR' else RandomForestClassifier(n_estimators=100,max_depth=2,min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
 return make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),e)
def metrics(y,p):
 tn,fp,fn,tp=confusion_matrix(y,p>=.5,labels=[0,1]).ravel();pr,re,_=precision_recall_curve(y,p);return dict(n=len(y),n_work=int(sum(y)),AUC=roc_auc_score(y,p) if len(set(y))==2 else np.nan,AP=average_precision_score(y,p),PR_AUC=auc(re,pr),balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2 if tp+fn and tn+fp else np.nan,precision=tp/(tp+fp) if tp+fp else 0,recall=tp/(tp+fn) if tp+fn else np.nan,TN=int(tn),FP=int(fp),FN=int(fn),TP=int(tp),brier=float(brier_score_loss(y,p)))
res=[];pred=[]
def record(family,validation,method,yy,pp,indices=None,scope='matched512'):
 if indices is None:indices=np.arange(len(a))
 res.append(dict(family=family,validation=validation,model=method,scope=scope,**metrics(yy,pp)))
 for i,v in zip(indices,pp):pred.append(dict(family=family,validation=validation,model=method,scope=scope,construct=names[i],C388_mean=a.C388_mean.iloc[i],work=int(y[i]),mutation_group=groups[i],score=float(v),prediction=int(v>=.5)))
# Frozen artifacts: zero use of new labels in fitting or thresholding.
for scope in ['matched512','full519']:
 b=features[features.scope==scope].set_index('construct').loc[names]
 for key,cols in [('CP_density3',cp3),('S12_only',[cp3[1]])]:
  for typ in ['RF','LR']:
   m=joblib.load(old/(key+'_'+typ+'.joblib'));p=m.predict_proba(b[cols].to_numpy())[:,1];record('frozen', 'transfer',key+'_'+typ,y,p,scope=scope)
# Original training ranges: target RNA changes 13->15nt, sequence prefix is explicit.
trainold=pd.read_csv(old/'main_features.csv');ood=[]
for _,r in a.iterrows():
 ood.append(dict(construct=r.construct,**{f+'_outside':bool(r[f]<trainold[f].min() or r[f]>trainold[f].max()) for f in cp3}))
pd.DataFrame(ood).to_csv(O/'frozen_range_check.csv',index=False)
pd.DataFrame(res).to_csv(O/'metrics.csv',index=False);print('FROZEN',pd.DataFrame(res)[['model','scope','AUC','balanced_accuracy','TN','FP','FN','TP']].to_string(index=False),flush=True)
def splits(validation,ix):
 gs=names[ix] if validation=='LOCO' else groups[ix]
 for g in sorted(set(gs)):
  yield ix[gs!=g],ix[gs==g]
imp=[]
for val in ['LOCO','mutation_group_out']:
 for key,cols in sets.items():
  X=a[cols].to_numpy()
  for typ in ['LR','RF']:
   pp=np.zeros(len(a))
   for tr,te in splits(val,np.arange(len(a))):
    m=model(typ).fit(X[tr],y[tr]);pp[te]=m.predict_proba(X[te])[:,1]
    if typ=='RF':
     for f,v in zip(cols,m[-1].feature_importances_):imp.append(dict(validation=val,feature_set=key,heldout=names[te[0]],feature=f,importance=v))
   record('retrained',val,key+'_'+typ,y,pp)
  print('CV',val,key,flush=True)
 # Position-only baseline: leave unseen group gives training prior.
 pp=[]
 for tr,te in splits(val,np.arange(len(a))):
  pass
 pp=np.zeros(len(a))
 for tr,te in splits(val,np.arange(len(a))):
  X=pd.get_dummies(a.mutation_group).to_numpy(float);m=model('LR').fit(X[tr],y[tr]);pp[te]=m.predict_proba(X[te])[:,1]
 record('position_baseline',val,'mutation_position_LR',y,pp)
pd.DataFrame(imp).to_csv(O/'RF_fold_importance.csv',index=False)
pd.DataFrame(res).to_csv(O/'metrics.csv',index=False);pd.DataFrame(pred).to_csv(O/'predictions.csv',index=False)
# Nested selection guards against choosing a feature/model combo on outer predictions.
candidates=[(k,t) for k in ['CP3','Structure6','CP_Structure9','RNA4','CP_RNA7'] for t in ['LR','RF']];choices=[]
for val in ['LOCO','mutation_group_out']:
 pp=np.zeros(len(a))
 for fi,(tr,te) in enumerate(splits(val,np.arange(len(a)))):
  inner=list(StratifiedKFold(3,shuffle=True,random_state=2026).split(tr,y[tr])) if val=='LOCO' else list(GroupKFold(3).split(tr,y[tr],groups[tr]))
  scores=[]
  for key,typ in candidates:
   X=a[sets[key]].to_numpy();ip=np.zeros(len(tr))
   for itr,ite in inner:
    train,test=tr[itr],tr[ite];assert len(set(y[train]))==2;m=model(typ).fit(X[train],y[train]);ip[ite]=m.predict_proba(X[test])[:,1]
   scores.append(roc_auc_score(y[tr],ip))
  best=int(np.argmax(scores));key,typ=candidates[best];X=a[sets[key]].to_numpy();m=model(typ).fit(X[tr],y[tr]);pp[te]=m.predict_proba(X[te])[:,1];choices.append(dict(validation=val,heldout=';'.join(names[te]),selected=key+'_'+typ,inner_AUC=scores[best]))
  if fi%5==0:print('NESTED',val,fi+1,flush=True)
 record('nested_selection',val,'inner_select_10_candidates',y,pp)
pd.DataFrame(choices).to_csv(O/'nested_selection.csv',index=False)
# Sensitivity to overlapping WT core and two zeros encoded as absent C>T variant records.
for subset,mask in [('exclude_WT',a.construct!='PUF12-9'),('exclude_undetected_zero',~a.zero_variant_record)]:
 ix=np.where(mask)[0]
 for key,typ in candidates:
  X=a[sets[key]].to_numpy();pp=np.zeros(len(a))
  for tr,te in splits('LOCO',ix):m=model(typ).fit(X[tr],y[tr]);pp[te]=m.predict_proba(X[te])[:,1]
  record('sensitivity',subset,key+'_'+typ,y[ix],pp[ix],ix)
# Final fits are for reuse only; never treated as performance.
# Keep all evaluation candidates above; export only retained research artifacts.
retention=json.loads((O/"model_retention.json").read_text())
for key,typ in candidates:
 if key+"_"+typ not in retention["retained_models"]:
  continue
 m=model(typ).fit(a[sets[key]].to_numpy(),y);joblib.dump(m,O/(key+'_'+typ+'.joblib'),compress=3)
pd.DataFrame(res).to_csv(O/'metrics.csv',index=False);pd.DataFrame(pred).to_csv(O/'predictions.csv',index=False);json.dump({'feature_sets':sets,'nested_candidates':candidates,'new_label':'mean C388 >= 0.5','preprocessing':'remove SGSETPG first7 residues for match to old 512aa range, known before label analysis','RF':dict(trees=100,depth=2,min_leaf=3,max_features=1.,class_weight='balanced'),'LR_C':.1,'threshold':.5,'primary_input':'C388 RNA context only; model->seed->construct equal means'},open(O/'manifest.json','w'),indent=2)
print('DONE',pd.DataFrame(res).query("family != 'sensitivity'")[['family','validation','model','AUC','balanced_accuracy','TN','FP','FN','TP']].to_string(index=False),flush=True)
