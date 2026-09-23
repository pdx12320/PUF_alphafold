import os,json
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score,average_precision_score,precision_recall_curve,auc,confusion_matrix
from xgboost import XGBClassifier
from stats import blocks,corr
O='architecture_validation';a=pd.read_csv(O+'/all_features.csv');z=np.load(O+'/all_arrays.npz');ids=[list(z['names']).index(n) for n in a.construct];D=z['density'][ids];CP=z['cp'][ids];ii,jj=np.triu_indices(432,12);C=CP[:,ii,jj];y=a.success.to_numpy();cp14=json.load(open('combined12/manifest.json'))['feature_sets']['CP_summary14'];cp3=['pp_nonlocal4_high_per_res','pp_nonlocal12_high_per_res','pp_nonlocal24_high_per_res'];struct=['plddt_core_mean','plddt_core_min','global_core_PAE','pae_contact_weighted','Rg_full_length_normalized','anisotropy'];repeat=['adj_CP_mean','adj_CP_min','nonadj_CP_mean'];feature_sets={'CP_summary14':cp14,'CP_density3':cp3,'structural_only':struct,'CP_structural':cp3+struct,'CP_interface':cp3+repeat,'full':cp3+struct+repeat,'S12_only':['pp_nonlocal12_high_per_res']}
def estimator(typ,yy):
 if typ=='LR':m=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026)
 elif typ=='RF':m=RandomForestClassifier(n_estimators=300,max_depth=2,min_samples_leaf=2,max_features='sqrt',class_weight='balanced',random_state=2026,n_jobs=1)
 else:m=XGBClassifier(n_estimators=100,max_depth=1,learning_rate=.05,min_child_weight=1,reg_lambda=10,reg_alpha=1,colsample_bytree=.8,subsample=1,scale_pos_weight=(len(yy)-sum(yy))/sum(yy),eval_metric='logloss',tree_method='hist',random_state=2026,n_jobs=1)
 return make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),m)
def metrics(yy,p):
 tn,fp,fn,tp=confusion_matrix(yy,p>=.5,labels=[0,1]).ravel();both=len(np.unique(yy))==2;pr,re,_=precision_recall_curve(yy,p) if both else (None,None,None)
 return dict(AUC=roc_auc_score(yy,p) if both else np.nan,AP=average_precision_score(yy,p) if both else np.nan,PR_AUC=auc(re,pr) if both else np.nan,balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2 if both else np.nan,precision=tp/(tp+fp) if tp+fp else 0,recall=tp/(tp+fn) if tp+fn else np.nan,TN=int(tn),FP=int(fp),FN=int(fn),TP=int(tp))
def local(tr):
 V=D[tr];valid=np.isfinite(V).all(0)&(np.nanstd(V,axis=0)>1e-6);rv=np.where(valid)[0];r=rv[np.argmax(abs(corr(V[:,rv],y[tr])))];groups=blocks(V[:,:396] if not valid[396:].all() else V);reg=max(groups,key=lambda g:abs(corr(V[:,g].mean(1)[:,None],y[tr])[0])) if groups else []
 va=np.isfinite(C[tr]).all(0)&(np.nanstd(C[tr],axis=0)>1e-6)&(np.nanmax(C[tr],axis=0)>=.1);pids=np.where(va)[0];p=pids[np.argmax(abs(corr(C[tr][:,pids],y[tr])))];values=np.c_[D[:,r],np.nanmean(D[:,reg],axis=1) if reg else np.zeros(len(y)),C[:,p]];return values,dict(residue_aligned=int(r+1),CCR_aligned=[int(x+1) for x in reg],pair_aligned=[int(ii[p]+1),int(jj[p]+1)])
records=[];imp=[];selection=[];foldmetrics=[];main=(a.repeat_count==12)&(a.batch!='TRM')
for subset,mask in [('PUF12_main',main),('plus_PUF11',(a.batch!='TRM')),('plus_TRM',(a.repeat_count==12))]:
 for validation in ['LOCO','architecture_out']:
  groups=a.construct if validation=='LOCO' else a.source_order
  for g in sorted(set(groups[mask])):
   te=np.where(mask&(groups==g))[0];tr=np.where(mask&(groups!=g))[0];assert len(set(y[tr]))==2
   L,sel=local(tr);selection.append(dict(subset=subset,validation=validation,heldout_group=g,**sel))
   for key,ff in feature_sets.items():
    if subset!='PUF12_main' and key not in ['CP_density3','S12_only','full']:continue
    X=a[ff].to_numpy();fn=ff[:]
    if key=='full':X=np.c_[X,L];fn+=['selected_residue','selected_CCR','selected_contact']
    for typ in (['RF'] if subset!='PUF12_main' else ['LR','RF','XGB']):
     m=estimator(typ,y[tr]).fit(X[tr],y[tr]);pred=m.predict_proba(X[te])[:,1];foldmetrics.append(dict(subset=subset,validation=validation,feature_set=key,model=typ,heldout_group=g,n_test=len(te),n_success=int(y[te].sum()),n_train_success=int(y[tr].sum()),**metrics(y[te],pred)))
     for i,v in zip(te,pred):records.append(dict(subset=subset,validation=validation,feature_set=key,model=typ,construct=a.construct.iloc[i],architecture=a.architecture.iloc[i],success=int(y[i]),score=float(v),prediction=int(v>=.5)))
     if subset=='PUF12_main' and typ=='RF':
      native=m[-1].feature_importances_;base=-(y[te]*np.log(np.clip(pred,1e-8,1))+(1-y[te])*np.log(np.clip(1-pred,1e-8,1)))
      for j,f in enumerate(fn):
       xm=X[te].copy();xm[:,j]=np.nanmedian(X[tr,j]);pp=m.predict_proba(xm)[:,1];loss=-(y[te]*np.log(np.clip(pp,1e-8,1))+(1-y[te])*np.log(np.clip(1-pp,1e-8,1)))
       imp.append(dict(validation=validation,feature_set=key,heldout_group=g,feature=f,MDI=float(native[j]),n_test=len(te),heldout_logloss_increase_sum=float((loss-base).sum())))
  print(subset,validation,'complete',flush=True)
pred=pd.DataFrame(records);pred.to_csv(O+'/heldout_predictions.csv',index=False);pd.DataFrame(foldmetrics).to_csv(O+'/per_group_metrics.csv',index=False);json.dump(selection,open(O+'/fold_local_selection.json','w'),indent=2);pd.DataFrame(imp).to_csv(O+'/fold_importance.csv',index=False)
res=[]
for (subset,val,key,typ),g in pred.groupby(['subset','validation','feature_set','model']):res.append(dict(subset=subset,validation=val,feature_set=key,model=typ,**metrics(g.success.to_numpy(),g.score.to_numpy())))
mt=pd.DataFrame(res);order=mt[mt.validation=='architecture_out'].copy();order['validation']='source_order_out';mt=pd.concat([mt,order]);mt.to_csv(O+'/validation_metrics.csv',index=False);porder=pred[pred.validation=='architecture_out'].copy();porder['validation']='source_order_out';pd.concat([pred,porder]).to_csv(O+'/heldout_predictions.csv',index=False)
# Fixed minimal deployable candidates; choose recommendation only after inspecting grouped validation.
for key in ['S12_only','CP_density3']:
 for typ in ['LR','RF']:
  ff=feature_sets[key];m=estimator(typ,y[main]).fit(a.loc[main,ff].to_numpy(),y[main]);path=O+'/'+key+'_'+typ+'.joblib';joblib.dump(m,path,compress=3);assert np.allclose(m.predict_proba(a.loc[main,ff].to_numpy()),joblib.load(path).predict_proba(a.loc[main,ff].to_numpy()))
json.dump({'feature_sets':feature_sets,'full_local_selection':'one residue + one coherent CCR + one contact, selected in training folds only','fixed_parameters':{'LR_C':.1,'RF_trees':300,'RF_depth':2,'RF_min_leaf':2,'XGB_depth':1,'XGB_lambda':10,'XGB_alpha':1},'threshold':.5},open(O+'/model_manifest.json','w'),indent=2)
print(mt[(mt.subset=='PUF12_main')&(mt.model=='RF')].to_string(index=False))
