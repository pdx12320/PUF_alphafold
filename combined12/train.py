import os,json,collections
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_auc_score,average_precision_score,balanced_accuracy_score,accuracy_score,confusion_matrix,precision_score,recall_score,f1_score
from xgboost import XGBClassifier
O='combined12';old=pd.read_csv('extended/all_construct_features.csv');old=old[~old.construct.str.startswith('puf_11')].copy();old['batch']='previous';old['design_id']=np.nan
new=pd.read_csv('new_batch/final_summary.csv');new=new[new.repeat_count==12].copy();new['success']=(new.design_id==1).astype(int);new['batch']='new';new['protein_length']=new.length
cp=[c for c in old if c.startswith(('pp_','pr_'))];struct=['plddt_protein_mean','pae_contact_weighted','Rg_full_length_normalized'];cols=['construct','batch','design_id','success','protein_length']+cp+struct
a=pd.concat([old[cols],new[cols]],ignore_index=True);assert len(a)==24 and a.success.sum()==4 and len(new)==12 and len(old)==12;a.to_csv(O+'/training_data.csv',index=False);y=a.success.values
sets={'CP_density3':['pp_nonlocal4_high_per_res','pp_nonlocal12_high_per_res','pp_nonlocal24_high_per_res'],'CP_summary14':cp,'CP_structure17':cp+struct}

def est(typ,yy):
 if typ=='LR':m=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026)
 elif typ=='RF':m=RandomForestClassifier(n_estimators=300,max_depth=2,min_samples_leaf=2,max_features='sqrt',class_weight='balanced',random_state=2026,n_jobs=1)
 elif typ=='XGB':m=XGBClassifier(n_estimators=100,max_depth=1,learning_rate=.05,min_child_weight=1,reg_lambda=10,reg_alpha=1,subsample=1,colsample_bytree=.8,scale_pos_weight=(len(yy)-sum(yy))/sum(yy),objective='binary:logistic',eval_metric='logloss',tree_method='hist',random_state=2026,n_jobs=1)
 else:m=DecisionTreeClassifier(max_depth=1,class_weight='balanced',random_state=2026)
 return make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),m)

def metrics(yy,p):
 tn,fp,fn,tp=confusion_matrix(yy,p>=.5,labels=[0,1]).ravel();return {'AUC':roc_auc_score(yy,p),'AP':average_precision_score(yy,p),'accuracy':accuracy_score(yy,p>=.5),'balanced_accuracy':balanced_accuracy_score(yy,p>=.5),'precision':precision_score(yy,p>=.5,zero_division=0),'recall':recall_score(yy,p>=.5),'F1':f1_score(yy,p>=.5),'TN':int(tn),'FP':int(fp),'FN':int(fn),'TP':int(tp)}
records=[];predictions=[];importance=[];coefs=[]
for key,features in sets.items():
 X=a[features].values;assert np.isfinite(X).all()
 for typ in ['LR','RF','XGB']+(['Stump'] if key=='CP_density3' else []):
  prob=np.zeros(len(a))
  for i in range(len(a)):
   tr=np.arange(len(a))!=i;m=est(typ,y[tr]).fit(X[tr],y[tr]);prob[i]=m.predict_proba(X[i:i+1])[0,1]
   if typ in ['RF','XGB']:
    for f,v in zip(features,m[-1].feature_importances_):importance.append({'model':key+'_'+typ,'heldout':a.construct.iloc[i],'feature':f,'importance':v})
  records.append({'validation':'combined_LOCO','model':key+'_'+typ,**metrics(y,prob)})
  for i,p in enumerate(prob):predictions.append({'validation':'combined_LOCO','model':key+'_'+typ,**a.iloc[i][['construct','batch','design_id','success']].to_dict(),'score':p,'prediction':int(p>=.5)})
  # Previous-only 12-repeat training is tested on the complete labeled new 12-repeat batch.
  tr=a.batch.values=='previous';te=~tr;m=est(typ,y[tr]).fit(X[tr],y[tr]);p=m.predict_proba(X[te])[:,1];records.append({'validation':'previous12_to_new12','model':key+'_'+typ,**metrics(y[te],p)})
  for i,v in zip(np.where(te)[0],p):predictions.append({'validation':'previous12_to_new12','model':key+'_'+typ,**a.iloc[i][['construct','batch','design_id','success']].to_dict(),'score':v,'prediction':int(v>=.5)})
  m=est(typ,y).fit(X,y);path=O+'/'+key+'_'+typ+'.joblib';joblib.dump(m,path,compress=3);assert np.max(abs(joblib.load(path).predict_proba(X)-m.predict_proba(X)))<1e-10
  if typ=='LR':
   for f,mu,scale,beta in zip(features,m[1].mean_,m[1].scale_,m[-1].coef_[0]):coefs.append({'model':key+'_'+typ,'feature':f,'mean':mu,'scale':scale,'standardized_coefficient':beta,'intercept':float(m[-1].intercept_[0])})
 # Manifest
json.dump({'feature_sets':sets,'labels':'previous original 3 successes; new ID1 only successful; only 12-repeat constructs','n_previous':12,'n_new':12,'n_success':4,'n_failure':20,'threshold':.5,'random_state':2026,'interpretation':'uncalibrated score; LOCO predictions are separate from final refitted model scores'},open(O+'/manifest.json','w'),indent=2)
pd.DataFrame(records).to_csv(O+'/model_metrics.csv',index=False);pd.DataFrame(predictions).to_csv(O+'/heldout_predictions.csv',index=False);pd.DataFrame(coefs).to_csv(O+'/LR_formulas.csv',index=False);imp=pd.DataFrame(importance);imp.to_csv(O+'/fold_feature_importance.csv',index=False);imp.groupby(['model','feature']).importance.agg(['mean','std']).reset_index().sort_values(['model','mean'],ascending=[True,False]).to_csv(O+'/feature_importance_summary.csv',index=False)
legacy=new[['construct','design_id','success','XGB_score','RF_score','LR_score']].copy();legacy.to_csv(O+'/previously_issued_scores_with_new_labels.csv',index=False);print(pd.DataFrame(records).to_string(index=False));print('Legacy previous14 frozen XGB:',metrics(new.success.values,new.XGB_score.values))
