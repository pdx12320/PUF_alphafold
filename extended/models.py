import os,json,collections
import numpy as np,pandas as pd,joblib
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score,average_precision_score,balanced_accuracy_score,accuracy_score,precision_score,recall_score,f1_score,matthews_corrcoef,confusion_matrix,log_loss
from xgboost import XGBClassifier
from common import build_selected
O=os.path.dirname(__file__);z=np.load(O+'/construct_arrays.npz');names=z['names'].tolist();arr={k:z[k] for k in z.files if k!='names'};old=pd.read_csv('previous/results/construct_features.csv').set_index('construct').loc[names];y=old.success.values;st=pd.read_csv(O+'/structural_construct_features.csv').set_index('construct').loc[names];st['Rg_full_length_normalized']=st.Rg_CA_A/old.protein_length**(1/3);basecols=[c for c in old if c.startswith(('pr_','pp_'))];B=old[basecols].values;E=np.c_[B,st.values];en=basecols+st.columns.tolist();F=old.join(st);formula_cols=['pp_nonlocal12_high_per_res','plddt_protein_mean','pae_contact_weighted','Rg_full_length_normalized'];Q=F[formula_cols].values;F.to_csv(O+'/all_construct_features.csv')

def estimator(typ,yy):
 if typ=='LR':m=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026)
 elif typ=='RF':m=RandomForestClassifier(n_estimators=300,max_depth=2,min_samples_leaf=2,max_features='sqrt',class_weight='balanced',random_state=2026,n_jobs=1)
 elif typ=='XGB':m=XGBClassifier(n_estimators=100,max_depth=1,learning_rate=.05,min_child_weight=1,reg_lambda=10,reg_alpha=1,subsample=1,colsample_bytree=.8,scale_pos_weight=(len(yy)-sum(yy))/sum(yy),objective='binary:logistic',eval_metric='logloss',tree_method='hist',random_state=2026,n_jobs=1)
 return make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),m)

def metric(yy,p):
 hard=p>=.5;tn,fp,fn,tp=confusion_matrix(yy,hard).ravel();return {'AUC':roc_auc_score(yy,p),'AP':average_precision_score(yy,p),'accuracy':accuracy_score(yy,hard),'balanced_accuracy':balanced_accuracy_score(yy,hard),'precision':precision_score(yy,hard,zero_division=0),'recall':recall_score(yy,hard),'F1':f1_score(yy,hard),'MCC':matthews_corrcoef(yy,hard),'TN':int(tn),'FP':int(fp),'FN':int(fn),'TP':int(tp)}

metrics=[];prediction=[];selection=[];importance=[];lossimportance=[]
for subset in ['all14','12repeat']:
 take=np.ones(len(y),bool) if subset=='all14' else np.array([not n.startswith('puf_11') for n in names]);indices=np.where(take)[0];out=collections.defaultdict(list)
 for i in indices:
  train=take.copy();train[i]=False;Xsel,selnames,meta,regions=build_selected(arr,y,train)
  if subset=='all14':selection.append({'heldout':names[i],'residues':[d['indices'][0]+1 for d in meta if d['kind']=='residue'],'CCRs':regions})
  feature_sets={'summary':(B,basecols),'structural':(E,en),'full':(np.c_[E,Xsel],en+selnames)}
  for mode,(X,fnames) in feature_sets.items():
   for typ in ['LR','RF','XGB']:
    model=estimator(typ,y[train]);model.fit(X[train],y[train]);prob=float(model.predict_proba(X[i:i+1])[0,1]);key=mode+'_'+typ;out[key].append(prob)
    if subset=='all14' and mode=='full' and typ=='RF':
     native_imp=model[-1].feature_importances_;pos=len(en)
     for j,d in enumerate(meta):
      for rid in d['indices']:importance.append({'heldout':names[i],'aligned_index':rid+1,'kind':d['kind'],'importance_allocated':float(native_imp[pos+j]/len(d['indices']))})
     # Held-out marginal replacement, uncalibrated score log-loss difference, no refit.
     ixres=[j for j,d in enumerate(meta) if d['kind']=='residue'];xp=np.repeat(X[i:i+1],len(ixres),axis=0)
     for k,j in enumerate(ixres):xp[k,pos+j]=np.nanmedian(X[train,pos+j])
     if len(ixres):
      probs=model.predict_proba(xp)[:,1];base_loss=-np.log(np.clip(prob if y[i] else 1-prob,1e-8,1));loss=-np.log(np.clip(probs if y[i] else 1-probs,1e-8,1))
      for k,j in enumerate(ixres):lossimportance.append({'heldout':names[i],'aligned_index':meta[j]['indices'][0]+1,'delta_heldout_logloss':float(loss[k]-base_loss)})
  model=estimator('LR',y[train]).fit(Q[train],y[train]);out['formula4_LR'].append(float(model.predict_proba(Q[i:i+1])[0,1]))
  print(subset,'fold',names[i],flush=True)
 for model,prob in out.items():
  metrics.append({'subset':subset,'model':model,**metric(y[indices],np.array(prob))})
  for i,v in zip(indices,prob):prediction.append({'subset':subset,'model':model,'construct':names[i],'actual':int(y[i]),'score':v,'predicted':int(v>=.5)})
pd.DataFrame(metrics).to_csv(O+'/LOCO_model_metrics.csv',index=False);pd.DataFrame(prediction).to_csv(O+'/LOCO_predictions.csv',index=False);json.dump(selection,open(O+'/fold_feature_selection.json','w'),indent=2)
imp=pd.DataFrame(importance);imp.to_csv(O+'/fold_RF_importance.csv',index=False);loss=pd.DataFrame(lossimportance);loss.to_csv(O+'/fold_residue_replacement_importance.csv',index=False)
ri=imp.groupby(['aligned_index','kind']).importance_allocated.sum().unstack(fill_value=0)/14;ri.columns=['RF_'+c+'_allocated_importance' for c in ri.columns];ri['selected_folds_direct_residue']=imp[imp.kind=='residue'].groupby('aligned_index').heldout.nunique();ri['mean_heldout_logloss_increase']=loss.groupby('aligned_index').delta_heldout_logloss.sum()/14;ri=ri.fillna(0).sort_values('RF_residue_allocated_importance',ascending=False);ri.to_csv(O+'/residue_model_importance.csv')
# Freeze final reference-cohort models, and export an explicit simple scoring formula.
Xsel,selnames,meta,regions=build_selected(arr,y,np.ones(len(y),bool));X=np.c_[E,Xsel];fnames=en+selnames
for typ in ['LR','RF','XGB']:
 model=estimator(typ,y).fit(X,y);joblib.dump(model,O+'/final_full_'+typ+'.joblib')
model=estimator('LR',y).fit(Q,y);joblib.dump(model,O+'/final_formula4_LR.joblib');sc=model[1];lr=model[-1];coef=lr.coef_[0];rawcoef=coef/sc.scale_;rawinter=float(lr.intercept_[0]-np.sum(coef*sc.mean_/sc.scale_));formula={'feature_names':formula_cols,'standardized_intercept':float(lr.intercept_[0]),'mean':sc.mean_.tolist(),'scale':sc.scale_.tolist(),'standardized_coefficients':coef.tolist(),'raw_intercept':rawinter,'raw_coefficients':rawcoef.tolist(),'threshold':.5,'interpretation':'Exploratory ranking score, not a calibrated probability; frozen on 14 constructs.'};json.dump(formula,open(O+'/scoring_formula.json','w'),indent=2)
json.dump({'feature_names':fnames,'base_features':en,'selected_feature_names':selnames,'selected_feature_metadata':meta,'CCRs':regions,'model_parameters':{'LR':{'C':.1,'class_weight':'balanced'},'RF':{'n_estimators':300,'max_depth':2,'min_samples_leaf':2},'XGB':{'n_estimators':100,'max_depth':1,'learning_rate':.05,'reg_lambda':10,'reg_alpha':1}},'training_names':names},open(O+'/final_model_manifest.json','w'),indent=2)
print(pd.DataFrame(metrics).to_string(index=False));print(json.dumps(formula,indent=2))
