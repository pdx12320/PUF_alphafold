import os,json,re
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,confusion_matrix
from xgboost import XGBClassifier
O=os.path.dirname(__file__);a=pd.read_csv(O+'/all_construct_features.csv');names=a.construct.tolist();y=a.success.values;cols=[c for c in a if c.startswith(('pr_','pp_'))];X=a[cols].values;families=json.load(open(O+'/architecture_groups.json'));groups=np.array([families[n] for n in names]);rows=[];pred=[]
def est(typ,yy):
 if typ=='LR':model=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=2026)
 elif typ=='RF':model=RandomForestClassifier(n_estimators=300,max_depth=2,min_samples_leaf=2,max_features='sqrt',class_weight='balanced',random_state=2026,n_jobs=1)
 else:model=XGBClassifier(n_estimators=100,max_depth=1,learning_rate=.05,min_child_weight=1,reg_lambda=10,reg_alpha=1,subsample=1,colsample_bytree=.8,scale_pos_weight=(len(yy)-sum(yy))/sum(yy),objective='binary:logistic',eval_metric='logloss',tree_method='hist',random_state=2026,n_jobs=1)
 return make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),model)
for typ in ['LR','RF','XGB']:
 p=np.empty(len(y))
 for group in np.unique(groups):
  test=groups==group;train=~test;assert len(set(y[train]))==2;p[test]=est(typ,y[train]).fit(X[train],y[train]).predict_proba(X[test])[:,1]
 rows.append({'model':'summary_'+typ,'validation':'leave_source_repeat_order_out','AUC':roc_auc_score(y,p),'balanced_accuracy':balanced_accuracy_score(y,p>=.5),'TN_FP_FN_TP':str(confusion_matrix(y,p>=.5).ravel().tolist())})
 for n,actual,score in zip(names,y,p):pred.append({'construct':n,'model':typ,'actual':actual,'score':score})
 joblib.dump(est(typ,y).fit(X,y),O+'/final_summary_'+typ+'.joblib')
json.dump({'columns':cols,'training_min':X.min(0).tolist(),'training_max':X.max(0).tolist()},open(O+'/summary_model_manifest.json','w'),indent=2);pd.DataFrame(rows).to_csv(O+'/architecture_holdout_metrics.csv',index=False);pd.DataFrame(pred).to_csv(O+'/architecture_holdout_predictions.csv',index=False)
# Descriptive seed robustness and source-matched associations.
z=np.load(O+'/seed_arrays.npz');ccr=json.load(open(O+'/CCR_definitions.json'));rstats=pd.read_csv(O+'/residue_associations.csv');imp=pd.read_csv(O+'/residue_model_importance.csv');rr=rstats.merge(imp,on='aligned_index',how='left').fillna({'selected_folds_direct_residue':0,'RF_residue_allocated_importance':0});rr=rr.merge(pd.read_csv(O+'/within_architecture_residue_correlation.csv'),on='aligned_index');rr.sort_values('RF_residue_allocated_importance',ascending=False).to_csv(O+'/top_residue_combined.csv',index=False)
seedrecords=[]
for reg in ccr:
 inds=reg['indices'];v=z['density'][:,inds].mean(1)
 for n,seed,value in zip(z['names'],z['seeds'],v):seedrecords.append({'CCR':reg['CCR'],'construct':n,'seed':seed,'density':float(value),'success':int(y[names.index(n)])})
pd.DataFrame(seedrecords).to_csv(O+'/CCR_seed_values.csv',index=False)
# Exact native positions for contact candidates/CCRs in every construct, no gap-to-zero conversion.
mp=json.load(open(O+'/core_mapping.json'));seqs=json.load(open('previous/results/sequences.json'));contacts=pd.read_csv(O+'/contact_tests_all14.csv').head(20);maprows=[]
for _,c in contacts.iterrows():
 for n in names:
  i,j=int(c.aligned_i)-1,int(c.aligned_j)-1;ni,nj=mp[n][i],mp[n][j];maprows.append({'construct':n,'aligned_i':i+1,'aligned_j':j+1,'native_i':ni+1,'native_j':nj+1,'aa_i':seqs[n]['protein'][ni],'aa_j':seqs[n]['protein'][nj]})
pd.DataFrame(maprows).to_csv(O+'/top_contact_native_mapping.csv',index=False)
looprows=[]
for n in names:
 m=mp[n];seq=seqs[n]['protein']
 for match in re.finditer('MNDGPHS',seq):
  before=[i for i,j in enumerate(m) if j>=0 and j<match.start()];after=[i for i,j in enumerate(m) if j>=match.end()];looprows.append({'construct':n,'motif':'MNDGPHS','native_start':match.start()+1,'native_end':match.end(),'left_repeat_slot':before[-1]//36+1 if before else None,'right_repeat_slot':after[0]//36+1 if after else None,'annotation':'sequence motif; native and engineered occurrences retained'})
pd.DataFrame(looprows).to_csv(O+'/loop_mapping.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
