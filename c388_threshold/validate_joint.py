import pathlib
exec(pathlib.Path('c388_threshold/joint.py').read_text().split('u=np.unique(rate);thresholds=')[0])
from sklearn.model_selection import StratifiedKFold,GroupKFold
rows=[];preds=[]
for c in [.3,.35,.4,.5,.55]:
 y=(rate>=c).astype(int)
 for key,typ in [('CP3','LR'),('CP_Structure9','LR')]:
  X=a[man['feature_sets'][key]].to_numpy()
  for val in ['LOCO','mutation_group_out']:
   groups=a.construct.to_numpy() if val=='LOCO' else a.mutation_group.to_numpy();p=np.zeros(n);tt=np.zeros(n)
   for g in sorted(set(groups)):
    tr=np.where(groups!=g)[0];te=np.where(groups==g)[0]
    folds=StratifiedKFold(3,shuffle=True,random_state=2026).split(tr,y[tr]) if val=='LOCO' else GroupKFold(3).split(tr,y[tr],a.mutation_group.to_numpy()[tr]);ip=np.zeros(len(tr))
    for it,iv in folds:
     m=model(typ).fit(X[tr[it]],y[tr[it]]);ip[iv]=m.predict_proba(X[tr[iv]])[:,1]
    s=sweep(y[tr],ip);s['distance']=(s.score_threshold-.5).abs();t=s.sort_values(['balanced_accuracy','distance'],ascending=[False,True]).iloc[0].score_threshold
    m=model(typ).fit(X[tr],y[tr]);p[te]=m.predict_proba(X[te])[:,1];tt[te]=t
   d=sweep(y,(p>=tt).astype(float));d=d[d.score_threshold==.5].iloc[0].to_dict();rows.append(dict(C388_threshold=c,model=key+'_'+typ,validation=val,median_training_threshold=np.median(tt),min_training_threshold=min(tt),max_training_threshold=max(tt),AUC=roc_auc_score(y,p),**d))
   for i in range(n):preds.append(dict(C388_threshold=c,model=key+'_'+typ,validation=val,construct=a.construct.iloc[i],y=y[i],score=p[i],training_threshold=tt[i],prediction=int(p[i]>=tt[i])))
pd.DataFrame(rows).to_csv(O/'joint_inner_threshold_validation.csv',index=False);pd.DataFrame(preds).to_csv(O/'joint_inner_threshold_predictions.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
