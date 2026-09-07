import os,json,re,zipfile,glob
import pandas as pd,numpy as np,joblib
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
O='new_batch';oldseq=json.load(open('previous/results/sequences.json'));old=pd.read_csv('extended/all_construct_features.csv').set_index('construct');manifest=json.load(open('extended/summary_model_manifest.json'));formula=json.load(open('extended/scoring_formula.json'));cols=manifest['columns']
# Reconstruct only unreadable RF artifact using original 14 training constructs and unchanged parameters.
rf=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),RandomForestClassifier(n_estimators=300,max_depth=2,min_samples_leaf=2,max_features='sqrt',class_weight='balanced',random_state=2026,n_jobs=1)).fit(old[cols].values,old.success.values);joblib.dump(rf,O+'/restored_frozen_RF.joblib',compress=3);joblib.load(O+'/restored_frozen_RF.joblib')
requests={}
for p in glob.glob('upload/*.zip'):
 z=zipfile.ZipFile(p)
 for n in z.namelist():
  if not n.endswith('job_request.json'):continue
  req=json.loads(z.read(n))[0];name=re.sub('_rna_.*','',n.split('/')[0]);sq=req['sequences'][0]['proteinChain']['sequence'];rna=req['sequences'][1]['rnaSequence']['sequence'];requests[name]={'protein':sq,'rna':rna,'repeat_count_TRM':len(re.findall('[SCN][YR][FV][IV][EQR]',sq)),'loop_positions':[m.start()+1 for m in re.finditer('MNDGPHS',sq)],'exact_old_matches':[k for k,v in oldseq.items() if v['protein']==sq and v['rna']==rna]}
raw=pd.read_csv(O+'/model_features.csv')
seed=raw.groupby(['construct','seed']).mean(numeric_only=True).drop(columns='model').reset_index();agg=seed.groupby('construct').mean(numeric_only=True).drop(columns='seed').reset_index();cols=manifest['columns']
for label,frame in [('seed',seed),('construct',agg)]:
 for typ in ['XGB','RF','LR']:
  model=(joblib.load(O+'/restored_frozen_RF.joblib') if typ=='RF' else joblib.load('extended/final_summary_'+typ+'.joblib'));frame[typ+'_score']=model.predict_proba(frame[cols].values)[:,1]
 model=joblib.load('extended/final_formula4_LR.joblib');frame['formula4_score']=model.predict_proba(frame[formula['feature_names']].values)[:,1]
 frame.to_csv(O+'/'+label+'_scores.csv',index=False)
# Numeric range audit is descriptive; novelty also assessed from architecture.
for i,row in agg.iterrows():
 n=row.construct;seq=requests[n]['protein'];nearest=max(oldseq,key=lambda k:__import__('difflib').SequenceMatcher(None,seq,oldseq[k]['protein'],autojunk=False).ratio());ratio=__import__('difflib').SequenceMatcher(None,seq,oldseq[nearest]['protein'],autojunk=False).ratio();outside=[c for c in cols if row[c]<old[c].min()-1e-6 or row[c]>old[c].max()+1e-6];ss=seed[seed.construct==n];agg.loc[i,'nearest_old']=nearest;agg.loc[i,'sequence_similarity_ratio']=ratio;agg.loc[i,'exact_old_match']=';'.join(requests[n]['exact_old_matches']);agg.loc[i,'n_seeds']=len(ss);agg.loc[i,'XGB_seed_positive_count']=int((ss.XGB_score>=.5).sum());agg.loc[i,'XGB_seed_min']=ss.XGB_score.min();agg.loc[i,'XGB_seed_max']=ss.XGB_score.max();agg.loc[i,'S12_seed_sd']=ss.pp_nonlocal12_high_per_res.std();agg.loc[i,'outside_summary_ranges']=';'.join(outside);agg.loc[i,'outside_summary_range_count']=len(outside)
agg=agg.sort_values(['XGB_score','pp_nonlocal12_high_per_res'],ascending=False);agg.to_csv(O+'/ranked_results.csv',index=False);json.dump(requests,open(O+'/sequences_and_audit.json','w'),indent=2)
print(agg[['design_id','construct','repeat_count','length','pp_nonlocal12_high_per_res','XGB_score','RF_score','LR_score','formula4_score','XGB_seed_positive_count','iptm','outside_summary_range_count']].to_string(index=False));print('models',len(raw),'seeds',len(seed),'constructs',len(agg),'missing ID6',6 not in agg.design_id.values,flush=True)
