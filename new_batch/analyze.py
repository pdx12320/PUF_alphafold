import os,sys,re,json,zipfile,glob,io,hashlib
import numpy as np,pandas as pd,joblib
from Bio.PDB import MMCIFParser
sys.path.insert(0,os.path.abspath('extended'))
from score_new import get_features
O='new_batch';rows=[];requests={};seen=set();parser=MMCIFParser(QUIET=True)
oldseq=json.load(open('previous/results/sequences.json'));old=pd.read_csv('extended/all_construct_features.csv').set_index('construct');manifest=json.load(open('extended/summary_model_manifest.json'));formula=json.load(open('extended/scoring_formula.json'))
for p in sorted(glob.glob('upload/folds_2026_09_07_02_*.zip')):
 z=zipfile.ZipFile(p)
 for member in z.namelist():
  if not re.search(r'_full_data_\d+\.json$',member):continue
  folder=member.split('/')[0];name,seed=re.match(r'(.*)_rna_.*_seed(\d+)$',folder).groups();seed=int(seed);model=int(re.search(r'_full_data_(\d+)',member)[1]);key=(name,seed,model);assert key not in seen;seen.add(key)
  req=json.loads(z.read(folder+'/fold_'+folder+'_job_request.json'))[0];sq=req['sequences'][0]['proteinChain']['sequence'];rna=req['sequences'][1]['rnaSequence']['sequence']
  if name in requests:assert requests[name]['protein']==sq and requests[name]['rna']==rna
  else:requests[name]={'protein':sq,'rna':rna,'repeat_count_TRM':len(re.findall('[SCN][YR][FV][IV][EQR]',sq)),'loop_positions':[m.start()+1 for m in re.finditer('MNDGPHS',sq)],'exact_old_matches':[n for n,s in oldseq.items() if s['protein']==sq and s['rna']==rna]}
  d=json.loads(z.read(member));chains=np.array(d['token_chain_ids']);a=np.where(chains=='A')[0];b=np.where(chains=='B')[0];full=np.array(d['contact_probs'],dtype=np.float32);assert len(a)==len(sq) and len(b)==len(rna);cp=full[np.ix_(a,a)];pr=full[np.ix_(a,b)];pae=np.array(d['pae'],dtype=np.float32)[np.ix_(a,a)];st=parser.get_structure('x',io.StringIO(z.read(member.replace('full_data','model').replace('.json','.cif')).decode()));res=list(st[0]['A'].get_residues());assert len(res)==len(a);coords=np.array([r['CA'].coord for r in res]);plddt=np.array([np.mean([atom.bfactor for atom in r]) for r in res]);f=get_features(cp,pae,coords,plddt,pr);conf=json.loads(z.read(member.replace('full_data','summary_confidences')))
  rows.append({'construct':name,'design_id':int(name.split('_')[0]),'seed':seed,'model':model,'length':len(sq),'repeat_count':requests[name]['repeat_count_TRM'],'loop_count':len(requests[name]['loop_positions']),'iptm':conf['iptm'],'ptm':conf['ptm'],'plddt_p10':np.quantile(plddt,.1),**f})
  if len(rows)%60==0:print('extracted',len(rows),flush=True)
raw=pd.DataFrame(rows);raw.to_csv(O+'/model_features.csv',index=False);seed=raw.groupby(['construct','seed']).mean(numeric_only=True).drop(columns='model').reset_index();agg=seed.groupby('construct').mean(numeric_only=True).drop(columns='seed').reset_index();cols=manifest['columns']
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
