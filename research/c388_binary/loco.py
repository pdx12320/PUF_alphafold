"""Nested identity-LOCO on the frozen development cohort; fixed five stay excluded."""
import argparse,os,sys,json,gzip,hashlib,time,importlib.util
from pathlib import Path
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
p=argparse.ArgumentParser();p.add_argument('--v4',type=Path,required=True);p.add_argument('--protocol',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
spec=importlib.util.spec_from_file_location('binary_runner',Path(__file__).with_name('train.py'));binary=importlib.util.module_from_spec(spec);spec.loader.exec_module(binary)
protocol=json.loads(args.protocol.read_text())
for rel,h in protocol['source_sha256'].items():assert binary.digest(args.v4/rel)==h
for rel,h in protocol['input_sha256'].items():assert binary.digest(args.v4/rel)==h
sys.path.insert(0,str(args.v4/'src'))
import core
binary.install_binary(core)
import train
import numpy as np,pandas as pd
from sklearn.metrics import accuracy_score,balanced_accuracy_score,f1_score,matthews_corrcoef,confusion_matrix
args.output.mkdir(parents=True,exist_ok=False)
ds=core.Dataset('C388','combined_PR');original=ds.d.copy(deep=True)
ids=np.flatnonzero(original.construct_id.isin(protocol['train']));fixed=np.flatnonzero(original.construct_id.isin(protocol['test']))
assert len(ids)==76 and len(fixed)==5
reps=protocol['representations'];models=protocol['refinement_models'];thresholds=protocol['candidates_a_pp']
groups=sorted(original.iloc[ids].sequence_identity.unique());rows=[];selections=[];started=time.monotonic()
binary.save(args.output,'protocol.json',dict(base_protocol_sha256=binary.digest(args.protocol),base_protocol=protocol,scope='76 development constructs only; five fixed tests excluded throughout',outer='leave one sequence identity out',inner='training-only stratified grouped folds, all thresholds/features/models reselected',pooled_AUC='not reported: different models and label thresholds across folds',runner_sha256=binary.digest(__file__),n_outer_folds=len(groups)))
for fold,group in enumerate(groups):
 te=ids[original.iloc[ids].sequence_identity.to_numpy()==group];tr=ids[original.iloc[ids].sequence_identity.to_numpy()!=group]
 ds.d=original.copy(deep=True);ds.d.loc[np.r_[te,fixed],[c for c in ds.d if c.endswith(('_fraction','_delta_pp'))]]=np.nan
 assert not set(original.iloc[tr].sequence_identity)&set(original.iloc[np.r_[te,fixed]].sequence_identity)
 seed=core.seed('C388-binary-nested-LOCO',20260923,group);records=[];audits=[];feasible={}
 for a in thresholds:
  y=binary.binary_labels(ds.d,'C388',a)
  if np.bincount(y[tr],minlength=2).min()<core.CFG['min_outer_class']['C388']:continue
  inner,kind=core.inner_splits(ds.d,tr,y,'LOCO',seed)
  if not inner:continue
  feasible[a]=inner
  for aa,bb,_ in inner:
   assert set(aa)|set(bb)==set(tr) and not set(aa)&set(bb)
   assert not set(original.iloc[aa].sequence_identity)&set(original.iloc[bb].sequence_identity)
  audits.append(dict(a=a,kind=kind,folds=[dict(train=original.iloc[aa].construct_id.tolist(),validation=original.iloc[bb].construct_id.tolist()) for aa,bb,_ in inner]))
  for c in reps:
   r=dict(a=a,b=0,representation=c,family=binary.family(c),model=protocol['screening_model'],stage='screen',status='invalid')
   try:
    mm,counts=train.evaluate(ds,c,r['model'],y,inner,seed);r.update(status='ok',metrics=mm,actual_counts=counts)
   except (ValueError,np.linalg.LinAlgError) as e:r['reason']=str(e)
   records.append(r)
  train.clear_fitted(ds)
 rank=lambda r:train.rank(r,'C388');short=[]
 for fam in sorted({binary.family(c) for c in reps}):short+=sorted([r for r in records if r['status']=='ok' and r['family']==fam],key=rank)[:2]
 for base in short:
  y=binary.binary_labels(ds.d,'C388',base['a'])
  for mc in models:
   r=dict(a=base['a'],b=0,representation=base['representation'],family=base['family'],model=mc,stage='refine',status='invalid')
   try:
    mm,counts=train.evaluate(ds,r['representation'],mc,y,feasible[r['a']],seed);r.update(status='ok',metrics=mm,actual_counts=counts)
   except (ValueError,np.linalg.LinAlgError) as e:r['reason']=str(e)
   records.append(r)
  train.clear_fitted(ds)
 good=[r for r in records if r['status']=='ok']
 if not good:raise RuntimeError(f'No estimable candidate in fold {fold}')
 best=min(good,key=rank);y=binary.binary_labels(ds.d,'C388',best['a']);majority=int(np.bincount(y[tr],minlength=2).argmax())
 # Matched-label family comparisons use the primary fold threshold.
 candidates={'primary':best}
 for fam in sorted({r['family'] for r in good}):
  pool=[r for r in good if r['family']==fam and r['a']==best['a']]
  if pool:candidates['matched_'+fam]=min(pool,key=rank)
 sealed=[]
 for version,r in candidates.items():
  m,st,z,emb=core.model_fit(ds,r['representation'],r['model'],tr,y,seed)
  if st['bad'][te].any():raise ValueError('Missing held-out coordinates')
  pred,scores,kind=core.model_scores(m,z[te])
  for i,q,ss in zip(te,pred,scores):sealed.append(dict(fold=fold,version=version,construct_id=original.iloc[i].construct_id,batch=original.iloc[i].batch,threshold_pp=-best['a'],prediction=int(q),predicted_class=binary.CLASSES[int(q)],scores=ss.tolist(),score_kind=kind,algorithm=r['model']['algorithm'],family=r['family'],actual_features=st['actual_CP_count'],majority_prediction=majority))
 binary.save(args.output,f'fold_{fold:03d}_sealed.json',sealed)
 # Only now evaluate this outer holdout.
 for r in sealed:
  actual=float(original.loc[original.construct_id==r['construct_id'],'C388_delta_pp'].iloc[0]);truth=int(actual>=r['threshold_pp']);rows.append(dict(**r,delta_pp=actual,true_label=truth,true_class=binary.CLASSES[truth],correct=r['prediction']==truth))
 binary.save(args.output,f'fold_{fold:03d}_search.json.gz',dict(train=original.iloc[tr].construct_id.tolist(),test=original.iloc[te].construct_id.tolist(),fixed_test_excluded=protocol['test'],seed=seed,inner=audits,records=records,selected=best,sealed_sha256=binary.digest(args.output/f'fold_{fold:03d}_sealed.json')))
 selections.append(dict(fold=fold,threshold_pp=-best['a'],family=best['family'],representation=best['representation']['id'],model=json.dumps(best['model']),inner_macro_F1=best['metrics']['macro_F1']))
 train.clear_fitted(ds)
 print(f'Completed {fold+1}/{len(groups)} folds; elapsed {time.monotonic()-started:.1f}s',flush=True)
frame=pd.DataFrame(rows);frame.to_csv(args.output/'predictions.csv',index=False);pd.DataFrame(selections).to_csv(args.output/'selections.csv',index=False)
metrics=[]
for version,g in frame.groupby('version'):
 y=g.true_label.to_numpy();pred=g.prediction.to_numpy()
 for name,q in [(version,pred)]+([('matched_majority',g.majority_prediction.to_numpy())] if version=='primary' else []):
  metrics.append(dict(version=name,n=len(y),correct=int(sum(y==q)),accuracy=accuracy_score(y,q),balanced_accuracy=balanced_accuracy_score(y,q),macro_F1=f1_score(y,q,average='macro'),MCC=matthews_corrcoef(y,q),confusion_matrix=confusion_matrix(y,q,labels=[0,1]).tolist()))
binary.save(args.output,'metrics.json',dict(results=metrics,elapsed_seconds=time.monotonic()-started,threshold_frequency=pd.DataFrame(selections).threshold_pp.value_counts().to_dict(),interpretation='Adaptive-threshold LOCO on 76 development constructs. Fold-specific label boundaries; not a fixed-threshold endpoint.'))
print(json.dumps(metrics),flush=True)
