#!/usr/bin/env python3
"""Joint five-construct test using the archived v4 training-only search.

Usage: python retrain_model2.py --v4 /path/to/retrain_classification_v4_dynamic_fullCP_20260921 --endpoint C295
The source archive is identified in ../results_20260922/source_archives.json.
"""
from pathlib import Path
import argparse,os,sys,json,hashlib,time,itertools,gzip
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
p=argparse.ArgumentParser();p.add_argument('--v4',type=Path,required=True);p.add_argument('--endpoint',choices=['C295','C388','C871'],required=True);args=p.parse_args()
sys.path.insert(0,str(args.v4/'src'))
import core,train
import numpy as np
import pandas as pd
from joblib import dump
OUT=Path(__file__).resolve().parent/'model2';OUT.mkdir(exist_ok=True)
TASK=args.endpoint
TARGETS={'P9-GNS':'25_p9_gns','P9-NPS':'30_p9_nps','P9-NTQ':'32_p9_ntq','P8-GVE':'p8_r6_gve','P4-R5-SNE+P7-R5-SNE':'p4_r5_sne_plus_p7_r5_sne'}
def save(name,obj):
 text=json.dumps(obj,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v))
 if name.endswith('_search.json'):(OUT/(name+'.gz')).write_bytes(gzip.compress(text.encode(),mtime=0))
 else:(OUT/name).write_text(text)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
ds=core.Dataset(TASK,'combined_PR');d=ds.d
held_sequences=set(d.loc[d.construct_id.isin(TARGETS.values()),'sequence_identity'])
te=np.flatnonzero(d.sequence_identity.isin(held_sequences));tr=np.flatnonzero(~d.sequence_identity.isin(held_sequences))
assert len(te)==5 and len(tr)==76 and set(d.iloc[te].construct_id)==set(TARGETS.values())
assert not set(d.iloc[tr].sequence_identity)&held_sequences
# Test labels are removed from the object seen by all training/search functions.
observed=d.iloc[te].copy()
for col in [c for c in d.columns if c.endswith(('_delta_pp','_fraction'))]:d.loc[te,col]=np.nan
s=core.seed('fixed-five-test',20260923,TASK,'combined_PR')
protocol=dict(endpoint=TASK,seed=s,cohort='combined_PR',train=d.iloc[tr].construct_id.tolist(),test=d.iloc[te].construct_id.tolist(),held_sequence_identities=sorted(held_sequences),selection='all five excluded together, including identical sequences; training-only inner CV',test_labels_masked=True,source_sha256={f: digest(args.v4/f) for f in ['src/core.py','src/train.py','config/protocol.json']},config=core.CFG)
save(TASK+'_protocol.json',protocol)
records=[];audits=[];feasible={};K=len(core.names(TASK))
grids=list(itertools.product(core.CFG['grids'][TASK],repeat=2)) if TASK=='C388' else [(a,0) for a in core.CFG['grids'][TASK]]
for a,b in grids:
 y=core.labels(d,TASK,a,b)
 if np.bincount(y[tr],minlength=K).min()<core.CFG['min_outer_class'][TASK]:continue
 inner,kind=core.inner_splits(d,tr,y,'LOCO',s)
 if not inner:continue
 for aa,bb,_ in inner:assert not set(aa)&set(te) and not set(bb)&set(te)
 feasible[a,b]=inner;audits.append(dict(a=a,b=b,kind=kind,folds=[dict(train=d.iloc[aa].construct_id.tolist(),validation=d.iloc[bb].construct_id.tolist()) for aa,bb,_ in inner]))
 for c in core.CFG['representations']:
  if not ds.available(c):continue
  rec=dict(a=a,b=b,representation=c,model=core.CFG['screening_classifier'],stage='screen',status='invalid')
  try:mm,counts=train.evaluate(ds,c,rec['model'],y,inner,s);rec.update(status='ok',metrics=mm,actual_counts=counts)
  except (ValueError,np.linalg.LinAlgError,MemoryError) as e:rec['reason']=str(e)
  records.append(rec)
 train.clear_fitted(ds);print(TASK,'screened',a,b,'candidates',len(records),flush=True)
shortlist=sorted([r for r in records if r['status']=='ok'],key=lambda r:train.rank(r,TASK))[:core.CFG['max_refine']]
for base in shortlist:
 y=core.labels(d,TASK,base['a'],base['b'])
 for mc in core.CFG['refinement_models']:
  rec=dict(a=base['a'],b=base['b'],representation=base['representation'],model=mc,stage='refine',status='invalid')
  try:mm,counts=train.evaluate(ds,rec['representation'],mc,y,feasible[rec['a'],rec['b']],s);rec.update(status='ok',metrics=mm,actual_counts=counts)
  except (ValueError,np.linalg.LinAlgError,MemoryError) as e:rec['reason']=str(e)
  records.append(rec)
 train.clear_fitted(ds)
save(TASK+'_inner_splits.json',audits)
save(TASK+'_search.json',records)
best=min([r for r in records if r['status']=='ok'],key=lambda r:train.rank(r,TASK));y=core.labels(d,TASK,best['a'],best['b'])
model,st,z,embedding=core.model_fit(ds,best['representation'],best['model'],tr,y,s)
assert not st['bad'][te].any()
pred,scores,kind=core.model_scores(model,z[te])
save(TASK+'_selected.json',best)
# Seal predictions before reintroducing any test outcome.
sealed=[dict(construct_id=d.iloc[i].construct_id,predicted_class=core.names(TASK)[int(q)],scores=sc.tolist(),score_kind=kind) for i,q,sc in zip(te,pred,scores)]
save(TASK+'_predictions_sealed.json',sealed)
modelpath=OUT/(TASK+'_model.joblib');dump(dict(classifier=model,state=st,representation=best['representation'],embedding=embedding,train_ids=protocol['train'],test_ids=protocol['test']),modelpath,compress=3)
actual=core.labels(observed,TASK,best['a'],best['b']);reverse={v:k for k,v in TARGETS.items()};rows=[]
for j,(i,q,obs) in enumerate(zip(te,pred,actual)):
 rows.append(dict(target=reverse[d.iloc[i].construct_id],construct_id=d.iloc[i].construct_id,endpoint=TASK,true_class=core.names(TASK)[int(obs)],predicted_class=core.names(TASK)[int(q)],correct=bool(q==obs),decrease_threshold_pp=best['a'],increase_threshold_pp=best['b'],measured_delta_pp=float(observed.iloc[j][TASK+'_delta_pp']),feature_family=best['representation']['family'],algorithm=best['model']['algorithm'],n_train=len(tr),n_test=len(te),validation='joint_fixed_five_test'))
pd.DataFrame(rows).to_csv(OUT/(TASK+'_predictions.csv'),index=False)
majority=int(np.bincount(y[tr],minlength=K).argmax())
save(TASK+'_evaluation.json',dict(correct=int(sum(pred==actual)),n_test=len(te),matched_majority_correct=int(sum(actual==majority)),majority_class=core.names(TASK)[majority],predictions_sha256=digest(OUT/(TASK+'_predictions_sealed.json')),model_sha256=digest(modelpath),evaluated_after_prediction_seal=True))
print(TASK,'DONE',int(sum(pred==actual)),'/5',best['a'],best['b'],best['model'],flush=True)
