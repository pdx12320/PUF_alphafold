from pathlib import Path
import json,gzip,hashlib,collections,sys
import numpy as np,pandas as pd,joblib
from sklearn.metrics import accuracy_score,balanced_accuracy_score,f1_score,confusion_matrix
root=Path.cwd();p=root/'outputs/c388_binary_fixed_five';v=root/'work/source_archive/PUF_alphafold/retrain_classification_v4_dynamic_fullCP_20260921'
proto=json.loads((p/'protocol.json').read_text());ev=json.loads((p/'evaluation.json').read_text());chosen=json.loads((p/'selected.json').read_text());splits=json.loads((p/'inner_splits.json').read_text());seal=json.loads((p/'predictions_sealed.json').read_text());d=pd.read_csv(v/'data/construct_registry.csv');m=joblib.load(p/'model.joblib')
checks={}
checks['prediction_seal_hash']=hashlib.sha256((p/'predictions_sealed.json').read_bytes()).hexdigest()==ev['sealed_sha256']
checks['input_hashes']=all(hashlib.sha256((v/f).read_bytes()).hexdigest()==h for f,h in proto['input_sha256'].items())
checks['runner_hash']=hashlib.sha256((root/'work/PUF_alphafold/research/c388_binary/train.py').read_bytes()).hexdigest()==proto['runner_sha256']
checks['source_hashes']=all(hashlib.sha256((v/f).read_bytes()).hexdigest()==h for f,h in proto['source_sha256'].items())
tr=d[d.construct_id.isin(proto['train'])];te=d[d.construct_id.isin(proto['test'])]
checks['split_counts']=len(tr)==76 and len(te)==5 and not set(tr.sequence_identity)&set(te.sequence_identity)
checks['inner_splits']=all(not set(f['train'])&set(f['validation']) and set(f['train'])|set(f['validation'])==set(proto['train']) and not set(d[d.construct_id.isin(f['train'])].sequence_identity)&set(d[d.construct_id.isin(f['validation'])].sequence_identity) for s in splits if s['status']=='ok' for f in s['folds'])
# Reconstruct inference using only saved artifact plus canonical CP input.
ids=m['state']['selected'];pick=m['state']['pick'];np_pp=121278
raw=np.load(v/'cp_features/C388_PR.npy');x=raw[:,ids-np_pp].astype(float)
if m['representation']['mode']=='delta':x=x-m['reference']
z=(x-m['state']['mean'])/m['state']['scale']
if m['embedding']:z=z@m['basis']
idx=d.index[d.construct_id.isin(proto['test'])].to_numpy();pred=m['classifier'].predict(z[idx]);margin=m['classifier'].decision_function(z[idx]);scores=np.c_[-margin,margin]
checks['artifact_prediction_roundtrip']=all(m['classes'][int(q)]==s['predicted_class'] for q,s in zip(pred,seal)) and np.allclose(scores,np.array([s['scores'] for s in seal]),atol=1e-10)
y=(te.C388_delta_pp.to_numpy()>=m['threshold_pp']).astype(int)
checks['metrics_recomputed']=np.isclose(accuracy_score(y,pred),ev['metrics']['accuracy']) and np.isclose(balanced_accuracy_score(y,pred),ev['metrics']['balanced_accuracy']) and np.isclose(f1_score(y,pred,average='macro'),ev['metrics']['macro_F1'])
checks['feature_dictionary_count']=len(pd.read_csv(p/'selected_CP_pairs.csv'))==m['state']['actual_CP_count']==876
records=json.loads(gzip.decompress((p/'search.json.gz').read_bytes()));ok=[r for r in records if r['status']=='ok'];simp={'Ridge':0,'LR':1,'LinearSVM':2,'LDA':3,'RBF':4,'RF':5,'ET':5}
def rank(r):return (-r['metrics']['macro_F1'],-r['metrics']['MCC'],simp[r['model']['algorithm']],np.mean(r['actual_counts']),r['representation']['id'],r['a'],r['b'])
checks['selection_matches_frozen_rank']=min(ok,key=rank)==chosen
checks['search_budget']=len(records)==388 and len(ok)==384
checks['high_variance_equivalent_to_all']=all(r['actual_counts']==chosen['actual_counts'] and r['metrics']==chosen['metrics'] for r in ok if r['a']==5 and r['stage']=='screen' and r['representation']['family']=='full_PR' and r['representation']['retain'] in [1000,'all'])
assert all(checks.values()),checks
summary=dict(checks=checks,n_checks=len(checks),actual_selected_features=m['state']['actual_CP_count'],eligible_features=len(m['state']['eligible']),exact_embedding_rank=m['state']['train_rank'],training_class_counts={'decrease':int(sum(tr.C388_delta_pp < m['threshold_pp'])),'unchanged':int(sum(tr.C388_delta_pp >= m['threshold_pp']))},training_batches=tr.batch.value_counts().to_dict(),test_batches=te.batch.value_counts().to_dict(),confusion_matrix=confusion_matrix(y,pred).tolist(),input_registry_sha256=hashlib.sha256((v/'data/construct_registry.csv').read_bytes()).hexdigest())
checks.update({k:bool(v) for k,v in checks.items()})
with (p/'validation_verified.json').open('x') as h:json.dump(summary,h,indent=2)
def summarize(r):return dict(family=r['family'],algorithm=r['model']['algorithm'],model=json.dumps(r['model']),threshold_pp=-r['a'],representation=r['representation']['id'],actual_inner_features=json.dumps(r['actual_counts']),inner_accuracy=r['metrics']['accuracy'],inner_macro_F1=r['metrics']['macro_F1'],inner_balanced_accuracy=r['metrics']['balanced_accuracy'])
families=[summarize(min([r for r in ok if r['family']==f],key=rank)) for f in sorted({r['family'] for r in ok})]
algorithms=[summarize(min([r for r in ok if r['model']['algorithm']==a],key=rank)) for a in sorted({r['model']['algorithm'] for r in ok})]
pd.DataFrame(families).to_csv(p/'family_comparison_internal.csv',index=False)
pd.DataFrame(algorithms).to_csv(p/'model_comparison_internal.csv',index=False)
with (p/'run.log').open('x') as h:h.write((root/'work/c388-full.log').read_text())
print(json.dumps(summary,indent=2));print(pd.DataFrame(algorithms)[['algorithm','threshold_pp','inner_macro_F1','inner_balanced_accuracy']].to_string(index=False))
