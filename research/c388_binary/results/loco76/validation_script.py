from pathlib import Path
import json,gzip,hashlib
import pandas as pd,numpy as np
from sklearn.metrics import accuracy_score,balanced_accuracy_score,f1_score,confusion_matrix
p=Path('outputs/c388_binary_loco76');proto=json.loads((p/'protocol.json').read_text());base=proto['base_protocol'];d=pd.read_csv('work/source_archive/PUF_alphafold/retrain_classification_v4_dynamic_fullCP_20260921/data/construct_registry.csv').set_index('construct_id');pred=pd.read_csv(p/'predictions.csv');metrics=json.loads((p/'metrics.json').read_text());folds=list(p.glob('fold_*_search.json.gz'));held=[];count=0;invalid=0
for path in folds:
 r=json.loads(gzip.decompress(path.read_bytes()));tr=set(r['train']);te=set(r['test']);held.extend(te)
 assert tr|te==set(base['train']) and not tr&te and not (tr|te)&set(base['test'])
 assert not set(d.loc[list(tr)].sequence_identity)&set(d.loc[list(te)].sequence_identity)
 for inner in r['inner']:
  for f in inner['folds']:
   a=set(f['train']);b=set(f['validation']);assert a|b==tr and not a&b
   assert not set(d.loc[list(a)].sequence_identity)&set(d.loc[list(b)].sequence_identity)
 seal=path.with_name(path.name.replace('_search.json.gz','_sealed.json'))
 assert hashlib.sha256(seal.read_bytes()).hexdigest()==r['sealed_sha256']
 for s in json.loads(seal.read_text()):
  q=pred[(pred.construct_id==s['construct_id'])&(pred.version==s['version'])].iloc[0]
  assert int(q.prediction)==s['prediction'] and q.threshold_pp==s['threshold_pp']
 count+=len(r['records']);invalid+=sum(x['status']!='ok' for x in r['records'])
assert len(held)==76 and len(set(held))==76 and set(held)==set(base['train']) and len(folds)==76
for version,g in pred.groupby('version'):
 assert len(g)==76 and g.construct_id.nunique()==76
 np.testing.assert_array_equal(g.true_label,(g.delta_pp>=g.threshold_pp).astype(int))
 m=next(x for x in metrics['results'] if x['version']==version)
 assert np.isclose(accuracy_score(g.true_label,g.prediction),m['accuracy'])
 assert np.isclose(balanced_accuracy_score(g.true_label,g.prediction),m['balanced_accuracy'])
 assert np.isclose(f1_score(g.true_label,g.prediction,average='macro'),m['macro_F1'])
 assert confusion_matrix(g.true_label,g.prediction,labels=[0,1]).tolist()==m['confusion_matrix']
with (p/'validation.json').open('x') as h:json.dump(dict(status='passed',outer_folds=76,unique_primary_predictions=76,fixed_five_excluded_everywhere=True,inner_identity_disjoint=True,sealed_predictions_verified=True,metrics_recomputed=True,candidates_evaluated=count,invalid_candidates=invalid),h,indent=2)
print('All fold isolation, prediction seal and metric checks passed',count,invalid)
