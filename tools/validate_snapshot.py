import hashlib,json
from pathlib import Path
import pandas as pd,numpy as np
root=Path(__file__).resolve().parents[1];manifest=json.loads((root/'docs/file_manifest.json').read_text());errors=[]
for row in manifest:
 p=root/row['path']
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append(row['path'])
if errors:raise SystemExit('Missing/changed snapshot files: '+', '.join(errors))
a=pd.read_csv(root/'architecture_validation/main_features.csv');assert len(a)==24 and a.success.sum()==4 and a.construct.nunique()==24
p=pd.read_csv(root/'architecture_validation/heldout_predictions.csv');p=p[p.subset=='PUF12_main'];assert p.groupby(['validation','feature_set','model']).construct.nunique().eq(24).all()
sort=['feature_set','model','construct'];a=p[p.validation=='architecture_out'].sort_values(sort);b=p[p.validation=='source_order_out'].sort_values(sort);assert np.array_equal(a.score.to_numpy(),b.score.to_numpy())
print(f'PASS: {len(manifest)} file hashes; 24 unique PUF12; 4 successes; complete held-out predictions.')
