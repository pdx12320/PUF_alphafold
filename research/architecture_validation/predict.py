import argparse,json,joblib
from pathlib import Path
import pandas as pd,numpy as np
p=argparse.ArgumentParser();p.add_argument('--features',required=True);p.add_argument('--out',default='PUF_scores.csv');args=p.parse_args();root=Path(__file__).resolve().parent;meta=json.load(open(root/'recommended_model.json'));a=pd.read_csv(args.features);cols=meta['features'];X=a[cols].to_numpy(float);assert np.isfinite(X).all(),'Missing/nonfinite features';m=joblib.load(root/meta['model']);a['RF_score']=m.predict_proba(X)[:,1];a['predicted_work']=a.RF_score>=meta['threshold'];a['out_of_training_range']=((X<np.array([meta['training_min'][c] for c in cols]))|(X>np.array([meta['training_max'][c] for c in cols]))).any(1);a.to_csv(args.out,index=False)
