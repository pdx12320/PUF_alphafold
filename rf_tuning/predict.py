"""Input: construct-level feature CSV matching the original extraction definitions."""
import argparse,json,joblib
from pathlib import Path
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--features',required=True);p.add_argument('--scope',default='joint',choices=['joint','CP_density3','CP_summary14','CP_structure17']);p.add_argument('--threshold',default='fixed',choices=['fixed','inner']);p.add_argument('--out',default='predictions.csv');args=p.parse_args()
root=Path(__file__).resolve().parent;spec=next(s for s in json.load(open(root/'final_model_manifest.json')) if s['scope']==args.scope);a=pd.read_csv(args.features);X=a[spec['features']].to_numpy();assert np.isfinite(X).all(),'Features must be finite';m=joblib.load(root/f'RF_tuned_{args.scope}.joblib');a['RF_score']=m.predict_proba(X)[:,1];t=.5 if args.threshold=='fixed' else spec['inner_selected_threshold'];a['threshold']=t;a['predicted_work']=(a.RF_score>=t).astype(int);a.to_csv(args.out,index=False)
