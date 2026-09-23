#!/usr/bin/env python3
"""Fit the fixed scaffold recipe after excluding all four test constructs."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
R=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent/'model1';OUT.mkdir(exist_ok=True)
p=json.loads((R/'results_20260922/scaffold_RF/audit/protocol.json').read_text())
d=pd.read_csv(R/'results_20260922/scaffold_RF/data/scaffolds24.csv');test=d.construct.isin(p['fixed4']);cols=p['CP_structure9']
assert test.sum()==4 and (~test).sum()==20
model=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),RandomForestClassifier(**p['RF']))
model.fit(d.loc[~test,cols],d.loc[~test,'success'])
z=d.loc[test,['construct','success']].copy();z['score']=model.predict_proba(d.loc[test,cols])[:,1];z['prediction']=(z.score>=.5).astype(int);z['correct']=z.prediction==z.success
z.to_csv(OUT/'predictions.csv',index=False)
d.loc[~test,['construct','success']+cols].to_csv(OUT/'training_data.csv',index=False)
(OUT/'split.json').write_text(json.dumps(dict(train=d.loc[~test,'construct'].tolist(),test=d.loc[test,'construct'].tolist(),test_excluded_from_fit_and_preprocessing=True,features=cols,RF=p['RF'],cutoff=.5,correct=int(z.correct.sum())),indent=2))
print(z.to_string(index=False))
