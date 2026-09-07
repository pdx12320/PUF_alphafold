"""Predict using the refitted 24-construct models. Default: CP_density3_RF.
python combined12/predict.py --zip new_folds.zip --out predictions.csv
python combined12/predict.py --features features.csv --model CP_structure17_RF --out predictions.csv
Scores are uncalibrated. ZIP mode supports the three-density feature set only.
"""
import os,argparse,zipfile,re,json
import numpy as np,pandas as pd,joblib
O=os.path.dirname(__file__)
a=argparse.ArgumentParser();a.add_argument('--zip',nargs='+');a.add_argument('--features');a.add_argument('--model',default='CP_density3_RF');a.add_argument('--out',default='predictions.csv');args=a.parse_args();manifest=json.load(open(O+'/manifest.json'));feature_set=args.model.rsplit('_',1)[0];columns=manifest['feature_sets'][feature_set]
if bool(args.zip)==bool(args.features):a.error('Use --zip or --features')
if args.features:df=pd.read_csv(args.features)
else:
 if feature_set!='CP_density3':a.error('ZIP mode supports CP_density3 models; other models need feature CSV')
 rows=[];seen=set()
 for path in args.zip:
  z=zipfile.ZipFile(path)
  for n in z.namelist():
   if not re.search(r'_full_data_\d+\.json$',n):continue
   folder=n.split('/')[0];match=re.match(r'(.*)_rna_.*_seed(\d+)$',folder)
   if not match:raise ValueError('Expected *_rna_*_seedN folder')
   name,seed=match[1],int(match[2]);model=int(re.search(r'_full_data_(\d+)',n)[1]);key=(name,seed,model)
   if key in seen:raise ValueError('Duplicate model: '+str(key))
   seen.add(key);req=json.loads(z.read(folder+'/fold_'+folder+'_job_request.json'))[0];sequence=req['sequences'][0]['proteinChain']['sequence'];count=len(re.findall('[SCN][YR][FV][IV][EQR]',sequence))
   if count!=12:continue
   d=json.loads(z.read(n));ids=np.where(np.array(d['token_chain_ids'])=='A')[0];cp=np.array(d['contact_probs'])[np.ix_(ids,ids)];L=len(ids);assert L==len(sequence);row={'construct':name,'seed':seed}
   for gap in [4,12,24]:row['pp_nonlocal%d_high_per_res'%gap]=float((cp[np.triu_indices(L,gap)]>.5).sum()/L)
   rows.append(row)
 if not rows:raise ValueError('No 12-repeat constructs identified; motif-based repeat counting may need an explicit mapping for other scaffolds')
 seeds=pd.DataFrame(rows).groupby(['construct','seed']).mean();df=seeds.groupby('construct').mean().reset_index()
X=df[columns].to_numpy(float)
if not np.isfinite(X).all():raise ValueError('Missing/nonfinite features')
m=joblib.load(O+'/'+args.model+'.joblib');df['model']=args.model;df['score']=m.predict_proba(X)[:,1];df['prediction']=np.where(df.score>=.5,'work_like','nonwork_like');df.to_csv(args.out,index=False);print(df[['construct','score','prediction']].to_string(index=False))
