"""Score new AlphaFold-server predictions with frozen CP-summary models or the four-feature logistic formula.
Usage: python score_new.py --zip folds.zip --out new_scores.csv
       python score_new.py --features all_construct_features.csv --out scores.csv
All seed/model files for a construct are averaged. Scores are not calibrated probabilities.
"""
import argparse,os,json,zipfile,re,io
import numpy as np,pandas as pd,joblib
from Bio.PDB import MMCIFParser
from scipy.special import expit
O=os.path.dirname(__file__)
def get_features(cp,pae,coords,plddt,pr):
 n=len(coords);sep=np.triu(np.ones((n,n),bool),12);p=(pae+pae.T)/2;rg=np.sqrt(np.mean(np.sum((coords-coords.mean(0))**2,axis=1)))
 features={'pp_nonlocal12_high_per_res':float((cp[sep]>.5).sum()/n),'plddt_protein_mean':float(np.mean(plddt)),'pae_contact_weighted':float((p[sep]*cp[sep]).sum()/cp[sep].sum()),'Rg_full_length_normalized':float(rg/n**(1/3))}
 for gap in [4,12,24]:
  mask=np.triu(np.ones((n,n),bool),gap);features['pp_nonlocal%d_per_res'%gap]=float(cp[mask].sum()/n);features['pp_nonlocal%d_high_per_res'%gap]=float((cp[mask]>.5).sum()/n)
 features.update(pr_cp_sum=float(pr.sum()),pr_cp_per_nt=float(pr.sum()/pr.shape[1]),pr_top20_mean=float(np.sort(pr.ravel())[-20:].mean()),pr_rna_max_mean=float(pr.max(0).mean()),pr_rna_max_min=float(pr.max(0).min()),pr_rna_coverage_05=float((pr.max(0)>.5).mean()),pr_protein_coverage_05=float((pr.max(1)>.5).mean()),pr_high_contacts_per_nt=float((pr>.5).sum()/pr.shape[1]))
 return features
def run():
 ap=argparse.ArgumentParser();ap.add_argument('--zip',nargs='+');ap.add_argument('--features');ap.add_argument('--construct');ap.add_argument('--model',choices=['formula4_LR','summary_LR','summary_RF','summary_XGB'],default='summary_XGB');ap.add_argument('--out',default='new_scores.csv');args=ap.parse_args();f=json.load(open(O+'/scoring_formula.json'))
 if bool(args.zip)==bool(args.features):ap.error('Choose exactly one of --zip or --features')
 if args.features:df=pd.read_csv(args.features)
 else:
  parser=MMCIFParser(QUIET=True);rows=[]
  for path in args.zip:
   z=zipfile.ZipFile(path)
   for member in z.namelist():
    if not re.search(r'_full_data_\d+\.json$',member):continue
    folder=member.split('/')[0];m=re.match(r'(.*)_rna_.*_seed(\d+)$',folder)
    if m:name,seed=m[1],int(m[2])
    else:raise ValueError('For ZIP input folders must identify constructs and seeds using *_rna_*_seedN; alternatively use --features.')
    if args.construct and name!=args.construct:continue
    d=json.loads(z.read(member));chains=np.array(d['token_chain_ids']);idx=np.where(chains=='A')[0];cp=np.array(d['contact_probs'],dtype=np.float32)[np.ix_(idx,idx)];pae=np.array(d['pae'],dtype=np.float32)[np.ix_(idx,idx)];text=z.read(member.replace('full_data','model').replace('.json','.cif')).decode();st=parser.get_structure('x',io.StringIO(text));res=list(st[0]['A'].get_residues());assert len(res)==len(idx);coords=np.array([r['CA'].coord for r in res]);plddt=np.array([np.mean([a.bfactor for a in r.get_atoms()]) for r in res]);rows.append({'construct':name,'seed':seed,**get_features(cp,pae,coords,plddt,np.array(d['contact_probs'],dtype=np.float32)[np.ix_(idx,np.where(chains=='B')[0])])})
  if not rows:raise ValueError('No matching predictions')
  sd=pd.DataFrame(rows).groupby(['construct','seed']).mean();df=sd.groupby('construct').mean().reset_index();df['n_seeds']=df.construct.map(sd.reset_index().groupby('construct').seed.nunique())
 X=df[f['feature_names']].to_numpy(float)
 if not np.isfinite(X).all():raise ValueError('Scoring requires all four finite features')
 standardized=(X-np.array(f['mean']))/np.array(f['scale']);score=expit(f['standardized_intercept']+standardized@np.array(f['standardized_coefficients']));
 if args.model!='formula4_LR':
  mf=json.load(open(O+'/summary_model_manifest.json'));features=df[mf['columns']].values
  if not np.isfinite(features).all():raise ValueError('All CP summary features must be finite')
  score=joblib.load(O+'/final_'+args.model+'.joblib').predict_proba(features)[:,1]
 df['model']=args.model;df['exploratory_score']=score;df['threshold_05_class']=np.where(score>=.5,'success_like','failure_like');df['structural4_any_abs_z_gt3']=(abs(standardized)>3).any(1);df.to_csv(args.out,index=False);print(df[['construct','exploratory_score','threshold_05_class','structural4_any_abs_z_gt3']].to_string(index=False))
if __name__=='__main__':run()
