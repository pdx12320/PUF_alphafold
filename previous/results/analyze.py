import json,zipfile,glob,re,os,itertools
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,confusion_matrix
O=os.path.dirname(__file__)
success={'puf12_r123_r567_r567_r678','puf12_r123_r567_r6r7_r5678','puf12_r123_r567_r5loop67_loopr678'}
rows=[]; matrices={}; seqs={}; provenance=[]
for path in sorted(glob.glob('upload/*.zip')):
 z=zipfile.ZipFile(path)
 for name in z.namelist():
  if not re.search(r'_full_data_\d.json$',name):continue
  folder=name.split('/')[0]; construct,seed=re.match(r'(.*)_rna_.*_seed(\d+)$',folder).groups(); seed=int(seed); model=int(re.search(r'_full_data_(\d)',name)[1])
  d=json.loads(z.read(name)); cp=np.array(d['contact_probs']); chains=np.array(d['token_chain_ids']); ids=np.array(d['token_res_ids']); a=np.where(chains=='A')[0]; b=np.where(chains=='B')[0]; p=cp[np.ix_(a,a)]; r=cp[np.ix_(a,b)]; n=len(a)
  assert np.isfinite(cp).all() and cp.min()>=0 and cp.max()<=1
  request=json.loads(z.read(folder+'/fold_'+folder+'_job_request.json'))[0]
  seq=request['sequences'][0]['proteinChain']['sequence']; rna=request['sequences'][1]['rnaSequence']['sequence']; assert n==len(seq) and len(b)==len(rna)==13
  if construct in seqs:assert seqs[construct]['protein']==seq
  seqs[construct]={'protein':seq,'rna':rna}
  summary=json.loads(z.read(name.replace('full_data','summary_confidences')))
  f={'construct':construct,'success':int(construct in success),'seed':seed,'model':model,'protein_length':n,'pr_cp_sum':r.sum(),'pr_cp_per_nt':r.sum()/len(b),'pr_top20_mean':np.sort(r.ravel())[-20:].mean(),'pr_rna_max_mean':r.max(axis=0).mean(),'pr_rna_max_min':r.max(axis=0).min(),'pr_rna_coverage_05':(r.max(axis=0)>.5).mean(),'pr_protein_coverage_05':(r.max(axis=1)>.5).mean(),'pr_high_contacts_per_nt':(r>.5).sum()/len(b)}
  for gap in [4,12,24]:
   mask=np.triu(np.ones((n,n),bool),k=gap);v=p[mask]; f['pp_nonlocal%d_per_res'%gap]=v.sum()/n; f['pp_nonlocal%d_high_per_res'%gap]=(v>.5).sum()/n
  for j in range(13):f['rna%02d_cp_sum'%(j+1)]=r[:,j].sum(); f['rna%02d_cp_max'%(j+1)]=r[:,j].max()
  for j,ind in enumerate(np.array_split(np.arange(n),12)):f['protein_bin%02d_pr_cp'%(j+1)]=r[ind].sum()/len(ind)
  f.update({'iptm':summary.get('iptm',np.nan),'ptm':summary.get('ptm',np.nan),'ranking_score':summary.get('ranking_score',np.nan),'mean_plddt':np.mean(d['atom_plddts']),'interface_pae':np.array(d['pae'])[np.ix_(a,b)].mean()})
  rows.append(f); matrices.setdefault(construct,[]).append(r); provenance.append({'archive':os.path.basename(path),'member':name,'construct':construct,'seed':seed,'model':model,'cp_symmetry_max_abs':float(abs(cp-cp.T).max())})
raw=pd.DataFrame(rows).sort_values(['construct','seed','model']); raw.to_csv(O+'/model_features.csv',index=False)
seed=raw.groupby(['construct','success','seed']).mean(numeric_only=True).drop(columns='model').reset_index(); seed.to_csv(O+'/seed_features.csv',index=False)
agg=seed.groupby(['construct','success']).mean(numeric_only=True).drop(columns='seed').reset_index();agg.to_csv(O+'/construct_features.csv',index=False)
json.dump(seqs,open(O+'/sequences.json','w'),indent=2);pd.DataFrame(provenance).to_csv(O+'/input_audit.csv',index=False)
np.savez_compressed(O+'/mean_protein_rna_CP.npz',**{k:np.mean(v,axis=0) for k,v in matrices.items()})
y=agg.success.values; combos=list(itertools.combinations(range(len(y)),3)); perms=np.zeros((len(combos),len(y)),int)
for i,c in enumerate(combos):perms[i,list(c)]=1
cpcols=[c for c in agg if c.startswith(('pr_','pp_','rna','protein_bin'))]
stats=[]
for c in cpcols+['protein_length','iptm','ptm','mean_plddt','interface_pae']:
 v=agg[c].values; auc=roc_auc_score(y,v); stats.append({'feature':c,'success_mean':v[y==1].mean(),'failure_mean':v[y==0].mean(),'auc_high_is_success':auc,'separation_auc':max(auc,1-auc),'MW_two_sided_p':mannwhitneyu(v[y==1],v[y==0],alternative='two-sided',method='exact').pvalue})
st=pd.DataFrame(stats).sort_values('separation_auc',ascending=False);st.to_csv(O+'/feature_comparison.csv',index=False)
# Family-wise exact permutation: strongest absolute class-mean rank separation among all CP features.
from scipy.stats import rankdata
ranks=np.stack([rankdata(agg[c]) for c in cpcols],axis=1)
aucs=(perms@ranks-6)/33; observed=max(st[st.feature.isin(cpcols)].separation_auc);maxnull=np.maximum(aucs,1-aucs).max(axis=1)
family_p=float(np.mean(maxnull>=observed-1e-10))
sets={'CP_summary':[c for c in cpcols if c.startswith(('pr_','pp_'))], 'CP_RNA_profile':[c for c in cpcols if c.startswith('rna')], 'length_only':['protein_length']}
pred=agg[['construct','success']].copy();metrics=[]
for label,cols in sets.items():
 X=agg[cols].values; folds=[]
 for i in range(len(y)):
  ix=np.arange(len(y))!=i; scaler=StandardScaler().fit(X[ix]); folds.append((ix,scaler.transform(X[ix]),scaler.transform(X[i:i+1])))
 def fitpred(labels):
  out=[]
  for ix,x,t in folds:
   m=LogisticRegression(C=.1,class_weight='balanced',solver='liblinear',random_state=0).fit(x,labels[ix]);out.append(m.predict_proba(t)[0,1])
  return np.array(out)
 prob=fitpred(y);obs=roc_auc_score(y,prob); null=np.array([roc_auc_score(p,fitpred(p)) for p in perms]);pred[label+'_prob']=prob
 metrics.append({'model':label,'n_features':len(cols),'LOCO_AUC':obs,'LOCO_balanced_accuracy':balanced_accuracy_score(y,prob>=.5),'confusion_matrix_TN_FP_FN_TP':confusion_matrix(y,prob>=.5).ravel().tolist(),'exact_permutation_p_AUC':float(np.mean(null>=obs-1e-10))})
pred.to_csv(O+'/leave_construct_out_predictions.csv',index=False)
json.dump({'n_constructs':len(agg),'n_seeds':len(seed),'n_models':len(raw),'models':metrics,'best_CP_feature_exploratory_AUC':observed,'CP_feature_search_familywise_exact_p':family_p},open(O+'/metrics.json','w'),indent=2)
print(agg[['construct','success','protein_length','pr_cp_sum','pr_rna_max_mean','pp_nonlocal12_per_res','iptm']].to_string(index=False));print(st.head(12).to_string(index=False));print(json.dumps(metrics,indent=2));print('family p',family_p)
