import os,sys,json,zipfile,re,io,difflib
import numpy as np,pandas as pd,joblib
from Bio.PDB import MMCIFParser
sys.path.insert(0,os.path.abspath('extended'))
from score_new import get_features
z=zipfile.ZipFile('upload/91-96.zip');O='trm_validation';requests={};rows=[];core_rows=[];parser=MMCIFParser(QUIET=True)
for n in z.namelist():
 if n.endswith('job_request.json'):
  q=json.loads(z.read(n))[0];prot=q['sequences'][0]['proteinChain']['sequence'];rna=q['sequences'][1]['rnaSequence']['sequence'];name=q['name'].rsplit('_',1)[0];requests[n.split('/')[0]]=dict(construct=name,target=q['name'].rsplit('_',1)[1],protein=prot,rna=rna,seed=q['modelSeeds'][0],success=1,scaffold_group='WT_PUF12-9_family')
for n in z.namelist():
 if not re.search('_full_data_\\d+.json$',n):continue
 q=requests[n.split('/')[0]];d=json.loads(z.read(n));ch=np.array(d['token_chain_ids']);ii=np.where(ch=='A')[0];jj=np.where(ch=='B')[0];cp=np.array(d['contact_probs'],dtype=np.float32);pae=np.array(d['pae'],dtype=np.float32)[np.ix_(ii,ii)];st=parser.get_structure('x',io.StringIO(z.read(n.replace('full_data','model').replace('.json','.cif')).decode()));res=list(st[0]['A'].get_residues());assert len(res)==len(q['protein'])==len(ii);assert len(jj)==len(q['rna']);coords=np.array([r['CA'].coord for r in res]);plddt=np.array([np.mean([at.bfactor for at in r]) for r in res]);f=get_features(cp[np.ix_(ii,ii)],pae,coords,plddt,cp[np.ix_(ii,jj)]);rows.append(dict(construct=q['construct'],target=q['target'],seed=q['seed'],model=int(re.search('full_data_(\\d+)',n)[1]),protein_length=len(ii),rna_length=len(jj),**f))
 # Sensitivity: remove extra N-terminal 7 residues from the matrices/coordinates; no refolding.
 core=get_features(cp[np.ix_(ii[7:],ii[7:])],pae[7:,7:],coords[7:],plddt[7:],cp[np.ix_(ii[7:],jj)])
 core_rows.append(dict(construct=q['construct'],target=q['target'],seed=q['seed'],model=int(re.search('full_data_(\\d+)',n)[1]),**core))
raw=pd.DataFrame(rows);raw.to_csv(O+'/model_features.csv',index=False);features=json.load(open('combined12/manifest.json'))['feature_sets'];agg=raw.groupby(['construct','target','seed']).mean(numeric_only=True).drop(columns='model').reset_index()
for s,ff in features.items():
 m=joblib.load('combined12/'+s+'_RF.joblib');agg[s+'_RF']=m.predict_proba(agg[ff].to_numpy())[:,1];raw[s+'_RF']=m.predict_proba(raw[ff].to_numpy())[:,1]
agg['success']=1;agg['scaffold_group']='WT_PUF12-9_family';agg.to_csv(O+'/frozen_RF_scores.csv',index=False);raw.to_csv(O+'/model_features_and_scores.csv',index=False)
seqs={q['construct']:q['protein'] for q in requests.values()};wt=seqs['WT_PUF12-9'];mut=seqs['P8-R6-GVE'];assert len(wt)==len(mut);diff=[dict(position_1based=i+1,WT=a,mutant=b,change=f'{a}{i+1}{b}') for i,(a,b) in enumerate(zip(wt,mut)) if a!=b];pd.DataFrame(diff).to_csv(O+'/mutations_vs_WT.csv',index=False)
old=json.load(open('previous/results/sequences.json'));audit=[]
for name,seq in seqs.items():
 sims=[]
 for k,v in old.items():
  s=v['protein'];mat=difflib.SequenceMatcher(None,seq,s,autojunk=False);sims.append((mat.ratio(),k,len(s),mat.find_longest_match().size))
 sims.sort(reverse=True);audit.append(dict(construct=name,protein_length=len(seq),nearest_previous_construct=sims[0][1],sequence_match_ratio=sims[0][0],nearest_length=sims[0][2],longest_exact_segment=sims[0][3],exact_previous_matches=[k for k,v in old.items() if seq==v['protein']]))
json.dump(dict(requests=requests,mutations_vs_WT=diff,previous_sequence_comparison=audit),open(O+'/sequence_audit.json','w'),indent=2)
print(agg[['construct','target','protein_length','rna_length']+[s+'_RF' for s in features]].to_string(index=False));print('DIFF',diff);print('AUDIT',audit)

core=pd.DataFrame(core_rows).groupby(['construct','target','seed']).mean(numeric_only=True).drop(columns='model').reset_index()
for key,ff in features.items():core[key+'_RF']=joblib.load('combined12/'+key+'_RF.joblib').predict_proba(core[ff].to_numpy())[:,1]
core.to_csv(O+'/Nterminal7_removed_sensitivity.csv',index=False)
print('CORE',core[['construct','target']+[s+'_RF' for s in features]].to_string(index=False))
