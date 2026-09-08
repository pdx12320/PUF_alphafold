import os,sys,json,re,pathlib,hashlib
import numpy as np,pandas as pd
import shlex
from cp_features import get_features

O=pathlib.Path('c388_analysis');u=json.load(open(O/'inventory.json'));wt=next(x['sequence'] for x in u if x['name']=='WT_PUF12-9_C388'); rows=[]; arrays={};mapping=[]
def name(n):return n.removesuffix('_C388').replace('_plus_','+').replace('WT_PUF12-9','PUF12-9')
for unit in u:
 if not unit['name'].endswith('_C388'):continue
 seq=unit['sequence'];assert len(seq)==519 and seq[:7]=='SGSETPG';n=name(unit['name']);mut=[i for i,(a,b) in enumerate(zip(wt,seq)) if a!=b];slots=sorted(set((i-40)//36+1 for i in mut)); assert all(1<=x<=12 for x in slots)
 mapping.append(dict(construct=n,length=519,analysis_length=512,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),mutated_positions=','.join(str(i+1) for i in mut),mutations=';'.join(wt[i]+str(i+1)+seq[i] for i in mut),mutation_group='+'.join('P'+str(i) for i in slots) or 'WT',seed=str(unit['seed'][0]),rna=unit['rna'],unit_key=unit['key']))
 for f in sorted((O/'units'/unit['key']).glob('*full_data*.json')):
  d=json.load(open(f));chains=np.array(d['token_chain_ids']);pi=np.where(chains=='A')[0];ri=np.where(chains=='B')[0];cp0=np.array(d['contact_probs'],np.float32);pa0=np.array(d['pae'],np.float32);cp0=cp0[np.ix_(pi,pi)];pa0=pa0[np.ix_(pi,pi)];pr0=np.array(d['contact_probs'],np.float32)[np.ix_(pi,ri)]
  cif=f.with_name(f.name.replace('full_data','model').replace('.json','.cif')).read_text().splitlines();cols=[l.strip().split('.')[1] for l in cif if l.startswith('_atom_site.')];ci={x:i for i,x in enumerate(cols)};ca={};pls={}
  for line in cif:
   if not line.startswith('ATOM '):continue
   v=shlex.split(line) if '\"' in line or "'" in line else line.split();assert len(v)==len(cols)
   if v[ci['label_asym_id']]!='A':continue
   rn=int(v[ci['label_seq_id']]);pls.setdefault(rn,[]).append(float(v[ci['B_iso_or_equiv']]))
   if v[ci['label_atom_id']]=='CA':ca[rn]=[float(v[ci[k]]) for k in ['Cartn_x','Cartn_y','Cartn_z']]
  assert sorted(ca)==list(range(1,520));coords0=np.array([ca[i] for i in range(1,520)]);pl0=np.array([np.mean(pls[i]) for i in range(1,520)]);assert len(pi)==519 and len(ri)==15
  for offset,scope in [(0,'full519'),(7,'matched512')]:
   cp=cp0[offset:,offset:];pa=(pa0[offset:,offset:]+pa0[offset:,offset:].T)/2;pr=pr0[offset:];coords=coords0[offset:];pl=pl0[offset:];L=len(pl);core=np.arange(40-offset,472-offset);ftr=get_features(cp,pa,coords,pl,pr);pcc=cp[np.ix_(core,core)];pac=pa[np.ix_(core,core)];tri=np.triu_indices(432,1);eig=np.linalg.eigvalsh(np.cov(coords[core].T));gap=abs(np.arange(L)[:,None]-np.arange(L)[None,:]);dens=(cp*(gap>=12)).sum(1)
   ftr.update(plddt_core_mean=float(pl[core].mean()),plddt_core_min_model=float(pl[core].min()),global_core_PAE=float(pac[tri].mean()),anisotropy=float(eig[-1]/eig[0]),RNA_meanCP=float(pr.mean()),RNA_high_per_res=float((pr>.5).sum()/L))
   adj=[];non=[]
   for a in range(12):
    ai=core[a*36:(a+1)*36];mask=np.triu(np.abs(ai[:,None]-ai[None,:])>=4,1);ftr['internal_P'+str(a+1)]=float(cp[np.ix_(ai,ai)][mask].sum()/36)
    for b in range(a+1,12):
     bi=core[b*36:(b+1)*36];mask=np.abs(ai[:,None]-bi[None,:])>=4;v=float(cp[np.ix_(ai,bi)][mask].mean());(adj if b==a+1 else non).append(v);ftr[f'P{a+1}_P{b+1}_CP']=v
   ftr.update(adj_CP_mean=float(np.mean(adj)),adj_CP_min=float(np.min(adj)),nonadj_CP_mean=float(np.mean(non)))
   rows.append(dict(construct=n,scope=scope,seed=str(unit['seed'][0]),model=int(re.search(r'_full_data_(\d+)',f.name)[1]),**ftr))
   if scope=='matched512':
    key=(n,str(unit['seed'][0])); arrays.setdefault(key,[]).append(dict(plddt=pl[core],density=dens[core],rna=pr[core].sum(1),cp=pcc,pae=pac))
 print('Extracted',n,flush=True)
raw=pd.DataFrame(rows);raw.to_csv(O/'model_features.csv',index=False);sf=raw.groupby(['construct','scope','seed']).mean(numeric_only=True).drop(columns='model');sf.to_csv(O/'seed_features.csv');agg=sf.groupby(['construct','scope']).mean().reset_index();ns=sorted(set(k[0] for k in arrays));arr={}
for k in ['plddt','density','rna','cp','pae']:
 arr[k]=np.stack([np.mean([np.mean([v[k] for v in vals],axis=0) for key,vals in arrays.items() if key[0]==n],axis=0) for n in ns])
for n,pl in zip(ns,arr['plddt']):agg.loc[agg.construct==n,'plddt_core_min']=pl.min()
agg.to_csv(O/'construct_features.csv',index=False);pd.DataFrame(mapping).to_csv(O/'sequence_mapping.csv',index=False);np.savez_compressed(O/'construct_arrays.npz',names=ns,**arr)
r=pd.DataFrame(json.load(open(O/'experimental_rows.json')));r['construct']=r.group.str.replace(r'^\d+\s+','',regex=True);r.to_csv(O/'experimental_replicates.csv',index=False);labels=r.groupby('construct').C388.agg(C388_mean='mean',C388_min='min',C388_max='max',n_replicates='size').reset_index();labels['work']=(labels.C388_mean>=.5).astype(int);labels['zero_variant_record']=labels.C388_max==0;labels['threshold_discordant']=(labels.C388_min<.5)&(labels.C388_max>=.5);labels.to_csv(O/'labels.csv',index=False);assert set(labels.construct)<=set(agg.construct)
print('DONE',len(labels),labels.work.sum(),len(agg.construct.unique()))
