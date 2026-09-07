import os,sys,json,re,zipfile,io,collections
import numpy as np,pandas as pd
from Bio.PDB import MMCIFParser
sys.path.insert(0,os.path.abspath('extended'));from score_new import get_features
O='architecture_validation';base=pd.read_csv('combined12/training_data.csv');old=pd.read_csv('extended/all_construct_features.csv');labels=dict(zip(old.construct,old.success));labels.update(dict(zip(base.construct,base.success)));labels.update({'WT_PUF12-9':1,'P8-R6-GVE':1})
t=pd.read_csv('extended/repeat_mapping.csv');templates={}
for _,r in t.sort_values('source_identity_excluding_TRM',ascending=False).iterrows():templates.setdefault(int(r.source_repeat_inferred),r.core_sequence)
acc={};raw=[];meta={};mappingrows=[];parser=MMCIFParser(QUIET=True)
paths=['upload/folds_2026_09_07_02_49.zip','upload/folds_2026_09_07_02_48.zip','upload/91-96.zip']
for path in paths:
 z=zipfile.ZipFile(path)
 for n in z.namelist():
  if not re.search(r'_full_data_\d+\.json$',n):continue
  folder=n.split('/')[0];q=json.loads(z.read(folder+'/fold_'+folder+'_job_request.json'))[0];match=re.match(r'(.*)_rna_.*_seed(\d+)$',folder);trm=match is None
  name=match[1] if match else q['name'].rsplit('_',1)[0]
  if name not in labels:continue
  sq=q['sequences'][0]['proteinChain']['sequence'];rna=q['sequences'][1]['rnaSequence']['sequence'];offset=7 if trm else 0;seq=sq[offset:];seed=str(q['modelSeeds'][0]) if trm else match[2];context=q['name'].rsplit('_',1)[1] if trm else 'RNA13';mid=int(re.search(r'_full_data_(\d+)',n)[1])
  if name not in meta:
   removed=set()
   for mt in re.finditer('MNDGPHS',seq):removed.update(range(mt.start(),mt.end()))
   keep=[i for i in range(len(seq)) if i not in removed];clean=''.join(seq[i] for i in keep)
   if trm:anchors=[44+36*i for i in range(12)]
   else:anchors=[mt.start() for mt in re.finditer('[SCN][YR][FV][IV][EQR]',clean)];assert anchors==[44+36*i for i in range(len(anchors))]
   mp=np.full(432,-1,int);order=[]
   for r,an in enumerate(anchors):
    inds=keep[an-11:an-11+36];s=''.join(seq[k] for k in inds);scores={k:sum(s[i]==v[i] for i in range(36) if i not in [11,12,15])/33 for k,v in templates.items()};best=max(scores,key=scores.get);order.append(best);mp[r*36:(r+1)*36]=inds
    for j,v in enumerate(inds):mappingrows.append(dict(construct=name,slot=r+1,position=j+1,aligned_index=r*36+j+1,native_residue=v+offset+1,analysis_residue=v+1,aa=seq[v],source=best,TRM=j in [11,12,15],near_loop_5=min([abs(v-k) for k in removed],default=999)<=5))
   meta[name]=dict(construct=name,success=int(labels[name]),protein_length=len(seq),uploaded_length=len(sq),repeat_count=len(anchors),source_order=','.join(map(str,order)),batch='TRM' if trm else ('previous' if path.startswith('inputs') else 'new'),mapping=mp.tolist(),rna=rna)
  mp=np.array(meta[name]['mapping']);valid=mp>=0;src=mp[valid];d=json.loads(z.read(n));ch=np.array(d['token_chain_ids']);ii=np.where(ch=='A')[0][offset:];jj=np.where(ch=='B')[0];cp0=np.array(d['contact_probs'],np.float32);cp=cp0[np.ix_(ii,ii)];pa0=np.array(d['pae'],np.float32);pa=pa0[np.ix_(ii,ii)];pa=(pa+pa.T)/2;pr=cp0[np.ix_(ii,jj)];st=parser.get_structure('x',io.StringIO(z.read(n.replace('full_data','model').replace('.json','.cif')).decode()));res=list(st[0]['A'].get_residues())[offset:];coords=np.array([r['CA'].coord for r in res]);pl=np.array([np.mean([at.bfactor for at in r]) for r in res]);L=len(coords);assert L==len(seq)
  f=get_features(cp,pa,coords,pl,pr);gap=abs(np.arange(L)[:,None]-np.arange(L)[None,:]);tri=np.triu(np.ones((L,L),bool),1);hc=np.triu((gap>=12)&(cp>.5)&(pa<5)&(pl[:,None]>=80)&(pl[None,:]>=80),1);rg=np.sqrt(np.mean(np.sum((coords-coords.mean(0))**2,axis=1)));eig=np.linalg.eigvalsh(np.cov(coords[src].T));f.update(plddt_min=float(pl.min()),global_PAE=float(pa[tri].mean()),Rg=float(rg),anisotropy=float(eig[-1]/max(eig[0],1e-8)),compactness=float(L/(4*np.pi*rg**3/3)),HC_nonlocal_density=float(hc.sum()/L))
  ac=np.full((432,432),np.nan,np.float32);ap=ac.copy();ac[np.ix_(valid,valid)]=cp[np.ix_(src,src)];ap[np.ix_(valid,valid)]=pa[np.ix_(src,src)];dens=np.full(432,np.nan);dens[valid]=(cp*(gap>=12)).sum(1)[src];pc=np.full(432,np.nan);pc[valid]=pl[src];rnares=np.full(432,np.nan);rnares[valid]=pr.max(1)[src]
  adj=[];non=[];adjpa=[];intern=[]
  for r in range(meta[name]['repeat_count']):
   ri=mp[r*36:(r+1)*36];mask=np.triu(abs(ri[:,None]-ri[None,:])>=4,1);f[f'internal_P{r+1}']=float(cp[np.ix_(ri,ri)][mask].sum()/36);intern.append(f[f'internal_P{r+1}'])
   for s in range(r+1,meta[name]['repeat_count']):
    sj=mp[s*36:(s+1)*36];mask=abs(ri[:,None]-sj[None,:])>=4;block=cp[np.ix_(ri,sj)];paeb=pa[np.ix_(ri,sj)];key=f'P{r+1}_P{s+1}';f[key+'_CP']=float(block[mask].mean());f[key+'_PAE']=float(paeb[mask].mean());(adj if s==r+1 else non).append(f[key+'_CP'])
    if s==r+1:adjpa.append(f[key+'_PAE'])
  f.update(adj_CP_mean=np.mean(adj),adj_CP_min=np.min(adj),nonadj_CP_mean=np.mean(non),adj_PAE_mean=np.mean(adjpa),internal_CP_mean=np.mean(intern))
  raw.append(dict(construct=name,context=context,seed=seed,model=mid,**f));key=(name,context,seed)
  if key not in acc:acc[key]=dict(n=0,cp=np.zeros_like(ac),pae=np.zeros_like(ac),density=np.zeros(432),plddt=np.zeros(432),rna=np.zeros(432))
  acc[key]['n']+=1
  for k,v in [('cp',ac),('pae',ap),('density',dens),('plddt',pc),('rna',rnares)]:acc[key][k]+=v
  if len(raw)%100==0:print('models',len(raw),flush=True)
raw=pd.DataFrame(raw);raw.to_csv(O+'/model_features.csv',index=False);seedf=raw.groupby(['construct','context','seed']).mean(numeric_only=True).drop(columns='model');seedf.to_csv(O+'/seed_features.csv');agg=seedf.groupby(['construct','context']).mean().groupby('construct').mean();names=sorted(meta);arrays={}
for k in ['cp','pae','density','plddt','rna']:
 arrays[k]=np.stack([np.mean([v[k]/v['n'] for key,v in acc.items() if key[0]==name],axis=0) for name in names])
np.savez_compressed(O+'/arrays.npz',names=names,**arrays);m=pd.DataFrame([{k:v for k,v in q.items() if k!='mapping'} for q in meta.values()]).set_index('construct');groups=sorted(m.source_order.unique());gmap={g:'G%02d'%(i+1) for i,g in enumerate(groups)};m['architecture']=m.source_order.map(gmap);m['n_seeds']=m.index.map(lambda n:sum(k[0]==n for k in acc));m['n_models']=m.index.map(raw.construct.value_counts());m.join(agg).to_csv(O+'/construct_features.csv');pd.DataFrame(mappingrows).to_csv(O+'/residue_mapping.csv',index=False);json.dump(meta,open(O+'/mapping.json','w'),indent=2);print('DONE',len(raw),len(meta),flush=True)
