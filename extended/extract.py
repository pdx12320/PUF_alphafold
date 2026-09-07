import os,json,re,zipfile,glob,io,collections
import numpy as np,pandas as pd
from Bio.PDB import MMCIFParser
from scipy.spatial.distance import pdist
O=os.path.dirname(__file__);seqs=json.load(open('previous/results/sequences.json'));names=sorted(seqs)
ref='puf12_r123_r567_r567_r678'; templates={};mapping={};ann=[]
for name in names:
 seq=seqs[name]['protein'];removed=set()
 for m in re.finditer('MNDGPHS',seq):removed.update(range(m.start(),m.end()))
 keep=[i for i in range(len(seq)) if i not in removed];clean=''.join(seq[i] for i in keep);anchors=[m.start() for m in re.finditer('[SCN][YR][FV][IV][EQR]',clean)]
 assert len(anchors) in [11,12] and anchors==[44+36*i for i in range(len(anchors))],(name,anchors)
 core=np.full(432,-1,int)
 for r,start in enumerate([a-11 for a in anchors]):core[r*36:(r+1)*36]=keep[start:start+36]
 mapping[name]=core
 for aidx,native in enumerate(core):
  if native<0:continue
  rep=aidx//36+1; pos=aidx%36+1
  dist=min([abs(native-j) for j in removed],default=999)
  ann.append({'construct':name,'aligned_index':aidx+1,'repeat_slot':rep,'repeat_position':pos,'native_residue':native+1,'aa':seq[native],'TRM':pos in [12,13,16],'loop_near_5aa':dist<=5,'loop_distance':dist})
json.dump({k:v.tolist() for k,v in mapping.items()},open(O+'/core_mapping.json','w'),indent=2)
# Source-repeat identities inferred by nearest 36-aa reference module, ignoring TRM positions.
refseq=seqs[ref]['protein'];refs=[ ''.join(refseq[j] for j in mapping[ref][36*r:36*(r+1)]) for r in range(12)]
for r,source in enumerate([1,2,3,5,6,7,5,6,7,6,7,8]):templates.setdefault(source,refs[r])
r4name='puf_12_r4r4';templates[4]=''.join(seqs[r4name]['protein'][j] for j in mapping[r4name][360:396])
module_rows=[]
for name,core in mapping.items():
 for r in range(12):
  inds=core[36*r:36*(r+1)]
  if min(inds)<0:continue
  sq=''.join(seqs[name]['protein'][j] for j in inds);scores={k:sum(sq[j]==t[j] for j in range(36) if j not in [11,12,15])/33 for k,t in templates.items()};best=max(scores,key=scores.get)
  module_rows.append({'construct':name,'repeat_slot':r+1,'native_start':int(inds[0]+1),'native_end':int(inds[-1]+1),'source_repeat_inferred':best,'source_identity_excluding_TRM':scores[best],'TRM':sq[11]+sq[12]+sq[15],'core_sequence':sq})
pd.DataFrame(module_rows).to_csv(O+'/repeat_mapping.csv',index=False)
acc={};features=[];parser=MMCIFParser(QUIET=True);cif_saved=set()
for path in sorted(glob.glob('inputs/folds*.zip')):
 z=zipfile.ZipFile(path)
 for member in z.namelist():
  if not re.search('_full_data_[0-4].json$',member):continue
  folder=member.split('/')[0];name,seed=re.match(r'(.*)_rna_.*_seed(\d+)$',folder).groups();seed=int(seed);model=int(re.search('_full_data_(\d)',member)[1]);d=json.loads(z.read(member));chain=np.array(d['token_chain_ids']);pa=np.where(chain=='A')[0];pb=np.where(chain=='B')[0];n=len(pa)
  cp=np.array(d['contact_probs'],dtype=np.float32);pae=np.array(d['pae'],dtype=np.float32);pp=cp[np.ix_(pa,pa)];pr=cp[np.ix_(pa,pb)];paepp=(pae[np.ix_(pa,pa)]+pae[np.ix_(pa,pa)].T)/2
  ci=member.replace('full_data','model').replace('.json','.cif');ct=z.read(ci).decode();st=parser.get_structure('x',io.StringIO(ct));res=list(st[0]['A'].get_residues());assert len(res)==n
  coords=np.array([r['CA'].coord for r in res]);plddt=np.array([np.mean([a.bfactor for a in r.get_atoms()]) for r in res]);assert all(r.id[1]==i+1 for i,r in enumerate(res))
  if name not in cif_saved:open(O+'/'+name+'_representative.cif','w').write(ct);cif_saved.add(name)
  inds=mapping[name];valid=np.where(inds>=0)[0];src=inds[valid];nonlocalmask=abs(np.arange(n)[:,None]-np.arange(n)[None,:])>=12
  a_cp=np.full((432,432),np.nan,np.float32);a_cp[np.ix_(valid,valid)]=pp[np.ix_(src,src)]
  a_pae=np.full((432,432),np.nan,np.float32);a_pae[np.ix_(valid,valid)]=paepp[np.ix_(src,src)]
  per={}
  for k,v in [('density',(pp*nonlocalmask).sum(1)),('high_density',((pp>.5)*nonlocalmask).sum(1)),('plddt',plddt),('rna',pr.sum(1)),('rna_max',pr.max(1)),('pae_contact',(paepp*pp*nonlocalmask).sum(1)/np.maximum((pp*nonlocalmask).sum(1),1e-8))]:
   tmp=np.full(432,np.nan,np.float32);tmp[valid]=v[src];per[k]=tmp
  per.update(cp=a_cp,pae=a_pae)
  key=(name,seed)
  if key not in acc:acc[key]={k:[] for k in per}
  for k,v in per.items():acc[key][k].append(v)
  rg=np.sqrt(np.mean(np.sum((coords-coords.mean(0))**2,axis=1)));cc=coords[src];rgcore=np.sqrt(np.mean(np.sum((cc-cc.mean(0))**2,axis=1)));eigs=np.linalg.eigvalsh(np.cov(cc.T));sep=np.triu(nonlocalmask,1)
  features.append({'construct':name,'seed':seed,'model':model,'plddt_protein_mean':plddt.mean(),'plddt_protein_p10':np.quantile(plddt,.1),'plddt_core_mean':plddt[src].mean(),'plddt_core_below70':np.mean(plddt[src]<70),'pae_protein_mean':paepp[sep].mean(),'pae_contact_weighted':(paepp[sep]*pp[sep]).sum()/pp[sep].sum(),'pae_RNA_symmetric':((pae[np.ix_(pa,pb)]+pae[np.ix_(pb,pa)].T)/2).mean(),'Rg_CA_A':rg,'Rg_core_CA_A':rgcore,'Rg_core_length_normalized':rgcore/len(src)**(1/3),'core_max_CA_distance_A':pdist(cc).max(),'core_anisotropy':eigs[-1]/max(eigs[0],1e-8),'core_density_mean':per['density'][valid].mean(),'core_density_sd':per['density'][valid].std(),'core_density_p10':np.quantile(per['density'][valid],.1),'core_high_density_mean':per['high_density'][valid].mean()})
  if len(features)%50==0:print('extracted',len(features),flush=True)
seedkeys=sorted(acc);data={k:np.stack([np.mean(acc[key][k],axis=0) for key in seedkeys]) for k in next(iter(acc.values()))};np.savez_compressed(O+'/seed_arrays.npz',**{k:v for k,v in data.items() if k not in ['cp','pae']},names=np.array([k[0] for k in seedkeys]),seeds=np.array([k[1] for k in seedkeys]))
means={k:np.stack([np.mean(data[k][np.array([key[0]==name for key in seedkeys])],axis=0) for name in names]) for k in data};np.savez_compressed(O+'/construct_arrays.npz',**means,names=np.array(names))
pd.DataFrame(features).to_csv(O+'/structural_model_features.csv',index=False);sf=pd.DataFrame(features).groupby(['construct','seed']).mean().drop(columns='model');sf.to_csv(O+'/structural_seed_features.csv');sf.groupby('construct').mean().to_csv(O+'/structural_construct_features.csv')
an=pd.DataFrame(ann);anno=[]
for _,row in an.iterrows():
 c=names.index(row.construct);i=int(row.aligned_index)-1;v=means['cp'][c,i];other=(np.arange(432)//36!=i//36)&(abs(np.arange(432)-i)>=12);annointerface=bool(np.nanmax(np.where(other,v,np.nan))>.5)
 rr=row.to_dict();rr.update(inter_repeat_interface=annointerface,rna_contact_region=bool(means['rna_max'][c,i]>.5),mean_plddt=float(means['plddt'][c,i]),nonlocal_CP_density=float(means['density'][c,i]));anno.append(rr)
pd.DataFrame(anno).to_csv(O+'/residue_mapping_all_constructs.csv',index=False)
print('DONE',len(names),len(seedkeys),len(features),flush=True)
