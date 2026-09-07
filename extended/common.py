import numpy as np
from scipy.stats import t as tdist

def pearson_cols(X,y):
 X=np.asarray(X,float);y=np.asarray(y,float);xc=X-X.mean(axis=0);yc=y-y.mean();den=np.sqrt((xc*xc).sum(axis=0)*(yc*yc).sum());return np.divide(yc@xc,den,out=np.zeros(X.shape[1]),where=den>1e-12)

def bh(p):
 p=np.array(p,float);order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty_like(q);out[order]=np.minimum(q,1);return out

def discover_ccr(D,y,max_regions=5):
 # Only columns present in every training construct; no imputation for discovery.
 D=np.asarray(D,float);valid=np.isfinite(D).all(0)&(np.nanstd(D,axis=0)>1e-6);groups=[]
 for rep in range(12):
  i=rep*36
  while i<(rep+1)*36:
   if not valid[i]:i+=1;continue
   group=[i];j=i+1
   while j<(rep+1)*36 and len(group)<12 and valid[j]:
    cor=pearson_cols(D[:,group],D[:,j])
    if np.min(cor)<.6:break
    group.append(j);j+=1
   if len(group)>=3:
    value=D[:,group].mean(1);r=float(pearson_cols(value[:,None],y)[0]);co=np.corrcoef(D[:,group],rowvar=False);groups.append({'indices':group,'r':r,'min_internal_r':float(co[np.triu_indices(len(group),1)].min()),'mean_internal_r':float(co[np.triu_indices(len(group),1)].mean())})
   i=max(i+1,j)
 groups=[g for g in groups if abs(g['r'])>=.4];groups.sort(key=lambda g:(-abs(g['r']),g['indices'][0]));return groups[:max_regions]

def build_selected(arr,y,train,max_regions=5,n_residues=8):
 D=arr['density'];groups=discover_ccr(D[train],y[train],max_regions)
 valid=np.isfinite(D[train]).all(0)&(np.nanstd(D[train],axis=0)>1e-6);corr=np.zeros(D.shape[1]);corr[valid]=pearson_cols(D[train][:,valid],y[train]);res=sorted(np.where(valid)[0],key=lambda i:(-abs(corr[i]),i))[:n_residues]
 vals=[];names=[];meta=[]
 for i in res:vals.append(D[:,i]);names.append('residue_density_%03d'%(i+1));meta.append({'kind':'residue','indices':[int(i)]})
 for n,g in enumerate(groups):
  inds=g['indices']
  for typ in ['density','plddt','pae_contact']:
   vals.append(np.nanmean(arr[typ][:,inds],axis=1));names.append('CCR%d_%s'%(n+1,typ));meta.append({'kind':'CCR','indices':inds,'metric':typ})
 return np.stack(vals,axis=1) if vals else np.empty((len(y),0)),names,meta,groups
