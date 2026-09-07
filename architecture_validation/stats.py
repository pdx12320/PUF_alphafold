import os,itertools,json
import numpy as np,pandas as pd
from scipy.stats import rankdata
O='architecture_validation';a=pd.read_csv(O+'/main_features.csv');z=np.load(O+'/all_arrays.npz');idx=[list(z['names']).index(n) for n in a.construct];D=z['density'][idx].astype(float);CP=z['cp'][idx].astype(float);y=a.success.to_numpy();n=len(y);k=int(y.sum());g=a.architecture.to_numpy();P=np.zeros((__import__('math').comb(n,k),n))
for i,c in enumerate(itertools.combinations(range(n),k)):P[i,list(c)]=1
P=P/k-(1-P)/(n-k)
choices=[list(itertools.combinations(np.where(g==v)[0],int(y[g==v].sum()))) for v in sorted(set(g))];RP=[]
for comb in itertools.product(*choices):
 yy=np.zeros(n)
 for c in comb:yy[list(c)]=1
 for v in set(g):yy[g==v]-=yy[g==v].mean()
 RP.append(yy)
RP=np.array(RP);yr=y.astype(float)
for v in set(g):yr[g==v]-=yr[g==v].mean()

def corr(X,yy):
 X=np.atleast_2d(X);xx=X-X.mean(0);yy=yy-yy.mean();den=np.sqrt((xx*xx).sum(0)*(yy*yy).sum());return np.divide(yy@xx,den,out=np.zeros(X.shape[1]),where=den>1e-12)
def bh(p):
 p=np.array(p);s=np.argsort(p);q=np.minimum.accumulate((p[s]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty(len(p));out[s]=np.minimum(q,1);return out

def stats(X,cols,family):
 X=np.asarray(X,float);out=[]
 for start in range(0,X.shape[1],128):
  V=X[:,start:start+128];de=V[y==1].mean(0)-V[y==0].mean(0);per=P@V;p=(np.abs(per)>=np.abs(de)[None,:]-1e-10).mean(0);r=corr(V,y);rho=corr(rankdata(V,axis=0),y);pooled=((k-1)*V[y==1].var(0,ddof=1)+(n-k-1)*V[y==0].var(0,ddof=1))/(n-2);effect=np.divide(de,np.sqrt(pooled),out=np.zeros(len(de)),where=pooled>1e-20)*(1-3/(4*(n-2)-1));vr=V.copy()
  for v in set(g):vr[g==v]-=vr[g==v].mean(0)
  adjusted=yr@vr/(yr@yr);obs=yr@vr;rp=(np.abs(RP@vr)>=np.abs(obs)[None,:]-1e-10).mean(0);ss=((V-V.mean(0))**2).sum(0);r2=1-np.divide((vr*vr).sum(0),ss,out=np.ones(len(ss)),where=ss>1e-20);partial=corr(vr,yr)
  for j,c in enumerate(cols[start:start+128]):out.append(dict(feature=c,family=family,minimum=V[:,j].min(),maximum=V[:,j].max(),success_mean=V[y==1,j].mean(),failure_mean=V[y==0,j].mean(),delta=de[j],hedges_g=effect[j],pearson_r=r[j],spearman_r=rho[j],p_exact=p[j],architecture_R2=r2[j],adjusted_success_beta=adjusted[j],partial_r=partial[j],p_within_architecture_exact=rp[j]))
 return pd.DataFrame(out)
# Label-independent coherent consecutive density blocks; correlations unchanged by subtracting a fixed reference.
def blocks(D):
 groups=[]
 for r in range(D.shape[1]//36):
  i=r*36
  while i<(r+1)*36:
   group=[i];j=i+1
   while j<(r+1)*36 and len(group)<8:
    if np.std(D[:,j])<1e-8 or np.min(corr(D[:,group],D[:,j]))<.6 or np.min(corr(rankdata(D[:,group],axis=0),rankdata(D[:,j])))<.6:break
    group.append(j);j+=1
   if len(group)>=3:groups.append(group)
   i=max(i+1,j)
 return groups
if __name__=='__main__':
 cpcols=json.load(open('combined12/manifest.json'))['feature_sets']['CP_summary14'];meta=['construct','success','protein_length','uploaded_length','repeat_count','source_order','batch','rna','n_seeds','n_models','architecture'];numeric=[c for c in a if c not in meta];parts=[]
 for family,cols in [('CP_summary',cpcols),('repeat',[c for c in numeric if c.startswith(('P','internal_P'))]),('structural',[c for c in numeric if c not in cpcols and not c.startswith(('P','internal_P'))])]:parts.append(stats(a[cols].to_numpy(),cols,family))
 parts.append(stats(D,[f'density_slot_{i+1}' for i in range(432)],'residue'))
 ii,jj=np.triu_indices(432,12);C=CP[:,ii,jj];take=(C.max(0)>=.1)&(C.std(0)>1e-6);ii=ii[take];jj=jj[take];C=C[:,take];np.savez_compressed(O+'/contact_candidates.npz',values=C,i=ii,j=jj,names=a.construct.to_numpy());parts.append(stats(C,[f'contact_{i+1}_{j+1}' for i,j in zip(ii,jj)],'contact'))
 gs=blocks(D);cc={f'CCR_new_{i+1:02d}':group for i,group in enumerate(gs)};json.dump(cc,open(O+'/new_CCR_definitions.json','w'),indent=2)
 if gs:parts.append(stats(np.column_stack([D[:,group].mean(1) for group in gs]),list(cc),'CCR'))
 allstats=pd.concat(parts,ignore_index=True);allstats['q_global']=bh(allstats.p_exact);allstats['q_within_arch_global']=bh(allstats.p_within_architecture_exact);allstats['q_family']=allstats.groupby('family').p_exact.transform(lambda p:bh(p));allstats.to_csv(O+'/all_univariate_statistics.csv',index=False)
 for fam,t in allstats.groupby('family'):t.sort_values(['p_exact','pearson_r'],ascending=[True,False]).to_csv(O+'/univariate_'+fam+'.csv',index=False)
 allval={**{c:a[c].to_numpy() for c in numeric},**{f'density_slot_{i+1}':D[:,i] for i in range(432)},**{name:D[:,v].mean(1) for name,v in cc.items()}};directions=[]
 for name,V in allval.items():
  for group in sorted(set(g)):
   mask=g==group;pos=mask&(y==1);neg=mask&(y==0);directions.append(dict(feature=name,architecture=group,n_success=int(pos.sum()),n_failure=int(neg.sum()),delta=V[pos].mean()-V[neg].mean() if pos.any() and neg.any() else np.nan))
 pd.DataFrame(directions).to_csv(O+'/within_architecture_directions.csv',index=False)
 # Numbered reference positions are positions of the original 512-aa success reference.
 rm=pd.read_csv(O+'/all_residue_mapping.csv');ref='puf12_r123_r567_r567_r678';refm=rm[rm.construct==ref].set_index('aligned_index');top=allstats[allstats.family=='contact'].sort_values(['p_exact','pearson_r'],ascending=[True,False]).head(100).copy();ri=[]
 for t in top.itertuples():
  i,j=map(int,t.feature.split('_')[1:]);r=t._asdict();r.update(slot_i=(i-1)//36+1,slot_j=(j-1)//36+1,position_i=(i-1)%36+1,position_j=(j-1)%36+1,reference_i=int(refm.loc[i,'native_residue']),reference_j=int(refm.loc[j,'native_residue']),reference_aa_i=refm.loc[i,'aa'],reference_aa_j=refm.loc[j,'aa']);ri.append(r)
 pd.DataFrame(ri).drop(columns='Index').to_csv(O+'/top_contact_pairs_mapped.csv',index=False)
 pairs=allstats[allstats.family=='repeat'];pairs[pairs.feature.str.endswith('_CP')].sort_values('p_exact').to_csv(O+'/top_repeat_interfaces.csv',index=False)
 res=allstats[allstats.family=='residue'].copy();res['aligned_index']=res.feature.str.split('_').str[-1].astype(int);res=res.merge(refm.reset_index()[['aligned_index','native_residue','aa','slot','position','TRM','near_loop_5']],on='aligned_index');res['same_reference_aa_fraction']=[(rm[(rm.aligned_index==i)&rm.construct.isin(a.construct)].aa==refm.loc[i,'aa']).mean() for i in res.aligned_index];res.sort_values('p_exact').to_csv(O+'/top_residues_mapped.csv',index=False)
 print('tested',len(allstats),'contacts',len(C.T),'CCRs',len(cc),'unrestricted permutations',len(P),'within-architecture permutations',len(RP));print(allstats[allstats.feature.isin(['pp_nonlocal12_high_per_res','P5_P6_CP','K204_slot_density','A235_slot_density','CCR01_fixed'])].to_string(index=False));json.dump({'n_permutations':len(P),'n_restricted_permutations':len(RP),'n_tests':len(allstats),'n_contacts':len(C.T),'n_CCR':len(cc)},open(O+'/stats_manifest.json','w'),indent=2)
