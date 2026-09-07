import os,json,itertools
import numpy as np,pandas as pd
from common import pearson_cols,bh,discover_ccr
O=os.path.dirname(__file__);z=np.load(O+'/construct_arrays.npz');names=z['names'].tolist();label=pd.read_csv('previous/results/construct_features.csv').set_index('construct');y=label.loc[names,'success'].values;ref='puf12_r123_r567_r567_r678';mapping=json.load(open(O+'/core_mapping.json'));ann=pd.read_csv(O+'/residue_mapping_all_constructs.csv');repeat=pd.read_csv(O+'/repeat_mapping.csv');refann=ann[ann.construct==ref].set_index('aligned_index')

def stats(X,y):
 X=np.asarray(X,float);a=X[y==1];b=X[y==0];delta=a.mean(0)-b.mean(0);r=pearson_cols(X,y);pool=np.sqrt(((len(a)-1)*a.var(0,ddof=1)+(len(b)-1)*b.var(0,ddof=1))/(len(y)-2));g=np.divide(delta,pool,out=np.full(len(delta),np.nan),where=pool>1e-10)*(1-3/(4*(len(y)-2)-1));comb=list(itertools.combinations(range(len(y)),int(sum(y))));perm=np.zeros((len(comb),len(y)))
 for k,c in enumerate(comb):perm[k,list(c)]=1
 # Equality tolerance preserves exact permutation ties.
 null=perm@X/sum(y)-(1-perm)@X/(len(y)-sum(y));p=(abs(null)>=abs(delta)[None,:]-1e-9).mean(0)
 return {'delta':delta,'pearson_r':r,'hedges_g':g,'p_exact':p,'q_BH':bh(p)}
# Primary tests only common aligned residue pairs with sequence separation >=12;
# fixed label-free contact filter max construct mean CP >=0.1 and nonzero variation.
for suffix,subset in [('all14',np.ones(len(y),bool)),('12repeat',np.array([not n.startswith('puf_11') for n in names]))]:
 cp=z['cp'][subset];yy=y[subset];i,j=np.triu_indices(432,k=12);X=cp[:,i,j];valid=np.isfinite(X).all(0)&(np.nanmax(X,axis=0)>=.1)&(np.nanstd(X,axis=0)>1e-6);i=i[valid];j=j[valid];X=X[:,valid];ss=stats(X,yy);df=pd.DataFrame({'aligned_i':i+1,'aligned_j':j+1,'reference_i':[mapping[ref][k]+1 for k in i],'reference_j':[mapping[ref][k]+1 for k in j],'repeat_i':i//36+1,'repeat_j':j//36+1,'repeat_pos_i':i%36+1,'repeat_pos_j':j%36+1,'success_mean_CP':X[yy==1].mean(0),'failure_mean_CP':X[yy==0].mean(0),**ss});df['direction']=np.where(df.delta>0,'enhanced','reduced');df['inter_repeat']=df.repeat_i!=df.repeat_j;df.sort_values(['q_BH','p_exact','delta'],ascending=[True,True,False]).to_csv(O+'/contact_tests_'+suffix+'.csv',index=False)
 print(suffix,'pairs',len(df),'BH significant',int((df.q_BH<.05).sum()),'min q',df.q_BH.min(),flush=True)
# Residue statistics on complete common slots; density is total full-protein nonlocal CP sum.
D=z['density'];valid=np.isfinite(D).all(0)&(np.nanstd(D,axis=0)>1e-6);idx=np.where(valid)[0];ss=stats(D[:,valid],y);rd=pd.DataFrame({'aligned_index':idx+1,'reference_residue':[mapping[ref][i]+1 for i in idx],'repeat_slot':idx//36+1,'repeat_position':idx%36+1,'success_density':D[y==1][:,idx].mean(0),'failure_density':D[y==0][:,idx].mean(0),**ss});rd=rd.merge(refann[['aa','TRM','loop_near_5aa','inter_repeat_interface','rna_contact_region','mean_plddt']],left_on='aligned_index',right_index=True);rd.sort_values('pearson_r',key=abs,ascending=False).to_csv(O+'/residue_associations.csv',index=False)
ccrs=discover_ccr(D,y,max_regions=100);out=[];definitions=[]
for k,g in enumerate(ccrs):
 inds=g['indices'];native=[mapping[ref][i]+1 for i in inds];v=D[:,inds].mean(1);ss=stats(v[:,None],y);ra=refann.loc[np.array(inds)+1];row={'CCR':'CCR%02d'%(k+1),'aligned_start':inds[0]+1,'aligned_end':inds[-1]+1,'reference_residues':','.join(map(str,native)),'repeat_slot':inds[0]//36+1,'repeat_positions':','.join(str(i%36+1) for i in inds),'n_residues':len(inds),'success_mean_density':v[y==1].mean(),'failure_mean_density':v[y==0].mean(),'pearson_r':g['r'],'min_internal_r':g['min_internal_r'],'mean_internal_r':g['mean_internal_r'],'hedges_g':ss['hedges_g'][0],'p_exact_descriptive':ss['p_exact'][0],'interface_fraction_reference':ra.inter_repeat_interface.mean(),'RNA_contact_fraction_reference':ra.rna_contact_region.mean(),'TRM_residues_reference':','.join(str(r) for r in ra[ra.TRM].native_residue),'loop_near_fraction_reference':ra.loop_near_5aa.mean()};out.append(row);definitions.append({'CCR':row['CCR'],**g})
pd.DataFrame(out).to_csv(O+'/CCR_discovery.csv',index=False);json.dump(definitions,open(O+'/CCR_definitions.json','w'),indent=2)
# Within-backbone centering: compare constructs sharing inferred source-repeat order.
family={n:','.join(map(str,repeat[repeat.construct==n].sort_values('repeat_slot').source_repeat_inferred)) for n in names};families=np.array([family[n] for n in names]);yc=y.astype(float).copy();dc=D.copy()
for f in np.unique(families):
 take=families==f;yc[take]-=yc[take].mean();dc[take]-=np.nanmean(dc[take],axis=0)
wr=pearson_cols(dc[:,valid],yc);pd.DataFrame({'aligned_index':idx+1,'within_architecture_pearson_r':wr}).to_csv(O+'/within_architecture_residue_correlation.csv',index=False)
json.dump(family,open(O+'/architecture_groups.json','w'),indent=2)
print(pd.DataFrame(out).head(10).to_string(index=False),flush=True)
