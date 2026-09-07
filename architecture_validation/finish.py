import json,os
import numpy as np,pandas as pd,joblib
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O='architecture_validation';a=pd.read_csv(O+'/main_features.csv');s=pd.read_csv(O+'/all_univariate_statistics.csv');m=pd.read_csv(O+'/validation_metrics.csv');p=pd.read_csv(O+'/heldout_predictions.csv');imp=pd.read_csv(O+'/fold_importance.csv');rm=pd.read_csv(O+'/all_residue_mapping.csv');z=np.load(O+'/all_arrays.npz');ix=[list(z['names']).index(n) for n in a.construct];D=z['density'][ix];ref=rm[rm.construct=='puf12_r123_r567_r567_r678'].set_index('aligned_index');cpcols=json.load(open('combined12/manifest.json'))['feature_sets']['CP_summary14'];defs={}
for gap in [4,12,24]:
 defs[f'pp_nonlocal{gap}_per_res']=f'sum CP_ij / L, protein i<j and j-i>={gap}'
 defs[f'pp_nonlocal{gap}_high_per_res']=f'count(CP_ij>0.5) / L, protein i<j and j-i>={gap}'
defs.update(pr_cp_sum='sum of all protein-RNA CP_ij',pr_cp_per_nt='protein-RNA CP sum / RNA nucleotide count',pr_top20_mean='mean of 20 largest protein-RNA CP values',pr_rna_max_mean='mean across RNA nucleotides of maximum CP over protein residues',pr_rna_max_min='minimum across RNA nucleotides of maximum CP over protein residues',pr_rna_coverage_05='fraction of RNA nucleotides with maximum protein CP>0.5',pr_protein_coverage_05='fraction of protein residues with maximum RNA CP>0.5',pr_high_contacts_per_nt='count protein-RNA CP>0.5 / RNA nucleotide count',plddt_protein_mean='mean full-protein per-residue pLDDT; per-residue mean of atom B-factors',plddt_core_mean='mean construct-averaged residue pLDDT over mapped repeat cores',plddt_core_min='minimum construct-averaged residue pLDDT over mapped repeat cores (not mean of per-model minima)',global_core_PAE='mean symmetric PAE over all unique mapped core residue pairs i<j',pae_contact_weighted='sum symmetric PAE_ij*CP_ij / sum CP_ij, full protein i<j and sequence gap>=12; average models then seeds',Rg='sqrt(mean squared CA distance to CA centroid), full protein; Angstrom',Rg_full_length_normalized='full-protein Rg / L^(1/3); Angstrom',anisotropy='largest/smallest eigenvalue of mapped core CA coordinate covariance; dimensionless',compactness='L/(4*pi*construct-mean Rg^3/3); geometry proxy, not physical volume occupancy',HC_nonlocal_density='count(core mean CP>0.5, symmetric mean PAE<5 A, both mean pLDDT>=80, native gap>=12) / mapped core length; thresholds applied after aggregation',adj_CP_mean='mean CP of adjacent repeat interfaces, equal weight per interface',adj_CP_min='minimum CP of adjacent repeat interfaces',nonadj_CP_mean='mean CP of nonadjacent repeat interfaces, equal weight per interface',adj_PAE_mean='mean symmetric PAE of adjacent interfaces',internal_CP_mean='mean of repeat internal densities',K204_slot_density='full-protein gap>=12 summed CP at aligned slot171 (reference K204; P5 position27)',A235_slot_density='full-protein gap>=12 summed CP at aligned slot202 (reference A235; P6 position22)',CCR01_fixed='mean full-protein nonlocal density at aligned145-147; original reference residues178-180, P5 positions1-3')
for c in s.feature:
 if c.startswith('internal_P'):defs[c]='sum intrarepeat CP for unique native-gap>=4 residue pairs / 36 mapped residues'
 elif c.startswith('P') and c.endswith('_CP'):defs[c]='mean cross-repeat CP over residue pairs with native sequence gap>=4'
 elif c.startswith('P') and c.endswith('_PAE'):defs[c]='mean symmetric cross-repeat PAE for same gap>=4 interface mask; Angstrom'
summary=s[s.feature.isin(defs)].copy();summary['definition']=summary.feature.map(defs);out=[]
for (val,fs,f),g in imp.groupby(['validation','feature_set','feature']):out.append(dict(validation=val,feature_set=fs,feature=f,mean_fold_MDI=g.MDI.mean(),sd_fold_MDI=g.MDI.std(),heldout_mean_logloss_increase=g.heldout_logloss_increase_sum.sum()/24))
i=pd.DataFrame(out);i.to_csv(O+'/feature_importance_summary.csv',index=False)
for val in ['LOCO','architecture_out']:
 t=i[(i.validation==val)&(i.feature_set=='CP_summary14')][['feature','mean_fold_MDI','heldout_mean_logloss_increase']].rename(columns={'mean_fold_MDI':val+'_RF_MDI','heldout_mean_logloss_increase':val+'_heldout_loss_increase'});summary=summary.merge(t,on='feature',how='left')
summary.to_csv(O+'/feature_dictionary_ranges_importance.csv',index=False);summary[summary.feature.isin(cpcols)].to_csv(O+'/CP_summary_dictionary.csv',index=False)
# CCR positions, coherence, and annotations.
from scipy.stats import spearmanr
cc=json.load(open(O+'/new_CCR_definitions.json'));cc={'CCR01_fixed':[144,145,146],**cc};reg=[]
for name,inds in cc.items():
 v=D[:,inds];R=np.corrcoef(v,rowvar=False);S=spearmanr(v,axis=0).statistic;mask=np.triu_indices(len(inds),1);r=dict(feature=name,aligned_residues=','.join(str(x+1) for x in inds),reference_residues=','.join(str(ref.loc[x+1,'native_residue']) for x in inds),repeat_slot=inds[0]//36+1,repeat_positions=','.join(str(x%36+1) for x in inds),min_internal_Pearson=float(R[mask].min()),min_internal_Spearman=float(S[mask].min()));reg.append(r)
pd.DataFrame(reg).merge(s,on='feature').to_csv(O+'/CCR_regions_annotated.csv',index=False)
top=pd.read_csv(O+'/top_contact_pairs_mapped.csv');native=[];directions=[];allcp=z['cp'][ix];lookup=rm[rm.construct.isin(a.construct)].set_index(['construct','aligned_index'])
for r in top.itertuples():
 ai,aj=map(int,r.feature.split('_')[1:]);v=allcp[:,ai-1,aj-1]
 for g in sorted(a.architecture.unique()):
  pos=(a.architecture==g)&(a.success==1);neg=(a.architecture==g)&(a.success==0);directions.append(dict(feature=r.feature,architecture=g,n_success=int(pos.sum()),n_failure=int(neg.sum()),delta=float(v[pos].mean()-v[neg].mean()) if pos.any() and neg.any() else np.nan))
 for j,t in a.iterrows():
  aa=lookup.loc[(t.construct,ai)];bb=lookup.loc[(t.construct,aj)];native.append(dict(feature=r.feature,construct=t.construct,architecture=t.architecture,native_i=int(aa.native_residue),native_j=int(bb.native_residue),aa_i=aa.aa,aa_j=bb.aa,slot_i=aa.slot,slot_j=bb.slot,source_i=aa.source,source_j=bb.source,CP=float(v[j])))
pd.DataFrame(native).to_csv(O+'/top_contact_native_mapping.csv',index=False);pd.DataFrame(directions).to_csv(O+'/top_contact_architecture_directions.csv',index=False)
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):fig.savefig(O+'/'+name+'.png',dpi=300,facecolor='white');fig.savefig(O+'/'+name+'.pdf',facecolor='white');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,5.4),layout='constrained')
for ax,fs in zip(axes,['CP_density3','CP_structural']):
 t=i[(i.validation=='architecture_out')&(i.feature_set==fs)].sort_values('mean_fold_MDI');ax.barh(t.feature,t.mean_fold_MDI,color='#0072B2');ax.set(title=fs+' RF',xlabel='Mean architecture-fold impurity importance');ax.set_xlim(left=0)
fig.suptitle('RF feature importance — training importance; correlated features share credit');save(fig,'feature_importance')
keys=['CP_summary14','CP_density3','structural_only','CP_structural','CP_interface','full','S12_only'];fig,axes=plt.subplots(1,2,figsize=(11,5.2),layout='constrained')
for ax,metric in zip(axes,['AUC','balanced_accuracy']):
 for j,val in enumerate(['LOCO','architecture_out']):
  t=m[(m.subset=='PUF12_main')&(m.model=='RF')&(m.validation==val)].set_index('feature_set').loc[keys];ax.plot(t[metric],np.arange(len(keys))+(.1 if j else -.1),'s' if j else 'o',color=['#0072B2','#D55E00'][j],label=val)
 ax.set(yticks=range(len(keys)),yticklabels=keys,xlim=(0,1.03),xlabel=metric);ax.axvline(.5,ls=':',color='gray');ax.legend(loc='lower left',fontsize=8)
fig.suptitle('PUF12 RF: construct-out vs architecture-out (source-order-out is identical)');save(fig,'validation_comparison')
idmap=pd.read_csv(O+'/plot_ID_mapping.csv');fig,axes=plt.subplots(1,3,figsize=(12,8.5),layout='constrained')
for ax,fs in zip(axes,['CP_density3','CP_structural','structural_only']):
 for j,val in enumerate(['LOCO','architecture_out']):
  t=p[(p.subset=='PUF12_main')&(p.model=='RF')&(p.feature_set==fs)&(p.validation==val)].set_index('construct').loc[idmap.construct];ax.scatter(t.score,np.arange(len(t))+(.13 if j else -.13),marker='s' if j else 'o',color=['#0072B2','#D55E00'][j],s=27,label=val)
 ax.set(xlim=(-.03,1.03),yticks=range(len(idmap)),yticklabels=[f'{r.ID} {r.architecture} '+('S' if r.success else 'F') for r in idmap.itertuples()],title=fs,xlabel='Held-out RF score');ax.axvline(.5,ls=':',color='black');ax.invert_yaxis();ax.legend(loc='lower right',fontsize=7)
fig.suptitle('Held-out predictions — S: experimental success; F: failure; threshold 0.5');save(fig,'heldout_probabilities')
# Simple frozen ranking formula and re-useable models.
model=joblib.load(O+'/S12_only_LR.joblib');sc=model[1];lr=model[-1];formula={'features':['pp_nonlocal12_high_per_res'],'mean':sc.mean_.tolist(),'scale':sc.scale_.tolist(),'coefficient':lr.coef_[0].tolist(),'intercept':float(lr.intercept_[0]),'threshold':.5,'role':'minimal ranking baseline, low precision at default threshold; not calibrated'};json.dump(formula,open(O+'/minimal_logistic_formula.json','w'),indent=2)
# Sensitivity of pooled heldout AUC to each individual positive; no retraining, not CI.
from sklearn.metrics import roc_auc_score
rows=[]
for fs in ['CP_density3','CP_structural','full','S12_only']:
 for typ in ['RF','LR']:
  for val in ['LOCO','architecture_out']:
   t=p[(p.subset=='PUF12_main')&(p.model==typ)&(p.feature_set==fs)&(p.validation==val)]
   for c in t[t.success==1].construct:
    q=t[t.construct!=c];rows.append(dict(feature_set=fs,model=typ,validation=val,removed_success=c,AUC=roc_auc_score(q.success,q.score)))
pd.DataFrame(rows).to_csv(O+'/AUC_remove_one_success_score_sensitivity.csv',index=False)
print('Finished outputs',formula)
