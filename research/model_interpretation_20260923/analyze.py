#!/usr/bin/env python3
"""Interpret frozen models and match exact substitutions to existing assays.
No fitting, threshold changes or test-panel selection is performed.
"""
from pathlib import Path
import argparse,sys,json,gzip,re
import numpy as np,pandas as pd
from joblib import load
parser=argparse.ArgumentParser();parser.add_argument('--v4',type=Path,required=True);args=parser.parse_args()
R=Path(__file__).resolve().parents[1];O=Path(__file__).resolve().parent
sys.path.insert(0,str(args.v4/'src'));import core
sys.path.insert(0,str(R/'aice_mpnn_20260922/code'));from puf12_config import PUF12_SEQ
pairs=pd.read_csv(args.v4/'cp_features/pair_dictionary.csv');mapping=pd.read_csv(args.v4/'data/residue_mapping.csv')
ref=mapping[['unified_reference_residue','WT_aa','P_position','repeat_internal_index']].drop_duplicates().sort_values('unified_reference_residue')
assert ref.unified_reference_residue.is_unique
ref=ref[(ref.unified_reference_residue>=17)&(ref.unified_reference_residue<=509)].copy();assert ''.join(ref.WT_aa)==PUF12_SEQ
ref['model3_position']=ref.unified_reference_residue-16;ref['P_position']=ref.P_position.fillna(0).astype(int);ref.to_csv(O/'reference_crosswalk.csv',index=False)
residue_rows=[];repeat_rows=[];trm_rows=[]
for task in ['C295','C388','C871']:
 pack=load(R/f'fixed_test_20260923/model2/{task}_model.joblib');ds=core.Dataset(task,'combined_PR');st=pack['state'];model=pack['classifier'];c=pack['representation'];n=len(st['selected'])
 if hasattr(model,'feature_importances_'):weights=model.feature_importances_;method='RF impurity importance'
 else:
  coeff=np.atleast_2d(model.coef_)
  if pack['embedding']:
   z=ds.original_z(st,c);basis=z[st['tr']].T@st['projection'];coeff=coeff@basis.T
   rebuilt=z@coeff.T+np.atleast_1d(model.intercept_)
   expected=np.asarray(model.decision_function(st['embedding']))
   assert np.allclose(rebuilt.ravel() if expected.ndim==1 else rebuilt,expected,atol=1e-8), "coefficient reconstruction mismatch"
  weights=np.abs(coeff).mean(axis=0);method='mean absolute coefficient on standardized original features'
 cp=weights[:n];sub=pairs.iloc[st['selected']].copy();assert (sub.block=='PR').all();sub['importance']=cp;sub['endpoint']=task;sub['method']=method;sub['cp_share']=cp/cp.sum();sub.to_csv(O/f'{task}_contact_importance.csv',index=False)
 for pos,g in sub.groupby('reference_i'):
  i=int(pos)-16;residue_rows.append(dict(endpoint=task,model2_reference_position=int(pos),model3_position=i,residue=PUF12_SEQ[i-1]+str(i),P=int(g.P_i.iloc[0]),n_selected_contacts=len(g),importance=float(g.importance.sum()),cp_share=float(g.cp_share.sum()),mean_contact_importance=float(g.importance.mean()),method=method))
 for pos,g in sub.groupby('P_i'):
  repeat_rows.append(dict(endpoint=task,P=int(pos),n_selected_contacts=len(g),importance=float(g.importance.sum()),cp_share=float(g.cp_share.sum()),mean_contact_importance=float(g.importance.mean()),method=method,cp_fraction_of_total_weight=float(cp.sum()/weights.sum())))
 if st['trmstate'] is not None:
  idx,mu,scale=st['trmstate'];cols=np.array(ds.trm.columns)[idx]
  for name,w in zip(cols,weights[n:]):trm_rows.append(dict(endpoint=task,feature=name,importance=w,total_weight_fraction=float(w/weights.sum())))
pd.DataFrame(residue_rows).sort_values(['endpoint','cp_share'],ascending=[True,False]).to_csv(O/'model2_residue_importance.csv',index=False)
pd.DataFrame(repeat_rows).sort_values(['endpoint','cp_share'],ascending=[True,False]).to_csv(O/'model2_repeat_importance.csv',index=False)
pd.DataFrame(trm_rows).to_csv(O/'model2_sequence_feature_importance.csv',index=False)
# Read all 10,000 generated sequences, excluding native by header identifier.
arrays={}
for model in ['proteinmpnn','ligandmpnn']:
 seqs=[];header='';parts=[]
 def flush():
  if parts and re.search(r'\bid=\d+|\bsample=\d+',header):seqs.append(''.join(parts))
 with gzip.open(R/f'aice_mpnn_20260922/results/mpnn_samples_{model}.fa.gz','rt') as h:
  for line in h:
   if line.startswith('>'):flush();header=line.strip();parts=[]
   else:parts.append(line.strip())
  flush()
 assert len(seqs)==10000 and all(len(s)==493 for s in seqs)
 arrays[model]=np.array([list(s) for s in seqs])
rep=[]
for P,g in ref[ref.P_position>0].groupby('P_position'):
 ix=g.model3_position.to_numpy(int)-1;trm=g[g.repeat_internal_index.isin([12,13,16])].model3_position.to_numpy(int)-1
 for model,arr in arrays.items():
  rep.append(dict(P=P,model=model,n_residues=len(ix),WT_retention=float((arr[:,ix]==np.array(list(PUF12_SEQ))[ix]).mean()),TRM_WT_retention=float((arr[:,trm]==np.array(list(PUF12_SEQ))[trm]).mean()),n_generated=10000))
pd.DataFrame(rep).to_csv(O/'model3_repeat_preferences.csv',index=False)
cons=pd.read_csv(R/'aice_mpnn_20260922/results/aice_single_ranked.csv');cons=cons[cons.consensus==True].copy();cons['mutation']=cons.wt+cons.pos_1based.astype(str)+cons.ligandmpnn_mut;cons=cons.merge(ref[['model3_position','P_position']],left_on='pos_1based',right_on='model3_position',validate='one_to_one');cons.to_csv(O/'model3_consensus_by_repeat.csv',index=False)
# Reconstruct exact experimentally tested substitutions from verified alignment.
d=core.D.copy();rows=[]
for row in d.itertuples():
 m=mapping[(mapping.construct_id==row.construct_id)&mapping.directly_mutated.astype(bool)].copy();m=m[(m.unified_reference_residue>=17)&(m.unified_reference_residue<=509)]
 idx=m.unified_reference_residue.to_numpy(int)-17;aa=m.mutant_aa.to_numpy();wt=m.WT_aa.to_numpy();assert len(idx)>0
 assert np.array_equal(wt,np.array(list(PUF12_SEQ))[idx])
 rec=dict(construct_id=row.construct_id,TRM=row.TRM,P=row.position_set,mutation=';'.join(f'{w}{i+1}{a}' for w,i,a in zip(wt,idx,aa)),n_substitutions=len(idx),all_substitutions_TRM=bool(m.TRM_slot.all()),exact_nontrm_consensus_overlap=';'.join(sorted(set(f'{w}{i+1}{a}' for w,i,a in zip(wt,idx,aa))&set(cons.mutation))))
 for task in ['C295','C388','C871']:
  rec[task+'_percent']=getattr(row,task+'_fraction')*100;rec[task+'_delta_pp']=getattr(row,task+'_delta_pp');rec[task+'_control_percent']=getattr(row,task+'_control_fraction')*100
 for model,arr in arrays.items():
  rec[model+'_exact_mutant_count']=int((arr[:,idx]==aa).all(axis=1).sum());rec[model+'_same_sites_WT_count']=int((arr[:,idx]==wt).all(axis=1).sum())
  rec[model+'_mutant_marginal_mean']=float((arr[:,idx]==aa).mean());rec[model+'_WT_marginal_mean']=float((arr[:,idx]==wt).mean())
 rec['positive_measured_example']=bool(rec['C388_delta_pp']>=-10 and rec['C295_delta_pp']<=0 and rec['C871_delta_pp']<=-10)
 rows.append(rec)
out=pd.DataFrame(rows);out.to_csv(O/'model3_exact_mutation_wetlab_comparison.csv',index=False)
print('MODEL2 P ranks\n',pd.DataFrame(repeat_rows).sort_values(['endpoint','cp_share'],ascending=[True,False]).groupby('endpoint').head(4)[['endpoint','P','cp_share']].to_string(index=False))
print('MODEL2 residues\n',pd.DataFrame(residue_rows).sort_values(['endpoint','cp_share'],ascending=[True,False]).groupby('endpoint').head(4)[['endpoint','residue','P','cp_share']].to_string(index=False))
print('MODEL3 consensus P\n',cons.groupby('P_position').mutation.apply(', '.join).to_string())
print('POSITIVE CASES\n',out[out.positive_measured_example][['construct_id','mutation','C388_delta_pp','C295_delta_pp','C871_delta_pp','proteinmpnn_exact_mutant_count','proteinmpnn_same_sites_WT_count','ligandmpnn_exact_mutant_count','ligandmpnn_same_sites_WT_count']].to_string(index=False))
print('Non-TRM exact overlaps',out.exact_nontrm_consensus_overlap.ne('').sum())
