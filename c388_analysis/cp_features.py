import numpy as np
def get_features(cp,pae,coords,plddt,pr):
 n=len(coords);sep=np.triu(np.ones((n,n),bool),12);p=(pae+pae.T)/2;rg=np.sqrt(np.mean(np.sum((coords-coords.mean(0))**2,axis=1)))
 features={'pp_nonlocal12_high_per_res':float((cp[sep]>.5).sum()/n),'plddt_protein_mean':float(np.mean(plddt)),'pae_contact_weighted':float((p[sep]*cp[sep]).sum()/cp[sep].sum()),'Rg_full_length_normalized':float(rg/n**(1/3))}
 for gap in [4,12,24]:
  mask=np.triu(np.ones((n,n),bool),gap);features['pp_nonlocal%d_per_res'%gap]=float(cp[mask].sum()/n);features['pp_nonlocal%d_high_per_res'%gap]=float((cp[mask]>.5).sum()/n)
 features.update(pr_cp_sum=float(pr.sum()),pr_cp_per_nt=float(pr.sum()/pr.shape[1]),pr_top20_mean=float(np.sort(pr.ravel())[-20:].mean()),pr_rna_max_mean=float(pr.max(0).mean()),pr_rna_max_min=float(pr.max(0).min()),pr_rna_coverage_05=float((pr.max(0)>.5).mean()),pr_protein_coverage_05=float((pr.max(1)>.5).mean()),pr_high_contacts_per_nt=float((pr>.5).sum()/pr.shape[1]))
 return features
