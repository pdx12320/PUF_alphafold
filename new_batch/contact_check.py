import json,zipfile,glob,re,numpy as np,pandas as pd
rows=[]
for p in glob.glob('upload/*.zip'):
 z=zipfile.ZipFile(p)
 for n in z.namelist():
  if not re.search(r'_full_data_\d+\.json$',n):continue
  folder=n.split('/')[0];name,seed=re.match(r'(.*)_rna_.*_seed(\d+)$',folder).groups();d=json.loads(z.read(n));cp=np.array(d['contact_probs']);N=d['token_chain_ids'].count('A');pp=cp[:N,:N];f={'construct':name,'seed':int(seed)}
  for i,j in [(204,239),(196,235)]:f['CP_%d_%d'%(i,j)]=pp[i-1,j-1]
  for label,inds in [('CCR01',[178,179,180]),('CCR02',[126,127,128])]:f[label+'_density']=np.mean([np.sum(pp[i-1,abs(np.arange(N)-(i-1))>=12]) for i in inds])
  rows.append(f)
s=pd.DataFrame(rows).groupby(['construct','seed']).mean();s.to_csv('new_batch/local_contact_seed_values.csv');s.groupby('construct').mean().to_csv('new_batch/local_contact_means.csv')
