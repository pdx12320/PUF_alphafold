from pathlib import Path
import sys,io,zipfile,json,re,shlex,collections
import numpy as np,pandas as pd
sys.path.insert(0,str(Path('sources/pdx12320/PUF_alphafold/c388_analysis').resolve()))
from cp_features import get_features
out=Path('inputs');out.mkdir(exist_ok=True)
old=Path('../work_cp_modeling');sel=pd.read_csv(old/'results/work_selection.csv');cons=sel.construct.tolist();wanted={c+'_'+s for c in cons+['wt_puf12_9'] for s in ['c295','c388','c871']}
rows=[];prseed=collections.defaultdict(list)
def parse_cif(raw):
    lines=raw.decode().splitlines();cols=[l.strip().split('.')[1] for l in lines if l.startswith('_atom_site.')];ci={x:i for i,x in enumerate(cols)};ca={};pl=collections.defaultdict(list)
    for line in lines:
        if not line.startswith('ATOM '):continue
        v=shlex.split(line) if '"' in line or "'" in line else line.split()
        if v[ci['label_asym_id']]!='A':continue
        i=int(v[ci['label_seq_id']]);pl[i].append(float(v[ci['B_iso_or_equiv']]))
        if v[ci['label_atom_id']]=='CA':ca[i]=[float(v[ci[k]]) for k in ['Cartn_x','Cartn_y','Cartn_z']]
    assert sorted(ca)==list(range(1,520))
    return np.array([ca[i] for i in range(1,520)]),np.array([np.mean(pl[i]) for i in range(1,520)])
def scan(z,source):
    for n in z.namelist():
        if '__MACOSX' in n:continue
        if n.endswith('.zip'):scan(zipfile.ZipFile(io.BytesIO(z.read(n))),source+'!'+n)
        elif n.endswith('_job_request.json'):
            req=json.loads(z.read(n))[0];name=re.sub(r'_2$','',req['name'].lower().replace('+','_plus_').replace('-','_'))
            if name not in wanted:continue
            prefix=n[:-len('job_request.json')];d=json.loads(z.read(prefix+'full_data_0.json'));cp=np.array(d['contact_probs']);pae=np.array(d['pae']);ch=np.array(d['token_chain_ids']);pi=np.where(ch=='A')[0][7:];ri=np.where(ch=='B')[0];assert len(pi)==512 and len(ri)==15
            pp=cp[np.ix_(pi,pi)];pr=cp[np.ix_(pi,ri)];pa=pae[np.ix_(pi,pi)];prseed[name].append(pr)
            # CP is identical across the five models in each task; geometry need not be.
            modelrows=[]
            for i in range(5):
                coord,pl=parse_cif(z.read(prefix+f'model_{i}.cif'));f=get_features(pp,pa,coord[7:],pl[7:],pr)
                core=np.arange(33,465);c=coord[7:][core];eig=np.linalg.eigvalsh(np.cov(c.T));f.update(plddt_core_mean=pl[7:][core].mean(),core_anisotropy=eig[-1]/eig[0]);modelrows.append(f)
            f=pd.DataFrame(modelrows).mean().to_dict();rows.append(dict(name=name,seed=str(req['modelSeeds']),source=source,**f))
for p in sorted(Path('../upload').glob('*.zip')):scan(zipfile.ZipFile(p),p.name)
assert set(prseed)==wanted
pd.DataFrame(rows).to_csv(out/'repository_features_by_seed.csv',index=False);pd.DataFrame(rows).groupby('name').mean(numeric_only=True).to_csv(out/'repository_features.csv')
np.savez_compressed(out/'pr_by_seed.npz',**{k:np.stack(v) for k,v in prseed.items()})
for n in ['work_selection.csv','mutation_position_groups.csv','model_performance.csv','out_of_fold_predictions.csv']:(out/('previous_'+n)).write_bytes((old/'results'/n).read_bytes())
for site in ['C388','C871','C295']:(out/(site+'_previous_summary.csv')).write_bytes((old/'results'/f'{site}_features_summary.csv').read_bytes())
(out/'construct_site_summary.csv').write_bytes((old/'inputs/construct_site_summary.csv').read_bytes())
print('Extracted',len(rows),'seed tasks for',len(prseed),'construct/site combinations',flush=True)
