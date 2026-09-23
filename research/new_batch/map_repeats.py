import zipfile,glob,json,re,pandas as pd,numpy as np
old=pd.read_csv('data/architecture_inputs/repeat_templates.csv');templates={}
for _,r in old.sort_values('source_identity_excluding_TRM',ascending=False).iterrows():templates.setdefault(int(r.source_repeat_inferred),r.core_sequence)
seen={};rows=[]
for p in glob.glob('upload/folds_2026_09_07_02_*.zip'):
 z=zipfile.ZipFile(p)
 for n in z.namelist():
  if not n.endswith('job_request.json'):continue
  name=re.sub('_rna_.*','',n.split('/')[0])
  if name in seen:continue
  seq=json.loads(z.read(n))[0]['sequences'][0]['proteinChain']['sequence'];removed=set()
  for m in re.finditer('MNDGPHS',seq):removed.update(range(m.start(),m.end()))
  keep=[i for i in range(len(seq)) if i not in removed];clean=''.join(seq[i] for i in keep);anchors=[m.start() for m in re.finditer('[SCN][YR][FV][IV][EQR]',clean)];assert anchors==[44+36*i for i in range(len(anchors))]
  sources=[]
  for j,anchor in enumerate(anchors):
   start=anchor-11;inds=keep[start:start+36];s=clean[start:start+36];sim={k:sum(s[i]==t[i] for i in range(36) if i not in [11,12,15])/33 for k,t in templates.items()};best=max(sim,key=sim.get);sources.append(best);rows.append({'construct':name,'design_id':int(name.split('_')[0]),'repeat_slot':j+1,'source_R_inferred':best,'identity_excluding_TRM':sim[best],'native_start':inds[0]+1,'native_end':inds[-1]+1,'TRM':s[11]+s[12]+s[15],'core_sequence':s})
  seen[name]={'repeat_order':sources,'first6_match_old_training':[1,2,3,5,6,7]==sources[:6]}
pd.DataFrame(rows).to_csv('new_batch/repeat_mapping.csv',index=False);json.dump(seen,open('new_batch/architecture_audit.json','w'),indent=2)
print({k:v['repeat_order'] for k,v in seen.items()})
