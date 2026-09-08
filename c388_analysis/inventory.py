import zipfile,io,json,re,pathlib,hashlib,openpyxl
O=pathlib.Path('c388_analysis');units=[]
def walk(z,label):
 for n in z.namelist():
  if n.startswith('__MACOSX/'):continue
  if n.endswith('.zip'):
   walk(zipfile.ZipFile(io.BytesIO(z.read(n))),label+'!'+n)
 for n in z.namelist():
  if n.endswith('_job_request.json'):
   d=json.loads(z.read(n));d=d[0] if isinstance(d,list) else d
   seq=d['sequences'][0]['proteinChain']['sequence'];rna=d['sequences'][1]['rnaSequence']['sequence'];name=d['name'];folder=n.rsplit('/',1)[0]+'/' if '/' in n else '';full=[x for x in z.namelist() if x.startswith(folder) and re.search('_full_data_\\d+\\.json$',x) and x.rsplit('/',1)[0]==n.rsplit('/',1)[0]] if folder else [x for x in z.namelist() if '/' not in x and re.search('_full_data_\\d+\\.json$',x)]
   key=hashlib.sha256((label+'!'+n).encode()).hexdigest()[:12]; dest=O/'units'/key;dest.mkdir(parents=True,exist_ok=True)
   (dest/'request.json').write_text(json.dumps(d));
   if 'c388' in name.lower():
    for f in full:
     for src in [f,f.replace('full_data','model').replace('.json','.cif')]: (dest/pathlib.Path(src).name).write_bytes(z.read(src))
   units.append(dict(key=key,name=name,sequence=seq,rna=rna,length=len(seq),n_models=len(full),seed=d['modelSeeds'],source=label+'!'+n))
for p in pathlib.Path('upload').glob('*.zip'):walk(zipfile.ZipFile(p),p.name)
json.dump(units,open(O/'inventory.json','w'),indent=2)
w=openpyxl.load_workbook(next(pathlib.Path('upload').glob('*.xlsx')),data_only=True)
rows=[]
for r in list(w['分组核对结果'].values)[1:]: rows.append(dict(group=r[0],C295=r[1],C388=r[2],C871=r[3],C388_note=r[6]))
json.dump(rows,open(O/'experimental_rows.json','w'),ensure_ascii=False,indent=2)
print('UNITS',len(units));print('\n'.join(str((u['name'],u['length'],u['n_models'])) for u in units if 'c388' in u['name'].lower()));print('GROUPS',sorted(set(r['group'] for r in rows)))
