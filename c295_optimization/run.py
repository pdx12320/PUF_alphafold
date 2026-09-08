import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,re,warnings,itertools
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold,LeaveOneOut
from sklearn.metrics import confusion_matrix,roc_auc_score,average_precision_score
from threadpoolctl import threadpool_limits
threadpool_limits(1);warnings.filterwarnings('ignore',category=UserWarning)
P=Path(__file__).resolve().parent;O=P/'results';O.mkdir(exist_ok=True)
df=pd.read_csv(P/'inputs/construct_site_summary.csv');df=df[df.site=='C295'].sort_values('construct').reset_index(drop=True)
cons=df.construct.tolist();y=(df.delta_editing_pp.values < -20).astype(int);assert len(y)==30 and y.sum()==10
seq=json.loads((P/'inputs/sequences.json').read_text());wt=seq['wt_puf12_9_c295']['protein'];rna=seq['wt_puf12_9_c295']['rna'];mats=np.load(P/'inputs/pr_by_seed.npz');w=mats['wt_puf12_9_c295'].mean(0)
aa='ACDEFGHIKLMNPQRSTVWY';hyd=dict(zip(aa,[1.8,2.5,-3.5,-3.5,2.8,-.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-.8,-.7,4.2,-.9,-1.3]));heavy=dict(zip(aa,[1,2,4,5,7,0,6,4,5,4,4,4,3,5,7,2,3,3,10,8]))
charge={a:1 if a in 'KR' else -1 if a in 'DE' else 0 for a in aa}
ps=[set(map(int,re.findall(r'(?:^|_plus_)p(\d+)_',c))) for c in cons]
seqrows=[];cprows=[];audit=[]
for i,c in enumerate(cons):
    protein=seq[c+'_c295']['protein'];assert len(protein)==len(wt)==519 and seq[c+'_c295']['rna']==rna
    idx=np.array([j for j,(a,b) in enumerate(zip(wt,protein)) if a!=b]);assert ';'.join(str(j+1) for j in idx)==df.mutation_positions.iloc[i]
    aold=[wt[j] for j in idx];anew=[protein[j] for j in idx]
    f={'n_aa_changes':len(idx),'n_P_positions':len(ps[i])}
    for p in range(1,13):f[f'P{p}']=int(p in ps[i])
    rs=list(map(int,re.findall(r'_r(\d+)_',c)))
    f['unparsed_R_count']=len(ps[i])-len(rs)
    for r in range(1,9):f[f'R{r}_count']=rs.count(r)
    for a in aa:f['aa_delta_count_'+a]=anew.count(a)-aold.count(a)
    for prop,vals in [('hydropathy',hyd),('charge_KR_DE',charge),('sidechain_heavy_atoms',heavy)]:
        delta=np.array([vals[b]-vals[a] for a,b in zip(aold,anew)])
        f[prop+'_delta_sum']=delta.sum();f[prop+'_abs_delta_sum']=np.abs(delta).sum()
    for a in aa:f['mutant_changed_aa_'+a]=anew.count(a)
    slot_indices={slot:[159+36*(p-4)+offset for p in sorted(ps[i])] for slot,offset in [(12,0),(13,1),(16,4)]}
    assert set(idx).issubset(set(j for jj in slot_indices.values() for j in jj))
    for slot,jj in slot_indices.items():
        for a in aa:f[f'TRM{slot}_aa_delta_'+a]=sum(protein[j]==a for j in jj)-sum(wt[j]==a for j in jj)
    
    seqrows.append(f)
    # Cached PR arrays omit input residues 1-7; retain original sequence coordinates here.
    coords=np.arange(8,520);dist=np.abs(coords[:,None]-(idx+1)[None,:]).min(1)
    masks={'global':np.ones(512,dtype=bool),'local4':dist<=4,'local18':dist<=18,'neighbor19_54':(dist>18)&(dist<=54),'distal55':dist>54}
    m=mats[c+'_c295'].mean(0);delta=m-w;f={}
    for region,mask in masks.items():
        dm=delta[mask];wm=w[mask];mm=m[mask]
        f[region+'_delta_sum']=dm.sum();f[region+'_delta_abs_sum']=np.abs(dm).sum();f[region+'_delta_rms']=np.sqrt((dm**2).mean());f[region+'_gain_sum']=np.maximum(dm,0).sum();f[region+'_loss_sum']=np.minimum(dm,0).sum();f[region+'_WT_interface_delta_sum']=dm[wm>=.1].sum()
        for nt in range(15):f[f'{region}_nt{nt+1}_{rna[nt]}_delta']=dm[:,nt].sum()
        f[region+'_delta_top3_per_residue']=np.sort(mm,axis=1)[:,-3:].sum()-np.sort(wm,axis=1)[:,-3:].sum()
    cprows.append(f)
    audit.append(dict(construct=c,positions=sorted(ps[i]),changes=[f'{wt[j]}{j+1}{protein[j]}' for j in idx],n_seeds=len(mats[c+'_c295']),window_sizes={k:int(v.sum()) for k,v in masks.items()}))
S=pd.DataFrame(seqrows,index=cons);C=pd.DataFrame(cprows,index=cons);X={'sequence':S,'local_CP':C,'sequence_plus_CP':pd.concat([S,C],axis=1)}
for k,a in X.items():a.to_csv(O/(k+'_features.csv'))
(O/'sequence_and_window_audit.json').write_text(json.dumps(audit,indent=2));df.assign(label_large_decrease=y).to_csv(O/'labels.csv',index=False)
specs=[dict(kind='L2',C=c) for c in [.01,.1,1,10]]+[dict(kind='EN',C=c,l1_ratio=l) for c in [.01,.1,1,10] for l in [.1,.5,.9]]+[dict(kind='LDA')]+[dict(kind='RF',max_depth=d,min_samples_leaf=l) for d in [1,2,3] for l in [3,5,8]]
def model(q):
    if q['kind']=='L2':m=LogisticRegression(C=q['C'],class_weight='balanced',solver='liblinear',random_state=2026)
    elif q['kind']=='EN':m=LogisticRegression(C=q['C'],l1_ratio=q['l1_ratio'],solver='saga',class_weight='balanced',max_iter=10000,random_state=2026)
    elif q['kind']=='LDA':m=LinearDiscriminantAnalysis(solver='lsqr',shrinkage='auto',priors=[.5,.5])
    else:m=RandomForestClassifier(n_estimators=100,max_depth=q['max_depth'],min_samples_leaf=q['min_samples_leaf'],max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
    return make_pipeline(VarianceThreshold(),StandardScaler(),m)
def fitpred(a,tr,te,q):
    if len(np.unique(y[tr]))<2:return np.full(len(te),y[tr].mean())
    return model(q).fit(a[tr],y[tr]).predict_proba(a[te])[:,1]
def metrics(yt,p,t,weights=None):
    yh=p>=t;tn,fp,fn,tp=confusion_matrix(yt,yh,labels=[0,1],sample_weight=weights).ravel()
    rec=tp/(tp+fn) if tp+fn else np.nan;spec=tn/(tn+fp) if tn+fp else np.nan
    return dict(BA=(rec+spec)/2,accuracy=(tn+tp)/(tn+fp+fn+tp),precision=tp/(tp+fp) if tp+fp else 0.,recall=rec,specificity=spec,TN=tn,FP=fp,FN=fn,TP=tp,AUC=roc_auc_score(yt,p,sample_weight=weights) if len(set(yt))==2 else np.nan,AP=average_precision_score(yt,p,sample_weight=weights) if sum(yt) else np.nan)
def threshold(yt,p,weights):
    values=np.unique(p);ts=np.unique(np.r_[0,.5,1.00000001,(values[:-1]+values[1:])/2]);best=None
    for t in ts:
        ba=metrics(yt,p,t,weights)['BA'];key=(-ba,abs(t-.5),t)
        if best is None or key<best[0]:best=(key,t,ba)
    return best[1:]
def purged(idx):
    folds=[]
    for pos in sorted(set.union(*(ps[j] for j in idx))):
        te=np.array([j for j in idx if pos in ps[j]]);tr=np.array([j for j in idx if pos not in ps[j]])
        assert not any(pos in ps[j] for j in tr)
        if len(tr) and len(te):folds.append((tr,te,str(pos)))
    return folds
allidx=np.arange(len(y));outer={'LOOCV':[(tr,te,str(te[0])) for tr,te in LeaveOneOut().split(y)],'Strict_position_out':purged(allidx)}
predictions=[];selected=[];splitaudit=[]
for validation,folds in outer.items():
 for fold,(tr,te,group) in enumerate(folds):
    if validation=='LOOCV':inner=[(tr[a],tr[b]) for a,b in StratifiedKFold(3,shuffle=True,random_state=2026).split(tr,y[tr])]
    else:inner=[(a,b) for a,b,g in purged(tr)]
    ids=np.concatenate([b for a,b in inner]);counts=np.bincount(ids,minlength=len(y));weights=1/counts[ids]
    candidates=[]
    for family,frame in X.items():
        a=frame.values
        choices=[]
        for qi,q in enumerate(specs):
            pp=np.concatenate([fitpred(a,ta,tb,q) for ta,tb in inner]);t,ba=threshold(y[ids],pp,weights)
            choices.append(dict(family=family,q=qi,threshold=float(t),inner_BA=float(ba)))
        winner=sorted(choices,key=lambda r:(-r['inner_BA'],r['q']))[0];candidates.append(winner)
    joint=sorted(candidates,key=lambda r:(-r['inner_BA'],list(X).index(r['family']),r['q']))[0]
    for winner in candidates+[dict(joint,selection_family='nested_family_selection')]:
        family=winner['family'];q=specs[winner['q']];p=fitpred(X[family].values,tr,te,q);label=winner.get('selection_family',family)
        selected.append(dict(validation=validation,fold=group,reported_family=label,**winner,parameters=json.dumps(q)))
        for i,pp in zip(te,p):predictions.append(dict(validation=validation,fold=group,family=label,construct=cons[i],y=int(y[i]),score=float(pp),threshold=winner['threshold'],prediction=int(pp>=winner['threshold']),delta_editing_pp=df.delta_editing_pp.iloc[i]))
    for label,pp in [('always_no_large_decrease',np.zeros(len(te))),('always_large_decrease',np.ones(len(te)))]:
        for i,p in zip(te,pp):predictions.append(dict(validation=validation,fold=group,family=label,construct=cons[i],y=int(y[i]),score=float(p),threshold=.5,prediction=int(p>=.5),delta_editing_pp=df.delta_editing_pp.iloc[i]))
    splitaudit.append(dict(validation=validation,fold=group,train=[cons[j] for j in tr],test=[cons[j] for j in te]))
    pd.DataFrame(predictions).to_csv(O/'outer_predictions.csv',index=False);pd.DataFrame(selected).to_csv(O/'fold_model_selection.csv',index=False)
    print('DONE',validation,group,[(z['family'],z['inner_BA'],specs[z['q']]) for z in candidates],flush=True)
pr=pd.DataFrame(predictions);rows=[]
for (validation,family),g in pr.groupby(['validation','family']):
    weights=1/g.construct.map(g.construct.value_counts()).values
    # Score thresholds vary by fold; calculate classification metrics from stored binary decisions.
    m=metrics(g.y.values,g.prediction.values,.5,weights)
    m['AUC']=roc_auc_score(g.y,g.score,sample_weight=weights);m['AP']=average_precision_score(g.y,g.score,sample_weight=weights)
    rows.append(dict(validation=validation,family=family,n_unique=g.construct.nunique(),n_predictions=len(g),**m))
pd.DataFrame(rows).to_csv(O/'performance.csv',index=False)
position=[]
for (fold,family),g in pr[pr.validation=='Strict_position_out'].groupby(['fold','family']):position.append(dict(position=fold,family=family,n=len(g),**metrics(g.y.values,g.prediction.values,.5)))
pd.DataFrame(position).to_csv(O/'per_position_performance.csv',index=False)
(O/'split_audit.json').write_text(json.dumps(splitaudit,indent=2));(O/'candidate_models.json').write_text(json.dumps(specs,indent=2))
print(pd.DataFrame(rows).to_string(index=False),flush=True)
