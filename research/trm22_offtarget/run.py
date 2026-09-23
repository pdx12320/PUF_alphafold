"""Frozen 22-construct off-target benchmark. Run from any directory."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import sys,json,re,io,contextlib,warnings,hashlib,subprocess
import numpy as np
import pandas as pd
import sklearn
import scipy.special
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier,ExtraTreesClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold,GroupKFold,LeaveOneOut
from sklearn.metrics import confusion_matrix,roc_auc_score,average_precision_score
from joblib import Parallel,delayed
from threadpoolctl import threadpool_limits
threadpool_limits(1)
P=Path(__file__).resolve().parent; O=P/'results'; O.mkdir(exist_ok=True)
sys.path.insert(0,str(P/'sources'))
from ContactSeek.ContactExtraction import query_cp_with_top_n
from ContactSeek.FindContactResidue import find_contact_residues
from ContactSeek.CCRegionFinding import find_consensus_contact_regions
MODELS=['LR_0.01','LR_0.1','LR_1','RF_depth1','RF_depth2','ExtraTrees_depth2','LDA_shrinkage','GaussianNB','SVM_RBF']
FAMILIES=['Total_CP','Interface_RMS','CP_summary','Repository','Local_deltaCP','Sequence','Sequence_Local','ContactSeek_CCR']
def classifier(name):
    if name.startswith('LR_'): m=LogisticRegression(C=float(name.split('_')[1]),class_weight='balanced',solver='liblinear',random_state=2026)
    elif name.startswith('RF_'): m=RandomForestClassifier(n_estimators=100,max_depth=int(name[-1]),min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
    elif name.startswith('Extra'): m=ExtraTreesClassifier(n_estimators=100,max_depth=2,min_samples_leaf=3,max_features=1.,class_weight='balanced',random_state=2026,n_jobs=1)
    elif name=='LDA_shrinkage': m=LinearDiscriminantAnalysis(solver='lsqr',shrinkage='auto',priors=[.5,.5])
    elif name=='GaussianNB': m=GaussianNB(priors=[.5,.5],var_smoothing=1e-5)
    else: m=SVC(C=1,kernel='rbf',gamma='scale',class_weight='balanced')
    return make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),m)
def fitpredict(x,y,tr,te,name):
    if len(np.unique(y[tr]))<2:return np.full(len(te),float(y[tr][0]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=RuntimeWarning)
        m=classifier(name).fit(x[tr],y[tr])
        return m.predict_proba(x[te])[:,1] if hasattr(m,'predict_proba') else scipy.special.expit(m.decision_function(x[te]))
def metrics(y,score,pred,w=None):
    tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1],sample_weight=w).ravel()
    rec=tp/(tp+fn) if tp+fn else np.nan; sp=tn/(tn+fp) if tn+fp else np.nan
    return dict(BA=(rec+sp)/2,accuracy=(tn+tp)/(tn+fp+fn+tp),precision=tp/(tp+fp) if tp+fp else 0.,recall=rec,specificity=sp,TP=tp,TN=tn,FP=fp,FN=fn,AUC=roc_auc_score(y,score,sample_weight=w) if len(set(y))==2 else np.nan,AP=average_precision_score(y,score,sample_weight=w) if sum(y) else np.nan)
def threshold(y,p,w):
    if len(set(y))<2:return .5,np.nan
    u=np.unique(p);ts=np.unique(np.r_[0,.5,np.nextafter(1.,2.),(u[:-1]+u[1:])/2])
    pred=p[:,None]>=ts; a=(w*y)[:,None]; b=(w*(1-y))[:,None]
    ba=.5*((a*pred).sum(0)/a.sum()+(b*~pred).sum(0)/b.sum())
    j=min(range(len(ts)),key=lambda j:(-ba[j],abs(ts[j]-.5),ts[j]))
    return float(ts[j]),float(ba[j])
def setup(site):
    cohort=pd.read_csv(P/'inputs/frozen_cohort.csv');cons=sorted(cohort.construct_key)
    d=pd.read_csv(P/'inputs/construct_site_summary.csv');d=d[d.site==site].set_index('construct').loc[cons].reset_index()
    assert len(d)==22 and d.construct.nunique()==22
    rep=pd.read_csv(P/'inputs/experimental_replicates.csv').groupby('construct')[site].mean()*100
    assert np.allclose(d.editing_pct,rep.loc[cons]) and np.allclose(d.delta_editing_pp,d.editing_pct-d.wt_editing_pct)
    mats=np.load(P/'inputs/pr_by_seed.npz');wt=mats['wt_puf12_9_'+site.lower()].mean(0)
    seq=json.loads((P/'inputs/sequences.json').read_text());wp=seq['wt_puf12_9_c295']['protein']
    raw=[];diff=[];local=[];sequence=[];audit=[];rms=[];tot=[]
    for c in cons:
        a=mats[c+'_'+site.lower()];assert a.shape[1:]==(512,15)
        qs=[query_cp_with_top_n(m.copy(),list(range(1,513)),[(0,5),(5,10),(10,15)],top_n=3,ref_array=wt.copy()) for m in a]
        raw.append(np.mean([q['top_n_input_prob'] for q in qs],0));diff.append(np.mean([q['diff_prob'] for q in qs],0))
        m=a.mean(0);dm=m-wt;rms.append(np.sqrt((dm**2).sum()/(519*15)));tot.append(dm.sum())
        protein=seq[c+'_c295']['protein'];idx=np.array([j+1 for j,(x,z) in enumerate(zip(wp,protein)) if x!=z]);assert len(protein)==519
        assert ';'.join(map(str,idx))==str(d.set_index('construct').loc[c,'mutation_positions'])
        distance=np.abs(np.arange(8,520)[:,None]-idx[None,:]).min(1)
        f={}
        for region,mask in {'global':distance>=0,'local4':distance<=4,'local18':distance<=18,'neighbor19_54':(distance>18)&(distance<=54),'distal55':distance>=55}.items():
            z=dm[mask];w=wt[mask]
            for name,val in {'sum':z.sum(),'abs_sum':abs(z).sum(),'rms':np.sqrt((z*z).mean()),'gain':np.maximum(z,0).sum(),'loss':np.minimum(z,0).sum(),'WT_interface_sum':z[w>=.1].sum(),'top3_delta':np.sort(m[mask],axis=1)[:,-3:].sum()-np.sort(w,axis=1)[:,-3:].sum()}.items():f[region+'_'+name]=val
            for nt in range(15):f[f'{region}_nt{nt+1}']=z[:,nt].sum()
        local.append(f)
        ps=set(map(int,re.findall(r'(?:^|_plus_)p(\d+)_',c)))
        sf={'n_mutations':len(idx),'n_positions':len(ps)}
        for pos in range(4,9):sf[f'P{pos}']=int(pos in ps)
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            sf['delta_'+aa]=sum(protein[j-1]==aa for j in idx)-sum(wp[j-1]==aa for j in idx)
            for slot,off in [(12,0),(13,1),(16,4)]:sf[f'TRM{slot}_{aa}']=sum((protein[159+36*(p-4)+off]==aa)-(wp[159+36*(p-4)+off]==aa) for p in ps)
        sequence.append(sf);audit.append(dict(construct=c,site=site,seeds=len(a),mutation_residues=idx.tolist(),positions=sorted(ps)))
    raw=np.array(raw);diff=np.array(diff);regional=np.concatenate([raw,diff],2)
    repo=pd.read_csv(P/'inputs/repository_features.csv',index_col=0);rx=repo.loc[[c+'_'+site.lower() for c in cons]].copy();rx.index=cons;rx-=repo.loc['wt_puf12_9_'+site.lower()].values
    loc=pd.DataFrame(local,index=cons);sq=pd.DataFrame(sequence,index=cons)
    frames={'Total_CP':pd.DataFrame({'delta_sum':tot},index=cons),'Interface_RMS':pd.DataFrame({'rms_519_denominator':rms},index=cons),'CP_summary':rx[[c for c in rx if c.startswith(('pr_','pp_nonlocal'))]],'Repository':rx,'Local_deltaCP':loc,'Sequence':sq,'Sequence_Local':pd.concat([sq,loc],axis=1)}
    # Source's interface statistic averages per-seed RMS; ours applies RMS after seed averaging.
    # This preserves the C388 formula for single seeds and is explicitly documented for multiseed cases.
    for k,v in frames.items():v.to_csv(O/f'{site}_{k}_features.csv')
    y=(d.delta_editing_pp.values<({'C871':-25,'C295':-20}[site])).astype(int)
    d['label_large_decrease']=y;d.to_csv(O/f'{site}_labels.csv',index=False)
    return cons,d,y,{k:v.values for k,v in frames.items()},raw,diff,regional,audit
def strict(idx,positions):
    return [(np.array([j for j in idx if p not in positions[j]],int),np.array([j for j in idx if p in positions[j]],int),f'P{p}') for p in sorted(set.union(*(positions[j] for j in idx))) if any(p not in positions[j] for j in idx)]
def worker(site,label_threshold,cons,y,static,raw,diff,regional,positions,groups,validation,tr,te,fold):
    threadpool_limits(1)
    if validation=='LOOCV':
        n=min(3,int(np.bincount(y[tr],minlength=2).min()))
        inner=[(tr[a],tr[b],str(k)) for k,(a,b) in enumerate(StratifiedKFold(n,shuffle=True,random_state=2026).split(tr,y[tr]))] if n>=2 else [(tr[a],tr[b],str(k)) for k,(a,b) in enumerate(LeaveOneOut().split(tr))]
    elif validation=='Position_group_out':inner=[(tr[a],tr[b],str(k)) for k,(a,b) in enumerate(GroupKFold(min(3,len(set(groups[tr])))).split(tr,y[tr],groups[tr]))]
    else:inner=strict(tr,positions)
    ids=np.concatenate([b for a,b,g in inner]);counts=np.bincount(ids,minlength=len(y));weights=1/counts[ids]
    fc={};regionlog=[]
    def features(train):
        key=tuple(train)
        if key not in fc:
            with warnings.catch_warnings(),contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter('ignore');keep=find_contact_residues(list(raw[train]),list(diff[train]),min_cp_threshold=.15,min_diff_threshold=.1,verbose=False)
                if keep.sum()>1:_,_,regions=find_consensus_contact_regions([dict(y_g3=np.zeros(len(train)),cp_raw_cas_nuc=list(raw[train]))],keep,correlation_threshold=.6,band_width=7,max_merge_iterations=10,verbose=False,protein_name='PUF')
                else:regions=[dict(positions=np.where(keep)[0])] if keep.any() else []
            ccr=np.concatenate([regional[:,r['positions']].mean(1) for r in regions],1) if regions else np.zeros((len(y),1))
            fc[key]=dict(static,ContactSeek_CCR=ccr)
            regionlog.append(dict(train=[cons[j] for j in train],regions=[list(map(int,r['positions'])) for r in regions]))
        return fc[key]
    for a,b,g in inner:features(a)
    features(tr);pred=[];choices=[]
    for family in FAMILIES:
        for name in MODELS:
            ip=np.concatenate([fitpredict(features(a)[family],y,a,b,name) for a,b,g in inner]);t,ba=threshold(y[ids],ip,weights)
            pp=fitpredict(features(tr)[family],y,tr,te,name)
            choices.append(dict(family=family,model=name,threshold=t,inner_BA=ba,score=pp))
            for policy,th in [('fixed_0.5',.5),('training_only_threshold',t)]:
                for j,p in zip(te,pp):pred.append(dict(site=site,validation=validation,fold=fold,family=family,model=name,policy=policy,construct=cons[j],y=int(y[j]),score=float(p),threshold=th,prediction=int(p>=th)))
    chosen=[]
    for family in FAMILIES+['Nested_all']:
        pool=choices if family=='Nested_all' else [c for c in choices if c['family']==family]
        win=min(pool,key=lambda c:(-c['inner_BA'],FAMILIES.index(c['family']),MODELS.index(c['model'])))
        chosen.append(dict(site=site,validation=validation,fold=fold,reported_family=family,**{k:v for k,v in win.items() if k!='score'}))
        for j,p in zip(te,win['score']):pred.append(dict(site=site,validation=validation,fold=fold,family=family,model='nested_model_selection',policy='training_only_threshold',construct=cons[j],y=int(y[j]),score=float(p),threshold=win['threshold'],prediction=int(p>=win['threshold'])))
    for value in [0,1]:
        for j in te:pred.append(dict(site=site,validation=validation,fold=fold,family='Baseline',model=f'always_{value}',policy='fixed_0.5',construct=cons[j],y=int(y[j]),score=float(value),threshold=.5,prediction=value))
    audit=dict(site=site,validation=validation,fold=fold,train=[cons[j] for j in tr],test=[cons[j] for j in te],inner=[dict(train=[cons[j] for j in a],test=[cons[j] for j in b],one_class_train=len(set(y[a]))<2) for a,b,g in inner],CCR=regionlog)
    for r in pred+chosen: r['label_threshold_pp']=label_threshold
    audit['label_threshold_pp']=label_threshold
    audit['train_positive']=int(y[tr].sum());audit['train_negative']=int(len(tr)-y[tr].sum())
    print('DONE',site,label_threshold,validation,fold,flush=True)
    return pred,chosen,audit
def main():
    outputs=[];inputaudit=[]
    for site in ['C871','C295']:
        cons,d,y,static,raw,diff,regional,audit=setup(site);inputaudit+=audit
        positions=[set(a['positions']) for a in audit];groups=np.array(['+'.join(map(str,sorted(p))) for p in positions]);idx=np.arange(22)
        jobs=[]
        for tr,te in LeaveOneOut().split(idx):jobs.append(('LOOCV',tr,te,cons[te[0]]))
        for g in sorted(set(groups)):jobs.append(('Position_group_out',idx[groups!=g],idx[groups==g],g))
        for tr,te,g in strict(idx,positions):jobs.append(('Strict_position_out',tr,te,g))
        thresholds=[t for t in range(-5,-61,-5) if min((d.delta_editing_pp.values<t).sum(),(d.delta_editing_pp.values>=t).sum())>=3]
        print('LABELS',site,thresholds,flush=True)
        results=Parallel(n_jobs=int(os.environ.get('TRM_JOBS','4')),verbose=5)(delayed(worker)(site,t,cons,(d.delta_editing_pp.values<t).astype(int),static,raw,diff,regional,positions,groups,*j) for t in thresholds for j in jobs)
        outputs+=results
        pd.DataFrame([r for a,b,c in outputs for r in a]).to_csv(O/'outer_predictions.csv',index=False)
        pd.DataFrame([r for a,b,c in outputs for r in b]).to_csv(O/'fold_selection.csv',index=False)
        (O/'split_audit.json').write_text(json.dumps([c for a,b,c in outputs],indent=2))
    pr=pd.read_csv(O/'outer_predictions.csv');rows=[];per=[]
    # Joint label/model/feature/score selection uses outer-training information only.
    sels=pd.read_csv(O/'fold_selection.csv'); joint=[];jointsel=[]
    audits={(a['site'],a['label_threshold_pp'],a['validation'],a['fold']):a for x,z,a in outputs}
    for (site,val,fold),g in sels[sels.reported_family=='Nested_all'].groupby(['site','validation','fold']):
        g=g[[min(audits[(site,r.label_threshold_pp,val,fold)]['train_positive'],audits[(site,r.label_threshold_pp,val,fold)]['train_negative'])>=3 for r in g.itertuples()]]
        assert len(g), 'No eligible training label'
        win=g.sort_values(['inner_BA','label_threshold_pp'],ascending=[False,False],kind='stable').iloc[0]
        q=pr[(pr.site==site)&(pr.validation==val)&(pr.fold==fold)&(pr.label_threshold_pp==win.label_threshold_pp)&(pr.family=='Nested_all')].copy()
        assert len(q)
        q['family']='Nested_label_feature_model';joint.extend(q.to_dict('records'));jointsel.append(win.to_dict())
    pd.DataFrame(jointsel).to_csv(O/'joint_label_selection.csv',index=False)
    jp=pd.DataFrame(joint);jp.to_csv(O/'joint_label_predictions.csv',index=False)
    jr=[]
    for key,g in jp.groupby(['site','validation']):
        w=1/g.construct.map(g.construct.value_counts()).values
        jr.append(dict(site=key[0],validation=key[1],n_unique=g.construct.nunique(),**metrics(g.y,g.score,g.prediction,w)))
    pd.DataFrame(jr).to_csv(O/'joint_label_performance.csv',index=False)
    keys=['site','label_threshold_pp','validation','family','model','policy']
    for key,g in pr.groupby(keys):
        weights=1/g.construct.map(g.construct.value_counts()).values
        rows.append(dict(zip(keys,key),n_unique=g.construct.nunique(),n_predictions=len(g),**metrics(g.y,g.score,g.prediction,weights)))
        for fold,h in g.groupby('fold'):per.append(dict(zip(keys,key),fold=fold,n=len(h),**metrics(h.y,h.score,h.prediction)))
    perf=pd.DataFrame(rows);perf.to_csv(O/'performance.csv',index=False);pd.DataFrame(per).to_csv(O/'per_group_performance.csv',index=False)
    (O/'input_audit.json').write_text(json.dumps(inputaudit,indent=2))
    provenance={'source_commit':json.loads((P/'SOURCE_MANIFEST.json').read_text())['commit'],'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__,'seed':2026,'models':MODELS,'families':FAMILIES,'inputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (P/'inputs').iterdir() if p.is_file()}}
    (O/'provenance.json').write_text(json.dumps(provenance,indent=2))
    print(perf[(perf.policy=='training_only_threshold')&(perf.validation!='LOOCV')].sort_values(['site','validation','BA'],ascending=[True,True,False]).groupby(['site','validation']).head(5).to_string(index=False))
if __name__=='__main__':main()
