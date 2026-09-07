import os,json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from Bio.PDB import MMCIFParser
O=os.path.dirname(__file__);z=np.load(O+'/construct_arrays.npz');names=z['names'].tolist();a=pd.read_csv(O+'/all_construct_features.csv').set_index('construct').loc[names];y=a.success.values;cp=z['cp'];D=z['density'];ref='puf12_r123_r567_r567_r678';mapping=json.load(open(O+'/core_mapping.json'));ccr=json.load(open(O+'/CCR_definitions.json'));contacts=pd.read_csv(O+'/contact_tests_all14.csv');ids=pd.read_csv('previous/results/construct_summary.csv').set_index('construct').ID.to_dict();success='#0072B2';failure='#B45B00'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def save(fig,name):
 fig.savefig(O+'/'+name+'.png',dpi=240,facecolor='white');fig.savefig(O+'/'+name+'.pdf',facecolor='white');plt.close(fig)
meanS=np.mean(cp[y==1,:396,:396],axis=0);meanF=np.mean(cp[y==0,:396,:396],axis=0);delta=meanS-meanF;delta[abs(np.arange(396)[:,None]-np.arange(396)[None,:])<12]=np.nan
np.savez_compressed(O+'/delta_CP_common11.npz',success=meanS,failure=meanF,delta=delta)
fig,ax=plt.subplots(2,2,figsize=(13,10),layout='constrained');im=ax[0,0].imshow(delta,cmap='RdBu_r',vmin=-1,vmax=1,origin='lower',interpolation='nearest');ticks=np.arange(11)*36+17.5;ax[0,0].set_xticks(ticks,range(1,12));ax[0,0].set_yticks(ticks,range(1,12));ax[0,0].set(title='A  Nonlocal ΔCP: success − failure',xlabel='Repeat slot (N → C)',ylabel='Repeat slot (N → C)');fig.colorbar(im,ax=ax[0,0],label='ΔCP')
sub=delta[144:180,180:216];im=ax[0,1].imshow(sub,cmap='RdBu_r',vmin=-.5,vmax=.5,origin='lower',interpolation='nearest',extent=(.5,36.5,.5,36.5));ax[0,1].set(title='B  Repeat-slot 5–6 interface',xlabel='Position in repeat slot 6',ylabel='Position in repeat slot 5');fig.colorbar(im,ax=ax[0,1],label='ΔCP (detail scale)')
for xp,yp,txt in [(26,27,'K204–E239'),(22,19,'L196–A235')]:ax[0,1].scatter(xp,yp,s=110,facecolors='none',edgecolors='#111111',lw=1);ax[0,1].annotate(txt,(xp,yp),(xp-14,yp+5),arrowprops={'arrowstyle':'-','lw':.7},fontsize=9)
x=np.arange(396)+1;dd=D[y==1,:396].mean(0)-D[y==0,:396].mean(0);ax[1,0].plot(x,dd,color='#333333',lw=.9);ax[1,0].axhline(0,color='grey',lw=.7);ax[1,0].set(title='C  Residue nonlocal-contact density',xlabel='Aligned core residue (11 common slots)',ylabel='Success − failure: sum of CP')
for reg in ccr:
 ix=reg['indices'];ax[1,0].axvspan(ix[0]+1,ix[-1]+1,color='#0072B2',alpha=.35);ax[1,0].annotate(reg['CCR'],(np.mean(ix)+1,dd[ix].mean()),xytext=(np.mean(ix)+1,1.4),arrowprops={'arrowstyle':'-','lw':.7},fontsize=9,rotation=45)
i,j=170,205;values=cp[:,i,j];rng=np.random.default_rng(2026)
for cl,label,col,mark in [(0,'Failure (n=11)',failure,'o'),(1,'Success (n=3)',success,'D')]:
 take=y==cl;xx=np.full(sum(take),cl)+rng.uniform(-.1,.1,sum(take));ax[1,1].scatter(xx,values[take],c=col,marker=mark,s=45,label=label);ax[1,1].plot([cl-.17,cl+.17],[values[take].mean()]*2,color='black',lw=2)
ax[1,1].set_xticks([0,1],['Failure','Success']);ax[1,1].set(ylim=(-.02,1.02),ylabel='Construct-mean CP',title='D  K204–E239: an exploratory contact');ax[1,1].text(.03,.97,'ΔCP = +0.271; r = 0.864\nExact p = 0.00275; BH q = 0.349',transform=ax[1,1].transAxes,va='top');fig.suptitle('Structural associations with construction outcome',fontsize=17);save(fig,'Structural_differences')
# CCR seed distributions: all constructs, no omission of unfavorable cases.
sv=pd.read_csv(O+'/CCR_seed_values.csv');order=sorted(names,key=lambda n:ids[n]);fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
for ax,reg in zip(axes,ccr):
 for k,n in enumerate(order):
  v=sv[(sv.CCR==reg['CCR'])&(sv.construct==n)].density.values;cl=int(a.loc[n,'success']);col=success if cl else failure;mark='D' if cl else 'o';ax.scatter(np.full(len(v),k),v,color=col,s=14,alpha=.45);ax.errorbar(k,v.mean(),yerr=v.std(ddof=1),fmt=mark,color=col,capsize=3)
 ax.set_xticks(range(14),[ids[n] for n in order],rotation=60);ax.set(xlabel='Construct ID',ylabel='Mean residue nonlocal CP density',title=reg['CCR']+' | reference residues '+','.join(str(mapping[ref][i]+1) for i in reg['indices']))
fig.suptitle('Candidate CCRs: construct means ± SD across seeds',fontsize=16);save(fig,'CCR_seed_profiles')
# Model ablation with common color limits and exact numbers.
m=pd.read_csv(O+'/LOCO_model_metrics.csv');fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
for ax,subset,title in zip(axes,['all14','12repeat'],['All 14 constructs','12-repeat constructs only']):
 mat=np.array([[m[(m.subset==subset)&(m.model==mode+'_'+typ)].AUC.iloc[0] for typ in ['LR','RF','XGB']] for mode in ['summary','structural','full']]);im=ax.imshow(mat,vmin=0,vmax=1,cmap='cividis',aspect='auto');ax.set_xticks(range(3),['Logistic','Random Forest','XGBoost']);ax.set_yticks(range(3),['CP summary','+ confidence / compactness','+ CCR / residue features']);ax.set_title(title)
 for i in range(3):
  for j in range(3):ax.text(j,i,f'{mat[i,j]:.3f}',ha='center',va='center',color='black' if mat[i,j]>.55 else 'white',fontsize=14)
fig.colorbar(im,ax=axes,label='LOCO ROC AUC');fig.suptitle('Feature additions do not consistently improve prediction',fontsize=16);save(fig,'Model_comparison')
# Reference-coordinate projection, no structural generation or smoothing.
st=MMCIFParser(QUIET=True).get_structure('x',O+'/'+ref+'_representative.cif');coords=np.array([r['CA'].coord for r in st[0]['A'].get_residues()]);cent=coords-coords.mean(0);_,_,v=np.linalg.svd(cent,full_matrices=False);xy=cent@v[:2].T;fig,ax=plt.subplots(figsize=(10,7),layout='constrained');ax.plot(xy[:,0],xy[:,1],color='#B8B8B8',lw=.8)
for rep in range(12):
 ind=np.array(mapping[ref][rep*36:(rep+1)*36]);cen=xy[ind].mean(0);ax.text(*cen,f'P{rep+1}',fontsize=9,color='#666666')
for reg,col in zip(ccr,['#0072B2','#009E73']):
 ind=[mapping[ref][i] for i in reg['indices']];ax.scatter(xy[ind,0],xy[ind,1],s=65,color=col,label=reg['CCR'],zorder=5)
for i,j,col in [(196,235,'#B45B00'),(204,239,'#B45B00'),(283,315,'#0072B2')]:
 ii,jj=i-1,j-1;ax.plot(xy[[ii,jj],0],xy[[ii,jj],1],color=col,lw=2,ls='--');ax.scatter(xy[[ii,jj],0],xy[[ii,jj],1],color=col,s=30,zorder=6)
for i,offset in [(196,(-30,20)),(235,(8,15)),(204,(-45,-25)),(239,(8,-25)),(283,(8,15)),(315,(-28,-24))]:
 ax.annotate(str(i),xy[i-1],textcoords='offset points',xytext=offset,fontsize=10,arrowprops={'arrowstyle':'-','color':'#444444','lw':.7})
ax.set_aspect('equal');ax.set(xlabel='Coordinate PC1 (Å)',ylabel='Coordinate PC2 (Å)',title='Reference PUF: CCRs and candidate contacts\nOrange = enhanced; blue dashed = reduced; P = ordinal repeat slot');ax.legend(loc='best');save(fig,'Reference_structure_map')
