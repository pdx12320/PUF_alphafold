import json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O='architecture_validation';a=pd.read_csv(O+'/main_features.csv').sort_values(['architecture','success','construct']);a['ID']=['C%02d'%(i+1) for i in range(len(a))];a[['ID','construct','architecture','success']].to_csv(O+'/plot_ID_mapping.csv',index=False);z=np.load(O+'/all_arrays.npz');ix=[list(z['names']).index(n) for n in a.construct];D=z['density'][ix];y=a.success.to_numpy();S=pd.read_csv(O+'/all_univariate_statistics.csv').set_index('feature');plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False});colors={0:'#0072B2',1:'#D55E00'};markers={0:'o',1:'^'};rng=np.random.default_rng(2026)
def save(fig,name):
 fig.savefig(O+'/'+name+'.png',dpi=300,facecolor='white');fig.savefig(O+'/'+name+'.pdf',facecolor='white');plt.close(fig)
features=['pp_nonlocal12_high_per_res','P5_P6_CP','K204_slot_density','CCR01_fixed'];labels=['Non-local contact density (S12)','P5–P6 mean interface CP','Reference K204 slot: CP density','Original CCR01 mean density']
fig,axes=plt.subplots(1,4,figsize=(13,3.7),layout='constrained')
for ax,f,label in zip(axes,features,labels):
 for cl in [0,1]:
  v=a.loc[y==cl,f].to_numpy();ax.scatter(np.full(len(v),cl)+rng.uniform(-.12,.12,len(v)),v,c=colors[cl],marker=markers[cl],s=35);ax.plot([cl-.2,cl+.2],[v.mean()]*2,color='black',lw=1.5)
 ax.set(xticks=[0,1],xticklabels=['Failure (20)','Success (4)'],title=label,ylabel='Construct-mean feature');ax.set_title(label+f"\nr={S.loc[f,'pearson_r']:.2f}; p={S.loc[f,'p_exact']:.3g}; family q={S.loc[f,'q_family']:.3g}",fontsize=9)
fig.suptitle('Success / failure distributions — one point per construct; black line = mean');save(fig,'success_failure_distributions')
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained');groups=sorted(a.architecture.unique())
for ax,f,label in zip(axes.flat,features,labels):
 for j,g in enumerate(groups):
  for cl in [0,1]:
   v=a[(a.architecture==g)&(a.success==cl)][f].to_numpy();ax.scatter(j+rng.uniform(-.12,.12,len(v)),v,c=colors[cl],marker=markers[cl],s=40,label=('Failure' if cl==0 else 'Success') if j==0 else None)
 ax.set(xticks=range(len(groups)),xticklabels=groups,title=label,ylabel='Construct-mean feature')
handles=[plt.Line2D([],[],color=colors[c],marker=markers[c],ls='',label=['Failure','Success'][c]) for c in [0,1]];fig.legend(handles=handles,loc='outside upper right',ncols=2);fig.suptitle('Architecture distributions — only G03 and G04 contain successes');save(fig,'architecture_feature_distributions')
# Equal architecture-level weighting is NOT used in these group means: each construct has equal weight.
M=np.zeros((len(a),12,12))
for i in range(12):
 for j in range(i+1,12):M[:,i,j]=M[:,j,i]=a[f'P{i+1}_P{j+1}_CP'].to_numpy()
M[:,np.arange(12),np.arange(12)]=np.nan;mean0=np.nanmean(M[y==0],axis=0);mean1=np.nanmean(M[y==1],axis=0);de=mean1-mean0;lim=np.nanmax(abs(de));fig,axes=plt.subplots(1,3,figsize=(12,4.4),layout='constrained')
for ax,v,title in zip(axes,[mean0,mean1,de],['Failure mean (n=20)','Success mean (n=4)','Success − failure']):
 im=ax.imshow(v,cmap=plt.get_cmap('RdBu_r' if title.startswith('Success −') else 'viridis').with_extremes(bad='#eeeeee'),vmin=-lim if title.startswith('Success −') else 0,vmax=lim if title.startswith('Success −') else np.nanmax(M),interpolation='nearest');ax.set(title=title,xticks=range(12),yticks=range(12),xticklabels=range(1,13),yticklabels=range(1,13),xlabel='Repeat position P',ylabel='Repeat position P');fig.colorbar(im,ax=ax,shrink=.75,label='Mean CP')
fig.suptitle('Repeat-pair interface CP (native sequence separation ≥4; diagonal excluded)');save(fig,'repeat_pair_heatmap');np.savez_compressed(O+'/repeat_heatmap_source.npz',construct=a.construct.to_numpy(),matrices=M)
# Requested per-construct adjacent interfaces.
v=a[[f'P{i}_P{i+1}_CP' for i in range(1,12)]].to_numpy();fig,ax=plt.subplots(figsize=(9,8.2),layout='constrained');im=ax.imshow(v,aspect='auto',cmap='viridis',interpolation='nearest');ax.set(yticks=range(len(a)),yticklabels=[f'{r.ID} {r.architecture} '+('S' if r.success else 'F') for r in a.itertuples()],xticks=range(11),xticklabels=[f'P{i}–P{i+1}' for i in range(1,12)],title='Adjacent interface CP by construct (S=success, F=failure)');plt.setp(ax.get_xticklabels(),rotation=45,ha='right');fig.colorbar(im,ax=ax,label='Mean CP');save(fig,'adjacent_interface_heatmap')
cc=json.load(open(O+'/new_CCR_definitions.json'));cc={'CCR01_original':[144,145,146],**cc};v=np.column_stack([D[:,indices].mean(1) for indices in cc.values()]);ref=list(a.construct).index('puf12_r123_r567_r567_r678');delta=v-v[ref];lim=np.max(abs(delta));fig,ax=plt.subplots(figsize=(9,8),layout='constrained');im=ax.imshow(delta,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim,interpolation='nearest');ax.set(yticks=range(len(a)),yticklabels=[f'{r.ID} {r.architecture} '+('S' if r.success else 'F') for r in a.itertuples()],xticks=range(len(cc)),xticklabels=list(cc),title='CCR density difference from original 512-aa success reference');plt.setp(ax.get_xticklabels(),rotation=40,ha='right');fig.colorbar(im,ax=ax,label='Δ summed non-local CP per residue');save(fig,'CCR_heatmap');pd.DataFrame(delta,index=a.construct,columns=cc).to_csv(O+'/CCR_heatmap_source.csv')
print('Static plots complete')
