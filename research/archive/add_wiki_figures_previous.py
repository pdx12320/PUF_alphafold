#!/usr/bin/env python3
"""Reproduce combined validation and three additional Wiki figures; no fitting."""
from make_figures import *
import os
if os.environ.get("FIGURE_QA_TOOLS"):
    sys.path.insert(0, os.environ["FIGURE_QA_TOOLS"])
    import make_figures
    from audit_panel_alignment import require_matplotlib_panel_alignment
    make_figures.require_matplotlib_panel_alignment = require_matplotlib_panel_alignment
from matplotlib.patches import Rectangle, Patch

# Figure contract: recorded construct-level classifications, descriptive S12
# distributions, paired endpoint classes, and coordinate-based candidate locations.
# No inferred temporal ordering, mechanistic effects or synthetic measurements.
DATA.mkdir(exist_ok=True); QA.mkdir(exist_ok=True)
p = pd.read_csv(ROOT / SOURCES['scaffold_predictions']).merge(pd.read_csv(ROOT / SOURCES['scaffold_labels']),on='id',validate='many_to_one')
fig, axes = plt.subplots(1,2,figsize=(183*MM,82*MM))
fig.subplots_adjust(left=.12,right=.97,bottom=.24,top=.78,wspace=.48)
for ax,val,label,n in zip(axes,['LOCO','fixed4'],['a','b'],[24,4]):
    d=p[(p.dataset=='PUF12_24')&(p.model=='RF_CP_structure9')&(p.validation==val)]
    c=matrix_counts(d.true_work,d.prediction,[0,1]); assert c.sum()==n
    confusion_panel(ax,c,['Non-work','Work'], 'LOCO · 24 constructs' if val=='LOCO' else 'Joint holdout · 4 constructs',19)
    panel_label(ax,label)
    ax.text(.5,-.34,f'{np.trace(c)}/{n} correct',transform=ax.transAxes,ha='center')
fig.text(.5,.94,'Scaffold classification · random forest',ha='center',fontsize=9)
save_figure(fig,'scaffold_validation_combined',axes)

s=pd.read_csv(ROOT/'results_20260922/scaffold_RF/data/scaffolds24.csv')
s[['construct','success','pp_nonlocal12_high_per_res']].to_csv(DATA/'scaffold_s12_distribution.csv',index=False)
fig,ax=plt.subplots(figsize=(120*MM,85*MM));fig.subplots_adjust(left=.20,right=.95,bottom=.22,top=.84)
for cls,color in [(0,BLUE),(1,ORANGE)]:
    v=s.loc[s.success==cls,'pp_nonlocal12_high_per_res'].sort_values().to_numpy()
    x=cls+np.linspace(-.12,.12,len(v))
    ax.scatter(x,v,s=25,color=color,edgecolor='white',linewidth=.5,zorder=3)
    ax.plot([cls-.2,cls+.2],[np.median(v)]*2,color=DARK,lw=1)
ax.set_xticks([0,1],['Non-work (n = 20)','Work (n = 4)']);ax.set_xlim(-.55,1.55)
ax.set_ylabel('S12 · high-confidence nonlocal contacts / residue')
ax.set_title('Observed scaffold feature distribution');ax.grid(axis='y',alpha=.18)
fig.text(.57,.04,'One point per construct; horizontal lines show medians.',ha='center',fontsize=6.5)
save_figure(fig,'scaffold_s12_distribution',[ax])

p=pd.read_csv(ROOT/SOURCES['five_construct_predictions'])
p.to_csv(DATA/'five_construct_class_pairs.csv',index=False)
names=['P9-GNS','P9-NPS','P9-NTQ','P8-GVE','P4-R5-SNE+P7-R5-SNE']
classes=['decrease','not_decrease','within_selected_interval','increase']
colors=['#287C9D','#E4ECEF','#DCCDAA','#C47742'];symbols=['D','ND','I','U']
fig,axes=plt.subplots(1,3,figsize=(183*MM,94*MM));fig.subplots_adjust(left=.29,right=.98,bottom=.27,top=.82,wspace=.32)
for i,(ax,site) in enumerate(zip(axes,['C295','C388','C871'])):
    d=p[p.endpoint==site].set_index('target').loc[names]
    vals=np.array([[classes.index(r.predicted_class),classes.index(r.true_class)] for r in d.itertuples()])
    from matplotlib.colors import ListedColormap
    ax.imshow(vals,cmap=ListedColormap(colors),vmin=-.5,vmax=3.5,aspect='auto')
    for y in range(5):
        for x in range(2):ax.text(x,y,symbols[vals[y,x]],ha='center',va='center',color='white' if vals[y,x] in [0,3] else DARK)
        if vals[y,0]!=vals[y,1]:ax.add_patch(Rectangle((-.49,y-.48),1.98,.96,fill=False,ec=DARK,lw=1.4))
    ax.set_xticks([0,1],['Predicted','Measured']);ax.set_yticks(range(5),names if i==0 else ['']*5)
    ax.tick_params(length=0);ax.set_title(site,pad=12);panel_label(ax,chr(97+i))
fig.legend(handles=[Patch(color=c,label=l) for c,l in zip(colors,['D: decrease','ND: no decrease','I: within interval','U: increase'])],loc='lower center',bbox_to_anchor=(.56,.06),ncol=2,frameon=False)
fig.text(.56,.025,'Outlined rows indicate disagreement; classes use fold-specific thresholds.',ha='center',fontsize=6.5)
save_figure(fig,'five_construct_class_pairs',axes)

m=pd.read_csv(ROOT/SOURCES['mpnn_site_scores']);m=m[m.consensus==True].sort_values('pos_1based')
assert len(m)==17 and not m.is_specificity_pos.any()
coords={};rna=[]
for line in (ROOT/'aice_mpnn_20260922/inputs/complex.pdb').read_text().splitlines():
    if not line.startswith('ATOM'):continue
    atom=line[12:16].strip();chain=line[21]; pos=int(line[22:26]);v=[float(line[j:j+8]) for j in [30,38,46]]
    if chain=='A' and atom=='CA':coords[pos]=v
    if chain=='R' and atom=='P':rna.append(v)
back=np.array(list(coords.values()));rna=np.array(rna);center=back.mean(0)
_,_,basis=np.linalg.svd(back-center,full_matrices=False)
b=(back-center)@basis.T;r=(rna-center)@basis.T
points=np.array([coords[int(x)] for x in m.pos_1based]);q=(points-center)@basis.T
src=m[['pos_1based','wt','ligandmpnn_mut']].copy();src['mutation']=src.wt+src.pos_1based.astype(str)+src.ligandmpnn_mut
src[['x_A','y_A','z_A']]=points;src.to_csv(DATA/'nontrm_structure_locations.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(183*MM,113*MM));fig.subplots_adjust(left=.05,right=.95,bottom=.37,top=.86,wspace=.14)
for i,(ax,dims) in enumerate(zip(axes,[(0,1),(0,2)])):
    x,y=dims;ax.plot(b[:,x],b[:,y],color='#9CAAB0',lw=.65,zorder=1)
    ax.plot(r[:,x],r[:,y],color=BLUE,lw=2,zorder=2)
    ax.scatter(q[:,x],q[:,y],s=45,color=ORANGE,edgecolor='white',lw=.5,zorder=3)
    if i==1:
        order=np.argsort(q[:,y]); groups=[(order[:9],-31),(order[9:],31)]
        for ids,rail in groups:
            ids=sorted(ids,key=lambda j:q[j,x])
            for j,labelx in zip(ids,np.linspace(-46,46,len(ids))):
                ax.annotate(str(j+1),xy=(q[j,x],q[j,y]),xytext=(labelx,rail),ha='center',va='center',fontsize=6,color=DARK,arrowprops=dict(arrowstyle='-',color=DARK,lw=.4,shrinkA=5,shrinkB=5),zorder=4)

    ax.set_xlim(-55,55);ax.set_ylim(-38,38);ax.set_aspect('equal');ax.axis('off')
    ax.set_title('Backbone view 1' if i==0 else 'Orthogonal backbone view',pad=10);panel_label(ax,chr(97+i))
for j,mut in enumerate(src.mutation):
    col=j//6;row=j%6;fig.text(.12+col*.29,.30-row*.038,f'{j+1:2d}  {mut}',fontsize=7)
fig.text(.5,.95,'17 consensus non-TRM substitutions on the input structure',ha='center',fontsize=9)
fig.text(.5,.04,'Grey: PUF Cα trace · blue: RNA P trace · orange: candidate positions',ha='center',fontsize=7)
fig.text(.5,.015,'Orthographic projections; residue numbering follows chain A of the 493-residue input.',ha='center',fontsize=6.5)
save_figure(fig,'nontrm_structure_locations',axes)
