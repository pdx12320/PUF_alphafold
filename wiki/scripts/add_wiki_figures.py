#!/usr/bin/env python3
"""Reproduce combined validation and S12 Wiki figures; no fitting."""
from make_figures import *
import os
if os.environ.get("FIGURE_QA_TOOLS"):
    sys.path.insert(0, os.environ["FIGURE_QA_TOOLS"])
    import make_figures
    from audit_panel_alignment import require_matplotlib_panel_alignment
    make_figures.require_matplotlib_panel_alignment = require_matplotlib_panel_alignment
from matplotlib.patches import Rectangle, Patch

# Figure contract: recorded construct-level classifications, descriptive S12
# distributions.
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

s=pd.read_csv(ROOT/'research/results_20260922/scaffold_RF/data/scaffolds24.csv')
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
