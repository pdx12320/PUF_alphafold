#!/usr/bin/env python3
"""Fixed-test Wiki figures from the 23 September run; no refitting.

Contract: Figure 1 shows fixed-four scores and their errors; Figure 2 describes
training-only S12; Figure 3 reports all 15 fixed-five endpoint calls.
Python; 600-dpi PNG and editable SVG/PDF; minimum text 6.5 pt. No inference
beyond this small test panel. Model 2 thresholds are common to its five rows.
"""
from make_figures import *
import make_figures as common
import os
plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["DejaVu Sans", "Arial", "Helvetica"], "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none"})
if os.environ.get('FIGURE_QA_TOOLS'):
 sys.path.insert(0,os.environ['FIGURE_QA_TOOLS'])
 from audit_panel_alignment import require_matplotlib_panel_alignment
 common.require_matplotlib_panel_alignment=require_matplotlib_panel_alignment
def save_figure(fig, stem, axes):
 if os.environ.get('FIGURE_QA_TOOLS'):
  fig.canvas.draw()
  require_matplotlib_panel_alignment(fig, axes=axes,
   json_out=QA / f'{stem}.alignment.json',
   overlay_svg=QA / f'{stem}.alignment.svg',
   tolerance_pt=1.5, gutter_tolerance_pt=1.5,
   require_panel_labels=len(axes)>1, strict=True)
 else:common.alignment_gate(fig, stem, axes)
 fig.savefig(OUT / f"{stem}.svg")
 fig.savefig(OUT / f"{stem}.pdf")
 fig.savefig(OUT / f"{stem}.png", dpi=600)
 plt.close(fig)

P=ROOT/'research/fixed_test_20260923'
d=pd.read_csv(P/'model1/predictions.csv');d['design']=d.construct.str.split('_').str[0];d=d.set_index('design').loc[['1','3','7','8']].reset_index()
d.to_csv(DATA/'fixed4_test_predictions.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(183*MM,85*MM));fig.subplots_adjust(left=.12,right=.98,bottom=.23,top=.78,wspace=.65)
ax=axes[0];ax.bar(range(4),d.score,color=[ORANGE if v else BLUE for v in d.success],width=.6)
ax.axhline(.5,color=DARK,lw=.8,ls='--');ax.set_ylim(0,1);ax.set_xticks(range(4),['Design '+x for x in d.design],rotation=20,rotation_mode='anchor',ha='right');ax.set_ylabel('Work score');ax.set_title('Four test predictions',pad=12)
for j,row in d.iterrows():
 inside=abs(row.score+.05-.5)<.08
 ax.text(j,row.score-.065 if inside else row.score+.035,f'{row.score:.3f}',ha='center',fontsize=6.5,color='white' if inside else DARK)
counts=matrix_counts(d.success,d.prediction,[0,1]);confusion_panel(axes[1],counts,['Non-work','Work'],'Experimental agreement · 3/4',3)
for ax,label in zip(axes,['a','b']):panel_label(ax,label)
fig.text(.5,.94,'20 training constructs → 4 fixed test constructs',ha='center',fontsize=9)
fig.text(.5,.03,'Orange: measured work · blue: measured non-work · dashed line: fixed 0.5 cutoff',ha='center',fontsize=6.5)
save_figure(fig,'fixed4_scaffold_test',axes)
s=pd.read_csv(P/'model1/training_data.csv');s[['construct','success','pp_nonlocal12_high_per_res']].to_csv(DATA/'training20_s12.csv',index=False)
fig,ax=plt.subplots(figsize=(120*MM,85*MM));fig.subplots_adjust(left=.20,right=.95,bottom=.22,top=.84)
for cls,color in [(0,BLUE),(1,ORANGE)]:
 v=s.loc[s.success==cls,'pp_nonlocal12_high_per_res'].sort_values().to_numpy();ax.scatter(cls+np.linspace(-.12,.12,len(v)),v,s=25,color=color,edgecolor='white',linewidth=.5,zorder=3);ax.plot([cls-.2,cls+.2],[np.median(v)]*2,color=DARK,lw=1)
ax.set_xticks([0,1],['Non-work (n = 17)','Work (n = 3)']);ax.set_xlim(-.55,1.55);ax.set_ylabel('S12 · high-confidence nonlocal contacts / residue');ax.set_title('Scaffold training-set feature distribution');ax.grid(axis='y',alpha=.18)
fig.text(.57,.04,'Training constructs only; horizontal lines show medians.',ha='center',fontsize=6.5)
save_figure(fig,'training20_s12',[ax])
p=pd.read_csv(P/'model2/predictions.csv');p.to_csv(DATA/'fixed5_endpoint_predictions.csv',index=False)
fig,axes=plt.subplots(1,3,figsize=(183*MM,84*MM));fig.subplots_adjust(left=.14,right=.985,bottom=.29,top=.75,wspace=.72)
for i,(site,labels,display) in enumerate([('C295',['decrease','not_decrease'],['Decrease','No decrease']),('C388',['decrease','within_selected_interval','increase'],['Decrease','Within\ninterval','Increase']),('C871',['decrease','not_decrease'],['Decrease','No decrease'])]):
 q=p[p.endpoint==site];c=matrix_counts(q.true_class,q.predicted_class,labels);confusion_panel(axes[i],c,display,f'{site} · {np.trace(c)}/5 correct',5);panel_label(axes[i],chr(97+i))
fig.text(.5,.93,'76 training constructs → 5 fixed test constructs',ha='center',fontsize=9)
fig.text(.5,.07,'All five excluded together · One training-selected threshold set per endpoint',ha='center',fontsize=7)
save_figure(fig,'fixed5_endpoint_test',axes)
