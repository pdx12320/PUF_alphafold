#!/usr/bin/env python3
"""Model 2 repeat-level contact importance; fixed fitted-model attribution.
All P1–P12 and outside-core contacts retained; each endpoint separately normalized.
"""
from make_figures import *
from pathlib import Path
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans','Arial'],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
p=ROOT/'research/model_interpretation_20260923/model2_repeat_importance.csv';d=pd.read_csv(p)
cols=list(range(1,13))+[0];values=d.pivot(index='endpoint',columns='P',values='cp_share').reindex(index=['C295','C388','C871'],columns=cols).fillna(0)*100
values.to_csv(DATA/'model2_repeat_importance_percent.csv',index_label='endpoint')
fig,ax=plt.subplots(figsize=(183/25.4,67/25.4));fig.subplots_adjust(left=.11,right=.88,bottom=.29,top=.79)
im=ax.imshow(values,cmap='Blues',vmin=0,vmax=20,aspect='auto')
ax.set_xticks(range(13),['P'+str(p) if p else 'Outside\ncores' for p in cols]);ax.set_yticks(range(3),values.index);ax.tick_params(length=0)
for i in range(3):
 for j in range(13):
  v=values.iloc[i,j];ax.text(j,i,f'{v:.1f}',ha='center',va='center',fontsize=6.5,color='white' if v>11 else DARK)
ax.set_title('Where do the fixed endpoint models assign contact importance?',pad=14)
cax=fig.add_axes([.91,.29,.015,.50]);cb=fig.colorbar(im,cax=cax);cb.set_label('Share of CP importance (%)')
fig.text(.49,.10,'Each row totals 100% · C388 sequence descriptors are summarized separately',ha='center',fontsize=6.5)
fig.text(.49,.04,'RF impurity importance (C295); absolute standardized coefficients (C388, C871)',ha='center',fontsize=6.5)
alignment_gate(fig,'model2_repeat_importance',[ax])
fig.savefig(OUT/'model2_repeat_importance.svg');fig.savefig(OUT/'model2_repeat_importance.pdf');fig.savefig(OUT/'model2_repeat_importance.png',dpi=600)
