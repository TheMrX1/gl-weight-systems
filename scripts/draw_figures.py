from pathlib import Path
import json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Arc,Rectangle
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'report/figures'
plt.rcParams.update({'font.family':'DejaVu Serif','mathtext.fontset':'dejavuserif','font.size':10,'pdf.fonttype':42})
BLUE='#215c85';ORANGE='#b76523'
def draw(ax,w,title,share,arcs,foot=None):
 n=len(w);angles=[math.pi/2-2*math.pi*i/n for i in range(n)];xy=[(math.cos(a),math.sin(a)) for a in angles]
 ax.add_patch(Circle((0,0),1,edgecolor='#7c8791',facecolor='none',lw=.8))
 for a in set(w):
  i,j=[i for i,b in enumerate(w) if a==b];ax.plot([xy[i][0],xy[j][0]],[xy[i][1],xy[j][1]],color=BLUE if a in share else ORANGE,lw=1,alpha=.8)
 for lo,hi,sel in arcs:
  a=90-360*(hi-.5)/n;b=90-360*(lo-.5)/n;ax.add_patch(Arc((0,0),2.07,2.07,theta1=a,theta2=b,color=BLUE if sel else ORANGE,lw=3))
 for i,(x,y) in enumerate(xy):
  ax.plot(x,y,'.',color=BLUE if w[i] in share else ORANGE,ms=4)
  ax.text(x*1.13,y*1.13,str(w[i]),ha='center',va='center',fontsize=7)
 ax.set_title(title,pad=9,fontsize=14);ax.set_aspect('equal');ax.set_xlim(-1.28,1.28);ax.set_ylim(-1.4 if foot else -1.23,1.24);ax.axis('off')
 if foot:ax.text(0,-1.38,foot,ha='center',fontsize=9)
def save(fig,name):
 fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight',pad_inches=.08);fig.savefig(OUT/(name+'.png'),bbox_inches='tight',dpi=150,pad_inches=.08);plt.close(fig)
w=json.loads((ROOT/'results/mutants11/words.json').read_text())
fig,axs=plt.subplots(1,2,figsize=(7,3.35))
draw(axs[0],w['D'],'$D$',set(range(4)),[(0,7,1),(7,17,0),(17,18,1),(18,22,0)])
draw(axs[1],w['D_prime'],"$D'$",set(range(4)),[(0,1,1),(1,11,0),(11,18,1),(18,22,0)])
fig.subplots_adjust(wspace=.2);save(fig,'mutants11')
d=json.loads((ROOT/'results/mutants14/certificate.json').read_text());fig,axs=plt.subplots(2,2,figsize=(7,7))
for ax,word,t,v in zip(axs.flat,d['words'],['$X$','$R(X)$','$S(X)$','$RS(X)$'],d['exact_integer_weights']):
 draw(ax,word,t,set(range(7)),[(0,10,1),(10,20,0),(20,24,1),(24,28,0)] if 'S' not in t else [(0,4,1),(4,14,0),(14,24,1),(24,28,0)],'$w_V='+str(v)+'$')
fig.subplots_adjust(hspace=.22,wspace=.2);save(fig,'mutants14')
fig,axs=plt.subplots(1,3,figsize=(7,2.2))
for ax,rows,title in zip(axs,[[3,3],[3,1,1],[2,1]],['$(3^2)$','$(3,1,1)$','$(2,1)$']):
 for j,r in enumerate(rows):
  for i in range(r):ax.add_patch(Rectangle((i,-j),1,-1,edgecolor=BLUE,facecolor='#edf3f7',lw=1))
 ax.set_aspect('equal');ax.set_xlim(-.2,4);ax.set_ylim(-3.5,.6);ax.set_title(title,fontsize=13);ax.axis('off')
save(fig,'young')
fig,ax=plt.subplots(figsize=(6.8,2.9));data=json.loads((ROOT/'results/graphs9/chromatic_coefficients_mod_p.json').read_text())['nonzero_coefficients_mod_p'];active={(d['genus'],d['x_degree']) for d in data}
for g in range(5):
 for k in range(1,10):
  ax.add_patch(Rectangle((k-.5,g-.5),1,1,facecolor=BLUE if (g,k) in active else '#f0f1f2',edgecolor='white'))
  if (g,k) in active:ax.text(k,g,'×',ha='center',va='center',color='white',fontsize=13)
ax.set_xlim(.5,9.5);ax.set_ylim(4.5,-.5);ax.set_xticks(range(1,10));ax.set_yticks(range(5));ax.set_xlabel('$k$: степень $x$');ax.set_ylabel('$g$: показатель в $N^{-2g}$');ax.spines[['top','right','bottom','left']].set_visible(False);save(fig,'coefficient_grid')
print('Four figures generated from certificate data')
