"""Render at final physical size; preview by resampling the actual 300-dpi PNG."""
import csv
import math
from collections import Counter
from pathlib import Path
from statistics import NormalDist

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = HERE / "data"
MASTER, PREVIEW = HERE / "masters_300dpi", HERE / "preview_72dpi"
for dest in MASTER, PREVIEW:
    dest.mkdir(exist_ok=True)
plt.rcParams.update({"font.family":"Arial","font.size":10, "mathtext.fontset":"custom",
                     "mathtext.rm":"Arial","mathtext.it":"Arial:italic","mathtext.bf":"Arial:bold",
                     "text.color":"#191d20","axes.labelcolor":"#191d20","savefig.facecolor":"white"})
def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f: return list(csv.DictReader(f))
def write(path, rows):
    with path.open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
studies={r["study_id"]:r for r in read(ROOT/"meta_analysis/studies.csv")}
comp=read(ROOT/"meta_analysis/comparisons.csv")
def name(s):
    r=studies[s]; suffix=s[-1] if s.endswith(("_a","_b")) else ""
    return f"{r['first_author']} ({r['year']}{suffix})"
rows=read(DATA/"forest_rows.csv"); models=read(DATA/"forest_models.csv")
bias=read(DATA/"publication_bias_summary.csv")
primary_models={r["analysis"].split()[0].lower():r for r in
                read(ROOT/"revision/analysis/results/revision_analysis_summary.csv")
                if r["analysis_family"]=="primary" and "RVE+CR2" in r["analysis"]}
manifest=[]; plotted=[]; W=492  # 173.57 mm; one drawing coordinate = one pt
def canvas(height):
    fig=plt.figure(figsize=(W/72,height/72),dpi=144)
    ax=fig.add_axes([0,0,1,1],xlim=(0,W),ylim=(height,0));ax.axis("off")
    return fig,ax
def label(ax,x,y,s,sz=10,**kw):
    return ax.text(x,y,s,fontsize=sz,va="center",**kw)
def rule(ax,x0,x1,y,color="#89939a",lw=.65):
    ax.plot([x0,x1],[y,y],color=color,lw=lw)
def save(fig,stem,n,low,study_pt,pooled=""):
    path=MASTER/(stem+".png"); fig.savefig(path,dpi=300); plt.close(fig)
    with Image.open(path) as im:
        width,height=im.size; preview=im.convert("RGB").resize(
            (round(width*72/300),round(height*72/300)),Image.Resampling.LANCZOS)
        preview.save(PREVIEW/path.name,dpi=(72,72))
        assert abs(im.info["dpi"][0]-300)<.1 and width*25.4/300<=174 and height*25.4/300<=234
    manifest.append(dict(file=path.name,pixels_300=f"{width}x{height}",
      pixels_72=f"{preview.width}x{preview.height}",width_mm=round(width*25.4/300,2),
      height_mm=round(height*25.4/300,2),study_rows_or_funnel_points=n,min_type_pt=low,
      study_label_pt=study_pt,pooled_kappa_95CI=pooled))
def qarel():
    active={r["study_id"] for r in comp if not r["exclude"].strip()}
    q=[r for r in read(ROOT/"meta_analysis/qarel_for_analysis.csv") if r["study_id"] in active]
    q.sort(key=lambda r:(studies[r["study_id"]]["year"],name(r["study_id"])))
    assert len(q)==len(active)==len({r["study_id"] for r in q})==79
    # Use both shape and colour. At reduction, retain a visible mark for every score.
    styles={"Yes":("o","#18815b"),"No":("X","#b93a31"),"Unclear":("D","#d39517"),"N/A":("s","#8a9297")}
    domains=read(ROOT/"revision/analysis/results/qarel_item_counts.csv")
    for d in domains:
        a=Counter(r[d["item"]] for r in q)
        assert a==Counter({k:int(d[col]) for k,col in zip(styles,["yes","no","unclear","not_applicable"])})
    fig,ax=canvas(625)
    for x,rs in [(6,q[:40]),(253,q[40:])]:
        label(ax,x,13,"Study",10,weight="bold")
        for j in range(7): label(ax,x+109+j*20,13,f"Q{j+1}",9,ha="center",weight="bold")
        rule(ax,x,x+237,23)
        for i,r in enumerate(rs):
            y=32+i*12.55
            if i%2==0: ax.add_patch(Rectangle((x-1,y-6),238,12.55,facecolor="#f4f6f7",edgecolor="none"))
            label(ax,x,y,name(r["study_id"]),9.5)
            for j in range(7):
                marker,c=styles[r[f"Q{j+1}"]]
                ax.plot(x+109+j*20,y,marker=marker,ms=5.8,mec=c,mfc=c,mew=.4)
    rule(ax,6,486,534)
    for x,(key,(m,c)) in zip([20,88,154,242],styles.items()):
        ax.plot(x,549,marker=m,ms=6,color=c);label(ax,x+9,549,key if key!="N/A" else "Not applicable",10)
    for i,d in enumerate(domains):
        col,ri=divmod(i,4);x=7+244*col;y=568+ri*14
        label(ax,x,y,d["item"],9.5,weight="bold");label(ax,x+21,y,d["domain"],9.5)
    save(fig,"Figure_2_QAREL",79,9,9.5)

def forest(stem,s,groups,last):
    selected=[r for r in rows if r["stratum"]==s and r["group"] in groups]
    primary=next(m for m in models if m["stratum"]==s and m["type"]=="RVE+CR2")
    opes=next(m for m in models if m["stratum"]==s and m["type"]=="OPES")
    items=[]
    for group in groups:
        sub=sorted([r for r in selected if r["group"]==group],key=lambda r:name(r["study_id"]))
        items.append(("group",f"{group}  ({len(sub)} studies)",None))
        items.extend(("study",name(r["study_id"]),r) for r in sub)
        m=next(m for m in models if m["stratum"]==s and m["type"]=="subtotal" and m["group"]==group)
        items.append(("subtotal","Subtotal",m))
    if last:
        clusters=primary_models[s]["clusters"]
        items.extend([("overall",f"Overall OPES ({opes['k']} studies)",opes),
                      ("primary",f"RVE+CR2 ({primary['k']} effects, {clusters} clusters)",primary)])
    im=(s=="inter-modality")
    dy=14.6; start=36; axis=start+len(items)*dy+4; h=axis+(58 if im else 31)
    fig,ax=canvas(h);lx,nx,left,right,end=6,166,190,357,485
    assert all(-.5<=float(r["lower"])<=float(r["upper"])<=1 for r in selected)
    xx=lambda v:left+(float(v)+.5)/1.5*(right-left)
    title="Study"; label(ax,lx,13,title,11,weight="bold")
    label(ax,nx,13,"n",11,ha="right",style="italic",weight="bold")
    label(ax,(left+right)/2,13,r"Agreement ($\kappa$)" if not im else r"Concordance ($\kappa$)",10.5,ha="center")
    label(ax,end,13,r"$\kappa$  [95% CI]",11,ha="right",weight="bold")
    rule(ax,lx,end,23); ax.plot([xx(0)]*2,[25,axis],color="#b1bbc1",lw=.65,zorder=0)
    def trip(r):return [float(r[k]) for k in ("estimate","lower","upper")]
    def fmt(v,nd):return f"{v:.{nd}f}".replace("-","−")
    for i,(kind,title,r) in enumerate(items):
        y=start+i*dy
        if kind=="group":
            ax.add_patch(Rectangle((lx-1,y-6.3),end-lx+2,12.6,color="#e9eff2",zorder=1))
            label(ax,lx+2,y,title,10.5,weight="bold");continue
        c="#91312f" if kind=="primary" else "#191d20"
        t=label(ax,lx,y,title,10 if kind=="study" else 9.8,weight="normal" if kind=="study" else "bold",color=c)
        e,lo,hi=trip(r)
        if kind=="study":
            if r["marker"]:
                fig.canvas.draw(); edge=t.get_window_extent().transformed(ax.transData.inverted()).x1
                label(ax,edge+1.4,y-3.3,r["marker"],8.2)
            label(ax,nx,y,r["n"],10,ha="right")
            rule(ax,xx(lo),xx(hi),y,"#566671",1.05)
            ax.plot(xx(e),y,"s",ms=2.8+1.35*math.sqrt(float(r["weight"])),color="#294451",mew=0)
            plotted.append(dict(panel=stem,study_id=r["study_id"],comparison_id=r["comparison_id"],
                                stratum=s,group=r["group"],estimate=r["estimate"]))
        else:
            ax.add_patch(Polygon([(xx(lo),y),(xx(e),y-3.9),(xx(hi),y),(xx(e),y+3.9)],color=c))
            if kind=="overall":rule(ax,lx,end,y-7)
        nd=2 if kind=="study" else 3
        label(ax,end,y,f"{fmt(e,nd)}  [{fmt(lo,nd)}, {fmt(hi,nd)}]",9.8,ha="right",color=c)
    rule(ax,left,right,axis)
    for v in (-.5,0,.5,1):
        x=xx(v);ax.plot([x,x],[axis,axis+3],lw=.65,color="black")
        label(ax,x,axis+12, f"{v:g}".replace("-","−"),9.5,ha="center")
    if im:
        label(ax,lx,axis+30,"a   Dual-energy CT–MRI comparison",9.5)
        label(ax,lx,axis+44,"b   Combined CT/MRI arm",9.5)
    summary="; ".join(f"{r['group'] if kind=='subtotal' else title}: {float(r['estimate']):.3f} [{float(r['lower']):.3f}, {float(r['upper']):.3f}]"
                      for kind,title,r in items if kind in ("overall","primary","subtotal"))
    save(fig,stem,len(selected),8.2 if im else 9.5,10,summary)

def funnel():
    d=[r for r in rows if r["stratum"]=="inter-reader"]; b=next(r for r in bias if r["stratum"]=="inter-reader")
    se=[math.sqrt(float(r["vi"])) for r in d]; x=[float(r["yi"]) for r in d]
    mean=math.atanh(float(b["opes_kappa"])); top=math.ceil(max(se)*10)/10
    rad=NormalDist().inv_cdf(.975)*top
    fig,ax=plt.subplots(figsize=(W/72,340/72))
    fig.subplots_adjust(left=.115,right=.976,top=.956,bottom=.154)
    ax.fill([mean,mean-rad,mean+rad],[0,top,top],color="#f1f4f5",zorder=0)
    ax.plot([mean-rad,mean,mean+rad],[top,0,top],"--",lw=1,color="#6b7880")
    ax.axvline(mean,color="#30495a",lw=1);ax.scatter(x,se,s=21,color="#2b4554",edgecolor="white",lw=.5)
    ax.set_ylim(top,0);ax.set_ylabel("Standard error",fontsize=12);ax.set_xlabel(r"Fisher $z$-transformed $\kappa$",fontsize=12)
    ax.spines[["top","right"]].set_visible(False);ax.tick_params(labelsize=10)
    save(fig,"Figure_S2A_funnel_inter_reader",len(x),10,"",f"OPES: {float(b['opes_kappa']):.3f} [{float(b['ci_lower']):.3f}, {float(b['ci_upper']):.3f}]")

qarel()
forest("Figure_3A_part1_CT","inter-reader",["CT"],False)
forest("Figure_3A_part2_other_modalities","inter-reader",["MRI","CEUS/US","Combined CT/MRI"],True)
forest("Figure_3B_inter_modality","inter-modality",["CT vs MRI","CT vs US-based"],True)
funnel()
# Verify actual rendered rows independently against raw OPES inclusion,
# variance availability, not just counts or the R-produced input table.
def number(s):
    try: return float(s)
    except (TypeError,ValueError): return float("nan")
def eligible(r):
    return (not r["exclude"].strip() and r["opes_include"]=="1" and math.isfinite(number(r["kappa"]))
            and ((number(r["ci_upper"])>number(r["ci_lower"])) or number(r["kappa_se"])>0 or number(r["n_cysts"])>3))
for s,k in [("inter-reader",46),("inter-modality",30)]:
    expected={(r["study_id"],r["comparison_id"]) for r in comp if r["analytic_stratum"]==s and eligible(r)}
    actual=[(r["study_id"],r["comparison_id"]) for r in plotted if r["stratum"]==s]
    assert len(actual)==len(expected)==k and set(actual)==expected and len(set(a[0] for a in actual))==k
write(DATA/"forest_rendered_row_audit.csv",plotted)
write(HERE/"figure_inventory.csv",manifest)
