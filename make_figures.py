#!/usr/bin/env python3
"""Render three main figures and five supplementary figures from archived outputs.
No structural contacts, evolutionary rates or statistical models are recomputed.
Run in the project's plotting environment: python make_figures.py.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import re
import shutil

HERE = Path(__file__).parent
OUT = HERE / "figures"
MM = 1 / 25.4
W2 = 170 * MM                      # double-column width
GREY = "#8a8a8a"                   # de-emphasised / excluded marks
BLUE = "#1f6fb4"                   # focal: datasets carried into the test
ORANGE = "#b4531f"                 # AFP Afp2
ORANGE_L = "#d9873f"               # AFP Afp3
BASE, SMALL, TICK = 9, 8, 7        # the three-step size ladder


def style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 600, "savefig.bbox": "tight",
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "legend.fontsize": SMALL, "xtick.labelsize": TICK, "ytick.labelsize": TICK,
        "axes.titlelocation": "left", "axes.titlepad": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "axes.labelcolor": "#222", "text.color": "#222",
        "pdf.fonttype": 42, "ps.fonttype": 42,     # editable text in the vector output
        "figure.constrained_layout.use": False,
    })


def save(fig: plt.Figure, stem: str) -> list[Path]:
    OUT.mkdir(exist_ok=True)
    paths = []
    for ext in ("png", "pdf"):
        p = OUT / f"{stem}.{ext}"
        fig.savefig(p)
        paths.append(p)
    return paths


def check(fig: plt.Figure, name: str) -> None:
    """Geometric legibility check: no two non-tick text boxes may overlap."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ticks = {t for ax in fig.axes for t in ax.get_xticklabels() + ax.get_yticklabels()}
    bad = [(a.get_text()[:24], b.get_text()[:24])
           for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:]
           if ba.overlaps(bb) and a not in ticks and b not in ticks]
    print(f"  {name}: {len(texts)} text objects, overlaps {bad if bad else 'none'}")


# ---------------------------------------------------------------- Figure 2
def figure_remodeling() -> None:
    R = pd.read_csv(HERE / "results/R_all_references.csv")
    prim = {("r2pyocin", "F1"): "R2 pyocin", ("afp", "F1"): "AFP Afp2",
            ("afp", "F2"): "AFP Afp3", ("salcis", "F1"): "SalCIS",
            ("diffocin", "F1"): "Diffocin"}
    R["label"] = [prim.get((s, f)) for s, f in zip(R.source, R.family)]
    R = R[R.label.notna()]
    order = R.groupby("label").R.median().sort_values(ascending=False).index.tolist()
    in_test = {"AFP Afp2", "AFP Afp3", "Diffocin"}
    rng = np.random.default_rng(7)

    fig, ax = plt.subplots(figsize=(W2, 92 * MM))
    for i, lab in enumerate(order):
        v = R[R.label == lab].R.to_numpy()
        col = BLUE if lab in in_test else GREY
        ax.scatter(v, i + rng.uniform(-0.26, 0.26, v.size), s=10, facecolor=col,
                   edgecolor="none", alpha=0.5, zorder=2)
        med = float(np.median(v))
        ax.plot([med, med], [i - 0.36, i + 0.36], color=col, lw=2.6, zorder=4,
                solid_capstyle="butt")
        ax.text(1.05, i, f"{v.size}", fontsize=TICK, va="center", ha="left", color=col)
        ax.text(1.13, i, f"{(v == 0).mean():.0%}", fontsize=TICK, va="center",
                ha="left", color=col)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    for t, l in zip(ax.get_yticklabels(), order):
        t.set_color("#222" if l in in_test else GREY)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.75, len(order) - 0.25)
    ax.set_xlabel("contact remodeling $R_i$   "
                  "(0 = contact set unchanged by contraction, 1 = fully replaced)")
    ax.set_title("Contact remodeling is heterogeneous within and between sheath proteins")
    ax.text(1.05, len(order) - 0.45, "sites", fontsize=TICK, ha="left", color=GREY)
    ax.text(1.13, len(order) - 0.45, "$R_i$=0", fontsize=TICK, ha="left", color=GREY)
    ax.tick_params(axis="y", length=0)
    for s in ("left",):
        ax.spines[s].set_visible(False)
    ax.text(0.0, -0.30, "blue = carried into the evolutionary test;  "
            "grey = structural analysis only.  Vertical bar = median.",
            transform=ax.transAxes, fontsize=TICK, color=GREY)
    check(fig, "Figure1")
    print("   ", save(fig, "Figure1"))
    plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def figure_datasets() -> None:
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    fig,ax=plt.subplots(figsize=(W2,105*MM))
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    ax.text(.5,.98,"Four structural systems evaluated",ha="center",va="top",fontsize=BASE+1)
    xs=[.125,.375,.625,.875]
    system=["R2 pyocin","SalCIS","AFP","Diffocin"]
    detail=["Excluded\nReference lineage\nambiguous","Excluded\nReference locus\npoorly recovered","Retained\nTwo positional\nprotein datasets","Retained\nAmbiguous-genome\nassignments removed"]
    for x,name,txt in zip(xs,system,detail):
        color=GREY if name in ["R2 pyocin","SalCIS"] else BLUE
        ax.add_patch(FancyBboxPatch((x-.112,.69),.224,.14,boxstyle="round,pad=.008",facecolor="#f2f5f8",edgecolor=color))
        ax.text(x,.76,name,ha="center",va="center",fontsize=BASE,fontweight="bold",color=color)
        ax.text(x,.61,txt,ha="center",va="top",fontsize=SMALL,color=color,linespacing=1.6)
        ax.annotate("",(x,.84),(.5,.91),arrowprops=dict(arrowstyle="->",lw=.7,color="#9aa5ad"))
    for x,y,txt,color in [(.49,.16,"AFP Afp2\n27 sequences · 84 sites",ORANGE),(.695,.16,"AFP Afp3\n31 sequences · 81 sites",ORANGE_L),(.90,.16,"Diffocin\n33 sequences · 58 sites",BLUE)]:
        ax.add_patch(FancyBboxPatch((x-.0925,y),.185,.135,boxstyle="round,pad=.007",facecolor="#f6f8fa",edgecolor=color))
        ax.text(x,y+.0675,txt,ha="center",va="center",fontsize=TICK,linespacing=1.7)
        parent=.625 if "AFP" in txt else .875
        ax.add_patch(FancyArrowPatch((parent,.415),(x,y+.145),arrowstyle="->",mutation_scale=9,lw=.9,color=color))
    ax.text(.05,.27,"Two systems retained\nThree fitted protein datasets\n223 fitted sites in total",ha="left",va="center",fontsize=SMALL,linespacing=1.8)
    ax.text(.05,.035,"Afp2 and Afp3 belong to one system; fitted sites are not system-level replicates.",fontsize=TICK,color=GREY)
    fig.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98)
    check(fig,"Figure2");print(save(fig,"Figure2"));plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def figure_association() -> None:
    pr=pd.read_csv(HERE/"results/partial_residuals.csv")
    P={r["dataset"]:r for r in json.loads((HERE/"results/hypothesis_test_primary.json").read_text())}
    order=["afp_Afp2","afp_Afp3","diffocin"]
    names=["AFP Afp2","AFP Afp3","Diffocin"];colors=[ORANGE,ORANGE_L,BLUE]
    fig,axs=plt.subplots(1,3,figsize=(W2,92*MM),sharex=True,sharey=True)
    fig.subplots_adjust(left=.10,right=.98,bottom=.24,top=.70,wspace=.15)
    for letter,ax,k,name,col in zip("abc",axs,order,names,colors):
        dd=pr[pr.protein==k];r=P[k];x=dd.resid_R.to_numpy();y=dd.resid_rate.to_numpy()
        # Algebraic identity check on stored residuals, without refitting the model.
        saved_slope=float(x@y/(x@x));assert np.isclose(saved_slope,r['beta']['R'],rtol=1e-10,atol=1e-10)
        assert abs(x.mean())<1e-10 and abs(y.mean())<1e-10
        lo,hi=r['block_bootstrap']['ci95']
        ax.axhline(0,color="#ddd",lw=.7);ax.axvline(0,color="#ddd",lw=.7)
        ax.scatter(x,y,s=15,facecolor=col,edgecolor="white",linewidth=.4,alpha=.8,zorder=3)
        xs=np.linspace(x.min(),x.max(),100);ax.plot(xs,r['beta']['R']*xs,color=col,lw=1.8,zorder=4)
        ax.set_title(f"{letter}  {name}",fontsize=BASE,fontweight="bold",pad=53)
        ax.text(0,1.05,f"$\\beta_R$ = {r['beta']['R']:+.3f}\nPartial $R^2$ = {r['partial_r2_R']:.3f}\n95% CI [{lo:+.3f}, {hi:+.3f}]",transform=ax.transAxes,fontsize=TICK,va="bottom",linespacing=1.5)
        ax.text(.04,.95,f"n = {r['sites_fitted']} sites",transform=ax.transAxes,va="top",fontsize=TICK,color=GREY)
        ax.set_xlim(-2.6,2.6);ax.set_ylim(-2.0,4.5)
    axs[0].set_ylabel("Residualized evolutionary rate")
    fig.supxlabel("Residualized contact remodeling",y=.115,fontsize=BASE)
    fig.text(.10,.02,"Axes adjusted for monomer RSA, WCN and GNM averaged across states.\nIntervals: bootstrap over ten-reference-position sequence windows.",fontsize=TICK,color=GREY,linespacing=1.5)
    check(fig,"Figure3");print(save(fig,"Figure3"));plt.close(fig)


def figure_metric() -> None:
    """Explanatory drawing; contact sets are hypothetical, not measured structures."""
    from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

    fig = plt.figure(figsize=(W2, 124 * MM))
    ax = fig.add_axes([0.02, 0.03, 0.96, 0.94])
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.015, 1.01); ax.axis("off")
    retained, lost, gained = BLUE, ORANGE, "#187c65"
    ax.text(0.0, 1.0, "From sheath contraction to a site-level evolutionary test",
            fontsize=BASE + 1, va="top")
    ax.text(0.0, 0.91, "a  Contraction of an assembled sheath", fontweight="bold", fontsize=SMALL)
    # Generic ring stacks surround a central tube; dimensions are illustrative.
    for x, w, ys, tag in [(0.11, 0.12, np.linspace(0.47, 0.79, 8), "Extended"),
                           (0.34, 0.19, np.linspace(0.53, 0.72, 8), "Contracted")]:
        ax.plot([x, x], [0.43, 0.81], color="#444", lw=3.2, zorder=1)
        for y in ys:
            ax.add_patch(Ellipse((x, y), w, 0.037, facecolor="#dce9f3",
                                 edgecolor=retained, lw=1.1, zorder=2))
        ax.text(x, 0.39, tag, ha="center", fontsize=SMALL)
    ax.add_patch(FancyArrowPatch((0.19, 0.64), (0.22, 0.64),
                                arrowstyle="-|>", mutation_scale=10, color="#444", lw=1.2))
    ax.text(0.0, 0.33, "Illustrative geometry; not to scale", fontsize=TICK, color=GREY)

    ax.text(0.51, 0.91, "b  Contact identity can change at equal count",
            fontweight="bold", fontsize=SMALL)
    for cx, state, labels in [(0.62, "$C_E$", ['a', 'b', 'c']),
                               (0.88, "$C_C$", ['a', 'd', 'c'])]:
        cy=0.67
        ax.add_patch(Circle((cx, cy), 0.024, facecolor="#333", edgecolor="none", zorder=4))
        ax.text(cx, cy, "$i$", color="white", ha="center", va="center", fontsize=SMALL, zorder=5)
        for dx,dy,lab in zip([-0.06,0.07,0.03],[0.09,0.05,-0.10],labels):
            col=retained if lab in ['a','c'] else lost if lab=='b' else gained
            ax.plot([cx,cx+dx],[cy,cy+dy],color=col,lw=1.4,zorder=2)
            ax.add_patch(Circle((cx+dx,cy+dy),0.021,facecolor=col,edgecolor="white",zorder=3))
            ax.text(cx+dx,cy+dy,lab,color="white",ha="center",va="center",fontsize=TICK)
        ax.text(cx,0.49,state,ha="center",fontsize=BASE)
    for x,label,col in [(0.51,"retained",retained),(0.67,"lost",lost),(0.79,"gained",gained)]:
        ax.plot([x,x+0.025],[0.43,0.43],color=col,lw=2)
        ax.text(x+0.035,0.43,label,va="center",fontsize=TICK)
    ax.text(0.51,0.35,"3 contacts per state; 2 shared; 4 in the union",fontsize=TICK)
    ax.text(0.51,0.29,r"$R_i=1-|C_E\cap C_C|/|C_E\cup C_C|=0.5$",fontsize=SMALL)

    ax.text(0.0, 0.23, "c  Structural comparison and evolutionary correspondence",
            fontweight="bold", fontsize=SMALL)
    boxes=[(0.0,0.22,"Paired structures\n4 systems / 5 proteins"),
           (0.265,0.22,"Positional orthology\n2 systems / 3 proteins"),
           (0.53,0.22,"Relative site rates\n+ monomer RSA/WCN/GNM"),
           (0.795,0.205,"Within-protein\nadjusted association")]
    for x,w,label in boxes:
        ax.add_patch(FancyBboxPatch((x,0.075),w,0.095,boxstyle="round,pad=0.009",
                                    facecolor="#f2f5f8",edgecolor="#bbc6cf",lw=0.8))
        ax.text(x+w/2,0.122,label,ha="center",va="center",fontsize=TICK)
    for x in [0.23,0.495,0.76]:
        ax.add_patch(FancyArrowPatch((x,0.122),(x+0.025,0.122),arrowstyle="-|>",
                                    mutation_scale=8,color="#555",lw=1))
    ax.text(0.0,0.0,"Contact example is hypothetical; no energies or transition trajectories are inferred.",
            fontsize=TICK,color=GREY)
    check(fig,"FigureS1a")
    print("   ",save(fig,"FigureS1a")); plt.close(fig)


def figure_uncertainty() -> None:
    """Describe existing per-site tables; do not estimate rates or refit models."""
    names={"diffocin":"Diffocin", "afp_Afp2":"AFP Afp2", "afp_Afp3":"AFP Afp3"}
    cols={"diffocin":BLUE,"afp_Afp2":ORANGE,"afp_Afp3":ORANGE_L}
    order=list(names)
    d=[pd.read_csv(HERE/f"results/sitetable_primary_{k}.csv") for k in order]
    Rs=[t.R.to_numpy() for t in d]
    uncertainty=[t.rate_std.to_numpy()/float(t.rate.std(ddof=1)) for t in d]
    fig,axs=plt.subplots(1,2,figsize=(W2,96*MM))
    fig.subplots_adjust(left=0.095,right=0.975,bottom=0.21,top=0.81,wspace=0.37)
    rng=np.random.default_rng(7)
    for ax,vals in zip(axs,[Rs,uncertainty]):
        bp=ax.boxplot(vals,positions=np.arange(3),widths=0.45,showfliers=False,patch_artist=True,
                      medianprops={"color":"#222","linewidth":1.3},
                      boxprops={"linewidth":0.7}, whiskerprops={"linewidth":0.7},
                      capprops={"linewidth":0.7})
        for j,(k,v) in enumerate(zip(order,vals)):
            bp['boxes'][j].set_facecolor(cols[k]); bp['boxes'][j].set_alpha(0.20)
            ax.scatter(j+rng.uniform(-0.16,0.16,len(v)),v,s=7,color=cols[k],alpha=0.5,zorder=3)
        ax.set_xticks(np.arange(3)); ax.set_xticklabels([names[k] for k in order],fontsize=TICK)
        ax.set_xlim(-0.6,2.6)
        ax.grid(axis='y',color='#e6e6e6',lw=0.5); ax.set_axisbelow(True)
    axs[0].set_ylim(-0.03,1.03); axs[0].set_ylabel("Contact remodeling $R_i$")
    axs[0].set_title("a  Remodeling at fitted sites",fontsize=SMALL)
    axs[1].set_ylim(bottom=0); axs[1].set_ylabel("Posterior SD / between-site rate SD")
    axs[1].set_title("b  Relative rate uncertainty",fontsize=SMALL)
    med=np.array([np.median(v) for v in uncertainty])
    axs[1].set_xticklabels([f"{names[k]}\nmedian {v:.0%}" for k,v in zip(order,med)],fontsize=TICK)
    fig.text(0.095,0.98,"Remodeling range and relative rate uncertainty in the fitted datasets",
             fontsize=BASE,va="top")
    fig.text(0.095,0.10,"Points: all fitted sites. Boxes: median and interquartile range.",fontsize=TICK,color=GREY)
    fig.text(0.095,0.055,"Descriptive summaries; uncertainty ratios were not used as regression weights.",
             fontsize=TICK,color=GREY)
    check(fig,"FigureS5")
    print("   ",save(fig,"FigureS5")); plt.close(fig)


def figure_structural_qc() -> None:
    text=(HERE/"sheath_pilot/qc_output/structural_qc_report.txt").read_text()
    patterns=[("Pseudomonas aeruginosa R2 pyocin",["R2 pyocin"]),("Salmonella salamae SalCIS",["SalCIS"]),("Clostridioides difficile diffocin",["Diffocin"]),("Serratia entomophila AFP",["AFP Afp2","AFP Afp3"])]
    values={}
    for tag,labs in patterns:
        section=text.split(tag,1)[1]
        # Each protein's first cutoff/domain line appears before the next system.
        cuts=re.findall(r'cutoff sweep : (.*)',section)[:len(labs)]
        domains=re.findall(r'domain sweep : (.*)',section)[:len(labs)]
        assert len(cuts)==len(labs)
        for lab,cut,domain in zip(labs,cuts,domains):
            c=[tuple(map(float,m)) for m in re.findall(r'([0-9.]+)A:n=(\d+),m=([0-9.]+)',cut)]
            d=[tuple(map(float,m)) for m in re.findall(r'>=(\d+):n=(\d+),m=([0-9.]+)',domain)]
            values[lab]=(c,d)
    fig,axs=plt.subplots(2,2,figsize=(W2,125*MM))
    colors=[GREY,"#766ba6",BLUE,ORANGE,ORANGE_L]
    for (lab,(c,d)),col in zip(values.items(),colors):
        for j,vals in enumerate([c,d]):
            arr=np.array(vals);axs[0,j].plot(arr[:,0],arr[:,1],"o-",label=lab,color=col,ms=3,lw=1)
            axs[1,j].plot(arr[:,0],arr[:,2],"o-",color=col,ms=3,lw=1)
    for i,ax in enumerate(axs.flat):
        ax.grid(color="#eee",lw=.5);ax.set_axisbelow(True)
        if i%2==0:ax.axvline(5,color="#ccc",lw=.8);ax.set_xticks([4,4.5,5,5.5,6])
        else:ax.axvline(3,color="#ccc",lw=.8);ax.set_xticks([2,3,4,5])
    axs[0,0].set_title("Contact-distance cutoff");axs[0,1].set_title("Minimum contact union")
    axs[0,0].set_ylabel("Eligible residues");axs[1,0].set_ylabel("Mean remodeling")
    axs[1,0].set_xlabel("Cutoff (Å)");axs[1,1].set_xlabel("Union threshold")
    axs[1,0].set_ylim(.2,.85);axs[1,1].set_ylim(.2,.85)
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc="lower center",ncol=5,frameon=False,fontsize=TICK)
    fig.subplots_adjust(left=.12,right=.97,top=.9,bottom=.15,hspace=.22,wspace=.27)
    check(fig,"FigureS1b");print(save(fig,"FigureS1b"));plt.close(fig)


def gene_track(ax,track,label,sheath_slots):
    from matplotlib.patches import FancyArrow
    tr=[g for g in track if -4<=g['slot']<=8]
    ax.set_xlim(-4.7,8.7);ax.set_ylim(-.7,1.2);ax.axis("off")
    ax.text(-4.6,1.05,label,fontsize=BASE,fontweight="bold",va="top")
    ax.plot([-4.5,8.5],[0,0],color="#ccc",lw=.8)
    for g in tr:
        slot=g['slot'];direction=1 if g.get('rel_strand',g.get('strand','+'))=='+' else -1
        col=ORANGE if slot in sheath_slots else '#a6afb7'
        ax.add_patch(FancyArrow(slot-.35*direction,0,.7*direction,0,width=.20,head_width=.34,head_length=.18,length_includes_head=True,color=col))
        lab=f"{slot}";ax.text(slot,-.26,lab,ha="center",va="top",fontsize=TICK)
    ax.text(-4.6,-.58,"Relative gene slot (arrow lengths schematic)",fontsize=TICK,color=GREY)


def figure_afp_locus() -> None:
    track=json.loads((HERE/"orthology/afp_reference_track.json").read_text())
    fig,ax=plt.subplots(figsize=(W2,43*MM));gene_track(ax,track,"AFP reference locus: three tandem sheath slots",[0,1,2])
    for x,name in [(0,'Afp2'),(1,'Afp3'),(2,'Third')]:ax.text(x,.40,name,ha="center",fontsize=TICK,color=ORANGE)
    fig.subplots_adjust(left=.025,right=.99,top=.98,bottom=.02)
    check(fig,"FigureS3a");print(save(fig,"FigureS3a"));plt.close(fig)


def figure_salcis_diffocin() -> None:
    tracks=json.loads((HERE/"orthology/salcis_diffocin_reference_tracks.json").read_text())
    tracks["salcis"]=json.loads((HERE/"orthology/salcis_reference_track.json").read_text())
    fig=plt.figure(figsize=(W2,119*MM));gs=fig.add_gridspec(3,1,height_ratios=[1,1,1.4],hspace=.24)
    for i,k in enumerate(['salcis','diffocin']):
        ax=fig.add_subplot(gs[i]);gene_track(ax,tracks[k]['track'],f"{'ab'[i]}  {'SalCIS' if k=='salcis' else 'Diffocin'} reference locus",tracks[k]['sheath_slots'])
    ax=fig.add_subplot(gs[2]);labs=[];yes=[];no=[]
    for path,name in [('salcis_loci.csv','SalCIS family'),('salcis_ent_loci.csv','SalCIS enterobacterial'),('diffocin_all_loci.csv','Diffocin family')]:
        t=pd.read_csv(HERE/'orthology'/path);n=int(t.anchored.astype(str).str.lower().eq('true').sum())
        labs.append(name);yes.append(n);no.append(len(t)-n)
    y=np.arange(3);ax.barh(y,yes,color=BLUE,label='Anchored');ax.barh(y,no,left=yes,color='#c9cfd4',label='Unanchored')
    ax.set_yticks(y);ax.set_yticklabels(labs,fontsize=TICK);ax.invert_yaxis();ax.set_xlabel('Screened loci');ax.set_title('c  Retained positional-anchoring records',fontsize=SMALL)
    for j,(a,b) in enumerate(zip(yes,no)):
        ax.text(a/2,j,str(a),ha='center',va='center',fontsize=TICK,color='white');ax.text(a+b/2,j,str(b),ha='center',va='center',fontsize=TICK)
    ax.legend(loc='upper right',frameon=False,fontsize=TICK)
    fig.subplots_adjust(left=.25,right=.97,bottom=.12,top=.97)
    check(fig,"FigureS4");print(save(fig,"FigureS4"));plt.close(fig)


if __name__ == "__main__":
    style()
    figure_remodeling();figure_datasets();figure_association()
    figure_metric();figure_structural_qc();figure_afp_locus();figure_salcis_diffocin();figure_uncertainty()
    for src,dest in [('orthology/r2pyocin_copy_number_audit.png','FigureS2a.png'),('orthology/r2pyocin_locus_comparison.png','FigureS2b.png'),('orthology/afp_joint_tree.png','FigureS3b.png')]:
        shutil.copy2(HERE/src,OUT/dest)
    print("Rendered three main figures and five supplementary figure groups from archived data.")
