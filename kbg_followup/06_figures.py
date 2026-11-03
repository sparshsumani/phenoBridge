"""Step 6: publication figures (SVG + PDF) -> fig/
 fig1_lollipop: variants along ANKRD11 by consequence class, pLDDT track, exon structure, NMD-escape boundary
 fig2_power:    simulated power of the primary per-term LRT
 fig3_forest:   (a) primary nested-model ORs for the fingerprint, (b) SV vs point-variant ORs (Firth)
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIG, OUT

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["FreeSerif"], "font.size": 7.5, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 300,
})
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
CLS = {  # validated categorical slots 1-5 (dataviz default palette) + distinct marker = secondary encoding
    "frameshift": ("#2a78d6", "o"), "nonsense": ("#eb6834", "s"), "splice": ("#1baf7a", "D"),
    "missense": ("#eda100", "^"), "inframe": ("#e87ba4", "v")}

F = pd.read_csv(OUT / "individuals.csv")
EX = pd.read_csv(OUT / "exons.csv")
FP = pd.read_csv(OUT / "fp_all_residues.csv")
meta = json.load(open(OUT / "meta.json"))
summ = json.load(open(OUT / "annotation_summary.json"))
L = 2663


def save(fig, name):
    for ext in ("svg", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ fig 1
pt = F[F.is_point == 1]
nmd_aa = (summ["last_junction_after_c"] - 50 + 1 + 2) / 3        # first codon whose PTC escapes by the 50-nt rule
fig, ax = plt.subplots(3, 1, figsize=(7.0, 3.5), sharex=True,
                       gridspec_kw=dict(height_ratios=[3.2, 1.0, 0.45], hspace=0.08))
a = ax[0]
for vc, (col, mk) in CLS.items():
    g = pt[pt.vclass == vc].groupby("aa_pos").size()
    if len(g) == 0:
        continue
    a.vlines(g.index, 0, g.values, color=GRID, lw=0.6, zorder=1)
    a.scatter(g.index, g.values, s=14, color=col, marker=mk, edgecolor="#fcfcfb", linewidth=0.5, zorder=3,
              label=f"{vc} ({int(g.sum())})")
a.axvline(nmd_aa, color=INK2, lw=0.7, ls="--", zorder=0)
a.text(nmd_aa - 15, a.get_ylim()[1] * 0.97 if a.get_ylim()[1] > 0 else 30, "NMD-escape\nboundary (50 nt)",
       ha="right", va="top", fontsize=6.5, color=INK2)
top = pt.groupby("aa_pos").size().sort_values(ascending=False).head(1)
hot = pt[pt.aa_pos == top.index[0]].hgvs_c.iloc[0].split(":")[-1]
a.annotate(f"{hot} ({int(top.iloc[0])})", (top.index[0], top.iloc[0]), xytext=(8, -2), textcoords="offset points",
           fontsize=6.5, color=INK2, va="center")
a.set_ylabel("Individuals")
a.set_ylim(0, top.iloc[0] * 1.12)
a.legend(frameon=False, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.13), handletextpad=0.3, columnspacing=1.2)
a.text(1.0, 0.72, f"+ {int(F.is_sv.sum())} individuals with 16q24.3 CNVs (not shown)", transform=a.transAxes,
       ha="right", fontsize=6.5, color=INK2)
# pLDDT
b = ax[1]
b.fill_between(FP.aa_pos, FP.seed_pLDDT, color="#2a78d6", alpha=0.18, lw=0)
b.plot(FP.aa_pos, FP.seed_pLDDT, color="#2a78d6", lw=0.5)
b.axhline(70, color=INK2, lw=0.5, ls=":")
b.set_ylim(0, 100); b.set_yticks([0, 50, 100]); b.set_ylabel("pLDDT")
for s, e in meta["segments"]:
    b.axvspan(s, e, color="#eda100", alpha=0.15, lw=0)
b.axvline(nmd_aa, color=INK2, lw=0.7, ls="--")
# exon / domain track
c = ax[2]
c.set_ylim(0, 1); c.set_yticks([])
c.spines["left"].set_visible(False)
for _, r in EX.dropna(subset=["aa_first"]).iterrows():
    lo, hi = r.aa_first, min(r.aa_last, L)
    c.add_patch(plt.Rectangle((lo, 0.55), hi - lo, 0.35, facecolor="#e6e5e1" if r.exon % 2 else "#c3c2b7",
                              edgecolor="#fcfcfb", lw=0.6))
    if hi - lo > 150:
        c.text((lo + hi) / 2, 0.72, f"exon {int(r.exon)}", ha="center", va="center", fontsize=6.3, color=INK)
c.text(L + 15, 0.72, "exons 10–13", ha="left", va="center", fontsize=6, color=INK2)
ank = meta["ank"]
c.add_patch(plt.Rectangle((ank[0][0], 0.08), ank[-1][1] - ank[0][0], 0.32, facecolor="#6250d6", alpha=0.8, lw=0))
c.text(ank[-1][1] + 20, 0.24, "ANK", fontsize=6.3, va="center", color=INK)
d0, d1 = meta["deg_region"]
c.add_patch(plt.Rectangle((d0, 0.08), d1 - d0, 0.32, facecolor="#e34948", alpha=0.75, lw=0))
c.text(d0 - 20, 0.24, "C-terminal degradation region (UniProt)", fontsize=6.3, va="center", ha="right", color=INK)
c.set_xlim(1, L)
c.set_xlabel("ANKRD11 residue (Q6UB99 / NP_037407.4)")
save(fig, "fig1_lollipop")

# ------------------------------------------------------------------ fig 2
G = pd.read_csv(OUT / "power_grid.csv")
fig, a = plt.subplots(figsize=(3.4, 2.3))
for alpha, ls, col in ((0.05, "-", "#2a78d6"), (0.005, "--", "#eb6834")):
    pv = G.pivot(index="OR_per_SD", columns="hpo", values=f"power_{alpha}")
    a.fill_between(pv.index, pv.quantile(0.1, axis=1), pv.quantile(0.9, axis=1), color=col, alpha=0.13, lw=0)
    a.plot(pv.index, pv.median(axis=1), color=col, lw=1.6, ls=ls, label=f"α = {alpha}")
a.axhline(0.8, color=INK2, lw=0.6, ls=":")
a.text(4.0, 0.82, "80%", ha="right", va="bottom", fontsize=6.5, color=INK2)
a.set_xlabel("True OR per SD of log$_{10}$ SemiCoreImpact")
a.set_ylabel("Power (per-term LRT)")
a.set_ylim(0, 1); a.set_xlim(1, 4)
a.legend(frameon=False, loc="lower right", title="median across 34 terms\n(band: 10th–90th pct.)",
         title_fontsize=6.3)
save(fig, "fig2_power")

# ------------------------------------------------------------------ fig 3
P = pd.read_csv(OUT / "primary.csv").sort_values("OR")
B1 = pd.read_csv(OUT / "baseline_B1_SV_vs_point.csv").sort_values("OR")
SHORT = {"Attention deficit hyperactivity disorder": "ADHD",
         "Delayed speech and language development": "Speech/language delay",
         "Intrauterine growth retardation": "Intrauterine growth restriction"}
fig, ax = plt.subplots(1, 2, figsize=(7.0, 4.6), gridspec_kw=dict(wspace=0.85))
for a, T, title, xl, col in (
        (ax[0], P, "(a) Fingerprint added to class + position (n = 256)", "OR per SD of log$_{10}$ SemiCoreImpact",
         "#2a78d6"),
        (ax[1], B1, "(b) 16q24.3 CNV vs point variant (n = 333)", "OR, CNV vs point variant (Firth)", "#eb6834")):
    y = np.arange(len(T))
    XL = (0.25, 5.0) if a is ax[0] else (0.05, 20.0)
    a.hlines(y, T.lo.clip(XL[0], XL[1]), T.hi.clip(XL[0], XL[1]), color=col, lw=0.9)
    for yy, (lo, hi) in enumerate(zip(T.lo, T.hi)):
        if lo < XL[0]:
            a.annotate("", (XL[0], yy), (XL[0] * 1.6, yy), arrowprops=dict(arrowstyle="->", color=col, lw=0.8))
        if hi > XL[1]:
            a.annotate("", (XL[1], yy), (XL[1] / 1.6, yy), arrowprops=dict(arrowstyle="->", color=col, lw=0.8))
    a.scatter(T.OR, y, s=12, color=col, zorder=3, edgecolor="#fcfcfb", linewidth=0.4)
    a.axvline(1, color=INK2, lw=0.6, ls="--")
    a.set_xscale("log")
    a.set_yticks(y)
    a.set_yticklabels([f"{SHORT.get(l, l)} ({n})" for l, n in zip(T.label, T.n)], fontsize=6.2)
    a.set_ylim(-0.7, len(T) - 0.3)
    a.set_xlabel(xl)
    a.set_title(title, loc="left", fontsize=7.5)
    for yy, (p, q) in enumerate(zip(T.p, T.q)):
        if p < 0.05:
            a.text(a.get_xlim()[1] if False else min(T.hi.iloc[yy], XL[1]) * 1.08, yy, f"p={p:.3f}", fontsize=5.8, va="center",
                   color=INK2)
    ticks = [0.1, 0.3, 1, 3, 10] if a is ax[1] else [0.3, 0.5, 1, 2, 4]
    a.set_xticks(ticks); a.set_xticklabels([str(t) for t in ticks]); a.set_xlim(*XL)
    a.minorticks_off()
save(fig, "fig3_forest")
print("figures written to", FIG)
