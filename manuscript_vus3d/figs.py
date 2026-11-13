"""Figures for the VUS-3D manuscript. Reads ../vus3d/out/ only."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent.parent / "vus3d" / "out"
FIG = Path(__file__).resolve().parent / "fig"
import os; PNG = Path(os.environ.get("PNGDIR", "/tmp"))
plt.rcParams.update({"font.family": "serif", "font.serif": ["FreeSerif", "Liberation Serif", "DejaVu Serif"],
                     "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e"})
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"M2-RWR": "#2a78d6", "M2-3D": "#eb6834", "M2-1D": "#1baf7a", "M1": "#52514e", "M0": "#a3a29c"}
MK = {"M2-RWR": "o", "M2-3D": "s", "M2-1D": "^"}
LAB = {"M2-RWR": "Network (RWR, PhenoBridge)", "M2-3D": "3D distance (Gaussian)", "M2-1D": "1D sequence window"}

C = pd.read_csv(OUT / "contrasts.csv"); C = C[C.metric == "strat"]
M = pd.read_csv(OUT / "model_metrics.csv")

# ---------------- Figure 2: forest plot of gains over M1
rows = [("primary", "Primary (n = 55,975; 3,153 genes)"), ("ge2_stars", "≥2 review stars in 2026"),
        ("prior_vus_only", "VUS in 2023 (not absent/conflicting)"), ("genes_ge20_labels", "Genes with ≥20 labelled residues"),
        ("strict_temporal", "Last evaluated after the 2023 cutoff"), ("evaluated_before_2024_09", "Evaluated before 2024-09"),
        ("evaluated_after_2024_09", "Evaluated 2024-09 or later"), ("include_same_residue", "Same-residue labels included"),
        ("non_vcep_pm1_genes", "Genes without a VCEP PM1 rule"), ("vcep_pm1_genes", "Genes with a VCEP PM1 rule (21)")]
fig, ax = plt.subplots(figsize=(6.9, 3.5))
XMAX = 0.03
for i, (tag, lab) in enumerate(rows):
    y0 = len(rows) - 1 - i
    if i % 2 == 0:
        ax.axhspan(y0 - 0.5, y0 + 0.5, color="#f3f2ee", zorder=0, lw=0)
    for k, m in enumerate(["M2-RWR", "M2-3D", "M2-1D"]):
        r = C[(C.analysis == tag) & (C.contrast == f"{m} - M1")].iloc[0]
        y = y0 + (0.22, 0, -0.22)[k]
        hi = min(r.ci_hi, XMAX)
        ax.plot([r.ci_lo, hi], [y, y], color=COL[m], lw=1.1, zorder=2, solid_capstyle="butt")
        if r.ci_hi > XMAX:
            ax.annotate("", xy=(XMAX, y), xytext=(XMAX - 0.0012, y),
                        arrowprops=dict(arrowstyle="-|>", color=COL[m], lw=1.0, mutation_scale=6))
        ax.plot(r.delta_auc, y, MK[m], ms=4.2, color=COL[m], mec="white", mew=0.5, zorder=3,
                label=LAB[m] if i == 0 else None)
ax.axvline(0, color=INK2, lw=0.7, zorder=1)
ax.axvline(0.01, color=INK, lw=0.8, ls=(0, (4, 2)), zorder=1)
ax.text(0.0103, len(rows) - 0.35, "pre-set minimum\nmeaningful gain (+0.01)", fontsize=6.8, color=INK, va="top")
ax.set_yticks(range(len(rows))); ax.set_yticklabels([l for _, l in rows][::-1], fontsize=7.4, color=INK)
ax.tick_params(axis="y", length=0)
ax.set_xlim(-0.004, XMAX); ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_xlabel("Gain in gene-stratified AUROC over AlphaMissense + gene prior (M1), 95% gene-bootstrap CI", color=INK)
ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
ax.spines["left"].set_visible(False)
ax.legend(loc="upper right", bbox_to_anchor=(1, 0.80), fontsize=7, frameon=True, framealpha=1, edgecolor=GRID, handletextpad=0.4)
fig.tight_layout(); fig.savefig(FIG / "fig2_forest.svg"); fig.savefig(PNG / "fig2_forest.png", dpi=170); plt.close(fig)

# ---------------- Figure 3: calibration of the propagation component
cal = pd.read_csv(OUT / "calibration_bins.csv"); amb = pd.read_csv(OUT / "explore_am_band_calibration.csv")
fig, axs = plt.subplots(1, 2, figsize=(6.9, 2.55), gridspec_kw=dict(width_ratios=[7, 5]))
def binlab(b):
    lo, hi = b.strip("(]").split(",")
    f = lambda v: v.replace("-", "−")
    if lo == "-inf":
        return f"≤ {f(hi)}"
    if hi == "inf":
        return f"> {f(lo)}"
    return f"({f(lo)}, {f(hi)}]"
for ax, d, title, err in ((axs[0], cal, "a  All test variants (pre-registered)", False),
                          (axs[1], amb, "b  AlphaMissense ambiguous band (exploratory)", True)):
    x = np.arange(len(d))
    for lv, ls, txt in ((2.406, (0, (4, 2)), "Supporting (2.41)"), (5.790, (0, (1, 1.5)), "Moderate (5.79)")):
        ax.axhline(lv, color=INK, lw=0.7, ls=ls); ax.axhline(1 / lv, color=INK, lw=0.7, ls=ls)
        ax.text(-0.5, lv * 1.07, txt, fontsize=6.3, ha="left", va="bottom", color=INK)
        ax.text(len(d) - 0.55, (1 / lv) * 0.93, txt.replace("Supporting", "benign Supporting").replace("Moderate", "benign Moderate"),
                fontsize=6.3, ha="right", va="top", color=INK)
    ax.axhline(1, color=GRID, lw=0.8, zorder=0)
    if err:
        ax.errorbar(x - 0.12, d.train_LR, yerr=[d.train_LR - d.train_LR_lo5, d.train_LR_hi95 - d.train_LR],
                    fmt="none", ecolor=COL["M2-RWR"], elinewidth=0.9, capsize=0)
    ax.plot(x - 0.12, d.train_LR, "o", ms=4.5, mfc="white", mec=COL["M2-RWR"], mew=1.1, label="TRAIN (out of fold)")
    ax.plot(x + 0.12, d.test_LR, "o", ms=4.5, color=COL["M2-RWR"], mec="white", mew=0.5, label="TEST (2026 labels)")
    ax.set_yscale("log"); ax.set_ylim(0.12, 9.5)
    ax.set_yticks([0.17, 0.25, 0.5, 1, 2, 4, 8]); ax.set_yticklabels(["0.17", "0.25", "0.5", "1", "2", "4", "8"])
    ax.minorticks_off()
    ax.set_xticks(x); ax.set_xticklabels([binlab(b) for b in d.bin], fontsize=6.6)
    ax.set_xlabel("Propagation component (log-odds added to M1)", color=INK)
    ax.set_title(title, loc="left", fontsize=8, color=INK, fontweight="bold")
    ax.set_xlim(-0.6, len(d) - 0.4)
axs[0].set_ylabel("Likelihood ratio, P vs B", color=INK)
axs[0].legend(loc="lower center", bbox_to_anchor=(0.5, 0.075), ncol=2, fontsize=6.6, frameon=False, handletextpad=0.2)
fig.tight_layout(w_pad=1.2); fig.savefig(FIG / "fig3_calibration.svg"); fig.savefig(PNG / "fig3_calibration.png", dpi=170); plt.close(fig)

# ---------------- Figure 1: pooled vs gene-stratified AUROC
m = M[M.analysis == "primary"].set_index("model")
order = ["M0", "M1", "M2-1D", "M2-3D", "M2-RWR"]
names = {"M0": "AlphaMissense", "M1": "+ gene prior", "M2-1D": "+ 1D labels", "M2-3D": "+ 3D labels", "M2-RWR": "+ network labels"}
fig, axs = plt.subplots(1, 2, figsize=(3.35, 2.15), sharey=True)
for ax, col, t in ((axs[0], "auroc", "Pooled AUROC"), (axs[1], "strat_auroc", "Gene-stratified AUROC")):
    for j, k in enumerate(order):
        c = COL.get(k, INK2)
        ax.plot([0.945, m.loc[k, col]], [j, j], color=GRID, lw=1.2, zorder=1, solid_capstyle="butt")
        ax.plot(m.loc[k, col], j, MK.get(k, "D"), ms=4.5, color=c, mec="white", mew=0.5, zorder=3)
        ax.text(m.loc[k, col] + 0.0012, j, f"{m.loc[k, col]:.3f}", va="center", fontsize=6.5, color=INK)
    ax.set_xlim(0.945, 0.985); ax.set_xticks([0.95, 0.96, 0.97, 0.98])
    ax.set_title(t, fontsize=7.6, color=INK); ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=6.6)
axs[0].set_yticks(range(len(order))); axs[0].set_yticklabels([names[k] for k in order], fontsize=7, color=INK)
axs[0].invert_yaxis(); axs[0].tick_params(axis="y", length=0); axs[1].tick_params(axis="y", length=0)
for a in axs: a.spines["left"].set_visible(False)
fig.tight_layout(w_pad=0.6); fig.savefig(FIG / "fig1_metric.svg"); fig.savefig(PNG / "fig1_metric.png", dpi=170); plt.close(fig)
print("ok")
