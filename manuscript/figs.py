"""Figures for the manuscript. Reads analysis_multigene/out/ and writes manuscript/fig/."""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "analysis_multigene" / "out"
FIG = Path(__file__).resolve().parent / "fig"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["FreeSerif"], "mathtext.fontset": "stix",
    "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6,
    "svg.fonttype": "path",
})
BLUE, RED, GREY, LGREY = "#2b5c8a", "#c0504d", "#8c8c8c", "#cfcfcf"

NAMES = {
    "S1_logSemi": "RWR semi-core mass (S1)",
    "S2_pLDDT": "Seed pLDDT (S2)",
    "S3_logCore": "RWR core mass (S3)",
    "S4_degree": "Contact degree (S4)",
    "S5_AlphaMissense": "AlphaMissense (S5)",
    "S6a_nbr1_pLDDT": "1-hop neighbour pLDDT (S6a)",
    "S6b_nbr2_pLDDT": "2-hop neighbour pLDDT (S6b)",
    "S7_in_domain": "In UniProt domain (S7)",
    "S8_logDisorder": "RWR disorder mass (S8)",
    "P0_rel_pos": "Relative position (P0)",
    "B_missense_vs_PTC": "Missense vs truncating (B)",
}


def panel(ax, t):
    ax.text(-0.02, 1.04, f"({t})", transform=ax.transAxes, fontweight="bold", va="bottom", ha="right")


def fig1():
    cpv = pd.read_csv(OUT / "X_carriers_per_variant.csv", index_col=0).carriers
    pg = pd.read_csv(OUT / "X_per_gene_calibration.csv")
    nul = np.load(OUT / "A_permnull_primary.npz")["S1_logSemi"]
    nul_i = np.load(OUT / "A_permnull_primary_individual_perm.npz")["S1_logSemi"]
    S = pd.read_csv(OUT / "summary.csv")
    r = S[(S.analysis == "primary") & (S.feature == "S1_logSemi")].iloc[0]

    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.05), gridspec_kw=dict(width_ratios=[1.12, 1.1, 1.1], wspace=.40))
    # (a) individuals by recurrence of their variant
    bins = [(1, 1, "1"), (2, 2, "2"), (3, 5, "3–5"), (6, 10, "6–10"), (11, 20, "11–20"), (21, 10**6, "21+")]
    vals = [int(cpv[(cpv >= lo) & (cpv <= hi)].sum()) for lo, hi, _ in bins]
    x = np.arange(len(bins))
    ax[0].bar(x, vals, color=[GREY] + [BLUE] * (len(bins) - 1), width=.75)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([b[2] for b in bins], fontsize=6.3)
    ax[0].set_xlabel("Carriers of the same variant")
    ax[0].set_ylabel("Individuals")
    share = cpv[cpv > 1].sum() / cpv.sum()
    ax[0].text(.97, .97, f"{share:.0%} of {cpv.sum():,} share\ntheir variant with\nanother individual",
               transform=ax[0].transAxes, ha="right", va="top", fontsize=7)
    panel(ax[0], "a")

    # (b) null distributions for S1
    lo = min(nul_i.min(), nul.min()) - 2
    hi = max(nul.max(), r.n_p05) + 4
    b = np.arange(lo, hi, 2)
    ax[1].hist(nul_i, bins=b, color=LGREY, label="Shuffle individuals")
    ax[1].hist(nul, bins=b, color=BLUE, alpha=.75, label="Shuffle variants")
    ax[1].axvline(r.n_p05, color=RED, lw=1.1)
    ax[1].text(r.n_p05 + 2, ax[1].get_ylim()[1] * .97, f"observed\n{int(r.n_p05)}", color=RED, fontsize=7, va="top")
    ax[1].axvline(r.expected_p05, color="k", lw=.6, ls=":")
    ax[1].set_xlabel("Tests with $p<0.05$ (of %d)" % r.n_tests)
    ax[1].set_ylabel("Permutations")
    ax[1].legend(frameon=False, fontsize=6.5, loc="center right", bbox_to_anchor=(1.02, .5), handlelength=1)
    ax[1].set_xlim(lo, hi + 30)
    panel(ax[1], "b")

    # (c) per-gene false-positive rate vs design effect
    G = pd.read_csv(OUT / "X_design_effect_genes.csv")
    from scipy.stats import norm
    ax[2].scatter(G.mean_de, G.null_rate_individual_S1_logSemi, s=9, color=LGREY, lw=0,
                  label="Shuffle individuals", zorder=2)
    ax[2].scatter(G.mean_de, G.null_rate_variant_S1_logSemi, s=11, color=BLUE, lw=0,
                  label="Shuffle variants", zorder=3)
    d = np.linspace(1, 9, 200)
    ax[2].plot(d, 2 * (1 - norm.cdf(1.96 / np.sqrt(d))), color="k", lw=.8, label="Predicted from DE", zorder=1)
    ax[2].axhline(.05, color="k", lw=.6, ls=":")
    ax[2].set_xscale("log")
    ax[2].set_xticks([1, 2, 4, 8])
    ax[2].set_xticklabels(["1", "2", "4", "8"])
    ax[2].minorticks_off()
    ax[2].set_ylim(-.02, 1.04)
    ax[2].set_xlabel("Design effect, mean over a gene's tests")
    ax[2].set_ylabel("Type I error of naive test")
    ax[2].legend(frameon=False, fontsize=6.5, loc="upper left", handlelength=1, handletextpad=.3)
    panel(ax[2], "c")
    fig.savefig(FIG / "fig1_recurrence.svg", bbox_inches="tight")
    fig.savefig(FIG / "fig1_recurrence.pdf", bbox_inches="tight")
    plt.close(fig)


def fig2():
    S = pd.read_csv(OUT / "summary.csv")
    nA = np.load(OUT / "A_permnull_primary.npz")
    nAi = np.load(OUT / "A_permnull_primary_individual_perm.npz")
    nB = np.load(OUT / "B_permnull.npz")
    order = ["P0_rel_pos", "B_missense_vs_PTC", "S5_AlphaMissense", "S1_logSemi", "S2_pLDDT", "S6a_nbr1_pLDDT",
             "S4_degree", "S7_in_domain", "S3_logCore", "S6b_nbr2_pLDDT", "S8_logDisorder"]
    fig, ax = plt.subplots(figsize=(3.4, 2.9))
    y = np.arange(len(order))[::-1]
    for yi, ft in zip(y, order):
        r = S[(S.analysis == "primary") & (S.feature == ft)].iloc[0]
        if ft.startswith("B_"):
            nv, ni = nB["variant"], nB["individual"]
        else:
            nv, ni = nA[ft], nAi[ft]
        n = r.n_tests
        for arr, col, dy, lw in ((ni, LGREY, -.14, 3.2), (nv, BLUE, .14, 3.2)):
            q = np.percentile(arr, [2.5, 97.5]) / n
            ax.plot(q, [yi + dy] * 2, color=col, lw=lw, solid_capstyle="butt")
            ax.plot(np.mean(arr) / n, yi + dy, "|", color="k", ms=4, mew=.6)
        ax.plot(r.n_p05 / n, yi, "o", color=RED, ms=3.6, zorder=4)
        ax.text(.262, yi, f"{r.perm_p:.3f}" if r.perm_p >= .001 else "<0.001", va="center", ha="right",
                fontsize=6.8, fontweight="bold" if r.perm_p < .05 else "normal")
    ax.axvline(.05, color="k", lw=.6, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([NAMES[o] for o in order], fontsize=7)
    ax.set_xlim(0, .265)
    ax.set_xticks([0, .05, .1, .15, .2])
    ax.set_xlabel("Fraction of tests with $p<0.05$")
    ax.text(.262, len(order) - .35, "perm. $p$\n(variants)", ha="right", va="bottom", fontsize=6.5)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=LGREY, lw=3.2), Line2D([], [], color=BLUE, lw=3.2),
         Line2D([], [], marker="o", color=RED, lw=0, ms=3.6)]
    ax.legend(h, ["95% null, shuffle individuals", "95% null, shuffle variants", "Observed"], frameon=False,
              fontsize=6.3, loc="lower center", bbox_to_anchor=(.42, 1.0), ncol=2, handlelength=1.3,
              columnspacing=.8)
    ax.set_ylim(-.6, len(order) - .4)
    fig.savefig(FIG / "fig2_features.svg", bbox_inches="tight")
    fig.savefig(FIG / "fig2_features.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig2()
    if (OUT / "X_per_gene_calibration.csv").exists():
        fig1()
