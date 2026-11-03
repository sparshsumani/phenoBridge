"""Step 7: markdown tables for RESULTS.md, generated from out/ (avoids hand-copied numbers).
Output: out/tables_md.txt"""
import json

import numpy as np
import pandas as pd

from common import OUT


def md(df, fmt=None):
    fmt = fmt or {}
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if c in fmt and pd.notna(v) and not isinstance(v, str):
                cells.append(fmt[c](v))
            elif isinstance(v, float):
                cells.append("" if np.isnan(v) else f"{v:.3g}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


p3 = lambda v: f"{v:.3f}" if v >= 0.001 else f"{v:.1e}"
f2 = lambda v: f"{v:.2f}"
pct = lambda v: f"{100 * v:.1f}%"
S = []

A = pd.read_csv(OUT / "variants_annotated.csv")
pt = A[A.vclass != "SV"]
ct = pd.crosstab(pt.vclass, pt.regex_trunc.map({0: "regex: not truncating", 1: "regex: truncating"}))
ct.insert(0, "VEP class (MANE)", ct.index)
S.append("### T1 regex flag vs VEP\n" + md(ct.reset_index(drop=True)))

P = pd.read_csv(OUT / "primary.csv")
t = P[["label", "n", "n_obs", "n_exc", "OR", "lo", "hi", "p", "q", "firth_p", "used"]].head(10).copy()
t["OR (95% CI) per SD"] = [f"{o:.2f} ({l:.2f}–{h:.2f})" for o, l, h in zip(t.OR, t.lo, t.hi)]
t = t[["label", "n", "n_obs", "n_exc", "OR (95% CI) per SD", "p", "q", "firth_p", "used"]]
S.append("### T2 primary (top 10 of %d by p)\n" % len(P) + md(t, {"p": p3, "q": p3, "firth_p": p3}))

sens = pd.read_csv(OUT / "primary_sens.csv")
S.append("### T3 sensitivity\n" + md(sens, {"min_p": p3, "min_q": p3, "top_OR": f2}))
om = json.load(open(OUT / "omnibus.json"))
S.append("### omnibus\n" + json.dumps(om))

rows = []
for k in ("B1_SV_vs_point", "B1s_SVdeletion_vs_point", "B2_nontrunc_vs_trunc", "B3i_PTC_position_per10pct",
          "B3ii_PTC_loses_Cterm_region"):
    T = pd.read_csv(OUT / f"baseline_{k}.csv")
    for _, r in T[T.p < 0.1].iterrows():
        d = dict(family=k, term=r.label, n=r.n, OR_CI=f"{r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f})", p=r.p, q=r.q)
        if "prev_exp" in T:
            d["prevalence exposed vs unexposed"] = f"{r.exp_pos}/{r.exp_n} ({100*r.prev_exp:.0f}%) vs {r.unexp_pos}/{r.unexp_n} ({100*r.prev_unexp:.0f}%)"
        rows.append(d)
S.append("### T4 baseline (terms with p<0.1)\n" + md(pd.DataFrame(rows), {"p": p3, "q": p3}))

R = pd.read_csv(OUT / "mc2023_reproduction.csv")
R["OR (95% CI)"] = [f"{o:.2f} ({l:.2f}–{h:.2f})" for o, l, h in zip(R.OR, R.lo, R.hi)]
R["this cohort"] = [f"{100*a:.1f}% vs {100*b:.1f}% (n={int(x)}/{int(y)})" if pd.notna(a) else "" for a, b, x, y in
                    zip(R.prev_exp, R.prev_unexp, R.exp_n.fillna(0), R.unexp_n.fillna(0))]
S.append("### T5 MC2023\n" + md(R[["claim", "reported", "this cohort", "n", "OR (95% CI)", "p", "fisher_p"]],
                                {"p": p3, "fisher_p": p3}))

MD = pd.read_csv(OUT / "power_mde.csv")
m = MD[["label", "n", "n_obs", "MDE_OR_perSD_0.05", "MDE_OR_perSD_0.005", "MDE_OR_per10fold_0.05"]].copy()
for c in ("MDE_OR_perSD_0.05", "MDE_OR_perSD_0.005", "MDE_OR_per10fold_0.05"):
    m[c] = m[c].map(lambda v: ">4 per SD" if pd.isna(v) else f"{v:.2f}")
S.append("### T6 power\n" + md(m))

E = pd.read_csv(OUT / "exploratory_disorder_ID.csv")
S.append("### T7 exploratory disorder-ID\n" + md(E[["analysis", "n", "n_exc", "OR_perSD", "p_ml", "p_firth"]],
                                                 {"p_ml": p3, "p_firth": p3, "OR_perSD": f2}))
MS = pd.read_csv(OUT / "missense_table.csv")
S.append("### T8 missense\n" + md(MS, {"SemiCoreImpact": f2, "CoreImpact": f2}))
(OUT / "tables_md.txt").write_text("\n\n".join(S) + "\n")
print("\n\n".join(S))
