"""Descriptive cohort structure for the 44 tested genes (added during manuscript revision).

Reports, per gene: missense carriers, distinct variants, source publications (PMIDs), and the
share of carriers contributed by the largest publication. Same carrier set as design_effect.py.
No association test is run.
"""
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_analysis as ra  # noqa: E402

OUT = HERE / "out"
DATA = "/home/user/data_ext/derived"


def main():
    genes, ind, var = ra.load(DATA)
    v = ra.build_carriers(genes, ind, var)
    v = v[v.has_structure_features.fillna(False).astype(bool) & (v.variant_class == "missense")].copy()
    v["nbr1_plddt"] = float("nan")
    v["nbr2_plddt"] = float("nan")
    f = ra.feature_table(v, "")
    f = f[f.S2_pLDDT.notna() & f.S1_logSemi.notna()]
    tested = pd.read_csv(OUT / "A_tests_primary.csv.gz")
    tested = tested[tested.feature == "S1_logSemi"].gene.unique()
    f = f[f.gene.isin(tested)].assign(pmid=lambda d: d.phenopacket_id.map(ind.set_index("phenopacket_id").pmid))
    rows = []
    for g, d in f.groupby("gene"):
        top = d.pmid.value_counts()
        rows.append(dict(gene=g, carriers=len(d), variants=d.vkey.nunique(), pmids=d.pmid.nunique(),
                         largest_pub_share=float(top.iloc[0] / len(d))))
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "X_cohort_structure_genes.csv", index=False)
    q = lambda s: [float(s.median()), float(s.quantile(.25)), float(s.quantile(.75))]
    out = dict(genes=len(G), carriers=int(G.carriers.sum()), pmids_total=int(f.pmid.nunique()),
               carriers_per_gene_median_iqr=q(G.carriers), variants_per_gene_median_iqr=q(G.variants),
               pmids_per_gene_median_iqr=q(G.pmids), largest_pub_share_median_iqr=q(G.largest_pub_share),
               genes_largest_pub_ge_half=int((G.largest_pub_share >= .5).sum()))
    json.dump(out, open(OUT / "X_cohort_structure.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
