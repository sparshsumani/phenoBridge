"""Step 3: duplicate / relatedness screen.

Rules (fixed before looking at any genotype-phenotype association):
  * source = publication token parsed from the individual ID ("Gnazzo, 2020_P10" -> Gnazzo2020; KBG## -> the
    67 new Martinez-Cayuelas 2023 patients, "MC2023").
  * Explicit family labels: Murray 2017 "(f.m)" family.member numbers; Sirmaci2011 "-F1"; trailing A/B suffix
    (KBG10A/B, KBG8A/B, KBG31A/B, Khalifa P1A/P1B). A labelled family whose members carry different ANKRD11
    alleles is flagged as an annotation inconsistency (an AD family must share the allele).
  * Same-source cluster: same source + identical allele -> presumed relatives -> one cluster, EXCEPT for recurrent
    alleles (seen in >=3 different sources, e.g. the c.1903_1907del hotspot), where sharing an allele within a
    source is not evidence of relatedness; those are clustered only through explicit family labels.
  * Cross-source candidate duplicate: identical allele, different source, compatible sex (equal or one U),
    >=5 jointly annotated HPO terms and 0 discordant calls  ("likely");  identical allele, different source,
    >=10 jointly annotated terms and <=1 discordant call regardless of sex ("possible").
  * Calibration: the same concordance statistics among pairs carrying *different* alleles (who cannot be the same
    person) give the false-positive rate of the phenotype criterion.
Outputs: out/dup_pairs.csv, out/clusters.csv (individualId, family_cluster, dup_cluster, keep flags),
         out/dup_summary.json
"""
import itertools
import json
import re

import numpy as np
import pandas as pd

from common import OUT

F = pd.read_csv(OUT / "individuals.csv")
meta = json.load(open(OUT / "meta.json"))
H = meta["hpo"]
P = {r.individualId: r[H].values.astype(float) for _, r in F.iterrows()}
sex = dict(zip(F.individualId, F.sex))
allele = dict(zip(F.individualId, F.alleleKey))
nann = dict(zip(F.individualId, F.n_annotated))


def source(i):
    if re.match(r"KBG\d", i):
        return "MC2023"
    m = re.match(r"([A-Za-z\- ]+?),? ?(\d{4})", i)
    return (m.group(1).strip().replace(" ", "") + m.group(2)) if m else i


src = {i: source(i) for i in F.individualId}


def family_label(i):
    m = re.match(r"Murray, 2017_P\d+ \((\d+)\.\d+\.?\)", i)
    if m:
        return f"Murray2017_fam{m.group(1)}"
    if "-F1" in i:
        return "Sirmaci2011_F1"
    m = re.match(r"(KBG\d+)[AB]$", i)
    if m:
        return m.group(1)
    m = re.match(r"(Khalifa, 2013_P\d+)[AB]$", i)
    if m:
        return m.group(1)
    return None


fam = {i: family_label(i) for i in F.individualId}

pairs = []
for a, b in itertools.combinations(F.individualId, 2):
    pa, pb = P[a], P[b]
    both = ~np.isnan(pa) & ~np.isnan(pb)
    pairs.append(dict(a=a, b=b, same_allele=allele[a] == allele[b], same_src=src[a] == src[b],
                      sex_ok=(sex[a] == sex[b]) or ("U" in (sex[a], sex[b])),
                      shared=int(both.sum()), disc=int((pa[both] != pb[both]).sum())))
R = pd.DataFrame(pairs)
R["likely_dup"] = R.same_allele & ~R.same_src & R.sex_ok & (R.shared >= 5) & (R.disc == 0)
R["possible_dup"] = R.same_allele & ~R.same_src & (R.shared >= 10) & (R.disc <= 1)

# calibration among pairs that cannot be the same person (different alleles)
null = R[~R.same_allele]
cal = {}
for k in (5, 10):
    m = null[(null.shared >= k) & null.sex_ok] if k == 5 else null[null.shared >= k]
    crit = (m.disc == 0) if k == 5 else (m.disc <= 1)
    cal[f"null_rate_{'likely' if k == 5 else 'possible'}_criterion"] = float(crit.mean())
    cal[f"null_pairs_{k}"] = int(len(m))

# expected number of false hits among same-allele cross-source pairs meeting the overlap requirement
sa = R[R.same_allele & ~R.same_src]
exp_likely = cal["null_rate_likely_criterion"] * int(((sa.shared >= 5) & sa.sex_ok).sum())
exp_possible = cal["null_rate_possible_criterion"] * int((sa.shared >= 10).sum())

# inconsistencies in labelled families
fam_rows = []
for f, g in F.assign(fam=F.individualId.map(fam)).dropna(subset=["fam"]).groupby("fam"):
    fam_rows.append(dict(family=f, members=list(g.individualId), alleles=sorted(set(g.alleleKey)),
                         consistent=g.alleleKey.nunique() == 1))

# union-find clusters
parent = {i: i for i in F.individualId}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra


# (1) family clusters: consistent labelled families + same-source same-allele
for fr in fam_rows:
    if fr["consistent"]:
        for m in fr["members"][1:]:
            union(fr["members"][0], m)
n_src_per_allele = F.assign(src=F.individualId.map(src)).groupby("alleleKey").src.nunique()
recurrent = set(n_src_per_allele[n_src_per_allele >= 3].index)
for _, r in R[R.same_allele & R.same_src].iterrows():
    if allele[r.a] not in recurrent:
        union(r.a, r.b)
famc = {i: find(i) for i in F.individualId}
# (2) add likely cross-source duplicates
for _, r in R[R.likely_dup].iterrows():
    union(r.a, r.b)
dupc = {i: find(i) for i in F.individualId}

C = F[["individualId", "alleleKey", "sex", "n_annotated"]].copy()
C["source"] = C.individualId.map(src)
C["family_label"] = C.individualId.map(fam)
C["family_cluster"] = C.individualId.map(famc)
C["family_dup_cluster"] = C.individualId.map(dupc)
# keep one per cluster: the most completely annotated member (ties: alphabetical)
for col, keep in (("family_cluster", "keep_one_per_family"), ("family_dup_cluster", "keep_one_per_family_or_dup")):
    best = C.sort_values(["n_annotated", "individualId"], ascending=[False, True]).groupby(col).individualId.first()
    C[keep] = C.individualId.isin(best.values).astype(int)
# exclusion of likely duplicates only (keep the better annotated member of each likely-dup pair)
drop = set()
for _, r in R[R.likely_dup].iterrows():
    drop.add(r.b if (nann[r.a], r.b) >= (nann[r.b], r.a) else r.a)
C["keep_excl_likely_dups"] = (~C.individualId.isin(drop)).astype(int)
C.to_csv(OUT / "clusters.csv", index=False)

cand = R[R.likely_dup | R.possible_dup].copy()
for c in ("a", "b"):
    cand[f"sex_{c}"] = cand[c].map(sex)
    cand[f"src_{c}"] = cand[c].map(src)
cand["allele"] = cand.a.map(allele)
cand.to_csv(OUT / "dup_pairs.csv", index=False)

summ = dict(
    n=len(F), n_sources=len(set(src.values())),
    recurrent_alleles={a: int(n_src_per_allele[a]) for a in sorted(recurrent)},
    labelled_families=fam_rows,
    n_family_clusters_multi=int((C.groupby("family_cluster").size() > 1).sum()),
    n_in_multi_family_clusters=int(C.groupby("family_cluster").individualId.transform("size").gt(1).sum()),
    n_keep_one_per_family=int(C.keep_one_per_family.sum()),
    n_keep_one_per_family_or_dup=int(C.keep_one_per_family_or_dup.sum()),
    n_likely_dup_pairs=int(R.likely_dup.sum()), n_possible_dup_pairs=int(R.possible_dup.sum()),
    expected_false_likely=exp_likely, expected_false_possible=exp_possible,
    n_keep_excl_likely_dups=int(C.keep_excl_likely_dups.sum()),
    calibration=cal,
    candidates=cand[["a", "b", "allele", "sex_a", "sex_b", "shared", "disc", "likely_dup", "possible_dup"]].to_dict("records"),
)
(OUT / "dup_summary.json").write_text(json.dumps(summ, indent=1, default=str))
print(json.dumps({k: v for k, v in summ.items() if k != "candidates"}, indent=1, default=str))
print(cand[["a", "b", "allele", "sex_a", "sex_b", "shared", "disc", "likely_dup", "possible_dup"]].to_string())
print(C[C.family_cluster.map(C.family_cluster.value_counts()) > 1].sort_values("family_cluster")[["individualId", "alleleKey", "family_cluster"]].to_string())
