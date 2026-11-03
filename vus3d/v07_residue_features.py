#!/usr/bin/env python3
"""Step 7: contact graphs and per-residue structural features for every feature protein
(v05 feature_protein: >=1 TEST variant and >=5 TRAIN-labelled residues; single-fragment AlphaFold v6).

Exactly the allelic/a09 per-protein computation (a09.one -> multigene/structfeat.compute_residue_features:
8 A C-alpha contact graph over all residues, pLDDT, degree_8A, Core/SemiCore/Disorder/CommBackbone RWR
masses at alpha 0.3 (no suffix), 0.15 (_a15), 0.5 (_a5) by exact sparse solve, SASA/RSA, rel_pos, region,
neighbourhood pLDDT, C-alpha count burial n_ca_10A, half-sphere exposure), a09.uniprot_features + s07.annotate
domains (UniProt REST / InterPro Pfam, cached in cache/) and per-residue AlphaMissense summary
(s07.am_mapping AF->AM numbering; am_mean / am_max / am_n_subst).
Proteins already computed by allelic/a09 (same AF model and code) are symlinked from derived_allelic/.
Resumable: proteins whose two output files exist are skipped (per-protein metadata in intermediate/structure_meta/).
Outputs residue_features/<UNIPROT>.parquet, contact_graphs/<UNIPROT>.parquet (i<j residue numbers, dist_A),
        intermediate/structure_meta.parquet, intermediate/uniprot_features.parquet, intermediate/features_log.json
"""
import os, sys, json, time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from common import INTER, OUT, RFDIR, GRAPHDIR, ALLELIC_RFDIR, ALLELIC_GRAPHDIR, ALLELIC_OUT, save_json, al_import, mg_import

a09 = al_import("a09_residue_features")
s07 = mg_import("s07_residue_features")
METADIR = INTER / "structure_meta"; METADIR.mkdir(exist_ok=True)


def add_annotations(df, acc, hgnc, amsub):
    ft, fmeta = a09.uniprot_features(acc, hgnc)
    ann = s07.annotate(df.pos.tolist(), ft)
    for k, v in ann.items():
        df[k] = v
    am_acc, posmap = s07.am_mapping(acc, "".join(df.aa_af), None, amsub) if amsub is not None else (None, None)
    fmeta["am_accession"] = am_acc
    fmeta["am_mapping"] = None if am_acc is None else ("identity" if posmap is None else "alignment")
    if am_acc is not None:
        s = amsub.groupby("pos").am_pathogenicity.agg(["mean", "max", "count"])
        ap = df.pos if posmap is None else df.pos.map(posmap)
        df["am_pos"] = pd.array(ap, dtype="Int64")
        df["am_mean"] = ap.map(s["mean"]).astype(float); df["am_max"] = ap.map(s["max"]).astype(float); df["am_n_subst"] = ap.map(s["count"])
    else:
        for c in ("am_pos", "am_mean", "am_max", "am_n_subst"):
            df[c] = np.nan
    return df, ft, fmeta


def main():
    g = pd.read_parquet(INTER / "genes.parquet")
    sel = g[g.feature_protein & g.af_pdb.notna()].sort_values("uniprot_length", ascending=False)
    am = pd.read_parquet(OUT / "alphamissense_feature_proteins.parquet")
    am = am.rename(columns={"uniprot": "uniprot_id"})
    am["uniprot_id"] = am.uniprot_id.astype(str)
    amg = {a: d.reset_index(drop=True) for a, d in am.groupby("uniprot_id", observed=True)}
    hg = dict(zip(sel.uniprot, sel.hgnc_id))
    # reuse allelic outputs (identical model + code)
    alle_meta = pd.read_parquet(ALLELIC_OUT / "intermediate/structure_meta.parquet") if (ALLELIC_OUT / "intermediate/structure_meta.parquet").exists() else None
    linked = []
    for r in sel.itertuples():
        a_rf, a_g = ALLELIC_RFDIR / f"{r.uniprot}.parquet", ALLELIC_GRAPHDIR / f"{r.uniprot}.parquet"
        if a_rf.exists() and a_g.exists() and alle_meta is not None and r.uniprot in set(alle_meta.uniprot):
            for src, dst in ((a_rf, RFDIR / f"{r.uniprot}.parquet"), (a_g, GRAPHDIR / f"{r.uniprot}.parquet")):
                if not dst.exists():
                    os.symlink(src, dst)
            mf = METADIR / f"{r.uniprot}.json"
            if not mf.exists():
                m = alle_meta[alle_meta.uniprot == r.uniprot].iloc[0].to_dict()
                m["source"] = "symlink_derived_allelic"
                mf.write_text(json.dumps(m, default=str))
            linked.append(r.uniprot)
    print(len(sel), "feature proteins;", len(linked), "linked from derived_allelic", flush=True)
    jobs = [(r.gene, r.uniprot, r.af_pdb) for r in sel.itertuples()
            if not ((RFDIR / f"{r.uniprot}.parquet").exists() and (GRAPHDIR / f"{r.uniprot}.parquet").exists() and (METADIR / f"{r.uniprot}.json").exists())]
    print(len(jobs), "to compute", flush=True)
    # prefetch domain annotations (REST, cached) in threads while the structure jobs run
    th = ThreadPoolExecutor(3)
    futs = {acc: th.submit(a09.uniprot_features, acc, hg[acc]) for _, acc, _ in jobs}
    t0 = time.time()
    with ProcessPoolExecutor(4, mp_context=mp.get_context("fork")) as ex:
        for k, (df, e, meta) in enumerate(ex.map(a09.one, jobs, chunksize=1)):
            acc = meta["uniprot"]
            futs[acc].result()
            df, ft, fmeta = add_annotations(df, acc, hg[acc], amg.get(acc))
            meta.update(fmeta); meta.pop("af_sequence", None); meta["source"] = "computed_vus3d"
            ft.to_parquet(METADIR / f"{acc}.features.parquet", index=False)
            df.to_parquet(RFDIR / f"{acc}.parquet", index=False)
            e.to_parquet(GRAPHDIR / f"{acc}.parquet", index=False)
            (METADIR / f"{acc}.json").write_text(json.dumps(meta, default=str))
            print(f"[{k + 1}/{len(jobs)} {time.time() - t0:.0f}s] {meta['gene']:10s} {acc} L={meta['L']} edges={meta['n_edges']} am={meta['am_mapping']} {meta['seconds']:.1f}s", flush=True)
    th.shutdown()
    metas = [json.loads((METADIR / f"{a}.json").read_text()) for a in sel.uniprot]
    md = pd.DataFrame(metas)
    for c in md.columns:
        if md[c].map(lambda x: isinstance(x, (list, tuple, dict))).any():
            md[c] = md[c].astype(str)
    md.to_parquet(INTER / "structure_meta.parquet", index=False)
    fts = [pd.read_parquet(METADIR / f"{a}.features.parquet") for a in sel.uniprot if (METADIR / f"{a}.features.parquet").exists()]
    if fts:
        pd.concat(fts, ignore_index=True).to_parquet(INTER / "uniprot_features.parquet", index=False)
    save_json(dict(proteins=len(md), linked_from_allelic=len(linked), computed=int((md.source == "computed_vus3d").sum()),
                   am_mapping={str(k): int(v) for k, v in md.am_mapping.value_counts(dropna=False).items()},
                   uniprot_entry_missing=sorted(md[md.uniprot_entry_found == False].uniprot.tolist()) if "uniprot_entry_found" in md else None),
              INTER / "features_log.json")


if __name__ == "__main__":
    main()
