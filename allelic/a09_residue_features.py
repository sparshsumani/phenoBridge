#!/usr/bin/env python3
"""Step 9: per-residue structural features for every residue of every feature protein
(genes eligible under any definition in a07), reusing multigene/structfeat.py and the s04/s07 logic.

Per protein (AlphaFold v6 single-fragment PDB, raw/alphafold/models/):
  structfeat.compute_residue_features (repo mechanistic_mvp graph/regions/backbone, exact RWR via
  sparse LU): pLDDT, degree_8A, Core/SemiCore/Disorder/CommBackbone(_strict)Impact at alpha 0.3
  (no suffix), 0.15 (_a15), 0.5 (_a5); sasa, rsa (Shrake-Rupley, Wilke max ASA), rel_pos, region,
  in_backbone, dist_to_plddt70_res.
  added here:
    ca_x, ca_y, ca_z           C-alpha coordinates (model numbering = UniProt canonical numbering)
    nbr1_mean_plddt, nbr1_n    mean pLDDT / count of residues at graph distance exactly 1 (8 A contact
                               graph over all residues, self excluded)
    nbr2_mean_plddt, nbr2_n    mean pLDDT / count of residues at graph distance 1 or 2 (self excluded)
    n_ca_10A                   burial proxy 1: number of OTHER C-alpha atoms within 10 A of this C-alpha
    hse_ca_up, hse_ca_down     burial proxy 2: C-alpha half-sphere exposure (Hamelryck 2005;
                               Bio.PDB.HSExposureCA, radius 13 A): number of C-alpha atoms in the half
                               sphere on the side-chain side (up) / opposite side (down); NaN for the
                               terminal residues (needs both neighbours)
  domain annotation (multigene/s04 query + s07.annotate): uniprot_domain, uniprot_region,
    uniprot_disordered_region, uniprot_topology, uniprot_motif, pfam_domain, in_domain,
    in_domain_broad, domain_name  (UniProt REST entry JSON and InterPro Pfam API, cached in cache/)
  AlphaMissense per residue (s07 am_mapping AF->AM numbering): am_pos, am_mean, am_max, am_n_subst
Outputs residue_features/<UNIPROT>.parquet, contact_graphs/<UNIPROT>.parquet (edge list i<j, dist_A),
        intermediate/structure_meta.parquet, intermediate/uniprot_features.parquet,
        alphamissense_feature_proteins.parquet (all AM substitutions of these proteins with af_pos)
"""
import sys, json, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree
from common import INTER, OUT, RFDIR, GRAPHDIR, mg_import, save_json

SF = mg_import("structfeat")
s04 = mg_import("s04_uniprot")
s06 = mg_import("s06_select_genes")
s07 = mg_import("s07_residue_features")


def uniprot_features(acc, hgnc_id):
    e = None
    if isinstance(hgnc_id, str):
        res = s04.uniprot_search(f"xref:hgnc-{hgnc_id.split(':')[1]} AND reviewed:true AND organism_id:9606")
        e = next((x for x in (res or {}).get("results", []) if x["primaryAccession"] == acc), None)
    if e is None:
        res = s04.uniprot_search(f"accession:{acc} AND reviewed:true")
        e = next((x for x in (res or {}).get("results", []) if x["primaryAccession"] == acc), None)
    feats = []
    if e:
        for f in e.get("features", []):
            if f["type"] in s04.FEATURE_TYPES:
                s, t = f["location"]["start"].get("value"), f["location"]["end"].get("value")
                if s is None or t is None:
                    continue
                feats.append(dict(uniprot=acc, source="UniProt", type=f["type"], description=f.get("description"),
                                  start=int(s), end=int(t), feature_id=None))
    pf = s04.interpro_pfam(acc)
    for r in (pf or {}).get("results", []):
        md = r["metadata"]
        for p in r.get("proteins", []):
            for loc in p.get("entry_protein_locations", []):
                for fr in loc.get("fragments", []):
                    feats.append(dict(uniprot=acc, source="Pfam", type="Pfam", description=md.get("name"),
                                      start=int(fr["start"]), end=int(fr["end"]), feature_id=md["accession"]))
    return pd.DataFrame(feats, columns=["uniprot", "source", "type", "description", "start", "end", "feature_id"]), \
        dict(uniprot_entry_found=e is not None, pfam_ok=pf is not None)


def one(args):
    gene, acc, path = args
    from Bio.PDB import HSExposureCA
    t0 = time.time()
    feats, meta, (G_all, nodes, W, core, semi, dis, bb) = SF.compute_residue_features(path)
    st = SF.read_structure(path)
    pos, xyz, pl, seq = SF.load_ca(st)
    assert np.array_equal(pos, feats["pos"]) and np.allclose(pl, feats["plddt"])
    df = pd.DataFrame(feats)
    df.insert(0, "uniprot", acc); df.insert(0, "gene", gene)
    df["ca_x"], df["ca_y"], df["ca_z"] = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    idx = {p: i for i, p in enumerate(pos)}
    n = len(pos)
    e = pd.DataFrame([(min(u, v_), max(u, v_), a["dist"]) for u, v_, a in G_all.edges(data=True)], columns=["i", "j", "dist_A"])
    e = e.sort_values(["i", "j"]).reset_index(drop=True)
    r_, c_ = e.i.map(idx).to_numpy(), e.j.map(idx).to_numpy()
    A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[r_, c_], np.r_[c_, r_])), shape=(n, n))
    A1 = (A > 0).astype(float)
    A2 = ((A1 + A1 @ A1) > 0).astype(float).tolil(); A2.setdiag(0); A2 = A2.tocsr()
    n1 = np.asarray(A1.sum(1)).ravel(); n2 = np.asarray(A2.sum(1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        df["nbr1_mean_plddt"] = np.where(n1 > 0, (A1 @ pl) / n1, np.nan)
        df["nbr2_mean_plddt"] = np.where(n2 > 0, (A2 @ pl) / n2, np.nan)
    df["nbr1_n"], df["nbr2_n"] = n1.astype(int), n2.astype(int)
    tr = cKDTree(xyz)
    df["n_ca_10A"] = np.array([len(x) - 1 for x in tr.query_ball_point(xyz, r=10.0)])
    model = next(st.get_models())
    HSExposureCA(model, radius=13.0)
    chain = next(model.get_chains())
    up, dn = {}, {}
    for res in chain.get_residues():
        if res.id[0].strip():
            continue
        x = res.xtra
        up[int(res.id[1])] = x.get("EXP_HSE_A_U", np.nan); dn[int(res.id[1])] = x.get("EXP_HSE_A_D", np.nan)
    df["hse_ca_up"] = df.pos.map(up).astype(float); df["hse_ca_down"] = df.pos.map(dn).astype(float)
    meta.update(gene=gene, uniprot=acc, af_pdb=path, seconds=time.time() - t0, n_edges=len(e))
    return df, e, meta


def main():
    g = pd.read_parquet(INTER / "genes.parquet")
    sel = g[g.eligible_any_definition & g.af_single_fragment & g.af_pdb.notna()]
    am = pd.read_parquet(INTER / "alphamissense_subset.parquet")
    jobs = [(r.gene, r.uniprot, r.af_pdb) for r in sel.sort_values("uniprot_length", ascending=False).itertuples()]
    print(len(jobs), "proteins", flush=True)
    metas, fts, ams = [], [], []
    hg = dict(zip(sel.uniprot, sel.hgnc_id))
    with ProcessPoolExecutor(4) as ex:
        for df, e, meta in ex.map(one, jobs):
            acc = meta["uniprot"]
            ft, fmeta = uniprot_features(acc, hg[acc])
            meta.update(fmeta)
            fts.append(ft)
            ann = s07.annotate(df.pos.tolist(), ft)
            for k, v in ann.items():
                df[k] = v
            am_acc, posmap = s07.am_mapping(acc, "".join(df.aa_af), None, am)
            meta["am_accession"] = am_acc
            meta["am_mapping"] = None if am_acc is None else ("identity" if posmap is None else "alignment")
            if am_acc is not None:
                a = am[am.uniprot_id == am_acc].copy()
                s = a.groupby("pos").am_pathogenicity.agg(["mean", "max", "count"])
                ap = df.pos if posmap is None else df.pos.map(posmap)
                df["am_pos"] = pd.array(ap, dtype="Int64")
                df["am_mean"] = ap.map(s["mean"]); df["am_max"] = ap.map(s["max"]); df["am_n_subst"] = ap.map(s["count"])
                inv = None if posmap is None else {b: a_ for a_, b in posmap.items()}
                a["af_pos"] = pd.array(a.pos if inv is None else a.pos.map(inv), dtype="Int64")
                a["uniprot"] = acc
                ams.append(a)
            else:
                for c in ("am_pos", "am_mean", "am_max", "am_n_subst"):
                    df[c] = np.nan
            df.to_parquet(RFDIR / f"{acc}.parquet", index=False)
            e.to_parquet(GRAPHDIR / f"{acc}.parquet", index=False)
            meta.pop("af_sequence", None)
            metas.append(meta)
            print(f"{meta['gene']:10s} {acc} L={meta['L']} edges={meta['n_edges']} am={meta['am_mapping']} {meta['seconds']:.1f}s", flush=True)
    pd.DataFrame(metas).to_parquet(INTER / "structure_meta.parquet", index=False)
    pd.concat(fts, ignore_index=True).to_parquet(INTER / "uniprot_features.parquet", index=False)
    if ams:
        pd.concat(ams, ignore_index=True).to_parquet(OUT / "alphamissense_feature_proteins.parquet", index=False)
    save_json(dict(proteins=len(metas), am_mapping={str(k): int(v) for k, v in pd.Series([m["am_mapping"] for m in metas]).value_counts(dropna=False).items()},
                   uniprot_entry_missing=[m["uniprot"] for m in metas if not m["uniprot_entry_found"]]), INTER / "features_log.json")


if __name__ == "__main__":
    main()
