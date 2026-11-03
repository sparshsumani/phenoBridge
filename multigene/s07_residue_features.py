#!/usr/bin/env python3
"""Step 7: per-residue structural features for every selected gene (or --all genes with a
single-fragment model) -> DERIVED/residue_features/<UNIPROT>.parquet, INTER/structure_meta.parquet.

Domain annotation columns:
  uniprot_domain   UniProt features Domain / Repeat / Zinc finger / DNA binding covering the residue
  uniprot_region   UniProt 'Region' features EXCLUDING the MobiDB-lite 'Disordered' predictions
  uniprot_disordered_region  residue inside a UniProt 'Region: Disordered' (sequence-based prediction)
  uniprot_topology Topological domain / Transmembrane / Intramembrane description
  uniprot_motif    Motif / Coiled coil / Compositional bias
  pfam_domain      Pfam (InterPro API) domain names
  in_domain        uniprot_domain or pfam_domain non-empty
  in_domain_broad  in_domain or uniprot_region or Motif non-empty
  domain_name      Pfam name if any, else UniProt domain name
AlphaMissense per-residue summary (am_mean, am_max over the 19 substitutions; am_pos = the
AlphaMissense residue index after mapping AF->AM sequence, see am_mapping) is merged if
INTER/alphamissense_subset.parquet exists (run s08 first, or rerun this step after s08).
"""
import sys, json, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from common import INTER, DERIVED

RFDIR = DERIVED / "residue_features"; RFDIR.mkdir(parents=True, exist_ok=True)


def annotate(pos, ft):
    L = len(pos)
    cols = {k: [[] for _ in range(L)] for k in
            ("uniprot_domain", "uniprot_region", "uniprot_disordered_region", "uniprot_topology", "uniprot_motif", "pfam_domain")}
    idx = {p: i for i, p in enumerate(pos)}
    for r in ft.itertuples():
        if r.source == "Pfam":
            key = "pfam_domain"
        elif r.type in ("Domain", "Repeat", "Zinc finger", "DNA binding"):
            key = "uniprot_domain"
        elif r.type == "Region":
            key = "uniprot_disordered_region" if (r.description or "").strip().lower() == "disordered" else "uniprot_region"
        elif r.type in ("Topological domain", "Transmembrane", "Intramembrane"):
            key = "uniprot_topology"
        else:
            key = "uniprot_motif"
        label = r.description or r.type
        if r.type in ("Transmembrane", "Intramembrane") and key == "uniprot_topology":
            label = f"{r.type}:{label}"
        for p in range(int(r.start), int(r.end) + 1):
            i = idx.get(p)
            if i is not None and label not in cols[key][i]:
                cols[key][i].append(label)
    out = {k: [";".join(x) if x else None for x in v] for k, v in cols.items()}
    out["uniprot_disordered_region"] = [x is not None for x in out["uniprot_disordered_region"]]
    out["in_domain"] = [bool(a or b) for a, b in zip(out["uniprot_domain"], out["pfam_domain"])]
    out["in_domain_broad"] = [bool(a or b or c) for a, b, c in zip(out["in_domain"], out["uniprot_region"], out["uniprot_motif"])]
    out["domain_name"] = [(b.split(";")[0] if b else (a.split(";")[0] if a else None)) for a, b in zip(out["uniprot_domain"], out["pfam_domain"])]
    return out


def am_mapping(acc, af_seq, other_hits, am):
    """AlphaMissense accession + {AF pos: AM pos} (None = identity) for one protein."""
    from s06_select_genes import align_map
    cands = [acc] + [x for x in (other_hits or "").split(",") if x]
    for c in cands:
        sub = am[am.uniprot_id == c].drop_duplicates("pos").sort_values("pos")
        if not len(sub):
            continue
        seq = ["X"] * int(sub.pos.max())
        for p_, r_ in zip(sub.pos, sub.ref):
            seq[p_ - 1] = r_
        seq = "".join(seq)
        if seq == af_seq:
            return c, None
        m, idn = align_map(af_seq, seq)
        return c, m
    return None, None


def one(args):
    gene, acc, path, useq = args
    import structfeat as SF
    t0 = time.time()
    feats, meta, _ = SF.compute_residue_features(path)
    df = pd.DataFrame(feats)
    df.insert(0, "uniprot", acc); df.insert(0, "gene", gene)
    df["aa_uniprot"] = [useq[p - 1] if isinstance(useq, str) and p <= len(useq) else None for p in df.pos]
    meta.update(gene=gene, uniprot=acc, af_pdb=path, seconds=time.time() - t0)
    return df, meta


def main():
    allgenes = "--all" in sys.argv
    g = pd.read_parquet(INTER / "gene_summary.parquet")
    sel = g[(g.af_n_fragments == 1) & g.af_pdb_exists] if allgenes else g[g.selected]
    ft = pd.read_parquet(INTER / "uniprot_features.parquet")
    amf = INTER / "alphamissense_subset.parquet"
    am = pd.read_parquet(amf) if amf.exists() else None
    jobs = [(r.gene, r.uniprot, r.af_pdb, r.sequence) for r in sel.itertuples()]
    other = dict(zip(sel.uniprot, sel.other_hits))
    metas = []
    with ProcessPoolExecutor(4) as ex:
        for df, meta in ex.map(one, jobs):
            acc = meta["uniprot"]
            ann = annotate(df.pos.tolist(), ft[ft.uniprot == acc])
            for k, v in ann.items():
                df[k] = v
            if am is not None:
                am_acc, posmap = am_mapping(acc, "".join(df.aa_af), other.get(acc), am)
                meta["am_accession"] = am_acc
                meta["am_mapping"] = None if am_acc is None else ("identity" if posmap is None else "alignment")
                if am_acc is not None:
                    a = am[am.uniprot_id == am_acc]
                    s = a.groupby("pos").am_pathogenicity.agg(["mean", "max", "count"])
                    ap = df.pos if posmap is None else df.pos.map(posmap)
                    df["am_pos"] = ap
                    df["am_mean"] = ap.map(s["mean"]); df["am_max"] = ap.map(s["max"]); df["am_n_subst"] = ap.map(s["count"])
            df.to_parquet(RFDIR / f"{acc}.parquet", index=False)
            meta.pop("af_sequence", None)
            metas.append(meta)
            print(f"{meta['gene']:10s} {acc} L={meta['L']} bb={meta['backbone_size']}/{meta['backbone_strict_size']} "
                  f"pairs_ok={meta['backbone_pairs_connected']}/{meta['backbone_pairs']} {meta['seconds']:.1f}s", flush=True)
    pd.DataFrame(metas).to_parquet(INTER / "structure_meta.parquet", index=False)


if __name__ == "__main__":
    main()
