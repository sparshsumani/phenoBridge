#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mechanistic_mvp.py

import argparse
import json
import math
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser

import gradio as gr
import matplotlib.pyplot as plt

import statsmodels.api as sm


# -------------------------
# Config
# -------------------------
DEFAULT_CUTOFF_A = 8.0          # contact threshold (Å)
DEFAULT_ALPHA = 0.3            # restart prob in RWR
DEFAULT_MAX_ITER = 200
DEFAULT_TOL = 1e-10

# pLDDT quantile definitions for core/semi/disorder
CORE_Q = 0.90                  # top 10% as "core"
SEMI_Q = 0.50                  # > median = semi-core; <= median = disorder

# backbone construction on G_core (high confidence edges)
CORE_MIN_PLD = 70.0            # "structured core residues only" graph
ANCHOR_N_RANGE = (1, 200)
ANCHOR_C_TAIL = 200            # last ~200 aa as C anchor
BACKBONE_TOPK = 200
BACKBONE_MAX_PAIRS = 2000

# phenotype evidence
DEFAULT_MIN_PREV = 0.05        # min prevalence to fit HPO model
SIG_P_THRESHOLD = 0.05         # badge threshold (p < 0.05)


# -------------------------
# PDB / AlphaFold parsing
# -------------------------
def load_alphafold_ca(pdb_path: str, chain_id: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns:
      unip_pos: (L,) int positions
      ca_xyz:   (L,3) float CA coordinates
      plddt:    (L,) float from B-factor (AlphaFold stores pLDDT there)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("AF", pdb_path)
    model = next(structure.get_models())

    if chain_id is None:
        chain = next(model.get_chains())
    else:
        chain = model[chain_id]

    unip_pos, ca_xyz, plddt = [], [], []
    for res in chain.get_residues():
        if res.id[0].strip():
            continue
        if "CA" not in res:
            continue
        ca = res["CA"]
        pos = int(res.id[1])
        unip_pos.append(pos)
        ca_xyz.append(ca.coord.astype(float))
        plddt.append(float(ca.bfactor))

    unip_pos = np.array(unip_pos, dtype=int)
    ca_xyz = np.array(ca_xyz, dtype=float)
    plddt = np.array(plddt, dtype=float)

    order = np.argsort(unip_pos)
    return unip_pos[order], ca_xyz[order], plddt[order]


def build_confidence_graph(
    ca_xyz: np.ndarray,
    unip_pos: np.ndarray,
    plddt: np.ndarray,
    cutoff_A: float = DEFAULT_CUTOFF_A,
    min_plddt_node: float = 0.0,
) -> nx.Graph:
    """
    Nodes are residue positions (unip_pos).
    Edge exists if CA distance <= cutoff_A.
    Nodes filtered by min_plddt_node.
    """
    keep = plddt >= float(min_plddt_node)
    xyz = ca_xyz[keep]
    pos = unip_pos[keep]
    conf = plddt[keep]

    G = nx.Graph()
    for p, c in zip(pos, conf):
        G.add_node(int(p), unip_pos=int(p), plddt=float(c))

    if len(pos) == 0:
        return G

    tree = cKDTree(xyz)
    pairs = tree.query_pairs(r=float(cutoff_A))
    for i_small, j_small in pairs:
        i = int(pos[i_small])
        j = int(pos[j_small])
        dist = float(np.linalg.norm(xyz[i_small] - xyz[j_small]))
        G.add_edge(i, j, dist=dist)
    return G


def region_nodes(G: nx.Graph, plddt_lo: Optional[float] = None, plddt_hi: Optional[float] = None) -> set:
    out = set()
    for n, a in G.nodes(data=True):
        c = float(a.get("plddt", 0.0))
        if plddt_lo is not None and c < plddt_lo:
            continue
        if plddt_hi is not None and c > plddt_hi:
            continue
        out.add(n)
    return out


# -------------------------
# Backbone (shortest-path hubs)
# -------------------------
def nodes_by_unip_range(G: nx.Graph, lo: int, hi: int) -> set:
    return {n for n, a in G.nodes(data=True) if lo <= int(a["unip_pos"]) <= hi}


def shortest_path_backbone(
    G: nx.Graph,
    sources: set,
    targets: set,
    max_pairs: int = BACKBONE_MAX_PAIRS,
    top_k: int = BACKBONE_TOPK,
    seed: int = 0,
) -> set:
    rng = random.Random(seed)
    sources = list(sources)
    targets = list(targets)
    if not sources or not targets:
        return set()

    pairs = [(s, t) for s in sources for t in targets]
    if len(pairs) > max_pairs:
        pairs = rng.sample(pairs, k=max_pairs)

    counts = {n: 0 for n in G.nodes()}
    for s, t in pairs:
        try:
            path = nx.shortest_path(G, s, t)
        except nx.NetworkXNoPath:
            continue
        for n in path[1:-1]:
            counts[n] += 1

    backbone = sorted(counts, key=lambda n: counts[n], reverse=True)[:top_k]
    return set(backbone)


# -------------------------
# RWR
# -------------------------
def rwr(
    G: nx.Graph,
    seed_node: int,
    alpha: float = DEFAULT_ALPHA,
    max_iter: int = DEFAULT_MAX_ITER,
    tol: float = DEFAULT_TOL
) -> Dict[int, float]:
    """Random Walk with Restart on an undirected graph. Returns node->prob."""
    if seed_node not in G:
        return {}

    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    nbrs = [list(G.neighbors(node)) for node in nodes]
    deg = np.array([max(1, len(nbrs[i])) for i in range(n)], dtype=float)

    p0 = np.zeros(n, dtype=float)
    p0[idx[seed_node]] = 1.0
    p = p0.copy()

    for _ in range(int(max_iter)):
        p_new = np.zeros(n, dtype=float)
        for i in range(n):
            if p[i] == 0:
                continue
            share = (1.0 - alpha) * p[i] / deg[i]
            for nb in nbrs[i]:
                j = idx[nb]
                p_new[j] += share
        p_new += alpha * p0

        if np.linalg.norm(p_new - p, ord=1) < tol:
            p = p_new
            break
        p = p_new

    s = float(p.sum())
    if s > 0:
        p = p / s
    return {nodes[i]: float(p[i]) for i in range(n)}


# -------------------------
# Fingerprint model
# -------------------------
@dataclass
class MechanisticModel:
    pdb_path: str
    chain_id: Optional[str]
    cutoff_A: float
    G_all: nx.Graph
    G_core: nx.Graph
    plddt_by_pos: Dict[int, float]
    core_nodes: set
    semi_nodes: set
    disorder_nodes: set
    comm_backbone: set
    alpha: float


def build_model(pdb_path: str, chain_id: Optional[str], cutoff_A: float, alpha: float) -> MechanisticModel:
    unip_pos, ca_xyz, plddt = load_alphafold_ca(pdb_path, chain_id=chain_id)

    G_all = build_confidence_graph(ca_xyz, unip_pos, plddt, cutoff_A=cutoff_A, min_plddt_node=0.0)
    G_core = build_confidence_graph(ca_xyz, unip_pos, plddt, cutoff_A=cutoff_A, min_plddt_node=CORE_MIN_PLD)

    plddt_by_pos = {int(p): float(c) for p, c in zip(unip_pos.tolist(), plddt.tolist())}

    q_core = float(np.quantile(plddt, CORE_Q))
    q_semi = float(np.quantile(plddt, SEMI_Q))

    core_nodes = region_nodes(G_all, plddt_lo=q_core)
    semi_nodes = region_nodes(G_all, plddt_lo=q_semi)
    disorder_nodes = region_nodes(G_all, plddt_hi=q_semi)

    n_lo, n_hi = ANCHOR_N_RANGE
    c_lo = int(unip_pos.max()) - int(ANCHOR_C_TAIL) + 1
    c_hi = int(unip_pos.max())
    N_anchor = nodes_by_unip_range(G_core, n_lo, n_hi)
    C_anchor = nodes_by_unip_range(G_core, c_lo, c_hi)
    comm_backbone = shortest_path_backbone(
        G_core, N_anchor, C_anchor,
        max_pairs=BACKBONE_MAX_PAIRS,
        top_k=BACKBONE_TOPK,
        seed=0
    )

    return MechanisticModel(
        pdb_path=pdb_path,
        chain_id=chain_id,
        cutoff_A=cutoff_A,
        G_all=G_all,
        G_core=G_core,
        plddt_by_pos=plddt_by_pos,
        core_nodes=core_nodes,
        semi_nodes=semi_nodes,
        disorder_nodes=disorder_nodes,
        comm_backbone=comm_backbone,
        alpha=alpha,
    )


def fingerprint_for_variant(model: MechanisticModel, aa_pos: int) -> Optional[Dict[str, float]]:
    aa_pos = int(aa_pos)
    if aa_pos not in model.G_all:
        return None

    p = rwr(model.G_all, seed_node=aa_pos, alpha=model.alpha)

    core = model.core_nodes
    semi_only = set(model.semi_nodes) - set(model.core_nodes)
    disorder = model.disorder_nodes

    core_imp = sum(p.get(n, 0.0) for n in core)
    semi_imp = sum(p.get(n, 0.0) for n in semi_only)
    dis_imp = sum(p.get(n, 0.0) for n in disorder)

    bb = model.comm_backbone
    bb_imp = sum(p.get(n, 0.0) for n in bb)

    seed_plddt = float(model.plddt_by_pos.get(aa_pos, np.nan))

    return {
        "aa_pos": float(aa_pos),
        "seed_pLDDT": seed_plddt,
        "CoreImpact": float(core_imp),
        "SemiCoreImpact": float(semi_imp),
        "DisorderImpact": float(dis_imp),
        "CommBackboneImpact": float(bb_imp),
        "alpha": float(model.alpha),
        "cutoff_A": float(model.cutoff_A),
    }


# -------------------------
# Phenotype add-on (p only, no q)
# -------------------------
TYPE_TO_VALUE = {"Observed": 1.0, "Excluded": 0.0, "Na": np.nan, "NA": np.nan, None: np.nan}

def load_kbgs_json(json_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, str]]:
    data = json.loads(Path(json_path).read_text())

    headers = pd.DataFrame(data.get("hpoHeaders", []))
    hpo_ids = headers["hpoId"].tolist()
    hpo_label_map = dict(zip(headers["hpoId"], headers["hpoLabel"]))

    pheno_rows = []
    for r in data.get("rows", []):
        ind_id = (r.get("individualData") or {}).get("individualId")
        hpo_data = r.get("hpoData", []) or []
        if len(hpo_data) < len(hpo_ids):
            hpo_data = hpo_data + [{"type": "Na"}] * (len(hpo_ids) - len(hpo_data))
        elif len(hpo_data) > len(hpo_ids):
            hpo_data = hpo_data[:len(hpo_ids)]

        row = {"individualId": ind_id}
        for i, cell in enumerate(hpo_data):
            t = (cell or {}).get("type")
            row[hpo_ids[i]] = TYPE_TO_VALUE.get(t, np.nan)
        pheno_rows.append(row)
    df_pheno = pd.DataFrame(pheno_rows)

    allele_rows = []
    for r in data.get("rows", []):
        ind_id = (r.get("individualData") or {}).get("individualId")
        for allele_key, count in (r.get("alleleCountMap") or {}).items():
            allele_rows.append({"individualId": ind_id, "alleleKey": str(allele_key), "alleleCount": int(count)})
    df_mut = pd.DataFrame(allele_rows).dropna(subset=["alleleKey"]).drop_duplicates()

    return df_pheno, df_mut, hpo_label_map


def allelekey_to_aa_pos(alleleKey: str) -> float:
    s = str(alleleKey)
    m = re.search(r"p\.?\(?[A-Za-z]{1,3}(\d+)[A-Za-z]{1,3}\)?", s)
    if m:
        return float(int(m.group(1)))
    m = re.search(r"c\.?(\d+)", s)
    if m:
        return float(int(math.ceil(int(m.group(1)) / 3.0)))
    return np.nan


def is_truncating(s: str) -> bool:
    return bool(re.search(r"(stop|ter|\*|nonsense|frameshift|fs|splice|del|dup|ins)", str(s), flags=re.I))


def parse_mutation_tokens(text: str) -> List[str]:
    if text is None:
        return []
    toks = re.split(r"[\s,;]+", text.strip())
    return [t for t in toks if t]


def fit_hpo_models_from_cohort(
    model: MechanisticModel,
    json_path: str,
    min_prev: float = DEFAULT_MIN_PREV,
    use_trunc_frac: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, float]], Dict[str, str], Dict[str, str]]:
    """
    Fit cohort evidence (p only): y_hpo ~ logSemiCore (+ TruncFrac).
    Returns:
      evidence_df: phenotype, label, p, n
      params: phenotype -> {const, logSemiCore, TruncFrac?}
      hpo_label_map
      diagnostics
    """
    df_pheno, df_mut, hpo_label_map = load_kbgs_json(json_path)

    dfv = df_mut.copy()
    dfv["aa_pos"] = dfv["alleleKey"].map(allelekey_to_aa_pos)
    dfv = dfv.dropna(subset=["aa_pos"]).copy()
    dfv["aa_pos"] = dfv["aa_pos"].astype(int)
    dfv["TRUNC_like"] = dfv["alleleKey"].map(is_truncating).astype(int)

    fps = dfv["aa_pos"].map(lambda p: fingerprint_for_variant(model, int(p)))
    dfv["SemiCoreImpact"] = fps.map(lambda d: 0.0 if d is None else float(d["SemiCoreImpact"]))

    mech_patient = dfv.groupby("individualId")[["SemiCoreImpact"]].max().reset_index()

    if use_trunc_frac:
        trunc_patient = dfv.groupby("individualId")["TRUNC_like"].mean().reset_index().rename(columns={"TRUNC_like": "TruncFrac"})
    else:
        trunc_patient = pd.DataFrame({"individualId": df_pheno["individualId"].unique(), "TruncFrac": 0.0})

    df_mech = (
        df_pheno.merge(mech_patient, on="individualId", how="left")
        .merge(trunc_patient, on="individualId", how="left")
        .fillna(0.0)
    )

    EPS = 1e-15
    df_mech["logSemiCore"] = np.log10(df_mech["SemiCoreImpact"] + EPS)

    hpo_cols = [c for c in df_pheno.columns if c != "individualId"]

    keep = []
    for h in hpo_cols:
        yv = df_mech[h].dropna()
        if len(yv) == 0:
            continue
        prev = float((yv == 1.0).mean())
        if min_prev <= prev <= 1.0 - min_prev:
            keep.append(h)

    rows = []
    params: Dict[str, Dict[str, float]] = {}

    for h in keep:
        y = df_mech[h]
        msk = y.isin([0.0, 1.0])
        y = y[msk].astype(int)

        X = df_mech.loc[msk, ["logSemiCore"]].copy()
        if use_trunc_frac:
            X["TruncFrac"] = df_mech.loc[msk, "TruncFrac"].astype(float)
        X = sm.add_constant(X, has_constant="add")

        try:
            res = sm.Logit(y, X).fit(disp=0)
        except Exception:
            continue

        pval = float(res.pvalues.get("logSemiCore", np.nan))
        rows.append({"phenotype": h, "label": hpo_label_map.get(h, ""), "p": pval, "n": int(len(y))})
        params[h] = {k: float(v) for k, v in res.params.items()}

    evidence = pd.DataFrame(rows).sort_values("p").reset_index(drop=True) if len(rows) else pd.DataFrame(rows)

    diagnostics = {
        "n_samples": str(df_pheno.shape[0]),
        "n_unique_aa_pos_mapped": str(int(dfv["aa_pos"].nunique())),
        "n_hpo_models": str(int(len(evidence))),
    }
    return evidence, params, hpo_label_map, diagnostics


def predict_top10_sig_table(
    params: Dict[str, Dict[str, float]],
    evidence: pd.DataFrame,
    hpo_label_map: Dict[str, str],
    logSemiCore: float,
    truncFrac: float,
    topk: int = 10,
) -> pd.DataFrame:
    out = []
    for h, p in params.items():
        z = p.get("const", 0.0) + p.get("logSemiCore", 0.0) * float(logSemiCore)
        if "TruncFrac" in p:
            z += p.get("TruncFrac", 0.0) * float(truncFrac)
        prob = 1.0 / (1.0 + math.exp(-z))
        out.append({"phenotype": h, "label": hpo_label_map.get(h, ""), "pred_prob": float(prob)})

    pred = pd.DataFrame(out).sort_values("pred_prob", ascending=False).reset_index(drop=True)

    if evidence is not None and len(evidence) > 0:
        ev = evidence[["phenotype", "p"]].copy()
        pred = pred.merge(ev, on="phenotype", how="left")
        pred["sig"] = pred["p"].map(lambda pv: "✅" if (pd.notna(pv) and pv < SIG_P_THRESHOLD) else "")
    else:
        pred["p"] = np.nan
        pred["sig"] = ""

    show = pred[["sig", "label", "phenotype", "pred_prob", "p"]].head(topk).copy()
    return show


# -------------------------
# Plot utilities
# -------------------------
def plot_fp(fp1: Optional[Dict[str, float]], fp2: Optional[Dict[str, float]] = None):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    keys = ["CoreImpact", "SemiCoreImpact", "DisorderImpact", "CommBackboneImpact"]
    x = np.arange(len(keys))

    if fp1 is not None:
        y1 = [fp1.get(k, 0.0) for k in keys]
        ax.bar(x - 0.2, y1, width=0.4, label=f"aa_pos={int(fp1['aa_pos'])}")

    if fp2 is not None:
        y2 = [fp2.get(k, 0.0) for k in keys]
        ax.bar(x + 0.2, y2, width=0.4, label=f"aa_pos={int(fp2['aa_pos'])}")

    ax.set_xticks(x)
    ax.set_xticklabels(keys, rotation=20, ha="right")
    ax.set_ylabel("Influence mass (RWR)")
    ax.legend()
    fig.tight_layout()
    return fig


# -------------------------
# Gradio app
# -------------------------
def make_app(
    model: MechanisticModel,
    evidence_df: Optional[pd.DataFrame] = None,
    hpo_params: Optional[Dict[str, Dict[str, float]]] = None,
    hpo_label_map: Optional[Dict[str, str]] = None,
):
    cache: Dict[int, Dict[str, float]] = {}

    def compute_fingerprint(aa_pos: int) -> Tuple[str, str]:
        aa_pos = int(aa_pos)
        if aa_pos not in cache:
            fp = fingerprint_for_variant(model, aa_pos)
            if fp is None:
                return "", f"aa_pos {aa_pos} not found in graph (check PDB residue numbering)."
            cache[aa_pos] = fp
        fp = cache[aa_pos]
        summary = (
            f"aa_pos={int(fp['aa_pos'])} | seed pLDDT={fp['seed_pLDDT']:.1f} | "
            f"alpha={fp['alpha']} | cutoff_A={fp['cutoff_A']}"
        )
        return json.dumps(fp, indent=2), summary

    def compare_fingerprints(pos1: int, pos2: int):
        fp1 = fingerprint_for_variant(model, int(pos1)) if pos1 is not None else None
        fp2 = fingerprint_for_variant(model, int(pos2)) if pos2 is not None else None

        msg = []
        if fp1 is None:
            msg.append(f"aa_pos {pos1} not found in graph.")
        if fp2 is None:
            msg.append(f"aa_pos {pos2} not found in graph.")
        note = "\n".join(msg) if msg else "OK"

        fig = plot_fp(fp1, fp2)
        return fig, note

    def phenotype_from_mutations(text: str):
        if hpo_params is None or hpo_label_map is None:
            return pd.DataFrame([{
                "sig": "",
                "label": "Provide --json to enable",
                "phenotype": "(disabled)",
                "pred_prob": np.nan,
                "p": np.nan
            }])

        toks = parse_mutation_tokens(text)
        aa_positions = []
        trunc_flags = []

        for t in toks:
            if re.fullmatch(r"\d+", t):
                aa_positions.append(int(t))
                trunc_flags.append(0)
            else:
                ap = allelekey_to_aa_pos(t)
                if not pd.isna(ap):
                    aa_positions.append(int(ap))
                trunc_flags.append(1 if is_truncating(t) else 0)

        semi_vals = []
        for p in aa_positions:
            fp = fingerprint_for_variant(model, int(p))
            if fp is not None:
                semi_vals.append(float(fp["SemiCoreImpact"]))
        semi = max(semi_vals) if len(semi_vals) else 0.0

        EPS = 1e-15
        logSemiCore = float(np.log10(semi + EPS))
        truncFrac = float(np.mean(trunc_flags)) if len(trunc_flags) else 0.0

        return predict_top10_sig_table(
            params=hpo_params,
            evidence=evidence_df,
            hpo_label_map=hpo_label_map,
            logSemiCore=logSemiCore,
            truncFrac=truncFrac,
            topk=10,
        )

    with gr.Blocks(title="Mechanistic MVP: fingerprint + phenotype (p evidence)") as demo:
        gr.Markdown("# Mechanistic MVP: 1) Fingerprint  2) Compare  3) Phenotype prediction")

        gr.Markdown("## 1) Compute fingerprint from one residue (aa_pos)")
        with gr.Row():
            aa_pos_in = gr.Number(value=600, label="Input aa_pos")
            btn_fp = gr.Button("Run fingerprint")

        fp_json = gr.Code(label="Fingerprint JSON", language="json")
        fp_summary = gr.Textbox(label="Fingerprint summary", lines=2)
        btn_fp.click(compute_fingerprint, inputs=[aa_pos_in], outputs=[fp_json, fp_summary])

        gr.Markdown("## 2) Compare fingerprints (two residues)")
        with gr.Row():
            aa1 = gr.Number(value=600, label="Input aa_pos #1")
            aa2 = gr.Number(value=2400, label="Input aa_pos #2")
            btn_cmp = gr.Button("Run compare")

        cmp_plot = gr.Plot(label="Fingerprint comparison plot")
        cmp_note = gr.Textbox(label="Compare notes", lines=2)
        btn_cmp.click(compare_fingerprints, inputs=[aa1, aa2], outputs=[cmp_plot, cmp_note])

        gr.Markdown("## 3) Predict phenotypes from mutation(s) (Top 10)")
        gr.Markdown(f"- Table uses **pred_prob** (per-input) and **p** (cohort-level evidence). Badge ✅ means p < {SIG_P_THRESHOLD}.")
        mut_text = gr.Textbox(
            label="Input mutation(s) (aa_pos or alleleKey tokens)",
            placeholder="Examples: 600 2400  OR  p.Arg345Cys  OR  c.1903_1907del ...",
            lines=1
        )
        btn_pheno = gr.Button("Run phenotype prediction")
        pheno_tbl = gr.Dataframe(label="sig | label | phenotype | pred_prob | p", interactive=False)
        btn_pheno.click(phenotype_from_mutations, inputs=[mut_text], outputs=[pheno_tbl])

    return demo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=str, required=True, help="Path to AlphaFold PDB (e.g., AF-Q6UB99-F1-model_v6.pdb)")
    ap.add_argument("--chain", type=str, default=None, help="Chain ID (default: first chain)")
    ap.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF_A, help="Contact cutoff in Å (default: 8.0)")
    ap.add_argument("--alpha", type=float, default=DEFAULT_ALPHA, help="RWR restart probability (default: 0.3)")
    ap.add_argument("--json", type=str, default=None, help="Optional cohort JSON (ANKRD11_KBGS_individuals.json) to enable phenotype table")
    ap.add_argument("--min_prev", type=float, default=DEFAULT_MIN_PREV, help="Min prevalence to fit HPO models (default: 0.05)")
    ap.add_argument("--share", action="store_true", help="Gradio share link")
    args = ap.parse_args()

    model = build_model(args.pdb, chain_id=args.chain, cutoff_A=args.cutoff, alpha=args.alpha)

    print("[INFO] Loaded PDB:", args.pdb)
    print("[INFO] G_all nodes/edges:", model.G_all.number_of_nodes(), model.G_all.number_of_edges())
    print("[INFO] G_core nodes/edges:", model.G_core.number_of_nodes(), model.G_core.number_of_edges())
    print("[INFO] Comm backbone size:", len(model.comm_backbone))

    evidence_df = None
    hpo_params = None
    hpo_label_map = None

    if args.json:
        evidence_df, hpo_params, hpo_label_map, diag = fit_hpo_models_from_cohort(
            model=model,
            json_path=args.json,
            min_prev=args.min_prev,
            use_trunc_frac=True
        )
        print("[INFO] Phenotype evidence diagnostics:", diag)
        if evidence_df is not None and len(evidence_df) > 0:
            print("[INFO] Evidence table head (p only):")
            print(evidence_df.head(10)[["phenotype", "label", "p", "n"]].to_string(index=False))

    app = make_app(model, evidence_df=evidence_df, hpo_params=hpo_params, hpo_label_map=hpo_label_map)
    app.launch(share=args.share)  # no show_api (gradio version compatibility)


if __name__ == "__main__":
    main()