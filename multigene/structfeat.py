"""Per-residue structural-context features from an AlphaFold model.

Fingerprint definitions follow mechanistic_mvp.py exactly (contact graph 8 A Calpha incl.
all residues; Core = pLDDT >= q0.9, SemiCore-only = q0.5 <= pLDDT < q0.9 (semi minus core),
Disorder = pLDDT <= q0.5; backbone = top-200 shortest-path hubs between N and C anchors on the
pLDDT>=70 subgraph, 2000 sampled pairs, seed 0).  The repo's networkx functions are used for the
graph/regions/backbone so results are identical; only RWR is re-implemented (exact linear solve).
"""
import gzip, io
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.DSSP import residue_max_acc
from Bio.PDB.Polypeptide import three_to_index, index_to_one

from common import import_mvp

M = import_mvp()
ALPHAS = (0.3, 0.15, 0.5)
MAXACC = residue_max_acc["Wilke"]  # Tien et al. 2013 theoretical max ASA


def read_structure(path):
    p = PDBParser(QUIET=True)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return p.get_structure("AF", io.StringIO(fh.read()))
    return p.get_structure("AF", str(path))


def load_ca(structure):
    """Same logic as mechanistic_mvp.load_alphafold_ca, from a parsed structure (first model, first chain)."""
    model = next(structure.get_models())
    chain = next(model.get_chains())
    pos, xyz, pl, aa = [], [], [], []
    for res in chain.get_residues():
        if res.id[0].strip() or "CA" not in res:
            continue
        pos.append(int(res.id[1])); xyz.append(res["CA"].coord.astype(float)); pl.append(float(res["CA"].bfactor))
        try:
            aa.append(index_to_one(three_to_index(res.get_resname())))
        except Exception:
            aa.append("X")
    o = np.argsort(pos)
    return np.array(pos)[o], np.array(xyz)[o], np.array(pl)[o], "".join(np.array(aa)[o])


def transition_matrix(G, nodes):
    """W = A D^-1 (column-stochastic for non-isolated nodes), deg = max(1, deg) as in repo rwr()."""
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    r, c = [], []
    for u, v in G.edges():
        r += [idx[u], idx[v]]; c += [idx[v], idx[u]]
    A = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(n, n))
    deg = np.maximum(1, np.asarray(A.sum(axis=0)).ravel())
    return A @ sp.diags(1.0 / deg), A


def rwr_iterative(W, s, alpha, max_iter=200, tol=1e-10):
    """Vectorised version of mechanistic_mvp.rwr (identical update / stopping / normalisation)."""
    n = W.shape[0]
    p0 = np.zeros(n); p0[s] = 1.0
    p = p0.copy()
    for _ in range(max_iter):
        p_new = (1 - alpha) * (W @ p) + alpha * p0
        if np.abs(p_new - p).sum() < tol:
            p = p_new; break
        p = p_new
    t = p.sum()
    return p / t if t > 0 else p


def rwr_solve_all(W, alpha):
    """Return LU of M^T where M = I - (1-alpha) W ; p_s = alpha M^-1 e_s."""
    n = W.shape[0]
    Mt = (sp.identity(n, format="csc") - (1 - alpha) * W).T.tocsc()
    return splu(Mt)


def region_impacts(lu, alpha, masks):
    """For every seed s: sum_{n in region} p_s[n] / sum_n p_s[n].  Uses M^T w = alpha*1_region."""
    n = next(iter(masks.values())).shape[0]
    tot = lu.solve(alpha * np.ones(n))
    return {k: lu.solve(alpha * m.astype(float)) / tot for k, m in masks.items()}


def anchors(pos):
    L = int(pos.max())
    if L >= 600:
        return (1, 200), (L - 200 + 1, L), "repo_1-200_lastL-199"
    k = L // 3
    return (1, k), (L - k + 1, L), "short_protein_first_last_third"


def compute_residue_features(path, uniprot_seq=None):
    st = read_structure(path)
    pos, xyz, pl, seq = load_ca(st)
    L = len(pos)
    G_all = M.build_confidence_graph(xyz, pos, pl, cutoff_A=8.0, min_plddt_node=0.0)
    G_core = M.build_confidence_graph(xyz, pos, pl, cutoff_A=8.0, min_plddt_node=M.CORE_MIN_PLD)
    q_core = float(np.quantile(pl, M.CORE_Q)); q_semi = float(np.quantile(pl, M.SEMI_Q))
    core = M.region_nodes(G_all, plddt_lo=q_core)
    semi = M.region_nodes(G_all, plddt_lo=q_semi)
    dis = M.region_nodes(G_all, plddt_hi=q_semi)
    (nlo, nhi), (clo, chi), anchor_rule = anchors(pos)
    NA = M.nodes_by_unip_range(G_core, nlo, nhi)
    CA = M.nodes_by_unip_range(G_core, clo, chi)
    bb = M.shortest_path_backbone(G_core, NA, CA, max_pairs=M.BACKBONE_MAX_PAIRS, top_k=M.BACKBONE_TOPK, seed=0)
    # path counts again (to know which backbone nodes have non-zero betweenness counts) - same rng/pairs
    import random, networkx as nx
    counts = {}
    if NA and CA:
        rng = random.Random(0)
        pairs = [(s, t) for s in list(NA) for t in list(CA)]
        if len(pairs) > M.BACKBONE_MAX_PAIRS:
            pairs = rng.sample(pairs, k=M.BACKBONE_MAX_PAIRS)
        n_ok = 0
        for s, t in pairs:
            try:
                path_ = nx.shortest_path(G_core, s, t); n_ok += 1
            except nx.NetworkXNoPath:
                continue
            for x in path_[1:-1]:
                counts[x] = counts.get(x, 0) + 1
        n_pairs = len(pairs)
    else:
        n_ok = n_pairs = 0
    bb_strict = {x for x in bb if counts.get(x, 0) > 0}

    nodes = list(G_all.nodes())
    assert nodes == list(pos), "node order must equal residue order"
    W, A = transition_matrix(G_all, nodes)
    degree = np.asarray(A.sum(axis=0)).ravel().astype(int)
    S = set
    masks = {
        "CoreImpact": np.isin(pos, list(core)),
        "SemiCoreImpact": np.isin(pos, list(S(semi) - S(core))),
        "DisorderImpact": np.isin(pos, list(dis)),
        "CommBackboneImpact": np.isin(pos, list(bb)),
        "CommBackboneImpact_strict": np.isin(pos, list(bb_strict)),
    }
    feats = {"pos": pos, "aa_af": list(seq), "plddt": pl, "degree_8A": degree}
    for a in ALPHAS:
        lu = rwr_solve_all(W, a)
        imp = region_impacts(lu, a, masks)
        sfx = "" if a == 0.3 else f"_a{str(a).replace('0.', '')}"
        for k, v in imp.items():
            feats[k + sfx] = v
    # RSA (Shrake-Rupley, residue level, isolated monomer)
    sr = ShrakeRupley()
    model = next(st.get_models())
    sr.compute(model, level="R")
    chain = next(model.get_chains())
    rsa_map = {}
    for res in chain.get_residues():
        if res.id[0].strip():
            continue
        mx = MAXACC.get(res.get_resname())
        rsa_map[int(res.id[1])] = (res.sasa, res.sasa / mx if mx else np.nan)
    feats["sasa"] = np.array([rsa_map[p][0] for p in pos])
    feats["rsa"] = np.array([rsa_map[p][1] for p in pos])
    feats["rel_pos"] = pos / float(L)
    hi = pos[pl >= 70]
    if len(hi):
        tr = cKDTree(hi.reshape(-1, 1).astype(float))
        feats["dist_to_plddt70_res"] = tr.query(pos.reshape(-1, 1).astype(float))[0].astype(int)
    else:
        feats["dist_to_plddt70_res"] = np.full(L, -1)
    feats["region"] = np.where(masks["CoreImpact"], "core", np.where(np.isin(pos, list(semi)), "semicore", "disorder"))
    # note: residues with pLDDT exactly == q0.5 are in both semi and disorder sets in the repo; region label = semicore
    feats["in_backbone"] = masks["CommBackboneImpact"]
    meta = dict(L=L, q_core=q_core, q_semi=q_semi, n_core=len(core), n_semi_only=len(S(semi) - S(core)),
                n_disorder=len(dis), G_all_edges=G_all.number_of_edges(), G_core_nodes=G_core.number_of_nodes(),
                G_core_edges=G_core.number_of_edges(), anchor_rule=anchor_rule, N_anchor=(nlo, nhi), C_anchor=(clo, chi),
                n_N_anchor_nodes=len(NA), n_C_anchor_nodes=len(CA), backbone_pairs=n_pairs, backbone_pairs_connected=n_ok,
                backbone_size=len(bb), backbone_strict_size=len(bb_strict), mean_plddt=float(pl.mean()),
                frac_plddt_lt50=float((pl < 50).mean()), frac_plddt_ge70=float((pl >= 70).mean()), af_sequence=seq)
    return feats, meta, (G_all, nodes, W, core, semi, dis, bb)
