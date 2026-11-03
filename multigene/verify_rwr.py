#!/usr/bin/env python3
"""Verify the fast RWR (vectorised iteration and exact sparse-LU solve) against the repo's
pure-Python mechanistic_mvp.rwr / fingerprint_for_variant on ANKRD11 (repo PDB)."""
import sys, json
import numpy as np
from common import REPO
import structfeat as SF

M = SF.M
pdb = REPO / "AF-Q6UB99-F1-model_v6.pdb"
feats, meta, (G_all, nodes, W, core, semi, dis, bb) = SF.compute_residue_features(pdb)
model = M.build_model(str(pdb), None, 8.0, 0.3)
assert model.comm_backbone == bb, "backbone differs from repo build_model"
assert model.core_nodes == core and model.semi_nodes == semi and model.disorder_nodes == dis
pos = feats["pos"]
seeds = [1, 150, 600, 1100, 1715, 2400, 2512, 2663]
out = {"max_abs_diff_prob_iter": 0.0, "max_abs_diff_prob_solve": 0.0, "max_abs_diff_fp_solve": 0.0}
for a in (0.3, 0.15, 0.5):
    lu = SF.rwr_solve_all(W, a)
    for s in seeds:
        i = int(np.where(pos == s)[0][0])
        ref = M.rwr(G_all, s, alpha=a)
        refv = np.array([ref[n] for n in nodes])
        it = SF.rwr_iterative(W, i, a)
        e = np.zeros(len(nodes)); e[i] = 1.0
        sol = a * SF.splu(SF.sp.identity(len(nodes), format="csc") - (1 - a) * W.tocsc()).solve(e)
        sol /= sol.sum()
        out["max_abs_diff_prob_iter"] = max(out["max_abs_diff_prob_iter"], float(np.abs(it - refv).max()))
        out["max_abs_diff_prob_solve"] = max(out["max_abs_diff_prob_solve"], float(np.abs(sol - refv).max()))
        if a == 0.3:
            fp = M.fingerprint_for_variant(model, s)
            for k in ("CoreImpact", "SemiCoreImpact", "DisorderImpact", "CommBackboneImpact"):
                out["max_abs_diff_fp_solve"] = max(out["max_abs_diff_fp_solve"], abs(fp[k] - feats[k][i]))
        else:
            sfx = f"_a{str(a).replace('0.', '')}"
            core_ref = sum(ref.get(n, 0) for n in core)
            out["max_abs_diff_fp_solve"] = max(out["max_abs_diff_fp_solve"], abs(core_ref - feats["CoreImpact" + sfx][i]))
# pLDDT equals repo loader
p2, x2, pl2 = M.load_alphafold_ca(str(pdb))
out["plddt_identical_to_repo_loader"] = bool(np.array_equal(pl2, feats["plddt"]) and np.array_equal(p2, pos))
out["meta"] = {k: v for k, v in meta.items() if k != "af_sequence"}
print(json.dumps(out, indent=1, default=str))
ok = max(out["max_abs_diff_prob_iter"], out["max_abs_diff_prob_solve"], out["max_abs_diff_fp_solve"]) < 1e-9
print("VERIFIED" if ok and out["plddt_identical_to_repo_loader"] else "MISMATCH")
sys.exit(0 if ok else 1)
