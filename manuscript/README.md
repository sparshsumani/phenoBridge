# Manuscript

**Variant recurrence and publication clustering in literature-derived Mendelian disease cohorts: implications for genotype–phenotype inference in rare-disease research**
Sparsh Sumani, Ashok Kumar Bhatia

| File | Role |
|---|---|
| `body.html` | **Editable source**: text, tables and figure placement. Citations are written as `{{key}}`. |
| `style.css` | Page layout: LaTeX-like two-column article style (the same as `../manuscript_vus3d/`), with numbered sections and bracketed citations. |
| `build.py` | Builds `paper.html` from `body.html` + `style.css`. It also holds the reference list (DOIs checked against Crossref), numbers citations in order of first use, and prints word counts. |
| `paper.html` | Generated. Do not edit by hand. |
| `Sumani_Bhatia_variant_recurrence.pdf` | Compiled manuscript. |
| `REVISION_NOTES.md` | Summary of changes, claims that need support, remaining weaknesses and the validation plan. |
| `figs.py` | Generates Figures 1–2 from `../analysis_multigene/out/`. Supplementary Figures S1–S2 are copied from `../kbg_followup/fig/`. |

## Rebuild

1. Run the data pipeline: `multigene/run_all.sh`. It writes its outputs to `$PHENOBRIDGE_DATA/derived`.
2. Run the analyses from inside `analysis_multigene/`:
   ```
   python3 run_analysis.py --perms 2000
   python3 extra_analyses.py
   python3 design_effect.py
   python3 cohort_structure.py   # descriptive publication summaries added in revision
   ```
3. Build the figures and the manuscript, then render:
   ```
   python3 manuscript/figs.py
   python3 manuscript/build.py
   node manuscript/render.js $PWD/manuscript/paper.html $PWD/manuscript/Sumani_Bhatia_variant_recurrence.pdf
   ```

The earlier 8-page version (Sumani and Fu, `Sumani_Fu_variant_recurrence.pdf`) remains in the git history. The 4-page KBG-only draft (`../paper/`) is also superseded; its errors are listed in `../kbg_followup/RESULTS.md`, §7.
