# Manuscript

**Variant recurrence inflates genotype–phenotype tests in curated patient cohorts: a calibrated evaluation of AlphaFold structural features in 44 Mendelian genes**
Sparsh Sumani, Yuzhang Fu

- `Sumani_Fu_variant_recurrence.pdf` is the built paper.
- `paper.html` is the source.
- `figs.py` generates Figures 1–2 from `../analysis_multigene/out/`.
  - Figures 3–4 are copied from `../kbg_followup/fig/`.

## Rebuild

1. Run the data pipeline: `multigene/run_all.sh`. It writes its outputs to `$PHENOBRIDGE_DATA/derived`.
2. Run the analyses, from inside `analysis_multigene/`:
   ```
   python3 run_analysis.py --perms 2000
   python3 extra_analyses.py
   python3 design_effect.py
   ```
3. Build the figures: `python3 manuscript/figs.py`.
4. Render the PDF (requires Playwright Chromium):
   ```
   node manuscript/render.js $PWD/manuscript/paper.html $PWD/manuscript/Sumani_Fu_variant_recurrence.pdf
   ```

The earlier 4-page KBG-only draft (`../paper/`) is superseded. Its errors are listed in `../kbg_followup/RESULTS.md`, §7.
