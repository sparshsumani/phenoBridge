# Manuscript: VUS-3D temporal benchmark

`Sumani_Fu_VUS3D_benchmark.pdf` reports the pre-registered VUS-3D study (`../vus3d/`).
All numbers come from `../vus3d/out/`; see `../vus3d/RESULTS.md` for the full tables.

Rebuild:
```bash
python3 figs.py        # fig/*.svg from ../vus3d/out/
python3 build.py       # body.html + _head.html -> paper.html (citations numbered by first use)
node render.js "$PWD/paper.html" "$PWD/Sumani_Fu_VUS3D_benchmark.pdf"
```
Affiliations, correspondence, author contributions and competing interests are placeholders.
