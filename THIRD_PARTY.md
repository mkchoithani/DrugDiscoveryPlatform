# Third-party components and data

## Bundled code

| Component | Version | Author | Licence |
|---|---|---|---|
| [3Dmol.js](https://3dmol.csb.pitt.edu/) (`ui/vendor/3Dmol-min.js`) | 3.0.0 | David Koes and contributors | BSD-3-Clause |

The upstream attribution header is preserved in the minified file.

## Runtime dependencies

RDKit (BSD-3-Clause), scikit-learn (BSD-3-Clause), NumPy (BSD-3-Clause),
FastAPI (MIT), uvicorn (BSD-3-Clause), pytest (MIT). Optional, for docking:
[smina](https://sourceforge.net/projects/smina/) (GPL-2.0) — invoked as an external
binary, never linked into this code.

## Data (not redistributed here)

This repository contains **no** third-party data. Everything is fetched on demand and
is excluded by `.gitignore`:

- **ChEMBL** bioactivity data — [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/),
  fetched via the public EBI API.
- **RCSB PDB** structures (e.g. 2RGP) — public domain, fetched from files.rcsb.org.

Regenerate with:

```bash
python3 -m dmta.cli chembl-fetch --target CHEMBL203
```

```bash
python3 -m dmta.cli structure --pdb 2RGP
```
