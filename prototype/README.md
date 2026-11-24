# JS prototype

The dependency-free version the platform grew out of — own SMILES parser, descriptors,
fingerprints, ridge ensemble and active learning in ~900 lines of JavaScript, no installs.

```bash
python3 -m http.server 8777   # then open http://localhost:8777
```

```bash
node run.mjs 8 14
```

Useful for reading the loop end to end in one sitting. The real platform is in `../dmta/`:
RDKit chemistry, reaction SMARTS, calibrated uncertainty, persistence and an API.
