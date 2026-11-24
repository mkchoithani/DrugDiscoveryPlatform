# Who else builds this, and what a custom build should copy

Research notes behind the architecture. Vendor claims are the vendors' own and are not
independently verified; where a page could not be read, that is stated.

## The commercial platforms

| Platform | Shape | Named components | What they claim |
|---|---|---|---|
| **PostEra Proton** | End-to-end DMT platform, delivered through partnerships | Not published as separate modules; "chemistry foundation models", synthesis search, active learning | Closes the Design-Make-Test cycle; no quantified results on the public page |
| **Iktos** | The most explicit DMTA loop, including its own robotic lab | Makya (generative design), Spaya (retrosynthesis), Ilaka (orchestration), Iktos Robotics | "six times more parallel projects with the same number of chemists"; candidates in under 24 months; a Chemspeed platform running 100 reactions/day across 55+ named reactions |
| **Schrödinger** | Physics-first, not a closed loop | Maestro, Glide, FEP+, LiveDesign | "Hit to development candidate in 10 months" (MALT1); 30+ years of R&D |
| **Insilico Chemistry42** | Generative chemistry engine | Not readable — the page renders client-side | Secondary sources describe multi-objective generative design; clinical assets are the differentiator |

The pattern worth noting: **everyone converged on the same loop**, and the differentiation is
in which arm they own. Schrödinger owns *predict accurately*. Iktos owns *make it physically*.
PostEra positions on *synthesis-aware design*. None of them differentiate on the loop itself —
which is exactly why a custom build is reasonable.

## What you would otherwise assemble from open source

| Need | Tool | Licence |
|---|---|---|
| Cheminformatics | **RDKit** | BSD |
| Property prediction | **Chemprop** (D-MPNN, with ensembling / mean-variance / evidential uncertainty) | MIT |
| Generative design | **REINVENT 4** (AstraZeneca; de novo, R-group, scaffold hopping, linker design) | Apache 2.0 |
| Retrosynthesis | **AiZynthFinder 4** (AstraZeneca, MCTS to purchasable precursors) or **ASKCOS** (MIT) | MIT |
| Make-on-demand space | **Enamine REAL** — 78B+ molecules in the navigator, reported >80% synthesis success in 3–4 weeks | commercial catalogue |

A custom platform is mostly *integration*, not invention. The parts that are genuinely yours
are the ones no vendor can supply: your assay data, your project's property profile, your
chemists' constraints, and the loop's bookkeeping.

## What the literature says about the loop itself

- Batch active learning beats chronological (human) ordering on real med-chem series, with
  reported reductions of ~58–62% in the training data needed on some targets.
- Selection strategy matters more than model architecture in several comparisons — which is
  why this build makes the acquisition function a first-class, swappable thing.
- Practitioner guidance worth heeding: **batch size should match reviewer capacity, not
  algorithmic efficiency.** If the system proposes more compounds than a chemist can review
  in a week, the batch is too large.

## Design decisions this drove

1. **Synthesisability by construction, not as a filter.** Compounds exist because a reaction
   SMARTS joined catalogue reagents, so a route is never missing — the Iktos/PostEra insight,
   achievable without a retrosynthesis engine at this scale.
2. **Uncertainty must be calibrated.** Raw ensemble spread is not a confidence interval; it is
   scaled against out-of-fold residuals before any acquisition uses it.
3. **Optimism applies to potency only.** Exploration should not buy a compound that is
   unmakeable or un-drug-like — so the optimistic potency is pushed back through the same
   multi-parameter profile rather than added to the final score.
4. **The lab is a CSV, not an API.** Every real engagement starts with a worklist out and
   results back. The simulated assay is a benchmark harness, not the product.
5. **The profile is data.** A project changes its potency goal, cLogP window and cost ceiling
   constantly; that must not be a code change.

## Sources

- [PostEra — about](https://postera.ai/about/)
- [Iktos — discovery platform](https://iktos.ai/discovery-platform) · [Iktos home](https://iktos.ai/)
- [Schrödinger — platform](https://www.schrodinger.com/platform/)
- [Insilico — Chemistry42](https://insilico.com/chemistry42) (page did not render server-side)
- [AiZynthFinder 4.0, *J. Cheminform.* 2024](https://link.springer.com/article/10.1186/s13321-024-00860-x)
- [ASKCOS: open source synthesis planning](https://arxiv.org/pdf/2501.01835)
- [REINVENT 4 — repository](https://github.com/MolecularAI/REINVENT4) · [paper](https://link.springer.com/article/10.1186/s13321-024-00812-5)
- [Chemprop — *JCIM* 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10777403/)
- [Enamine REAL Space Navigator](https://enamine.net/library-synthesis/real-compounds/real-space-navigator) · [BioSolveIT REAL Space](https://www.biosolveit.de/chemical-spaces/real-space/)
- [Deep batch active learning for drug discovery (preprint)](https://www.biorxiv.org/content/10.1101/2023.07.26.550653v3.full.pdf)
- [Query matters: selection strategies in active learning, *JCIM* 2026](https://pubmed.ncbi.nlm.nih.gov/41748299/)
