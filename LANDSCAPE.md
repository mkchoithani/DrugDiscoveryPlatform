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

## The real competitors, module by module

Updated 2026-11-26 from the vendors' own product pages. These three are the companies we
would actually meet in a procurement process — unlike the funded cohorts in
[COMPETITORS.md](COMPETITORS.md), none of which sells software. Vendor claims are the
vendors' own and are quoted, not endorsed.

### Iktos — owns MAKE

| Module | What it does |
|---|---|
| **Makya** | Generative design, sold as a SaaS platform. Modes: growing, linking, cyclization, macrocyclization, fine tuning. Scoring via 3D ligand- and structure-based modelling, SAR modelling, ADMET. External models connect over an API |
| **Spaya** / **Spaya API** | Retrosynthesis to commercially available starting materials. Reinforcement learning plus Monte-Carlo tree search over a template network "trained on 6M reactions", "+2000 name reactions", +35 supplier partnerships |
| ↳ **RScore** · **RSpred Score** | Proprietary synthetic-accessibility score, and a fast estimate of it *without running the retrosynthesis* |
| **Ilaka** | Orchestration AI — raw-material ordering through to directing robots, with a cluster synthesis strategy to prioritise routes |
| **Iktos Robotics** | Chemspeed synthesis platform; **MT bench®** cell-based phenotypic assay; purification and analysis |

Claims as stated: *"six times more parallel projects with the same number of chemists"*,
*"Reduce the timeline to under 24 months"*, *"100 reactions per day in parallel across more
than 55 different named reactions"*. Sold as collaborations; Spaya deploys **on AWS or in
the customer's own VPC**, with Knime, Pipeline Pilot and Jupyter integrations, and can be
retrained on a customer's reaction data "in a couple of weeks".

### Schrödinger — owns PREDICT

| Module | What it does |
|---|---|
| **Maestro** | The modelling workbench everything else runs inside |
| **Glide** | Docking |
| **FEP+** | Free-energy perturbation — the accurate, expensive binding-affinity calculation |
| **LiveDesign** | Collaborative design and data platform; a biologics edition exists |
| **BioLuminate** | Biologics modelling |
| **Bunsen** | An agentic AI "co-scientist", early access from August 2026 |
| Solution areas | Structure prediction and target enablement, hit discovery, hit-to-lead and lead optimisation, drug formulation, antibody design, peptide discovery, enzyme engineering, macrocycle modelling |
| **Modeling Services** | They run the calculations for you |

Claims as stated: *"Built upon more than 30 years of R&D"*; a case study headed
*"Hit to development candidate in 10 months"*. No public pricing; Q3 2026 ACV was reported
at $45M, up 17%.

### XtalPi — owns BOTH, plus the lab

| Module | What it does |
|---|---|
| **ID4** | "Inclusive Digital Drug Discovery and Development" — ideation, drug-like property evaluation and optimisation, ADMET, synthesis, biological functional studies |
| **ID4Inno** | The small-molecule offering: AI models paired with automated experiments |
| The engine | Reported 300+ AI models plus a physics-based computational chemistry platform, and an AI-and-robotics R&D infrastructure |

Reported 100+ discovery projects supported; H1 2026 revenue RMB 393.6M, down from
RMB 517.1M a year earlier; a collaboration with potential value over US$400M.

## Where we actually stand against them

Honestly, and feature by feature, **we lose nearly everywhere.**

| Capability | Iktos | Schrödinger | XtalPi | Us |
|---|---|---|---|---|
| Generative design | Makya | — | ID4 | **parked, incomplete** |
| Retrosynthesis | Spaya, 6M reactions | — | ID4 | **none** — we enumerate instead |
| Synthetic accessibility | RScore / RSpred | — | ID4 | RDKit SA score |
| Docking | — | Glide | ID4 | built, **measured, rejected** |
| Free-energy affinity | — | **FEP+** | physics platform | **none** |
| ADMET | Makya scoring | property prediction | ID4 | **logD only; solubility failed** |
| Multi-parameter scoring | Makya | LiveDesign | ID4 | **yes, with class-specific profiles** |
| Active learning / acquisition | implicit | implicit | implicit | **yes, first-class and swappable** |
| Calibrated uncertainty | not stated | not stated | not stated | **yes, scaled on out-of-fold residuals** |
| Robotic synthesis | **Chemspeed, 100 rxn/day** | — | **automated lab** | **none — a CSV** |
| Phenotypic assay | **MT bench®** | — | yes | none |
| Biologics / peptides | — | BioLuminate, peptide, antibody | yes | **none** |
| Lab-software integration | Knime, Pipeline Pilot, Jupyter | LiveDesign | yes | **none** |
| Customer-premises deployment | **Spaya on your VPC** | not stated | not stated | **yes, fully local** |
| Agentic assistant | Ilaka (orchestration) | **Bunsen** | — | none |
| **Published negative results** | **none found** | **none found** | **none found** | **three rejections, three ties** |
| **Paired confidence intervals** | **none found** | **none found** | **none found** | **on every comparison** |

### Three things that table tells us

**1. Our feature gap is real and large.** Anyone comparing line by line will pick one of
them. We should stop pretending otherwise internally — the Tier 3 recommendation in
[FEATURES.md](FEATURES.md) against building modalities is sound precisely *because* we
cannot win this table by filling it in.

**2. The bottom two rows are the entire business.** Every claim on those pages is an
aggregate or a case study: six times more projects, under 24 months, hit to candidate in
ten months. **Not one of them states a control arm, a split, or an interval.** That is the
same finding as the funded cohorts, now confirmed against the three companies who actually
sell software. It is a consistent, structural gap in how this industry reports, and it is
the one row where we are not merely competitive but alone.

**3. A correction to our own positioning.** We wrote in FEATURES.md that keeping customer
data on the customer's own machine was a property nobody else claimed. **Iktos already
offers Spaya on a customer VPC.** Our claim has to narrow to what is actually true: we run
*entirely* locally with no model provider in the loop at all, which a SaaS generative
platform deployed into a VPC does not match — but "nobody else offers on-premises" was
wrong, and was exactly the kind of unchecked claim we audit other people for.

### What this does not change

Iktos owns making molecules. Schrödinger owns predicting them accurately. XtalPi owns
both and the building. We own the loop's bookkeeping and the examiner — and the examiner
is the only thing on that table that none of the three has.

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
