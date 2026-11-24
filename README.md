# DMTA platform

A custom Design–Make–Test platform for medicinal chemistry: enumerate a library that is
synthesisable by construction, score it against competing project objectives, plan batches that
share intermediates, and let a calibrated, uncertainty-aware model choose what to make next.

Real RDKit chemistry, **real ChEMBL bioactivity data**, a REST API, a SQLite campaign store, a
browser UI, a CLI, and benchmarks that check the loop actually beats screening.

- [RESULTS.md](RESULTS.md) — what it does on 8,978 real EGFR compounds, including the validation
  that fails and why that matters.
- [LANDSCAPE.md](LANDSCAPE.md) — how this compares to Proton, Iktos, Schrödinger and the
  open-source stack, and why it is built this way.

## Run it

```bash
pip install -r requirements.txt
```

```bash
python3 -m uvicorn dmta.api:app --port 8000
```

Then open http://localhost:8000 — the UI is served by the same process.

Pull a real target's measured activity data and check whether the model can predict it:

```bash
python3 -m dmta.cli chembl-fetch --target CHEMBL203
```

```bash
python3 -m dmta.cli chembl-validate --target CHEMBL203
```

Run the loop retrospectively over that real data — a measurement is revealed only when the loop
decides to make the compound:

```bash
python3 chembl_benchmark.py CHEMBL203 10 24
```

Prepare a protein structure for docking and see exactly what was stripped from it:

```bash
python3 -m dmta.cli structure --pdb 2RGP
```

Check whether docking into that structure actually ranks compounds whose potency you
already know — run this **before** trusting a docking score for anything:

```bash
python3 dock_validate.py 2RGP 250 CHEMBL203
```

Headless, which is how a computational chemist would really drive it:

```bash
python3 -m dmta.cli simulate --campaign kinaseX --cycles 8 --batch 16
```

```bash
python3 -m dmta.cli propose --campaign kinaseX --size 16 --out worklist.csv
```

```bash
python3 -m dmta.cli submit --campaign kinaseX --results worklist.csv
```

Compare acquisition strategies on identical landscapes:

```bash
python3 benchmark.py 8 16 42 7 101
```

```bash
python3 -m pytest tests/ -q
```

## Architecture

```
                 ┌──────────────┐
   catalogue ───▶│   library    │  reaction SMARTS × building blocks
   reagents      │ enumeration  │  → every product carries its route
                 └──────┬───────┘
                        ▼
  project profile ─▶ ┌────────┐      ┌──────────┐      ┌─────────────┐
  (data, not code)   │ scoring│─────▶│ selection│─────▶│  worklist   │──▶ lab / CRO
                     │  MPO   │      │acquisition      │   (CSV)     │
                     └────▲───┘      │+ constraints    └─────────────┘
                          │          └──────────┘             │
                     ┌────┴─────┐                             │ results
                     │  model   │◀────────────────────────────┘
                     │ RF + σ   │     SQLite: campaigns, batches,
                     └──────────┘     measurements, model runs
```

| Module | Responsibility |
|---|---|
| [`dmta/chem.py`](dmta/chem.py) | RDKit descriptors, Morgan fingerprints, QED, Ertl SA score, structural alerts, feature vectors |
| [`dmta/library.py`](dmta/library.py) | Reaction SMARTS, reagent catalogue, two-step enumeration, route construction, intermediate grouping |
| [`dmta/scoring.py`](dmta/scoring.py) | Desirability curves and the weighted multi-parameter profile |
| [`dmta/models.py`](dmta/models.py) | RF ensemble, per-tree σ, calibration against out-of-fold residuals, CV metrics |
| [`dmta/select.py`](dmta/select.py) | UCB / EI / Thompson / greedy / uncertainty / random, diversity penalty, budget and scaffold caps |
| [`dmta/assay.py`](dmta/assay.py) | Assay providers: simulated landscape, retrospective replay of real data, CSV from a real lab, worklist writer |
| [`dmta/chembl.py`](dmta/chembl.py) | ChEMBL fetch, cleaning (salts, replicates, censored values), Murcko scaffolds |
| [`dmta/validate.py`](dmta/validate.py) | Random / scaffold / time splits, enrichment, does σ rank error, docking agreement |
| [`dmta/structure.py`](dmta/structure.py) | PDB fetch, receptor cleaning, binding box taken from the crystal ligand, manifest |
| [`dmta/docking.py`](dmta/docking.py) | 3D embedding, smina docking with caching, loud failure when the engine is missing |
| [`dmta/store.py`](dmta/store.py) | SQLite schema and queries |
| [`dmta/campaign.py`](dmta/campaign.py) | The orchestrator — one cycle, replay from the database |
| [`dmta/api.py`](dmta/api.py) | FastAPI service, serves the UI |
| [`dmta/cli.py`](dmta/cli.py) | `list` / `propose` / `submit` / `simulate` / `export` |
| [`ui/`](ui/) | Browser client — no chemistry in the browser, it only calls `/api` |
| [`prototype/`](prototype/) | The dependency-free JS version this grew out of; useful for understanding the loop in one file |

### Four decisions worth knowing about

**Synthesisability is structural, not a filter.** A product exists only because
`[C:1](=[O:2])[OX2H1].[N;H2,H1:3]>>[C:1](=[O:2])[N:3]` and friends actually fired on catalogue
reagents. 924 compounds from 5 scaffolds × 14 amines × 14 boronic acids, deduplicated by InChIKey
and screened for structural alerts. No compound can be proposed without a route, because a
compound without a route was never created.

**Uncertainty is calibrated.** Raw ensemble spread is not a confidence interval. The model
cross-validates, compares predicted σ against out-of-fold residuals, and scales. Without this
the acquisition function chases noise.

**Optimism applies to potency, then re-enters the profile.** UCB does not add σ to the final
score. It computes `μ + βσ` as an optimistic *potency*, then runs the full multi-parameter
profile again with that value. A compound the model is unsure about gets made only if an
optimistic potency would still clear MW, cLogP, route confidence and cost. Adding σ to the final
score instead — which is the obvious implementation — made UCB **lose to random selection** in
the first benchmark run.

**The lab boundary is a CSV.** `propose` writes a worklist with an empty result column;
`submit` reads it back. The simulated assay is a benchmark harness, and the UI says so.

## Does it work on real data?

Validated on four real ChEMBL datasets spanning three protein classes — EGFR and VEGFR2
(kinases), carbonic anhydrase II (zinc metalloenzyme) and the μ-opioid receptor (GPCR) —
making 240 compounds from each. Actives found, with the multiple of each library's base rate:

| strategy | EGFR | VEGFR2 | CA II | μ-opioid |
|---|---|---|---|---|
| **greedy** | **142.5 (2.37×)** | **72 (1.74×)** | **154.5 (2.18×)** | 134.5 (1.97×) |
| expected improvement | 121.5 (2.02×) | 67 (1.62×) | 136.5 (1.92×) | **136 (1.99×)** |
| random | 61 (1.02×) | 36 (0.87×) | 69.5 (0.98×) | 67.5 (0.99×) |

Identical ordering on all four, at roughly twice the base rate, and the chosen compounds are
more drug-like besides (mean QED 0.76–0.81 against random's 0.52–0.57). Random recovers each
library's base rate within noise, which is the check that the harness itself is unbiased.

And the model's generalisation, which is the number that actually decides whether to trust it:

Spearman ρ, four targets across three protein classes:

| split | EGFR | VEGFR2 | CA II | μ-opioid |
|---|---|---|---|---|
| random | 0.87 | 0.84 | 0.78 | 0.85 |
| scaffold | 0.81 | 0.78 | 0.56 | 0.73 |
| **time** | **0.11** | **0.44** | **0.42** | **0.36** |

**The pipeline generalises beyond kinases** — a GPCR and a zinc metalloenzyme come through
the same code with scaffold ρ of 0.73 and 0.56 — **and EGFR turns out to be the outlier.**
It collapses on the time split; the other three retain signal (ρ 0.36–0.44, 1.9–3.1×
enrichment). Both have ~90% of test compounds on unseen
scaffolds, so structural novelty is not the cause — EGFR's post-2021 chemistry simply got
0.84 log units more potent while VEGFR2's did not. It is label drift, not scaffold novelty.
(It is also not covalent inhibitors; that hypothesis was tested and disproved.)

A second measured finding changed the product: **σ ranks error at ρ≈0.44 on a random split,
≈0.36 on scaffolds and ≈0.12 on future compounds — inverting on μ-opioid.** The uncertainty
every acquisition function spends batch on is least trustworthy in production, so
`chembl-validate` now grades σ per split and recommends a starting β from the worst one.

The rule that survives: **trust the model to rank within a series, never to compare across
eras or compound classes** — which is exactly the regime a closed loop operates in, retrained
every cycle and asked only to order the next batch. Full diagnosis in [RESULTS.md](RESULTS.md).

## Does it work on the simulated landscape?

8 cycles × 16 compounds (128 of 924 made), averaged over three independent landscapes:

| strategy | best found | top 5% found | prospective ρ | reactions |
|---|---|---|---|---|
| **UCB** | 8.22 | **18.3 / 46** | 0.49 | 208 |
| greedy | 8.33 | 19.3 / 46 | 0.53 | 208 |
| expected improvement | 8.22 | 19.0 / 46 | 0.52 | 207 |
| uncertainty (pure explore) | 8.43 | 7.7 / 46 | 0.64 | 212 |
| random screening | 8.23 | 6.7 / 46 | 0.69 | 227 |

**2.75× enrichment** over random screening for UCB, and shared intermediates cut the batch from
256 reactions to 208 (19% fewer). A single longer run reached 3.96× with prospective Spearman ρ
climbing to 0.78 by cycle 4.

Two honest readings of that table:

- **Greedy matches UCB here.** On a 924-compound library with a decent model, exploitation is
  hard to beat; exploration earns its keep on larger spaces and rougher landscapes. The platform
  does not hide this — it ships the benchmark that shows it.
- **Prospective R² is low (0.1–0.4) while ρ is 0.5–0.8.** That is expected and it is the right
  metric: each batch is *selected*, not random, so the model is always extrapolating into a
  narrow slice. Ranking is what the loop consumes, not absolute accuracy.

## What is real and what is simulated

Two modes, and the UI and API both say which one a campaign is in:

- **Real data** (`--chembl CHEMBL203`, or a `CHEMBL…` target over the API): measured IC50 values
  from the literature, replayed. A compound's value is hidden until the loop selects it. There
  are no routes or costs — these compounds came from papers, not from our building blocks — so
  the ChEMBL profile leans on the Ertl SA score instead.
- **Enumerated library**: the make-on-demand space with full routes, yields, costs and lead
  times, against a simulated assay.

Real in both: the chemistry (RDKit), the reactions (SMARTS that genuinely fire), the descriptors,
the models, the calibration, the acquisition functions, the persistence, the API, the CSV handoff.

Simulated in the second mode only: **the assay.** `SimulatedAssay` holds a hidden landscape — a random weighting over
fingerprint bits, a lipophilicity term, a size penalty, and bonuses for two substructures the
"target" likes. It exists to benchmark selection strategies, and nothing here predicts real
potency. Also simplified: reagent prices and lead times are plausible, not quoted; yields are
per-reaction medians, not condition-specific.

## Taking it to production

1. **Point it at your data.** `CsvAssay` reads results from a lab; `chembl.py` shows the full
   fetch-and-clean path for public data. Run `chembl-validate` on your own target first — the
   scaffold-split ρ is what tells you whether a loop is worth starting.
2. **Swap the model.** `PropertyModel` is one class with `fit` / `predict` / `apply`. Chemprop
   (D-MPNN with mean-variance estimation) drops in here; keep the calibration step.
3. **Widen the library.** Add reactions and reagents to `library.py`, or enumerate against an
   Enamine REAL subset. Past ~10⁵ compounds, replace the in-memory list with a vector index and
   score a candidate pool rather than the full space.
4. **Add the arms this does not have.** Docking or FEP as another objective; a retrosynthesis
   engine (AiZynthFinder, ASKCOS) if you need routes to molecules you did not enumerate; a
   generative model (REINVENT 4) if enumeration stops covering the space you care about.
5. **Operational fit matters more than any of the above.** Set batch size by what your chemists
   can review in a week, not by what the algorithm would like.
