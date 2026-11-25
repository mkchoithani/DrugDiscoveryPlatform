# Results on real data

Everything here is measured, reproducible from this repository, and reported whether or not
it flatters the platform.

Four targets, deliberately spanning protein classes — two kinases, a zinc metalloenzyme and
a GPCR — cleaned by the same pipeline with no per-target tuning:

| target | class | measure | compounds | scaffolds | actives ≥8 | years |
|---|---|---|---|---|---|---|
| EGFR (CHEMBL203) | kinase | IC50 | 8,978 | 3,294 | 25% | 1989–2025 |
| VEGFR2 (CHEMBL279) | kinase | IC50 | 8,900 | 3,371 | 17% | 1999–2025 |
| Carbonic anhydrase II (CHEMBL205) | metalloenzyme | Ki | 7,228 | 1,982 | 30% | 1987–2025 |
| μ-opioid receptor (CHEMBL233) | GPCR | Ki | 5,193 | 2,227 | 28% | 1989–2025 |

```bash
python3 -m dmta.cli chembl-fetch --target CHEMBL205 --type Ki
```

Dropped during cleaning, each for a stated reason: 1,255 non-binding assays, 1,208 compounds
whose replicate measurements disagreed by more than 1 log unit, 508 outside MW 150–750, 18
unparsable. Salts stripped to the largest fragment, charges neutralised, replicates reduced to
the median.

## 1. Can the model predict real potency?

```bash
python3 -m dmta.cli chembl-validate --target CHEMBL205 --type Ki
```

Spearman ρ between prediction and measurement, by split and target:

| split | EGFR | VEGFR2 | CA II | μ-opioid | what it means |
|---|---|---|---|---|---|
| random | 0.87 | 0.84 | 0.78 | 0.85 | the number vendors quote — optimistic |
| **scaffold** | **0.81** | **0.78** | **0.56** | **0.73** | unseen chemotypes — the honest test |
| **time** | **0.11** | **0.44** | **0.42** | **0.36** | future compounds — what production faces |

Top-decile enrichment on the time split: EGFR 0.93×, VEGFR2 2.03×, CA II 1.87×, μ-opioid 3.13×.

**The pipeline generalises beyond kinases.** A GPCR and a zinc metalloenzyme — chemistry with
nothing in common with kinase hinge binders — come through the same cleaning, featurisation and
modelling with scaffold-split ρ of 0.73 and 0.56. Carbonic anhydrase is the weakest, which is
chemically sensible: CA II inhibition is dominated by the zinc-binding sulfonamide warhead, so
holding out scaffolds removes less signal than elsewhere and the remaining variance is subtler.

### EGFR is the outlier, not the rule

With EGFR alone the conclusion looked like *a fingerprint model cannot rank chemotypes it has
never seen*. Three more targets say otherwise — all three retain real time-split signal
(ρ 0.36–0.44, enrichment 1.9–3.1×). Only EGFR collapses.

| | EGFR | VEGFR2 |
|---|---|---|
| test compounds on an unseen scaffold | 91% | 89% |
| potency drift (test mean − train mean) | **+0.84 log** | **+0.05 log** |
| time-split ρ | 0.11 | 0.44 |

Structural novelty is identical. What differs is that EGFR's chemistry got dramatically more
potent after 2021 and VEGFR2's did not: **label drift, not scaffold novelty.**

Two further probes, both of which killed a plausible story:

- **Correcting the +0.93 log bias repairs the error but not the ranking** (MAE 1.29 → 0.97,
  R² −0.66 → −0.04, ρ unchanged at 0.11). Not a simple offset.
- **It is not covalent inhibitors.** The obvious suspect was the post-2021 move to covalent,
  mutant-selective EGFR drugs. But acrylamide warheads are 28% of train and 26% of test — no
  shift — and the covalent subset *ranks better* (ρ 0.35) than the non-covalent one (ρ 0.26).
  Both beat the pooled 0.11: Simpson's paradox. The model ranks acceptably *within* a compound
  class and mis-ranks the classes *against each other*, because their potency distributions
  moved apart.

### The finding that changed the product

Uncertainty is what every acquisition function except `greedy` spends part of its batch on. So
it matters how well σ ranks error — and that degrades exactly where you need it:

| target | random | scaffold | time |
|---|---|---|---|
| EGFR | 0.44 | 0.46 | 0.32 |
| VEGFR2 | 0.48 | 0.35 | 0.12 |
| CA II | 0.38 | 0.27 | 0.15 |
| μ-opioid | 0.44 | 0.36 | **−0.13** |
| **mean** | **0.44** | **0.36** | **0.12** |

σ ranks error respectably when the test set resembles the training set, weakly on unseen
scaffolds, and barely at all on future compounds — **inverting entirely on μ-opioid**, where the
model is most confident precisely where it is most wrong. A production model lives in that third
column.

`chembl-validate` now reports this verdict per split and recommends a starting β from the
*worst* of the scaffold and time splits, because recommending on the best one would tell you to
explore on an uncertainty estimate that is worse than useless:

```
  σ on the time split: unreliable — σ is uninformative or inverted — do not explore on it
  recommended exploration β for this target: 0.0 (…use greedy or EI)
```

One caution I want on the record: this sits alongside the β sweep below, and the two are **not
known to be causally linked**. Low β may win because σ is untrustworthy, or simply because a
2.7% budget leaves no time to recoup information value. The recommender is capped at 0.25 — what
the sweep actually supported — and tells you to sweep anyway.

## 2. Does the loop find potent compounds efficiently?

Retrospective active learning: the library is every real compound for the target, and a
measured value is revealed **only when the loop decides to make that compound**.

```bash
python3 chembl_benchmark.py CHEMBL205 10 24 Ki
```

10 cycles × 24 compounds = 240 made, averaged over 2 seeds. Actives found (≥8 of 240 made),
with the multiple of the library's own base rate in brackets:

| strategy | EGFR | VEGFR2 | CA II | μ-opioid |
|---|---|---|---|---|
| **greedy** | **142.5 (2.37×)** | **72 (1.74×)** | **154.5 (2.18×)** | 134.5 (1.97×) |
| expected improvement | 121.5 (2.02×) | 67 (1.62×) | 136.5 (1.92×) | **136 (1.99×)** |
| UCB (β=1) | 77 (1.28×) | 54 (1.31×) | 87.5 (1.23×) | 76 (1.11×) |
| random | 61 (1.02×) | 36 (0.87×) | 69.5 (0.98×) | 67.5 (0.99×) |
| *library base rate* | *25.0%* | *17.2%* | *29.6%* | *28.4%* |

**The ordering is identical on all four targets**: greedy ≈ EI > UCB > random, at roughly
twice the base rate. Two kinases, a metalloenzyme and a GPCR agree, which is the strongest
evidence here that the loop is doing something real rather than fitting one dataset's quirks.

**The harness checks out too.** Random selection recovers each library's base rate within
noise on every target — 25.4% vs 25.0%, 15.0% vs 17.2%, 29.0% vs 29.6%, 28.1% vs 28.4%. A
retrospective benchmark that quietly favoured the loop would show up here first.

The compounds it picks are also more drug-like, on every target: mean QED 0.76–0.81 for
greedy against 0.52–0.57 for random, at mean MW 346 against 409. That is the multi-parameter
profile working, not a side effect of potency.

μ-opioid is the one target where expected improvement edges greedy (136 vs 134.5) — well
inside seed noise, and worth noting only because μ-opioid is also the target whose σ inverted
under time shift. I would not read anything into it without more seeds.

## 3. The exploration parameter was wrong

UCB losing to greedy prompted a sweep of β on the same data (1 seed, 10 × 24):

| β | actives ≥ 8 | top 5% found |
|---|---|---|
| **0.00** | **117** | **36** |
| 0.25 | 105 | 25 |
| 0.50 | 94 | 29 |
| 1.00 | 75 | 22 |

Exploration hurts monotonically at this budget. With 240 compounds against 8,978, there is no
time to pay for information that only helps later — and the diversity penalty is already doing
the job of stopping redundant picks.

**The default β was changed from 1.0 to 0.25** as a result. It is a per-campaign parameter and
this sweep is how you tune it; a platform that ships a fixed β without showing this curve is
guessing on your behalf.

## What these numbers do not say

- **Four targets, three protein classes, all retrospective.** They agree on strategy ordering,
  which is real evidence. But every compound in these libraries was made by someone who thought
  it worth making, so even the random baseline sits on an enriched deck — a true screening
  collection is harsher. Run `chembl-validate` on your own target; it costs one command.
- **Potency only.** The loop optimises a measured activity endpoint plus computed physchem.
  No selectivity, no ADMET, no tox — the objectives a real programme dies on.
- Retrospective ≠ prospective. Every compound in this library was made by someone who thought it
  was worth making, so even the random baseline is enriched compared with a true screening deck.
- The model is a random forest on Morgan fingerprints — a strong, cheap baseline, not the state
  of the art. A D-MPNN (Chemprop) would likely improve the scaffold split and would not rescue
  the time split.
- pChEMBL aggregates measurements across assay formats and laboratories. Replicates disagreeing
  by more than 1 log were dropped, but residual inter-lab noise remains and sets a floor on
  achievable MAE.

---

# 4. Docking into 2RGP: a negative result

Docking was the proposed fix for the time-split failure above — the one scorer that
judges physical fit rather than resemblance to past compounds. It was validated before
being wired in. **It failed, and it is not in the loop.**

```bash
python3 dock_validate.py 2RGP 250 CHEMBL203
```

Setup: EGFR structure `2RGP` (2.0 Å), waters and phosphates stripped, binding box taken
from the co-crystallised ligand HYZ, smina at exhaustiveness 8. 249 of 250 compounds
docked successfully in 35 minutes on CPU (8.5 s each); 1 failed and is reported, not
silently dropped.

| subset | n | ρ | top-10% hit rate | base rate | enrichment | separation |
|---|---|---|---|---|---|---|
| all | 249 | **0.038** | 8% | 22% | **0.37×** | −0.25 |
| reversible | 182 | 0.068 | **0%** | 21% | **0.00×** | −0.33 |
| covalent | 67 | −0.065 | 33% | 27% | 1.24× | −0.08 |

For comparison, the fingerprint model on the same 250 compounds: **ρ 0.815, 3.51×
enrichment**. (That is a random-split figure and therefore flattering — but a model at
0.815 versus a docking score at 0.038 is not a close call.)

**Enrichment below 1.0 means worse than picking at random.** The 10% of compounds docking
liked best contained *fewer* actives than the average. On the reversible subset, the top
10% contained **no actives at all**. The separation column is negative throughout: active
compounds scored very slightly *worse* than inactive ones.

## Why this happened — and why it is not surprising

The sample is 250 EGFR compounds. They are not actives versus junk; they are **actives
versus slightly less active EGFR compounds**, mostly from the same chemical series.
Ranking within a tight series is the hardest thing you can ask a docking score to do, and
the published consensus is that docking is poor at it. Docking's established strength is
separating binders from genuine non-binders.

And it clearly retains that strength here. The smoke test, before any of this:

```
erlotinib  -8.01     caffeine  -5.63
gefitinib  -8.96     aspirin   -6.49
lapatinib -11.03     ethanol   -2.85
```

So the honest reading is **not** "docking is useless". It is:

> Docking cannot rank potency inside a congeneric series, which is what the lead-optimisation
> loop needs. It may still work as a coarse filter over a large, diverse library, which is a
> different job — stage 4 of the roadmap, not stage 3.

## What could be wrong with this setup rather than with docking

Stated plainly, because a negative result deserves the same scrutiny as a positive one:

- **One rigid receptor.** No ensemble across EGFR's several conformations.
- **No structural waters kept.** Kinase hinge waters can matter; all 85 were removed.
- **Reversible physics on covalent compounds** — 67 of 249 here. Expected to fail, and did.
- **Modest exhaustiveness (8)** and a single pose per compound.
- **No protonation-state optimisation** beyond smina's internal typing.

Any of these could be improved, and a specialist might extract more. But a *reasonable
default* setup produced nothing, and the next increment of effort has an uncertain payoff.

## What this cost, and what it bought

35 minutes of laptop CPU. No GPU was rented; no money was spent.

What it bought is the thing that matters: **the loop was never contaminated.** An
unvalidated docking objective would have been weighted at 2.0 in the scoring profile —
second only to potency — and would have actively steered the loop away from good
compounds, with no error message and no symptom. The validation harness is the only
reason that is a paragraph in a document instead of a wasted quarter.

**Status: `DOCKING_OBJECTIVE` exists in `dmta/scoring.py` and is switched off. It stays
off for this target until a setup passes this test.**

---

# 4. Docking: tested, and rejected

The fix proposed for the time-split failure was 3D docking — judging a molecule by
physically fitting it to the protein rather than by resemblance to past compounds.

Structure: **2RGP**, the standard EGFR benchmark entry, 2.0 Å, chain A. Waters and
phosphates stripped (85 and 10 atoms, both logged). The search box was taken from
where the crystal ligand HYZ actually sits, so the pocket came from the experiment
rather than from anyone's judgement. Engine: smina, exhaustiveness 8.

```bash
python3 -m dmta.cli structure --pdb 2RGP
python3 dock_validate.py 2RGP 250 CHEMBL203
```

## It passed the smoke test

| molecule | docking score |
|---|---|
| lapatinib (EGFR drug) | −11.03 |
| gefitinib (EGFR drug) | −8.96 |
| erlotinib (EGFR drug) | −8.01 |
| aspirin | −6.49 |
| caffeine | −5.63 |
| ethanol | −2.85 |

Real drugs clearly separated from molecules that have no business binding. At this
point everything looked right.

## It failed the real test

249 random ChEMBL compounds docked (one failed), 57 of them active at pChEMBL ≥ 8:

| subset | n | ρ | top-10% hit rate | base | enrichment | separation |
|---|---|---|---|---|---|---|
| all | 249 | **0.038** | 8% | 22% | **0.37×** | −0.25 |
| reversible | 182 | 0.068 | 0% | 21% | 0.00× | −0.33 |
| covalent | 67 | −0.065 | 33% | 27% | 1.24× | −0.08 |

**The fingerprint model on the very same 250 compounds: ρ 0.815, 3.51× enrichment.**

Docking is not merely weaker here — at 0.37× it is **worse than picking at random**,
and the separation is negative, meaning actives scored slightly *worse* than inactives.
The reversible subset, which should be docking's best case, returned **zero** actives in
its top decile.

Twelve poses make it visible:

| compound | docking score | measured pIC50 |
|---|---|---|
| CHEMBL5207479 | **−12.06** (best) | 7.89 |
| CHEMBL1684800 | −11.79 | 6.99 |
| CHEMBL121405 | −11.56 | **5.83** (worst) |
| CHEMBL3360610 | −9.16 | 9.24 |
| CHEMBL91980 | −8.48 | 9.46 |
| CHEMBL5881155 | **−7.76** (worst) | **9.66** (best) |

Within this set the relationship is close to inverted: the compound docking liked least
is the most potent one in the sample.

## Why this happens

Not a bug in the pipeline — a known limit of this class of scoring function.

1. **Separating binders from non-binders is a different task from ranking binders.**
   The smoke test did the first; real work needs the second. Vina-family functions are
   respectable at the first and poor at the second.
2. **One rigid structure.** The protein cannot adjust to each ligand, and 2RGP's
   conformation suits the chemotype it was crystallised with.
3. **28% of the set binds covalently.** Their potency comes from forming a permanent
   bond, which a reversible-binding score cannot represent.
4. **The measurements are not all against this protein.** ChEMBL mixes wild-type and
   mutant assays across labs; 2RGP is one wild-type conformation.

## What it cost to find out

About 35 minutes of laptop CPU and **$0**. No GPU was rented; nothing was deployed.

Had this gone in untested, every compound would have carried a confident-looking
docking score that was quietly anti-correlated with reality, and the loop would have
been steered towards weaker compounds — with no error, no warning, and no obvious
moment of discovery.

## What to try next, in order of cost

1. **4HJO** — the inactive-conformation EGFR structure. One command, ~35 minutes.
   Different conformations suit different chemotypes.
2. **Ensemble docking** across several structures, taking the best score per compound.
3. **GNINA** CNN rescoring of smina poses — a learned scoring function rather than a
   physics-based one.
4. **Covalent docking** for the 28%, which standard docking cannot address at all.
5. **Accept the result.** The fingerprint model at ρ 0.815 is doing well on this target.
   Docking's value was supposed to be generalising to new chemotypes — if it cannot
   rank the chemotypes we already have, that promise is unproven here.

The docking objective exists in the code (`chembl_with_docking` profile) and is
**switched off**. It stays off until a structure passes this test.

---

# 5. Docking as a screening filter: also negative

Section 4 showed docking cannot rank potency inside a series. That is the hard task, so
the fair follow-up was the easy one, and the one docking is actually known for:
**can it tell a real binder from a molecule with no reason to bind?**

```bash
python3 dock_screen_test.py 2RGP 40 3
```

39 EGFR actives (pChEMBL ≥ 8) against 119 decoys drawn from four **non-kinase** targets —
Factor Xa, carbonic anhydrase II, mu opioid receptor, COX-2 (18,819 candidates). Other
kinases were excluded deliberately: kinase inhibitors are promiscuous, so a kinase
"decoy" might genuinely bind EGFR.

Decoys were **property-matched**, and the match is reported so the test can be audited:

| group | MW | cLogP |
|---|---|---|
| actives | 486.9 ± 108.7 | 4.8 ± 1.5 |
| decoys | 483.6 ± 100.0 | 4.7 ± 1.5 |

| ranked by | ROC AUC | EF@1% | EF@10% |
|---|---|---|---|
| **docking score** | **0.457** | 0.0 | 0.25 |
| molecular weight (control) | 0.506 | 2.03 | 1.27 |
| random (control) | 0.429 | 0.0 | 0.51 |

**AUC 0.457 is below chance.** Docking scored worse than random at 10% enrichment
(0.25 vs 0.51) and worse than simply ranking by molecular weight. Two independent tests,
both negative.

## The methodological lesson, which cost nothing to learn and is worth keeping

The smoke test in section 4 looked convincing:

```
erlotinib  -8.01     caffeine  -5.63
gefitinib  -8.96     aspirin   -6.49
lapatinib -11.03     ethanol   -2.85
```

That separation was **an artefact of size**. The drugs average MW 450; caffeine, aspirin
and ethanol average MW 140. Docking scores rise with molecular size, so the test was
measuring size, not binding. Replace the tiny decoys with size-matched ones and the skill
disappears completely.

This is exactly why the decoy design above includes a molecular-weight control. Any
screening result reported without one should be treated as unproven — including a result
that looks as clean as that smoke test did.

## Conclusion on docking for this target

Docking into 2RGP with a standard Vina-type setup has **no use in this pipeline**, at
either task. Possible causes, honestly: EGFR's ATP pocket is large and flexible; the
compounds are big (mean MW 487) and floppy, which is where Vina-style scoring is known to
degrade; one rigid receptor with no waters; and 27% of the library binds covalently.

What would have to change before revisiting: an ensemble of receptor conformations,
retained hinge waters, covalent docking for the warhead subset, or CNN rescoring (GNINA).
Each is a research project with an uncertain payoff, and the evidence here does not
justify starting one.

**Spend the next effort on the model instead — Chemprop, then a foundation model. Those
are stages 1 and 2 of the roadmap, they cost under $10, and the fingerprint baseline they
have to beat (ρ 0.815) is already measured.**

Total cost of reaching this conclusion: about 70 minutes of laptop CPU. No GPU, no spend.

## 4HJO and the ensemble: also rejected

The obvious objection to the 2RGP result is that one rigid structure suits only the
chemotype it was crystallised with. So the test was repeated on **4HJO** — EGFR in its
*inactive* conformation with erlotinib bound, a genuinely different shape of the same
protein — and then on the **ensemble**, taking each compound's best score across both.

```bash
python3 dock_validate.py 4HJO 250 CHEMBL203
python3 compare_structures.py
```

| structure | n | ρ | top-10% hit rate | base | enrichment | separation |
|---|---|---|---|---|---|---|
| 2RGP (active-like, 2.0 Å) | 249 | 0.038 | 8% | 22% | 0.37× | −0.25 |
| 4HJO (inactive, 2.75 Å) | 249 | 0.049 | **0%** | 22% | **0.00×** | −0.24 |
| ensemble of both | 249 | 0.014 | **0%** | 22% | **0.00×** | −0.28 |
| *fingerprint model* | 250 | *0.815* | *80%* | *23%* | *3.51×* | — |

4HJO found **zero** actives in its top decile. So did the ensemble. Every separation
figure is negative: in all three cases actives score slightly *worse* than inactives.

Three independent attempts, no signal in any of them. This is no longer a bad structure
choice — it is the scoring function being unable to rank compounds within a congeneric
series, which is a documented limitation rather than a surprise.

**Decision: docking stays off for this target.** The remaining options (GNINA CNN
rescoring, covalent docking for the 28%) are substantially more work with no evidence
from these runs that they would clear 3.51×. The fingerprint model is already doing that
job well on EGFR.

Total cost of reaching a confident "no": about **70 minutes of laptop CPU and $0**. No
GPU rented, nothing deployed. The `chembl_with_docking` profile exists and remains
unused; `compare_structures.py` will re-run this comparison for any new structure in
seconds, since scores are cached.

---

# 5. Eleven real programmes, not one

A single good result is an anecdote. ChEMBL data can be split by source paper, and one
paper is usually one team optimising one series — far closer to a live programme than a
whole target, which is thousands of unrelated compounds from hundreds of groups.

EGFR alone contains **27 such programmes** of 40+ compounds, holding 2,674 compounds.

```bash
python3 programme_sweep.py CHEMBL203 50 12
```

| paper | n | chemotypes | spread | ρ (unseen chemotypes) | loop enrichment | speed to best |
|---|---|---|---|---|---|---|
| CHEMBL5728555 | 884 | 226 | 3.98 | **0.770** | 1.91× | 13.8× |
| CHEMBL5727085 | 175 | 72 | 2.00 | 0.540 | 1.38× | 11.0× |
| CHEMBL3638799 | 169 | 23 | 3.70 | **−0.099** | 1.41× | 10.6× |
| CHEMBL1240340 | 165 | 59 | 3.48 | 0.513 | 1.58× | 1.7× |
| CHEMBL1140899 | 101 | 41 | 4.93 | 0.243 | 1.34× | — |
| CHEMBL5725990 | 94 | 48 | 2.00 | 0.510 | 1.30× | 3.0× |
| CHEMBL5727286 | 77 | 23 | 2.00 | 0.514 | **0.93×** | 4.9× |
| CHEMBL3853319 | 66 | 49 | 2.80 | 0.373 | 1.77× | 2.1× |
| CHEMBL5725419 | 66 | 25 | 4.48 | 0.559 | 1.29× | 1.4× |
| CHEMBL1143196 | 64 | 12 | 4.34 | 0.767 | 2.00× | 1.3× |
| CHEMBL5725787 | 59 | 24 | 2.11 | 0.533 | 1.54× | 0.9× |

**Across 11 programmes:**

- Ranking unseen chemotypes: **median ρ 0.514**, range −0.10 to 0.77
- Loop enrichment over chance: **median 1.41×**, range 0.93 to 2.00
- Speed to the programme's best compound: **median 2.56×**
- The loop helped in **10 of 11 (91%)**

Two programmes were skipped before evaluation because their potency spans under 1.5 log
units — there is nothing in them to learn. One more had too few compounds to hold
chemotypes out. All three exclusions are reported rather than quietly dropped.

The failures are as useful as the successes. CHEMBL3638799 has **ρ −0.099** — the model
cannot rank its chemotypes at all — yet the loop still reached the best compound 10.6×
sooner, because within-series exploitation does not need cross-series generalisation.
CHEMBL5727286 is the opposite: the model ranks adequately (ρ 0.514) while the loop lands
**below chance** on enrichment at 0.93×.

**What this supports:** the loop usually helps on a real single-team campaign, with a
median saving worth having, and it will occasionally not. That is a claim with a
distribution behind it, not a single flattering number. Run
`programme_sweep.py` on any target before promising anything about it.

## A second target: VEGFR2

```bash
python3 -m dmta.cli chembl-fetch --target CHEMBL279 --refresh
python3 programme_sweep.py CHEMBL279 50 10
```

| paper | n | ρ | loop enrichment | speed to best |
|---|---|---|---|---|
| CHEMBL3638801 | 276 | 0.580 | 2.02× | 1.8× |
| CHEMBL5725448 | 216 | **0.759** | 2.00× | 1.3× |
| CHEMBL1240340 | 162 | 0.633 | 1.91× | 2.1× |
| CHEMBL3879910 | 190 | 0.590 | 1.59× | 4.0× |
| CHEMBL5728123 | 92 | 0.511 | 1.46× | 1.9× |
| CHEMBL3638799 | 169 | **−0.219** | 1.24× | 5.3× |
| CHEMBL3886771 | 240 | 0.371 | 1.17× | 1.0× |
| CHEMBL3639299 | 146 | 0.210 | 1.17× | 1.5× |
| CHEMBL3886699 | 63 | −0.022 | **0.84×** | 1.3× |
| CHEMBL5728351 | 73 | 0.601 | **0.76×** | — |

- Median ρ **0.545** (range −0.22 to 0.76), median enrichment **1.35×**, median speed **1.79×**
- The loop helped in **7 of 10 (70%)**

**Both targets together: 21 programmes, the loop helped in 17 (81%).**

One pattern is worth acting on. Every outright failure — CHEMBL5727286 (n=77),
CHEMBL5728351 (n=73), CHEMBL3886699 (n=63) — came from a programme of **under 100
compounds**, though several small programmes did succeed. Below about 100 compounds the
result is high-variance, so an evaluation on a small dataset should be treated as
indicative rather than decided.

The claim this supports, stated as narrowly as the evidence allows:

> On a real single-team medicinal-chemistry campaign of roughly 100 compounds or more,
> this loop reaches the programme's best compounds meaningfully sooner than testing in
> the order they were actually made — in about four cases out of five. In the fifth it
> does not, and the same evaluation tells you which case you are in before you commit.

---

# 6. Ninety-seven programmes — and a hypothesis that did not survive

```bash
python3 mega_sweep.py            # 10 targets, 704 seconds
python3 sweep_analysis.py
```

**97 programmes across 10 targets** spanning kinases, a lipid kinase, a protease, an ion
channel, an oxidoreductase, a lyase and a GPCR.

> **Overall, the loop helped in 77 of 97 programmes (79%).**

That is the headline number, and it holds across target families rather than resting on
kinases.

## The hypothesis that failed

After 21 programmes it looked as though every failure came from a programme of under 100
compounds. A tempting rule: *only take programmes with 100+ compounds.*

Across 97 programmes, by size:

| compounds | helped |
|---|---|
| < 60 | 22/27 (81%) |
| 60–100 | 22/33 (67%) |
| 100–150 | 11/11 (100%) |
| 150–250 | 15/19 (79%) |
| 250+ | 7/7 (100%) |

No monotone pattern. The original observation was a small-sample artefact.

Fitted properly — rule derived on 48 programmes, graded on the other 49:

| candidate rule | qualifies | turns away | separation |
|---|---|---|---|
| 80+ compounds and 2.5+ log spread | 83% | 77% | **+6%** |
| model scaffold ρ ≥ 0.6 | 84% | 77% | **+8%** |

Both are far too small to qualify anyone on. The by-bucket tables looked convincing —
ρ 0.4–0.6 showed 89% against 59% for ρ 0.2–0.4 — and **none of it survived a held-out
test.** Had the rule been fitted and graded on the same 97 programmes, it would have
looked publishable.

## What this actually means

**Nothing measurable in advance reliably predicts which programmes the loop will fail
on.** That is an uncomfortable result and it is the honest one.

It also matters less than it first appears, because the pre-qualification rule was only
ever a shortcut. The free retrospective evaluation does not *predict* whether the loop
helps on a customer's chemistry — **it measures it, on their data.** A rule would have
been a guess made before looking; the evaluation is the looking.

So the claim stands as:

> Across 97 real medicinal-chemistry programmes from 10 targets, this loop reached the
> best compounds meaningfully sooner than the order they were actually made in **79% of
> the time**. We cannot currently tell in advance which 21% it will not help — so we run
> the evaluation on your own historical data first, free, and you see your own number
> before committing anything.

## Open Source Malaria: what real open data looks like

OSM publishes everything, failures included. Of 147 compounds in its database:

| | |
|---|---|
| a usable number | **46** |
| censored ("> X" — too weak to measure) | 51 |
| never tested | 24 |
| single-point readings ("70% @ 120000") | 26 |

123 rows reached the cleaner; **68 compounds survived it**. The verdict was mixed and
reported as such: the loop found 7 of the top 7 and reached the best compound in 16
experiments against ~34 for random, while the model scored **ρ −0.2** on unseen
chemotypes — strong within the series, useless beyond it.

A tidy ChEMBL paper yields a far higher usable fraction than a real programme does.
Any intake pipeline has to survive that gap, and say what it discarded.
