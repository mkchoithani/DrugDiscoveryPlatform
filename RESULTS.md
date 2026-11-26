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

---

# 7. Chemprop: tested, and not adopted

A D-MPNN graph neural network, wrapped behind the same interface as the random forest
and given exactly the same scaffold splits of exactly the same compounds, so the only
difference between the two columns is the model.

```bash
python3 model_benchmark.py 18
```

**18 programmes across 4 targets, 179 seconds:**

| | random forest | chemprop |
|---|---|---|
| median ρ | **0.510** | 0.425 |
| median MAE | **0.403** | 0.438 |
| programmes ranked better | 11/18 | 7/18 (39%) |
| typical fit time | **0.7 s** | 29 s |

Median paired difference in ρ: **−0.026**, with a spread of 0.135 across programmes.

> **Verdict: no reliable difference.** The gap is far smaller than the variation between
> programmes. Keep the random forest — it trains in seconds, needs no GPU, and the
> neural network buys nothing measurable on this data.

Stated fairly to Chemprop: it ran with default hyperparameters, 40 epochs, no validation
split and no early stopping, while the forest had its uncertainty cross-validated. A
tuned D-MPNN would likely close the 0.026 gap. What the result rules out is a *large*
win — nobody is leaving 0.2 of ρ on the table by keeping the simpler model.

**Two bugs found while building this**, both of the silent kind:

1. **Segfault, exit 139, no traceback.** PyTorch and RDKit each ship an OpenMP runtime
   and loading both crashes on macOS. Fixed with three environment variables that must
   be set before torch is imported.
2. **Predictions came back one row short.** Chemprop's dataloader drops a trailing batch
   of one sample, because batch norm cannot normalise a single row — so a 129-compound
   test set returned 128 predictions and the arrays silently disagreed. `drop_last` is
   now forced off, the batch size is nudged when a remainder of one would occur, and a
   length check raises rather than returning a short answer.

The second one is the more instructive: had the shapes happened to broadcast, it would
have scored the wrong compounds against the wrong truths and reported a number.

---

# 8. RMT-denoised linear model: tested, and not adopted

Lee et al., PNAS 2019, 116, 3373 — random matrix theory to separate real chemical
correlations from undersampling noise. It was the better-motivated next model after
Chemprop lost 7: built for the "features ≫ molecules" regime we are actually in,
prospectively validated, and interpretable. Implemented in `dmta/rmt_model.py` and run
on the same 24 programmes, same scaffold splits, same protocol.

```bash
python3 model_benchmark.py 24 forest,rmt-project,rmt-clip
```

```
24 programmes in 1002s

  model         median rho  median MAE    wins   fit seconds
  forest             0.562       0.401      15          17.3
  rmt-project        0.389       0.447       5         457.9
  rmt-clip           0.258       1.709       4         527.0

  paired against the forest, programme by programme:
    rmt-project median -0.098   spread 0.226   2 s.e. 0.092   -> WORSE than the forest
    rmt-clip   median -0.127   spread 0.295   2 s.e. 0.121   -> WORSE than the forest
```

**Verdict: keep the random forest.** The better RMT variant is worse by 0.098 of ρ, which
is just outside twice its standard error, and it takes **27× longer to fit**. It does win
outright on 5 of 24 programmes, so it is not useless — but nothing here justifies
replacing the incumbent.

## Two variants, because "clean the correlations" is ambiguous

The method keeps the eigenvalues of the feature correlation matrix that lie above the
Marchenko–Pastur edge, σ²(1 + √(p/n))², and treats the rest as noise. What you do next is
not stated by the eigenvalue clipping itself, and the two readings behave very differently:

| variant | what it does | median ρ |
|---|---|---|
| `rmt-clip` | invert the cleaned matrix literally — generalised least squares | 0.258 |
| `rmt-project` | solve inside the signal subspace, discard the noise complement | 0.389 |

`clip` is the textbook formula and it is the wrong tool for ranking. Inverting *whitens*:
a signal direction with eigenvalue 90 gets damped by 1/90, while a noise direction sitting
at the flat level keeps weight ≈ 1. It systematically throws away the structure RMT just
found. Its median MAE of 1.709 — against the forest's 0.401 — is that effect showing up
in absolute terms.

## The bug this nearly shipped with

The first version scored ρ between −0.44 and +0.37, near zero on almost every programme,
and the obvious conclusion was "the paper's method does not transfer". That conclusion
would have been wrong. Inspecting the fitted model showed:

```
n_signal 325   edge 0.0   noise_level 0.0
noise directions get weight 1/floor = 100000000.0
```

The textbook recipe fits the bulk variance self-consistently: set the edge from a
variance, re-estimate the variance from whatever falls below the edge, repeat. With
p = 2060 features and n = 325 compounds there are 1,735 **zero** eigenvalues, and they
count towards the trace — so every pass pulls the estimate down, which lowers the edge,
which admits more eigenvalues as signal, which lowers it further. It converged to
edge = 0, declared all 325 eigenvalues signal, and clamped the noise floor to 1e-8,
giving the noise directions a weight of 10⁸.

For a *correlation* matrix the fit is also unnecessary: the diagonal is all ones, so the
trace is exactly p and the null variance is exactly 1. Removing the iteration took one
programme from ρ −0.074 to **0.714**.

Two tests now guard this: pure random data must yield **zero** signal eigenvalues, and a
fitted model must report an edge above 1.0. Either would have caught it immediately.

**The lesson is the one this repository keeps relearning.** A plausible negative result
is not self-validating. "The paper's method does not work here" and "I implemented it
wrongly" produce the same table, and only reading the fitted internals tells them apart.

---

# 9. Auditing our own benchmark against PostEra's DEL Simulator

Every number above section 1 comes from real measured compounds. The simulated benchmark
in the README does not: its **2.75× enrichment** was measured on a landscape we wrote
ourselves (`dmta/assay.py`'s `SimulatedAssay`). A selection strategy can look good on one
landscape for reasons that have nothing to do with the method, and nothing in this
repository could have revealed it.

PostEra's DEL Simulator (*Digital Discovery* 2026, 5, 4976, `postera-ai/del-simulator`)
exists for the same reason ours does, and gives us a second ground-truth landscape that
we did not author — with affinities anchored on **2,660 real soluble-epoxide-hydrolase
binders** rather than on a random weighting we chose. We generated its example library
(105,197 distinct molecules) and ran our own acquisition comparison inside it, holding
everything else fixed.

```bash
python3 del_benchmark.py 8 16 <del-data-dir> 42 7 101
```

```
                                  UCB over random
  our landscape (924 cmpds / 5 scaffolds)        2.75x
  their landscape, uniform sample                0.71x
  their landscape, dense sample                  0.70x
```

On their landscape the loop is **worse than random screening**, in both sample shapes.
The 2.75× reproduced exactly on ours, so the protocol was identical.

## But the mechanism is not what the headline suggests

Prospective ρ on their landscape was **−0.04 to +0.11** — the model had no predictive
signal at all. Active learning cannot beat random selection with a blind model; that is
close to a tautology, and it makes the comparison a poor test of the *acquisition
function*. So the real question is whether their landscape is learnable at all:

| landscape | n_train | ρ | MAE | R² |
|---|---|---|---|---|
| ours (simulated) | 100 | **0.892** | 0.298 | 0.707 |
| ours (simulated) | 700 | **0.965** | 0.180 | 0.891 |
| theirs (5k pool) | 700 | 0.121 | 0.476 | 0.161 |
| theirs (20k pool) | 2,000 | −0.033 | 0.444 | 0.185 |

**Their landscape is essentially unlearnable from Morgan fingerprints** — ρ ≈ 0 even with
2,000 training compounds. Their affinity is a kernel on Tversky similarity to the
reference actives, so it is driven by rare high-similarity events: a needle landscape, not
a smooth one. (R² stays mildly positive while ρ is ~0, meaning the model learns the mean
and not the ranking.)

So the finding is not "UCB is bad". It is this:

> **Our simulated landscape is far too easy.** ρ 0.89 from 100 training compounds, where
> the median across 24 real ChEMBL programmes is 0.56 (section 8) and the best real
> programme reaches 0.78. A landscape that a model can nail from 100 compounds will make
> any sensible acquisition function look good.

## What changes, and what does not

**Changed.** The README's 2.75× is now reported as what it is: the enrichment on a
landscape we wrote, which flatters the loop, and an upper bound rather than a forecast.

**Not changed — and this is the important part.** The loop's real-data evidence does not
rest on the simulator. Section 2 measures enrichment on **four real targets** with real
measured values revealed only on selection: 2.37×, 1.74×, 2.18× and 1.97×, with random
selection recovering each library's base rate within noise. Section 6 finds the loop
helps in 77 of 97 real programmes. None of that is affected by what our simulator does,
because none of it uses our simulator.

**The precondition is now explicit, and it is measurable.** The loop beats screening when
its model has ranking signal on the target, and fails when it does not. We already report
that gate — the scaffold-split ρ — and on real medicinal-chemistry programmes it is met
(median 0.56) while on a uniform draw from a 105k DEL it is not (≈0). That is a statement
about *when to run a loop at all*, which is more useful than the 2.75× ever was.

## An honest weakness in this audit

The "dense" arm was meant to separate "their landscape is different" from "their sample is
sparse", by taking products that share two of three reaction cycles. It barely moved the
variable: 754 distinct Murcko scaffolds in 924 compounds, against 903 for the uniform
sample and **5** for our own library. So that control did not cleanly test what it was
built to test, and the learnability table above — not the dense arm — is what actually
settles the question.

---

# 10. Predicted pharmacophore maps: tested, and not adopted

Boe, Robinson, Lee et al., *Digital Discovery* 2026, 5, 5001 — a model that takes a
rigid pocket and predicts, over a 3D grid, the probability that each of six interaction
feature types belongs at each point. It was the most promising unimplemented idea we had,
because it is the paper's own answer to why docking disappoints, and docking failed here
three times (sections 4–5).

**Their release and the half that is missing.** `postera-ai/phc4-field-generator` gives
featurization, model and inference — the *field*. It does not include the screening side:
taking that field and ranking a library against it. `dmta/pharmacophore.py` is that half.
It reduces the field to sparse typed hot spots, annotates each ligand conformer with
RDKit's `BaseFeatures` factory (whose families are the model's six channels, exactly), and
finds the largest geometrically consistent correspondence between ligand features and hot
spots — a maximum-weight clique where every pairwise distance must agree within a
tolerance. No pose search, so none of the pose-dependence the paper argues against.

Their model ran on our cleaned 2RGP receptor and crystal ligand — the same inputs the
failed docking used — in about four CPU-minutes for 79,507 grid points.

```bash
python3 pharmacophore_validate.py /tmp/phc4/out_2rgp/phc4_field.npz 1500 CHEMBL203
```

1,500 real EGFR compounds, measured IC50, "active" = this dataset's top decile
(pIC50 ≥ 8.70). Scoring ran at ~10 compounds/s; all 1,500 embedded.

| method | ρ | top decile | base | enrichment |
|---|---|---|---|---|
| pharmacophore score | 0.174 | 0.133 | 0.107 | **1.24×** |
| pharmacophore match count | 0.070 | 0.167 | 0.107 | 1.55× |
| **heavy atom count** | **0.340** | 0.147 | 0.107 | **1.37×** |
| QED | −0.293 | 0.073 | 0.107 | 0.68× |
| random forest (scaffold split) | 0.766 | 0.459 | 0.104 | 4.42× |

**Verdict: fails, and Phase B should not be built.** 1.24× against a bar of 1.5×, and —
the part that settles it — **it loses to ranking compounds by molecular size.** Heavy atom
count gets 1.37× and has twice the rank correlation, for no compute at all. Its
correlation with the pharmacophore score is only 0.412, so the score is not *merely*
size; but whatever else it contains is not helping.

## The time-split result, and why it is not a win

The second bar was signal on the time split, where the 2D model collapses. On the surface
it looks passed:

| method | ρ | enrichment |
|---|---|---|
| pharmacophore (no training) | 0.179 | 1.65× |
| random forest (trained on earlier compounds) | 0.133 | 1.37× |

**That difference is one compound.** The late period has 411 compounds, so the top decile
is 41, and a base rate of 0.107 makes those enrichments **7 actives versus 6**. Simulating
41 draws at a common hit rate, a gap that large or larger appears **65% of the time by
chance**. It is not evidence of anything, and reporting 1.65× against 1.37× as a win would
be exactly the kind of number this platform exists to catch.

One real observation survives: the pharmacophore score is **stable** across splits — ρ
0.174 over the whole set, 0.179 on unseen-chemotype compounds — while the forest falls
from 0.766 to 0.133. That stability is structural, since the score requires no training
and so cannot degrade when the chemotype shifts. It is the one property that would justify
revisiting this. But a stable 0.18 is not useful, and stability at a low level is not the
same as signal.

## Where this test is unfair to the paper, stated plainly

The paper claims efficient **hit finding** — picking actives out of a large library, tested
prospectively on AKR1C1. Our test is **potency ranking within a mature congeneric series**:
1,500 elaborated EGFR inhibitors that are nearly all active, where the question is which
is more potent. Those are different tasks, and the second is harder for a coarse
representation by construction — a pharmacophore says "this molecule can reach the right
contacts", which barely distinguishes two compounds that both can.

So the honest claim is narrow: **predicted pharmacophore maps do not help rank the next
batch in a programme**, which is what this platform does. It is not a finding that the
method fails at the job its authors claim for it. The fair test of that would be actives
against property-matched decoys (DUD-E or LIT-PCBA style), which we do not currently run.

Four further caveats belong to my implementation, not the model: only 3 conformers per
compound; matches capped at 8 points of a 19-point query; no steric or clash term, so a
molecule can satisfy the query without physically fitting; and clique matching rather than
the alignment-and-vector-database screening the paper describes.

## What it cost, and the running tally

About 20 CPU-minutes and no GPU. Including this, the score on PostEra's published methods
is now:

| capability tried | outcome |
|---|---|
| Docking, 3 structures | rejected (0.0–0.37×) |
| Chemprop D-MPNN | rejected (ρ 0.425 vs 0.510) |
| RMT linear model | rejected (ρ 0.389 vs 0.562) |
| Pharmacophore maps | rejected (1.24×, loses to heavy-atom count) |
| **Random forest on 2D fingerprints** | **still the incumbent, 4.42×** |

Four attempts to beat a random forest on Morgan fingerprints have now failed, three of
them methods that are more sophisticated by any reasonable measure. That is worth stating
as a finding rather than a series of disappointments: **on this task, at this data scale,
the simple model is hard to beat, and the remaining value is not in the model.**

---

# 11. Compounds too weak to measure: kept, and worth keeping — but not for the reason expected

A potency reported as `> 10 µM` is not missing data. It is a measurement: this compound
is definitely bad. `dmta/ingest.py` discarded every one of them, which threw away the
most confident statement the assay made and biased the training set towards actives.

This is not a corner case. In the Open Source Malaria file, **50 of 123 compounds** are
reported as `> X`. Ingest now keeps them:

```
sample_data/osm_series.csv:   68 observed, 50 censored (50 too weak, 0 too potent)
sample_data/acme_egfr_programme.csv: 436 observed, 28 censored
```

Open Source Malaria goes from **68 usable compounds to 118**.

## The direction trap, which gets its own test

| as reported | means | on the p-scale |
|---|---|---|
| `IC50 > 10 µM` | weaker than 10 µM | pIC50 **<** 5.0 — left-censored |
| `IC50 < 1 nM` | more potent than 1 nM | pIC50 **>** 9.0 — right-censored |
| `pIC50 > 8` | already log | pIC50 **>** 8 — right-censored, **no flip** |

pIC50 = −log10(molar) *decreases* as concentration rises, so the inequality flips —
except when the column is already on a log scale, where it does not. Reading `IC50 > 10
µM` as "potent" would invert the dataset silently, so `censoring_on_p_scale` is the one
place that decision is made and all four cases are a test.

## Does using them actually improve the model?

Known-answer test: real programmes where every potency *is* known, with the weakest 40%
of each training set artificially censored — roughly Open Source Malaria's real rate.
Four handlings on identical scaffold splits, scored against held-out true values.

```bash
python3 censored_benchmark.py 12 0.4
```

| arm | median ρ | **mean ρ** | worst | 10th pct | programmes with ρ < 0.3 |
|---|---|---|---|---|---|
| oracle (true values) | 0.581 | **0.511** | −0.070 | 0.286 | 2 |
| drop (previous behaviour) | **0.553** | 0.420 | −0.143 | **−0.053** | 4 |
| bound (naive substitution) | 0.515 | 0.462 | −0.127 | 0.126 | 3 |
| censored EM | 0.508 | **0.472** | −0.132 | **0.253** | 3 |

**The median says dropping is fine. The median is the wrong statistic here.** Paired
against dropping, programme by programme, the median difference is −0.004 with twice its
standard error at 0.093 — no reliable difference, and that is what the benchmark's own
verdict reported.

But the mean moves by **+0.052** and the tenth percentile moves by **+0.306**. Censored
compounds do not improve the typical programme. **They prevent the disasters.**

Where dropping does worst, the gap is large:

| paper | drop | censored | oracle |
|---|---|---|---|
| CHEMBL3886771 | **−0.143** | 0.253 | 0.367 |
| CHEMBL3638801 | **0.185** | 0.386 | 0.581 |
| CHEMBL5728058 | **0.188** | 0.439 | 0.404 |

Those are programmes where discarding 40% of the training data leaves the model with too
narrow a range of potency to learn anything — in one case anti-correlated with reality.
The censored compounds restore the range. (This three-programme slice is a post-hoc
grouping and is shown as illustration; the mean and tenth percentile above are the
evidence, and neither requires choosing a subgroup after the fact.)

**Verdict: keep them, as insurance rather than as an improvement.** The honest claim is
not "the model gets better" — on a typical programme it does not. It is "the model stops
failing badly on the programmes where potent compounds are the minority", which is a
property worth having and which the median conceals.

## Two further honest notes

**The EM machinery barely beats naive substitution.** Censored EM averages +0.052 ρ over
dropping; simply asserting the bound as a value averages +0.042. The sophisticated
treatment buys 0.010 ρ over the crude one. EM is kept because it is principled and costs
nothing at this scale, but it is not where the value came from — **keeping the compounds
at all is**.

**Even the oracle only gains +0.091.** Having the *true* values for the weak compounds
improves the mean by less than 0.1 ρ. So the ceiling on this whole line of work is low,
and censored EM already captures about 57% of it. That is a useful thing to know before
anyone invests further in censored-data modelling.

---

# 12. The fair test for docking and pharmacophores: hit finding, not potency ranking

Docking and pharmacophore maps were both rejected here (sections 4–5 and 10), and both
rejections carried the same caveat: **neither method claims to do the job we tested it
on.** We measured potency ranking inside a mature EGFR series, where nearly every
compound is already active. Both are hit-finding tools — pull the few actives out of a
library of mostly-inactive molecules. That is a coarser question, and a coarse
representation can answer it while being useless at the fine one.

So this is the test they are entitled to.

**Benchmark.** DUD-E (Mysinger et al., *J. Med. Chem.* 2012, 55, 6582), EGFR target:
542 known actives and 35,050 decoys **matched to the actives on physical properties**
(MW, logP, rotatable bonds, H-bond counts, charge) while topologically dissimilar.
Decoys subsampled to 5,000; actives split by scaffold, decoys split at random, and every
arm scored on the same held-out half — 271 actives against 2,500 decoys, 9.78% active.

**Metrics changed with the question.** Hit finding is judged on whether actives land at
the *top*, because that is all anyone screens: enrichment factor at 1%, and BEDROC
(α = 20, weighting roughly the top 8%). ROC AUC is reported but is the least informative
of the three here.

```bash
python3 hitfinding_benchmark.py /tmp/dude /tmp/phc4/out_2rgp/phc4_field.npz 5000
```

| method | AUC | **EF 1%** | EF 5% | **BEDROC** |
|---|---|---|---|---|
| pharmacophore (2RGP field) | 0.580 | **1.46** | 1.33 | **0.163** |
| heavy atom count | 0.516 | **6.57** | 3.31 | **0.346** |
| QED | 0.372 | 0.00 | 0.07 | 0.026 |
| 2D forest *(see below)* | 0.996 | 10.23 | 10.23 | 0.991 |

**Verdict: the earlier rejection stands. It was not an artefact of the wrong task.**

The pharmacophore field does beat chance on AUC, and it beats the size baseline there
(0.580 against 0.516). But it is **4× worse on enrichment at 1%** (1.46 against 6.57) and
2× worse on BEDROC. It sorts the whole list a little better than size does, while putting
far less in the top 1% — and the top 1% is the only part a screen ever looks at.

## My own verdict function got this wrong first

The benchmark's automated verdict compared **AUC alone**, and printed:

> "on hit finding the pharmacophore reaches AUC 0.580, beating the size baseline (0.516).
> It was rejected on the wrong test. Revisit it."

That conclusion was wrong, and it was wrong in the direction I was hoping for — I had
written the caveat about unfair testing myself, so a result vindicating it was the one I
expected. A method can hold a respectable AUC by sorting the bottom 90% well. The verdict
now requires EF at 1% **and** BEDROC to beat the baseline, and it prints the opposite
conclusion on the same numbers.

## The 2D forest row is a measurement of the benchmark, not of the model

AUC **0.996** and EF at 1% of **10.23** — which is the arithmetic ceiling at a 9.78%
active rate. A fingerprint model put essentially every active in the top 1%.

That is DUD-E's documented bias, measured. Its actives are literature compounds from
medicinal-chemistry papers; its decoys are ZINC catalogue molecules. A 2D model can
separate "looks like a published kinase inhibitor" from "looks like a catalogue compound"
almost perfectly **without knowing anything about EGFR**. Property matching controls for
size and lipophilicity; it cannot control for the fact that the two sets come from
different universes. This is precisely the critique that led to LIT-PCBA, and the
docstring in `dmta/hitfinding.py` predicted it before the run.

**Consequence: this benchmark cannot adjudicate between a structure-based method and a
ligand-based one.** It can only tell us whether a structure-based method clears the cheap
descriptor baselines, which the pharmacophore field does not.

One more sign the benchmark is imperfectly matched: **heavy atom count alone reaches EF
1% of 6.57.** If the decoys were truly property-matched, molecular size should carry no
early-enrichment signal at all.

## What this costs the earlier conclusions

Nothing is overturned, and one caveat is retired. Sections 4–5 and 10 each said the fair
test had not been run. It has now, on the standard benchmark for the purpose, and the
answer is the same: **neither 3D method beats a cheap descriptor here.** The remaining
honest caveat is narrower — a cleaner benchmark (LIT-PCBA) might rank them differently,
and our implementation of the screening half is our own.

---

# 13. LIT-PCBA: the unbiased benchmark, and the strongest version of the test

Section 12 ran the hit-finding test on DUD-E and then had to discount its own result:
a 2D fingerprint model reached **AUC 0.996** and the arithmetic ceiling on enrichment,
because DUD-E's actives are literature compounds while its decoys are untested ZINC
molecules. A fingerprint tells those two populations apart without knowing anything
about the target, so that benchmark cannot compare a structure-based method with a
ligand-based one.

LIT-PCBA (Tran-Nguyen, Jacquemard and Rognan, *J. Chem. Inf. Model.* 2020, 60, 4263)
exists for exactly that reason. Actives **and** inactives come from the same PubChem
confirmatory bioassays, so an inactive is a *measured* inactive. Its splits are
AVE-unbiased, removing the train/validation analogue redundancy. And the class balance
is real.

**Target: MAPK1 (ERK2).** A kinase, so closest to the EGFR work. Pharmacophore field
from **4QTE** — ERK2 with VTX-11e at 1.5 Å. LIT-PCBA designates no single reference
structure, so that choice is a judgement call and is recorded as one.

```
train: 231 actives / 46,317 inactives
test:   77 actives / 15,250 inactives     test hit rate 0.50%
```

Inactives were deliberately **not** capped. The 0.5% hit rate is the point — DUD-E's
10% is twenty times easier — and capping would have inflated every enrichment figure.

```bash
python3 litpcba_benchmark.py /tmp/litpcba/ave/MAPK1 /tmp/phc4/out_4qte/phc4_field.npz
```

| method | AUC | EF 1% | EF 5% | BEDROC |
|---|---|---|---|---|
| pharmacophore (4QTE field) | 0.550 | 2.60 | 1.04 | 0.082 |
| heavy atom count | **0.581** | 0.00 | 0.26 | 0.041 |
| QED | 0.391 | 0.00 | 0.26 | 0.021 |
| 2D forest | **0.673** | **6.50** | **4.16** | **0.202** |

## First result: the unbiasing works

| | DUD-E | LIT-PCBA |
|---|---|---|
| 2D forest AUC | 0.996 | **0.673** |
| 2D forest EF 1% | 10.23 (the ceiling) | **6.50** |

The fingerprint model's near-perfect separation collapses. **So this benchmark can
adjudicate between a structure-based method and a ligand-based one**, and DUD-E could
not. That is the most useful thing in this section, and it retires the open caveat left
at the end of section 12.

## Second result: count the molecules before believing anything

The benchmark's automated verdict said the pharmacophore "beats molecular size on early
enrichment (EF1% 2.6 vs 0.0) — on an unbiased benchmark that is a real result". It is
not, and the script printed the reason immediately above its own conclusion: the top 1%
is 153 compounds, holding **0.77 actives at chance**.

Converting every enrichment into molecules, with a hypergeometric tail probability:

| method | EF 1% | actives in the top 153 | P(≥ that) by chance |
|---|---|---|---|
| pharmacophore | 2.60 | **2** | **0.18** |
| heavy atom count | 0.00 | 0 | 1.00 |
| 2D forest | 6.50 | **5** | **0.001** |

**The pharmacophore found two actives where chance gives 0.8.** That happens 18% of the
time with no method at all. Its apparent 2× advantage over molecular size on BEDROC
(0.082 against 0.041) is the same two compounds seen from another angle, and it *loses*
to molecular size on AUC (0.550 against 0.581).

The 2D forest's five actives, at p = 0.001, is real signal.

## Verdict

**The rejection stands, and this is now the strongest form of the test.** On the one
benchmark that is unbiased by construction, judged on the metric a screen actually uses:

- the pharmacophore field is **statistically indistinguishable from chance** (p = 0.18);
- the plain 2D random forest finds **2.5× its early enrichment** and is significant;
- and the structure-based method still loses to a fingerprint, on a benchmark
  specifically designed not to favour fingerprints.

Three benchmarks have now been run against this method — potency ranking in a
congeneric series, DUD-E, and LIT-PCBA — and all three say the same thing.

## What I got wrong twice, in the same direction

Section 12's verdict function compared AUC alone and declared the pharmacophore a
winner. I fixed it to require early enrichment. **This section's verdict function then
declared it a winner again**, because it compared enrichment without asking how many
molecules the enrichment was made of.

Both errors ran in the direction I had reason to want: I wrote the "we tested it on the
wrong task" caveat myself, so a result vindicating it was the one I was looking for.
A verdict function is not a safeguard against that — it encodes the same expectations as
the person who wrote it. **What caught it both times was converting a ratio back into a
count of actual molecules**, which is now printed next to every enrichment figure this
benchmark reports.

The honest remaining caveat is narrow: one target, 77 actives, and our own
implementation of the screening half. A method that is chance-level on MAPK1 might not
be on a target whose pocket is more distinctive. But nothing here supports building
Phase B.

---

# 14. Docking on the hit-finding test: the first positive signal in this project

Docking was rejected three times on potency ranking — 0.37× on 2RGP, 0.00× on 4HJO and
on the ensemble (sections 4–5). Section 12 gave the pharmacophore field its fair test on
DUD-E; this is docking's. 1,021 held-out compounds (decoys capped at 1,500, so 26.5%
active), smina against 2RGP, **6,884 seconds — 6.8 s per compound**, 6 embed failures
and 1 dock failure out of 1,021.

| method | AUC | EF 1% | EF 5% | BEDROC |
|---|---|---|---|---|
| **docking (smina, 2RGP)** | **0.629** | 3.74 † | 2.42 | **0.606** |
| heavy atom count | 0.527 | 3.77 † | 2.81 | 0.592 |
| pharmacophore (2RGP field) | 0.573 | 1.88 | 1.55 | 0.382 |
| QED | 0.357 | 0.00 | 0.22 | 0.087 |
| 2D forest *(benchmark artefact)* | 0.994 | 3.77 † | 3.77 | 0.999 |

† **saturated — see below.**

**Docking is the best structure-based method here, and the first one in this project to
show real signal.** AUC 0.629 against the size baseline's 0.527, and BEDROC 0.606 against
0.592. It also beats the pharmacophore field on both (0.629/0.606 against 0.573/0.382).

That is a genuine change from sections 4–5. Docking cannot rank potency within a
congeneric series, and it can tell an EGFR-like molecule from a non-binder at better
than chance. Those are different jobs and it is competent at one of them.

## The metric that had no resolution left, and the verdict that used it anyway

Enrichment at 1% is capped at 1/base_rate. At a 26.5% active rate the cap is **3.77**,
and on 1,021 compounds the top 1% is **ten molecules**. So:

| method | EF 1% | actives in the top 10 |
|---|---|---|
| heavy atom count | 3.77 | 10.0 |
| 2D forest | 3.77 | 10.0 |
| docking | 3.74 | 9.9 |

Three methods pinned at the ceiling. The benchmark's verdict compared docking's 3.74
against molecular size's 3.77 and printed *"does NOT beat molecular size"* — a comparison
of 9.9 actives against 10.0, out of 10. Meaningless, and in this case it gave the wrong
answer.

**That is the third verdict-function error in three sections**, and notably the first in
the *conservative* direction — the previous two flattered the method I had reason to
favour. The common cause is the same: a ratio was compared without asking how many
molecules it stood for. `hit_finding_report` now returns `ef_ceiling` and
`ef_1pct_saturated`, and the verdict refuses to decide on a saturated column, falling
back to AUC and BEDROC.

Why the active rate was 26.5% at all: decoys were capped at 1,500 to keep a 6.8 s/compound
docking run inside a few hours. That cap is what saturated the metric. The uncapped DUD-E
run in section 12 was 9.78% active and did not saturate.

## This result is on the benchmark we already disqualified

Section 13 established that DUD-E cannot adjudicate structure-based against ligand-based
methods: its 2D forest reaches AUC 0.994 here, essentially perfect, because its decoys are
untested ZINC molecules. **So docking's AUC 0.629 is measured on a set whose separability
is partly artificial**, and it is not yet evidence that docking would help on real
screening data.

**The decisive test is now running:** the same docking arm on LIT-PCBA MAPK1, where
inactives are experimentally measured and the AVE unbiasing works (section 13). Inactives
capped at 2,000 for a 3.71% hit rate — stated because the cap raises the apparent rate,
and chosen so EF 1% has headroom rather than saturating. At 6.8 s per compound that is
about four hours for 2,077 compounds. **Until it reports, the honest summary is: docking
shows its first real signal, on the one benchmark we have already shown cannot be
trusted for this comparison.**

**It has since reported, and it went against this section — see section 15.** The signal
did not survive: docking's edge over molecular size falls to +0.045 AUC with a confidence
interval that crosses zero, and its early enrichment is one molecule in the top 21.

---

# 15. The decisive test: docking's signal did not survive an unbiased benchmark

Section 14 ended on an open question. Docking had beaten molecular size on DUD-E
(AUC 0.629 against 0.527) — the first positive structure-based result in the project — but
DUD-E was already disqualified for exactly this comparison. The same docking arm has now
run on LIT-PCBA MAPK1, where every inactive was experimentally measured.

```bash
python3 -u litpcba_benchmark.py /tmp/litpcba/ave/MAPK1 \
        /tmp/phc4/out_4qte/phc4_field.npz 2000 --dock
```

2,077 test compounds (77 actives, inactives capped at 2,000 for a 3.71% hit rate),
docked against 4QTE. 2,066 of 2,077 produced a pose, in 87 minutes at 2.8 s/compound.

| method | n | AUC | EF 1% | EF 5% | BEDROC |
|---|---|---|---|---|---|
| pharmacophore (MAPK1 field) | 2075 | 0.551 | 2.57 | 1.04 | 0.105 |
| heavy atom count | 2077 | 0.577 | 0.00 | 0.26 | 0.053 |
| QED | 2077 | 0.384 | 1.28 | 0.26 | 0.031 |
| **2D forest** | 2077 | **0.723** | **6.42** | **3.11** | **0.219** |
| docking (smina, 4QTE) | 2066 | 0.622 | 1.28 | 0.52 | 0.100 |

## Count the molecules first

The discipline that caught all three verdict errors, applied before reading anything into
the table. The top 1% is **21 compounds**; at chance it holds 0.8 of the 77 actives.

| method | EF 1% | actives in the top 21 | P by chance |
|---|---|---|---|
| pharmacophore | 2.57 | 2 | 0.18 |
| heavy atom count | 0.00 | 0 | 1.00 |
| QED | 1.28 | 1 | 0.55 |
| **2D forest** | **6.42** | **5** | **0.0008** |
| **docking** | 1.28 | **1** | **0.55** |

**Docking's early enrichment is one molecule.** On the metric that decides whether a
screen is worth running — what you find in the fraction you can afford to test — docking
is indistinguishable from drawing at random. Only the 2D forest clears chance.

## The AUC gap, tested properly

AUC 0.622 against molecular size's 0.577 still looks like a win. It is not one. Because
every method ranks the *same* compounds, the two AUCs are correlated, and an unpaired
error bar would overstate the uncertainty of their difference. The pose cache makes the
docking scores free to re-read, so the difference was bootstrapped directly:

```bash
python3 litpcba_paired_auc.py /tmp/litpcba/ave/MAPK1 4QTE 2000
```

4,000 resamples, actives and inactives resampled separately to hold the class ratio.

| comparison | mean ΔAUC | 95% CI | P(Δ ≤ 0) |
|---|---|---|---|
| docking − size | **+0.045** | **[−0.013, +0.099]** | 0.063 |
| docking − 2D forest | −0.100 | [−0.174, −0.023] | 0.995 |
| 2D forest − size | +0.145 | [+0.070, +0.217] | <0.001 |

**Docking's advantage over counting heavy atoms includes zero.** The forest's advantage
over both does not. And docking is *significantly worse* than the 2D fingerprint model it
was supposed to beat on the one task where 3D should win.

On DUD-E the docking-minus-size gap was +0.102. Here it is +0.045 with a CI that crosses
zero. **Roughly half the DUD-E gap was decoy bias**, and what remains cannot be
distinguished from noise.

## Where this test is weak, stated plainly

The confidence interval's upper bound is +0.099. With 77 actives this test **cannot
exclude** a true gap as large as the one DUD-E showed. "Not significant" is not "shown to
be zero", and I am not claiming it is. Two things keep the verdict where it is anyway:

- The enrichment result is not a near miss. One active in the top 21 is what chance
  produces 55% of the time. A method with real early signal does not look like that.
- This is now the **fourth** independent test docking has failed here — three on potency
  ranking (sections 4–6) and one on hit finding. The DUD-E arm is the only positive, and
  the benchmark it ran on was disqualified before the arm existed.

A fifth test worth running one day: more actives (a second LIT-PCBA target, or MAPK1
uncapped at 61,567 inactives), and GNINA's CNN rescoring rather than smina's empirical
function. Neither is likely to overturn an enrichment of one molecule.

## What this costs the previous section

Section 14's conclusion — *"docking shows its first real signal"* — stands as written,
because it was explicitly conditioned on a benchmark it also called untrustworthy. But the
condition has now resolved against it. The correct summary of docking in this project:

> **Docking failed potency ranking on three structures and an ensemble. On hit finding it
> beat molecular size on a benchmark whose decoys were never tested, and failed to beat it
> on a benchmark whose inactives were. It is not in the product.**

## The running tally, with "rejected" split into what it actually means

The word was doing three different jobs, and collapsing them overstated the evidence.
A capability that failed a fair test of its own claim is not the same as one that showed
no gain *on our data*, and neither is the same as one where the thing we tested was our
own code.

| kind | what it licenses us to say |
|---|---|
| **A · failed a fair test** | The method's own claim was tested and did not hold. An independent published negative, where one exists, makes this stronger. |
| **B · no evidence of gain here** | A tie. The paired interval contains zero. Keeping the incumbent is a *cost* decision, not a measured win — and it says nothing about the method in the regime its authors validated. |
| **C · our implementation failed** | We built the part the authors do not release. The method may be fine; our version of it is not. |

| capability | tests | kind | outcome |
|---|---|---|---|
| the DMTA loop itself | 4 targets, 97 programmes | — | **passed** — 2.37× base rate, 77/97 won |
| censored values | OSM, 68 → 118 compounds | — | **kept**, for the stated reason only |
| docking (smina) | 3 structures + ensemble + DUD-E + LIT-PCBA | **A** | rejected ×4. Sections 4–6 and 15. Agrees with the LIT-PCBA authors' own finding, and sections 4–6 re-confirmed a negative the docking literature had already established |
| RMT, as our potency predictor | 2 modes, 24 programmes | **A** | rejected — −0.173 median ρ (section 8), and the decision survives relocation: on the authors' own task it still loses to the forest, −0.036 AUC, CI [−0.053, −0.014] (section 17) |
| RMT cleaning, the paper's own claim | 9 LIT-PCBA targets, cleaned vs uncleaned | **B** | **not a rejection** — paired CI [−0.028, +0.006] contains zero. No detectable gain over a cross-validated ridge on the raw matrix, and no evidence against the method either (section 17) |
| predicted pharmacophore maps | potency, DUD-E, LIT-PCBA | **C** | rejected ×3 — PostEra release the field generator, not the scorer. The scorer is ours, so three agreeing benchmarks reject *our* screening half (section 10) |
| Chemprop (D-MPNN) | 24 programmes + a 250→7,173 learning curve | **B** | **not a rejection** — paired CI [−0.097, +0.024] contains zero. The gap plateaus at −0.034 from n≈2,000 upward, so the forest is level with it at every size this data reaches, at 96× less compute (section 16) |
| de novo generation (BRICS) | smoke test | — | incomplete, parked |
| **the model on logD 7.4** | 4,198 compounds, random + scaffold | — | **passed** — ρ 0.757 on held-out scaffolds against cLogP's 0.219, and the control that strips cLogP out of the model's own features still wins by +0.478 (section 18) |
| **the model on solubility** | 884 compounds, random + scaffold | **B** | **the forest is not shipped for it** — ties a 12-descriptor ridge on ranking, CI [−0.033, +0.032], and loses on error, MAE 1.053 vs 0.731. `choose_model` now measures the two with scaffolds held out and ships the ridge, at the cost of a worse σ (section 18) |

**Eight capabilities examined, across ten claims. Three kept — the loop, censored values,
and now logD. Two real rejections with independent support, one rejection that is about
what we built rather than what was published, three ties — one we had been calling a
rejection, one we had never tested at all, and one where the tie on ranking hides a loss on
error. Five of our own published claims corrected, five verdict-function errors found and
pinned with tests, and one defect in a published metric (section 18) found by a benchmark
arm that existed only to be beaten.**

Section 17 is the one that cuts both ways and is worth reading as a pair: the rejection of
RMT *as our model* got stronger, and the rejection of *the paper's claim* dissolved into a
tie — because until section 17 we had never tested the paper's claim at all.

The thing that keeps earning its place is not a model. It is the habit of converting every
ratio back into a count of molecules, and every gap into a confidence interval, before
believing it — and then applying that same standard to the conclusions we have already
written down.

---

# 16. Applying section 15's standard to our own model rejections

Section 15 bootstrapped a paired difference and found docking's apparent edge over
molecular size contained zero. The model comparisons in sections 7 and 8 never got that
treatment — they were decided on a median gap and a 2 s.e. screen, then summarised in the
tally with the word **rejected**. That word claims more than a median supports, so the
same test was applied to them.

```bash
python3 model_benchmark.py 24 forest,chemprop,rmt-project --save /tmp/model_scores.json
python3 model_paired_test.py /tmp/model_scores.json
```

24 programmes from 4 targets, identical scaffold splits, 616 s. The **programme** is the
unit of resampling, not the compound: a model that wins on large programmes and loses on
small ones has not won, and resampling compounds would hide that.

| | forest | chemprop | rmt-project |
|---|---|---|---|
| median ρ | **0.562** | 0.425 | 0.389 |
| median MAE | **0.401** | 0.435 | 0.447 |
| programmes won | 12 | 8 | 4 |
| total fit seconds | **10.6** | 195.3 | 409.9 |

## The two verdicts disagree, and that is the finding

| | mean paired Δρ | 95% CI (bootstrap) | Wilcoxon | sign test | verdict |
|---|---|---|---|---|---|
| **rmt-project** | −0.121 | **[−0.214, −0.032]** | P = 0.006 | 5/24, P = 0.007 | **a real rejection** |
| **chemprop** | −0.036 | **[−0.097, +0.024]** | P = 0.281 | 9/24, P = 0.308 | **a tie** |

**The RMT rejection survives.** All three tests agree, the interval excludes zero, and the
effect is large. Section 8's conclusion stands as measured — with the separate caveat,
already in section 8, that we used RMT as a *predictor* while the paper is about cleaning
a covariance matrix.

**The Chemprop rejection does not survive.** The interval contains zero. Nine wins out of
24 is what a coin gives 31% of the time. Section 7's own text said *"no reliable
difference"* and was right; the **tally's label was wrong**, and the deck repeated it.

## What the interval's width says, which is the part worth keeping

The honest reading is not "the gap is zero". It is:

> The interval runs to **−0.097**. A D-MPNN worse than the forest by up to a tenth of ρ is
> still consistent with this data — and so is one **better by 0.024**. Twenty-four
> programmes of 60–884 compounds cannot resolve a difference this size.

So the defensible statement about Chemprop is a **cost** decision, not a measured win:

> No detectable gain, **18× the fit time** (195 s against 10.6 s across 24 programmes), a
> torch dependency, and two silent bugs found in one afternoon of integrating it. We do
> not pay for that. **This is not evidence that Chemprop is worse.**

## The sign test, which is what gets quoted

Section 7 reported "11 of 18 programmes". That reads as decisive and is not: **P = 0.48**
two-sided. A win record discards the size of every difference, and is the weakest of the
three tests here precisely because it is the easiest to quote. `model_paired_test.py`
prints it alongside the other two for that reason, and a test pins its values.

## The correction

| | before | after |
|---|---|---|
| Chemprop | "rejected" | **no evidence of a difference** — a tie; the forest is kept on cost |
| RMT | "rejected" | rejected, and now with a CI that excludes zero |
| the tally | one word, "rejected" | three kinds: failed a fair test / no gain here / our implementation |

**That is the fourth of our own published claims to be corrected**, and the first where the
error was in a *label* rather than in a number. The number in section 7 was right all
along; the summary of it was not. A tally is where overclaiming hides, because it is the
part people read.

## The learning curve: Chemprop in the regime its paper validated

Section 7 ran Chemprop on programmes of 175–884 compounds. That is the small-data regime,
where fingerprints plus a forest are documented to be competitive, and it is *not* the
regime the D-MPNN papers validate. So the test we owed the paper: one target's entire
cleaned set, one **fixed** scaffold test set, training set subsampled upward.

```bash
python3 chemprop_scale.py CHEMBL203 IC50 250,500,1000,2000,4000,0
```

EGFR, 8,978 compounds → 7,173 train / 1,805 test. 17 minutes.

| n train | forest ρ | chemprop ρ | paired Δ | forest MAE | cp MAE | forest s | cp s |
|---|---|---|---|---|---|---|---|
| 250 | 0.701 | 0.468 | **−0.233** | 0.861 | 1.019 | 0.9 | 14.5 |
| 500 | 0.751 | 0.637 | −0.114 | 0.726 | 0.871 | 0.6 | 26.8 |
| 1,000 | 0.768 | 0.686 | −0.082 | 0.712 | 0.810 | 0.9 | 52.5 |
| 2,000 | 0.776 | 0.750 | −0.026 | 0.700 | 0.723 | 1.3 | 103.7 |
| 4,000 | 0.796 | 0.754 | −0.042 | 0.668 | 0.695 | 2.3 | 207.5 |
| **7,173** | **0.810** | 0.777 | **−0.033** | 0.646 | 0.662 | **3.8** | 365.7 |

**Two regimes.** The gap collapses from −0.233 to −0.026 as data goes from 250 to 2,000,
and then stops: over the last three sizes the slope is **−0.014 per decade** and the level
is **−0.034 ± 0.008**. Chemprop is not catching up any more; both models improve together.

So the data-size hypothesis is **half right, and the half that matters is wrong**. Size
explains the large gap at small n. It does not explain the residual, and the residual
(−0.034) is inside the tie interval from the paired test — so at every size this target can
reach, the two models are level and the forest costs 96× less (3.8 s against 365.7 s).

> **There is no n in reach where switching pays.** That is a stronger conclusion than
> section 7's, and it is the first time measuring a rejection properly has made the
> original decision *more* defensible rather than less.

## My verdict function got this wrong too — the fourth time, and the same cause

The script first fitted **one log-linear slope** through all six points:

```
  trend: +0.125 of rho per 10x more training data
  -> The gap is closing as data grows (+0.125 per decade) but has
     not crossed by n=7173. Size is part of the explanation.
     Re-test if a programme ever reaches several times this much data.
```

That slope extrapolates to a crossover at **n = 7,111**. The measured difference at
n = 7,173 is **−0.033** — nowhere near zero. The slope is carried almost entirely by the
n = 250 point, and the verdict it produced would have sent us looking for a crossover that
the last three points rule out.

**Same class of error as the saturated enrichment factor in section 14:** a summary
statistic used where it is not valid. Fit a straight line through a curve with a knee in
it and the line describes the knee, not where you are standing.

`read_trend()` now fits the recent regime separately, refuses to extrapolate a whole-curve
slope that disagrees with it, and distinguishes *plateaued* from *flat from the start* from
*still closing* from *crossed*. A test pins all four against this real curve.

**Four verdict-function errors now, in four sections.** Two flattered a method I had reason
to favour, one was conservative, and this one invented a future crossover. The only thing
that has ever caught them is refusing to read a summary number without asking what it is a
summary *of* — molecules, in sections 13–15; a two-regime curve, here.

---

# 17. RMT on the task its authors actually validated

Section 8 rejected `dmta/rmt_model.py` on 24 ChEMBL programmes: forest median ρ 0.562
against RMT's 0.389. The effect was large and properly measured, and the rejection stood
for a day. Then [PAPERS.md](PAPERS.md) recorded, against our own result:

> We used it as a predictor. The paper is about cleaning a covariance matrix, and never
> claimed to beat a random forest at QSAR. The effect we measured is large, so rejecting
> it *as our predictor* is safe — but **its actual contribution is untested here.**

Lee *et al.*, PNAS 2019, **116**, 3373 clean the correlations between chemical features
computed from known actives *and* known inactives, and use the cleaned matrix to pick
candidates out of a library. Their reported success is a virtual screen — prospectively,
four confirmed M1 muscarinic agonists. Not a Spearman correlation against continuous
potency. **Section 8 graded the method on a task it was never offered for.**

This section runs the task it was offered for. `rmt_screen_test.py`, 15 LIT-PCBA targets,
9 with enough actives in the validation split to measure anything.

## What is compared, and why it is not "RMT versus our forest"

A method can be worth less than the incumbent and still have a real contribution. Those
are different questions, and section 16 is the record of what conflating them costs. So
the headline here is:

> **rmt-cleaned versus rmt-uncleaned**

Both arms take the identical path through `RmtModel` — same features, same
standardisation, same generalised-least-squares solve, same prediction. The only
difference is whether the Marchenko–Pastur clipping happens. Whatever separates them **is
the cleaning**, and nothing else. `clean=False` is an experimental control, not a second
method, and a test pins it to be exactly `(C_raw + λI)⁻¹` computed densely — because if
the control is quietly a second cleaned model, the whole comparison measures nothing.

The control needs that ridge to exist at all: with p = 2,060 features and n in the
hundreds, the raw correlation matrix is singular. That ridge is precisely the ad-hoc
regularisation RMT is proposed as a principled replacement for, so choosing it badly would
hand the comparison to RMT by default. **It is chosen by cross-validation on the training
split only**, per target, and the chosen value is printed.

## The regime is the paper's, and the benchmark is unbiased

LIT-PCBA's inactives are experimentally measured, not untested catalogue molecules, and
the splits are AVE-unbiased — so a 2D fingerprint cannot win by recognising ZINC the way
it did on DUD-E (section 12). Hit rates run 1.2%–7.6% before the inactive cap. Features
outnumber training actives by 5× to 30× on eight of the nine targets. That is the
"features ≫ molecules" regime the paper addresses.

## AUC per target

| target | n | actives | rmt **clip** | rmt **project** | **uncleaned** | ols | forest | size |
|---|---|---|---|---|---|---|---|---|
| ALDH1 | 3343 | 1343 | 0.679 | 0.683 | 0.708 | 0.671 | **0.737** | 0.572 |
| ESR1_ant | 819 | 25 | 0.493 | 0.529 | 0.493 | 0.497 | 0.490 | **0.557** |
| FEN1 | 2091 | 91 | 0.833 | 0.893 | 0.918 | 0.625 | **0.924** | 0.499 |
| GBA | 2041 | 41 | 0.615 | 0.785 | 0.757 | 0.566 | **0.807** | 0.639 |
| KAT2A | 2048 | 48 | 0.596 | 0.670 | 0.677 | 0.539 | **0.703** | 0.521 |
| MAPK1 | 2077 | 77 | 0.645 | 0.649 | 0.694 | 0.603 | **0.723** | 0.577 |
| MTORC1 | 2024 | 24 | 0.558 | 0.624 | 0.644 | 0.507 | **0.683** | 0.524 |
| PKM2 | 2136 | 136 | 0.619 | 0.705 | 0.707 | 0.510 | **0.753** | 0.514 |
| VDR | 2165 | 165 | 0.731 | 0.773 | **0.815** | 0.618 | 0.809 | 0.382 |

**Plain least squares collapses exactly as the paper says it does** — 0.510 on PKM2, 0.507
on MTORC1, barely above chance where both RMT variants and the ridged control reach
0.62–0.71. *The paper's diagnosis of the problem is correct.* The question is only whether
its specific remedy beats the obvious one.

## The four verdicts

Bootstrapped over **targets**, not compounds: a method that wins on one huge target and
loses on eight small ones has not won, and resampling compounds would let ALDH1's 3,343
carry the result.

| comparison | mean ΔAUC | 95% CI | wins | Wilcoxon | verdict |
|---|---|---|---|---|---|
| cleaning (**clip**) vs uncleaned | −0.071 | [−0.096, −0.046] | 1/9 | 0.008 | loses |
| **cleaning (project) vs uncleaned** | **−0.011** | **[−0.028, +0.006]** | 2/9 | 0.300 | **TIE** |
| RMT (clip) vs forest | −0.096 | [−0.129, −0.062] | 1/9 | 0.009 | loses |
| **RMT (project) vs forest** | **−0.036** | **[−0.053, −0.014]** | 1/9 | 0.039 | **loses** |

## What this licenses us to say

**The paper is not refuted.** In its better mode the cleaning is a **tie** against a
cross-validated ridge on the raw matrix. An effect up to 0.028 of AUC, either way, is
still inside that interval — a statement about what nine targets can resolve, not a
measurement of no effect. By section 16's taxonomy this is **kind B**, the same verdict
Chemprop got, for the same reason.

**Section 8's decision survives relocation.** Moved onto the task the authors validated,
in the regime they addressed, against experimentally measured inactives, RMT still loses
to the random forest: −0.036, CI [−0.053, −0.014], one win in nine. The hole PAPERS.md
flagged is closed, and it closed in favour of the original call — not because we had
graded it on the wrong task, but on its own task.

So the tally entry splits in two, and neither half is the one we had:

| claim | kind | verdict |
|---|---|---|
| RMT cleaning beats no cleaning | **B** | tie — CI [−0.028, +0.006] |
| RMT beats our forest, as a screen | **A** | loses — CI [−0.053, −0.014] |

## My verdict function got this wrong too — the fifth time

I pre-specified `clip` as the headline arm, on a theoretical argument: a screen is a
discriminant, clipping inverts the cleaned matrix literally, therefore clipping is the
discriminant solve. The argument is not obviously wrong. It is simply false here:

> **`project` beat `clip` on 9 of 9 targets**, mean +0.060, CI [+0.031, +0.094].

And the choice decides the verdict. Against the control, `clip` gives **−0.071 (loses)**
and `project` gives **−0.011 (a tie)**. Had I reported only my pre-specified arm — the
orthodox thing to do — I would have published a rejection of the paper's central claim on
the strength of a modelling choice I had reasoned my way into and never checked. **In the
direction that flatters our own earlier decision**, which is the direction these errors
have gone four times out of five.

The fix is not "pre-specify better". It is that **the pre-specification must be reported
alongside the result, not instead of it**. `rmt_screen_test.py` now carries both modes
against both rivals and prints all four, labelled *pre-specified* and *better*. A reader
can see that the prediction failed and judge the result anyway.

## Where this test is weak, stated plainly

- **Two targets may be handicapped.** On ESR1_ant and KAT2A the control's ridge ran to the
  top of the search grid (1000). There the control may be limited by the grid rather than
  genuinely beaten. The script warns by name rather than leaving it to be noticed.
- **ESR1_ant is noise.** Every arm sits at 0.49–0.56 and heavy-atom count "wins" it. It is
  kept in the tally because dropping the targets where nothing works is how a method gets
  a reputation it has not earned.
- **The top-1% counts are inflated.** The inactive cap keeps every active but samples the
  inactives, so ALDH1's hit rate reaches 40% and a "top 1%" there is barely a screen. AUC
  is a rank statistic over active/inactive pairs and is unaffected by the class ratio,
  which is why every verdict above is read off AUC and the counts are context only.
- **Nine targets is few.** The tie interval is ±0.028 wide. A real cleaning effect of
  0.02 AUC would not have been detected here.
- **Cost is not the argument.** All arms fit in 1.3–1.4 s. Unlike Chemprop, RMT was not
  kept out on compute; it was kept out on the measurement.

---

# 18. ADME: the first endpoint that is not potency

Every number above this line is a potency number, and [§1](#1-can-the-model-predict-real-potency)
already says potency is not what kills programmes: *"No selectivity, no ADMET, no tox — the
objectives a real programme dies on."* This is the first test of the platform's model on an
endpoint a real project would actually optimise against.

```bash
python3 adme_test.py both
```

Two public sets, through the platform's own `load_csv`, `PropertyModel`, `random_split` and
`scaffold_split`:

| set | endpoint | rows | usable | source |
|---|---|---|---|---|
| `logd` | logD 7.4 | 4,200 | 4,198 | MoleculeNet Lipophilicity, curated from ChEMBL |
| `solubility` | log solubility (mol/L) | 1,128 | 884 | MoleculeNet ESOL / Delaney |

## The trap this test had to be built around

`dmta.chem.feature_vector` is a Morgan fingerprint concatenated with twelve physchem
descriptors, and the second of those is **Crippen cLogP**. For potency that is one weak
feature among 2,060. For these two endpoints it is nearly the answer: logD *is* lipophilicity
measured, cLogP is lipophilicity computed, and solubility is largely the same quantity with
the sign reversed.

So "our model beats cLogP" would be a claim about a model that was handed cLogP as an input.
The comparison is only worth anything with a fourth arm: **the identical forest with the
physchem block stripped**, which separates having learned chemistry from having read the
baseline off its own feature vector. A test asserts that the strip really removes cLogP,
because if it did not, the control would silently contain the thing it exists to exclude and
the verdict would look fine while measuring nothing.

There is **no time split**: neither public set carries dates. That matters more here than
anywhere else in this document, because the time split is the column where this platform's
σ collapsed (0.44 → 0.36 → **0.12**, inverting to −0.13 on μ-opioid) and it is the column a
production model lives in. §18 covers the two easier columns only.

## logD: a real win, and the control confirms it

Spearman ρ, mean over 5 split seeds:

| arm | random | scaffold |
|---|---|---|
| **forest** (the platform's model) | **0.799** | **0.757** |
| forest, fingerprint only | 0.724 | 0.703 |
| ridge on 12 physchem descriptors | 0.528 | 0.380 |
| cLogP alone | 0.400 | 0.219 |
| heavy-atom count | 0.169 | 0.090 |
| constant predictor | — | — |

On the scaffold split (seed 0, 3,353 train / 845 test), paired bootstrap over test compounds:

| comparison | Δρ | 95% CI | verdict |
|---|---|---|---|
| forest − cLogP | +0.534 | [+0.469, +0.597] | **wins** |
| forest − descriptor ridge | +0.386 | [+0.333, +0.441] | **wins** |
| forest − heavy atoms | +0.641 | [+0.577, +0.706] | **wins** |
| **fingerprint-only forest − cLogP** | **+0.478** | **[+0.407, +0.549]** | **wins** |
| forest − fingerprint-only forest | +0.056 | [+0.039, +0.072] | the physchem block is worth this much |

**The control is the result.** Strip every descriptor, including cLogP, and the forest still
beats cLogP by +0.478 on held-out scaffolds. The whole physchem block is worth +0.056. The
model is not reading the answer off its own input; the fingerprint genuinely carries
lipophilicity information. MAE 0.654 against cLogP's 0.952, R² 0.483 against 0.020.

σ ranks error at +0.324 on the scaffold split — comparable to the potency targets, and usable.

## Solubility: a tie with a twelve-descriptor linear model, and a loss on error

| arm | random | scaffold |
|---|---|---|
| **forest** | **0.922** | **0.917** |
| forest, fingerprint only | 0.785 | 0.618 |
| ridge on 12 physchem descriptors | 0.877 | **0.914** |
| cLogP alone | 0.830 | 0.892 |
| heavy-atom count | 0.529 | 0.519 |

Scaffold split, seed 0 (707 train / 177 test):

| comparison | Δρ | 95% CI | verdict |
|---|---|---|---|
| forest − cLogP | +0.024 | [−0.008, +0.058] | **tie** |
| forest − descriptor ridge | **−0.001** | **[−0.033, +0.032]** | **tie** |
| forest − heavy atoms | +0.413 | [+0.319, +0.514] | wins |
| fingerprint-only forest − cLogP | **−0.282** | [−0.385, −0.190] | **loses** |

And on absolute error the ranking tie becomes an outright loss:

| arm | MAE | R² |
|---|---|---|
| forest | 1.053 | 0.612 |
| **ridge on 12 descriptors** | **0.731** | **0.821** |
| cLogP alone | 0.828 | 0.777 |
| forest, fingerprint only | 1.504 | 0.238 |

**A twelve-number linear model beats the full platform on this endpoint** — better MAE,
better R², indistinguishable ρ — and plain cLogP beats the forest on MAE too. The 2,048
fingerprint bits are doing the work of the 1968-vintage Yalkowsky relationship, worse.

### Why, and it is not a tuning problem

Solubility spans −11.6 to +1.14 log units. A forest cannot predict outside the range of its
training leaves, so held-out scaffolds in the tails get squashed toward the mean; a linear
model extrapolates. Ranking survives this (ρ ties) and absolute error does not (MAE loses by
0.32 log units). This is a structural property of tree ensembles on a wide-range regression
endpoint, not a hyperparameter we failed to find.

σ ranks error at only **+0.157** on the scaffold split — weak or uninformative. Every
acquisition function except `greedy` spends part of its batch on σ, so on this endpoint the
loop's exploration would be spending budget on noise.

## The bug this test found in our own metric

`dmta.models._rank` assigned ranks in `argsort` order with **no tie handling at all**. For
ordinary predictions that is a rounding-level error. For a *constant* predictor it is not:
the arbitrary tie order gives the constant arm real rank variance, the "is this constant"
guard never fires, and Spearman comes back as a small random number instead of undefined.

It showed itself as an impossibility in the first run — a point estimate outside its own
confidence interval:

```
  forest vs mean
    d-rho +0.665, 95% CI [+0.666, +0.812]
```

`_rank` now returns mid-ranks and agrees with `scipy.stats.spearmanr` exactly, ties included;
`regression_metrics` returns `None` for a constant prediction; and the constant arm is scored
on MAE only, where a floor actually means something.

**No published number moved.** Re-running `chembl-validate` on EGFR after the fix reproduces
σ~err of 0.437 / 0.464 / 0.321 — the 0.44 / 0.46 / 0.32 reported in §1 — because a 300-tree
forest on distinct compounds produces almost no exact ties. The defect was real, it was in a
metric the split evaluation publishes, and it was only reachable from a benchmark that
carried a degenerate floor. That is an argument for always carrying one.

## What this licenses us to say

- **logD: the platform's model works on an ADME endpoint**, and the claim survives the
  control that strips the trivial baseline out of its own input. This is the first
  non-potency capability to pass.
- **Solubility: it does not beat a twelve-descriptor linear model**, ties it on ranking and
  loses to it on error. We should ship the ridge for this endpoint and say so.
- The two answers disagree, which is the useful part: "the model does ADME" is not a claim
  this platform can make. It does *logD*. Each endpoint needs its own examination, and the
  cheap linear baseline has to be in every one of them.

## Where this test is weak, stated plainly

- **No time split, on either set.** The column that broke potency prediction is untested
  here, and there is no reason to expect ADME to be kinder.
- **Two endpoints, both public, both on MoleculeNet.** No clearance, no permeability, no
  hERG, no plasma protein binding. logD and solubility are the two most physics-like ADME
  endpoints and plausibly the two most favourable.
- **884 of 1,128 solubility rows are usable.** `standardize` drops anything under 8 heavy
  atoms, which removes the small solutes. That is what the platform really does to an
  uploaded CSV, and it makes the task harder rather than easier — but a fifth of the file
  is gone and the reader should know.
- **The bootstrap resamples compounds, not targets.** Unlike §16 and §17, there is one
  dataset per endpoint, so the interval describes uncertainty *within* that endpoint and
  says nothing about the next one. The 5-seed spread is reported alongside for that reason.
- **The solubility test set is 177 compounds.** The ±0.033 tie interval is as tight as that
  allows, and a real effect of 0.03 ρ would not have been detected.

---

# 19. The property profile rejects drugs that already exist

Everything before this section judged a compound's properties with one profile,
`DEFAULT_PROFILE`, and presented it as *the* property profile. It is not neutral. It
encodes **oral small-molecule** drug-likeness, and it was applied to every campaign
without anyone — me included — asking what happens when the programme is not an oral
small-molecule programme.

[COMPETITORS.md](COMPETITORS.md) turned up a company for which that is not hypothetical:
Bugworks, developing a broad-spectrum antibacterial now in Phase 1. Gram-negative
antibacterials are famously not drug-like by the oral criteria.

```bash
python3 f9_profile_test.py
```

## The test: marketed drugs are a known answer

A drug that is approved and sold is, by definition, good enough. If our profile marks
down drugs that already exist, the profile is wrong and **the number of them is the size
of the error**.

| arm | what it is | n |
|---|---|---|
| antibacterial | ATC **J01**, approved, small molecules | 138 |
| control | ATC **A, C, L, N**, approved, not also J01 | 1,032 |

Both arms are marketed drugs. Neither is a set of failures. A profile that genuinely
measures "drug-likeness" should not separate them.

Only the **descriptor** objectives are scored. Potency, route confidence, cost and lead
time are unknown for a marketed drug, and `Objective.value` returns 0.5 for a missing
value — including them would add the same constant to both arms and halve every
difference for no reason. The rejection line is pre-specified as the control arm's 25th
percentile.

## Result: the profile we ship marks down real antibiotics

| | antibacterials | control | difference in rejection rate |
|---|---|---|---|
| whole class | 47.1% rejected, median 0.708 | 25.0%, median 0.855 | **+22.1%**, CI [+13.3%, +30.9%] |
| **≤600 Da only** | **40.7%** rejected, median 0.770 | 21.0%, median 0.865 | **+19.6%**, CI [+10.8%, +28.7%] |

The second row is the one that matters. Antibacterials include glycopeptides of 1,400+ Da,
and a small-molecule profile marking *those* down is arithmetic, not a finding. Below
600 Da, where molecular weight cannot be the explanation, the bias is almost unchanged.

The worst-scoring marketed antibacterials under our own profile:

| score | drug | MW | cLogP | TPSA | QED |
|---|---|---|---|---|---|
| **0.000** | daptomycin | 1620.7 | −5.62 | 702.0 | 0.01 |
| 0.051 | plazomicin | 592.7 | −5.36 | 269.3 | 0.09 |
| 0.105 | amikacin | 585.6 | −8.42 | 331.9 | 0.11 |
| 0.113 | streptomycin | 581.6 | −8.16 | 336.4 | 0.07 |
| 0.153 | vancomycin | 1449.3 | 0.11 | 530.5 | 0.07 |
| 0.166 | cefiderocol | 752.2 | −0.18 | 256.9 | 0.05 |

Daptomycin scores **exactly zero**. Three of these are on the WHO Essential Medicines
list. A loop running this profile would never propose anything resembling them.

## The repair, and the two ways I got it wrong first

### Attempt 1: widen the windows by hand — caught by a pre-written check

I widened cLogP, widened TPSA, raised the HBD ceiling and dropped QED. It "worked": the
gap inverted to −13.4%. It also lifted the **control** arm's median from 0.855 to 0.906
and scored the median antibacterial exactly **1.000**.

That is not a class-specific profile. It is a lower bar wearing one. The degenerate-fix
check existed in the script before the first run precisely because this is the obvious
way to fake the result, and it caught it.

### Attempt 2: fit the windows to the class — and pick the width on the wrong measurement

Fitting the windows to the class's own percentile spread is defensible, but circular: a
profile fitted to a set will score that set well. So the fit uses **half** the marketed
antibacterials and every number is reported on the **held-out half**.

The width then had to be chosen. I swept it and required three things: don't over-reject
the class, don't lift the control arm, and **don't saturate** — a profile that gives a
whole class full marks cannot order it, and ordering is the only thing a profile does.

| window | rejected, abx | rejected, control | difference | below 0.5 | control median |
|---|---|---|---|---|---|
| **10–90** | 10.1% | 25.0% | −14.9% | **7.2%** | 0.805 |
| 20–80 | 15.9% | 25.0% | −9.1% | 13.0% | 0.660 |
| 25–75 | 18.8% | 21.2% | −2.4% | 18.8% | 0.600 |
| 30–70 | 26.1% | 21.5% | +4.6% | 26.1% | 0.540 |
| 35–65 | 30.4% | 25.0% | +5.4% | 31.9% | 0.500 |

I first chose **25–75**, because at 10–90 my saturation column read 68.8% and failed.

Then a test I had written failed on **ciprofloxacin** — one of the most-prescribed
antibacterials in the world — which scored 0.836 under that window, because its TPSA of
74.6 sits below the fitted lower edge. Fitting to the middle 50% of a class excludes the
edges of that same class.

## The sixth verdict-function error, and this one had a victim

**I measured saturation on the property terms alone.** But the property terms are not
what ranks compounds — **potency is**, at weight 3.0 against the guards' 5.0. With the
potency term restored and a realistic spread of activity, the identical profile at 10–90
saturates at **7.2%**, not 68.8%.

The number was not wrong. It was a measurement of the wrong quantity, and it pointed the
decision in a specific direction: toward a narrower profile that rejects **18.8%** of
marketed antibacterials instead of 7.2%, ciprofloxacin among them.

This is the sixth time in this project that my own checking logic, rather than a model,
produced the error — and the first time the mistake would have shipped as a *product
default* rather than as a sentence in a document. The pattern across all six is the same:
**the check measured something adjacent to the thing it claimed to measure.** The property
score is adjacent to the full score. A tie is adjacent to a rejection. A plateau is
adjacent to a trend.

What saved it was a unit test written from a real drug rather than from the fitted
distribution. The distribution cannot tell you that ciprofloxacin matters; only chemistry
can.

## What ships

`PROFILES["antibacterial"]`, fitted at the 10th–90th percentile of the marketed class,
with:

- **QED dropped, not down-weighted.** QED is fitted to oral drugs. Scoring a class that is
  deliberately not oral-drug-like by how much it resembles oral drugs assumes the
  conclusion.
- **Rotatable bonds tightened to ≤5**, the eNTRy accumulation guideline for entry into
  Gram-negative cells — the one term taken from the literature rather than the data, and
  the only place this profile is *stricter* than the default.
- Held-out antibacterials ≤600 Da: rejected **39.3% → 4.9%**.

### What this profile is not

It over-corrects. On the held-out half it now rejects antibacterials **14.9% less** often
than control drugs, CI [−21.8%, −6.6%], and that interval excludes zero. Parity would be
zero. **This is a profile for antibacterial programmes, not a better universal profile** —
which is the honest conclusion, and the one the architecture already supported: LANDSCAPE.md
records "the profile is data" as a design decision, and this is the first time that
decision has been load-bearing rather than decorative.

### The eNTRy criteria, as context rather than verdict

| arm | primary amine | ionizable N | rotB ≤ 5 | median spherocity |
|---|---|---|---|---|
| antibacterial (≤600 Da) | 17.1% | 43.1% | 66.7% | 0.204 |
| control (≤600 Da) | 6.0% | 29.5% | 65.5% | 0.218 |

The nitrogen criteria separate the arms; rotatable bonds do not. Globularity is the
published third criterion and RDKit does not ship the authors' definition — spherocity is
reported as the nearest available measure, and **their 0.25 threshold is deliberately not
applied to it**, because a cut-off defined on one measure does not transfer to another.

## What this cost and what it bought

Two ChEMBL queries, a few minutes of CPU, no GPU. It found a defect that had been shipping
in every campaign since the profile was written, and it is the first finding in this
project that came from a named prospective customer rather than from a paper.

**It does not make the loop work prospectively.** That remains unaddressed, and a
discovery organisation of Bugworks' size will ask about it sooner than a large pharma
would, because it has fewer compounds to waste.

---

# 20. The molecules never leave the machine, and now there is a test

Every number in this document is about whether a prediction is any good. This section is
about a different kind of claim, and it is the first one here that was true before anyone
checked it.

A customer's compound structures are the asset. Most AI discovery tools send them to
somebody else's computer. Ours does not — not by design, but because a random forest on
CPU over a local SQLite store never needed to. Palo Alto Networks paid **$4.5B** in
November 2026 for a company whose central pitch is substantially that property, so the
demand for it is not speculative. We had it and had never said so.

The problem with saying so is that "your data never leaves" is unfalsifiable from the
outside. This section makes it falsifiable.

```bash
python3 egress_test.py
```

## The result

```
Static scan of dmta/ for anything that can reach a network
  chembl.py      opt-in fetcher    import urllib.request; Request() at line 60
  structure.py   opt-in fetcher    import urllib.request; Request() at line 65

  2 of 26 modules touch the network, 0 of them unexpectedly.

Controls
  guard raises on a remote address:      yes
  loopback and AF_UNIX still permitted:  yes

The cycle itself, with every non-loopback address blocked
  3 cycles completed in 6.8s
  compounds made and tested: 36 of 924
  remote hosts contacted: 0
```

**The full Design–Make–Test cycle — enumerate, score, select, plan, measure, retrain —
runs to completion with every non-loopback connection blocked.** Two of twenty-six
modules can reach a network at all, and both are explicit one-off fetch steps the user
invokes by name (`chembl-fetch`, `structure`).

## Why a passing run is evidence rather than decoration

A test that cuts the network and passes proves nothing on its own: it passes identically
if the guard is installed, if the guard is a no-op, and if the cycle silently did no
work. So the run carries a **positive control** — a deliberate connection to a remote
address, which must raise. If it does not, the script prints that the result is void and
exits non-zero before reporting anything else.

The second control is the opposite error. A guard that blocked *all* sockets would sound
stronger and be useless: it would block loopback, so the customer could not run the UI,
and `AF_UNIX`, so scikit-learn's own worker processes (`n_jobs=-1`) would not start. It
would then fail for reasons with nothing to do with data leaving, and whoever hit that
would switch it off rather than debug it. So the line is drawn at **any connection to a
non-loopback address over AF_INET or AF_INET6**, which is exactly the set capable of
carrying a molecule off the machine.

## My first control was wrong, and it failed in the direction that looks like success

The first version of the local-traffic control asserted on `is_local_address("/tmp/x.sock")`
and reported a failure. The guard was fine. **The control was exercising a code path that
does not run** — for an `AF_UNIX` socket the guard checks the socket family first and
returns, so `is_local_address` is never called with a path at all.

It is worth recording because of which way it broke. It reported `CLAIM NOT SUPPORTED`
against a claim that was true, which is the harmless direction. The same mistake made one
step later — asserting on a helper the guard *did* call, but with an input the real code
never produces — would have reported support for a claim nobody had tested. A control has
to use the thing the customer uses, or it is testing a helper nobody called.

## The test that matters is the static one

The run-time result is a fact about the code as it stands today. One `import requests` in
`select.py` would quietly end it, and the run-time test would keep passing as long as
that import was not reached on the happy path.

So `test_only_the_two_opt_in_fetchers_touch_the_network` scans every module in `dmta/`
for network imports and calls, and fails if the set is anything other than those two. It
is crude and over-sensitive on purpose — a false positive costs someone thirty seconds of
reading, a false negative costs the claim. I verified it is not vacuous by adding
`import requests` to `select.py`; it failed and named the module.

## What this does not show, stated so nobody has to find it

- **One process, run time only.** Anything a dependency does at *import* time happens
  before the guard is installed and is invisible here.
- **Python's `socket` module only.** A C extension opening a raw socket would not be seen.
- **Nothing about disk.** A molecule written to a file that some other process uploads has
  still left, and no test here would notice.
- **The two fetchers really do reach the network** when invoked. The claim is about the
  cycle, not about the whole program.
- **It is not a security audit.** It is a check that one specific, checkable thing is
  true, which is the only kind of claim this document makes.

## What it licenses us to say

> The Design–Make–Test cycle completes with every non-loopback connection blocked, in a
> run that also proves the block works. Two modules of twenty-six can reach a network,
> both are opt-in fetches of public data, and a test fails if a third appears.

Not "we are secure". Not "your IP is safe". Those are claims about a system, and we have
tested a process. **100 tests pass.**

This cost an afternoon, needed no new science, and is the only item on the feature list
in [FEATURES.md](FEATURES.md) that was already true before it was built — which is also
the reason it had never been claimed. It is worth asking what else is in that category.

---

# 22. What an answer cost, reported next to the answer

*Section 21 is reserved for the Cell Painting work running in a sibling session.*

Every number in this file has a price, and until now that price appeared only in prose:
"a few hours of CPU", "96× less compute to fit", "two ChEMBL queries and a few minutes".
Those are the strongest claims this project has for a customer who cannot afford a physics
stack — and they were footnotes, written by hand, and never checked.

[COMPETITORS.md](COMPETITORS.md) made the case that in the market we can actually reach,
cost per outcome is the differentiator being funded rather than an apology for a simple
model. This section makes it a measured, reported number.

```bash
python3 -m dmta.cli evaluate --csv data/CHEMBL203-IC50.csv --cost-rate a100_80
```

## The result, on the real EGFR dataset

| step | wall seconds | CPU seconds |
|---|---|---|
| three splits, model fitted and scored on each | 22.4 | 65.3 |
| retrospective loop replay | 33.4 | 121.7 |
| **total** | **55.8** | **187.0** |

One verdict on 8,978 real compounds — three splits, a calibration check and an
eleven-cycle loop replay — for **$0.01–0.02** of rented machine, or **nothing** on a
laptop that was already bought.

Against the published cost of making and testing a single compound, **$565–$1,000**, the
entire analysis costs **under one thousandth of one compound**. That is the sentence worth
carrying into a conversation, and it is now generated rather than asserted.

## Four ways a cost number lies, and what stops each one here

This was a short build. Almost all of it is the guards, because a cost figure is trivial
to produce and trivial to make dishonest, and every failure below still renders a
plausible-looking page.

### 1. Time is measured. Money is assumed. They must not look alike.

A wall-clock second is a fact: a clock was read twice. A dollar is that fact multiplied by
a rental rate this process never observed. So `measured` always carries seconds; `assumed`
carries the rate **and its source**; `estimate` exists only when a rate was supplied. The
HTML prints the rate and its provenance inline — *"assuming $1.19/hour (A100 80GB PCIe);
that rate is an assumption, not a measurement: Runpod catalogue read, ROADMAP.md,
2026-11-24"* — so the page cannot be read as having observed a price.

### 2. No rate must mean no money, never a silent zero

**A zero that means "unknown" is the most expensive wrong number available**, because it
reads as *free* and gets repeated by whoever saw it. With no rate configured the report
says so in words and prints no figure at all.

The `local` rate is a genuine zero and a different thing: hardware already owned has a
marginal cash cost of nil. Its label says exactly that, so the two zeroes cannot be
confused.

### 3. The flattering denominator must not become the headline

Cost per compound scored is a fraction of a cent, always, for any platform, because
scoring is vectorised and libraries are large. It is true and it decides nothing.

**The headline denominator is the verdict: one run answers one question, and that is what
is being bought.** Per-compound is still printed — hiding it would be its own dishonesty —
with a line saying why it is not the number that matters. A test asserts the per-compound
figure is the cheaper of the two, so it cannot quietly be promoted.

### 4. A GPU rate on a run that never used a GPU

The platform's own finding is that a CPU random forest ties Chemprop at a fraction of the
compute. Pricing that run at an A100 rate overstates the bill — **in the direction that
makes our own cost argument look worse**, which is the safe direction, but a reader should
not have to work it out. The report now says so where it happens.

Cost is charged on **wall-clock**, not CPU seconds, because a rented machine bills by the
clock whether every core is busy or not. CPU time is reported beside it so a parallel run
is visible as parallel: above, 187 CPU seconds inside 56 wall seconds.

## The anchor is a range, and stays a range

| source | per compound |
|---|---|
| ChemRxiv 2025-01, AI-driven DMTA system — synthesis plus assay, one project | **€565** |
| arXiv 2026, D2ML direct-to-biology — custom synthesis | **~$1,000** |

Both are **single-project self-reports, not industry benchmarks**, and the report says so
under the figure. Collapsing them into one number would manufacture a precision the
literature does not have — and that invented number would then be quoted back at us. A
test asserts the constant stays a pair with two sources, and that the equivalence in
compounds reports the modest end first.

## What this does not do

It prices **compute**, which is the cheap half. It does not price the chemist reading the
batch, the assay, or the twelve months of calendar time a real programme spends. A
platform that costs two cents to run and proposes a compound nobody makes has cost the
programme a cycle, and no ledger here will show that.

So this is an honest number about the thing it measures, and the thing it measures is not
the whole cost of a decision. The loop still has no prospective result, and compute was
never what was blocking that.

## Shipping the ridge, and how the first attempt failed

§18 said we should ship the ridge for solubility. Writing `if endpoint == "solubility"`
would have been a guess dressed as a decision, and the next endpoint would arrive with no
rule at all — so `choose_model` fits both and keeps whichever measures better, on the
training split only.

**The first version picked the forest on solubility, five seeds out of five.** It was not
broken. It cross-validated on *random* folds, and on random folds the forest really is
better: MAE 0.624 against the ridge's 0.729. The forest's advantage on this endpoint is an
artefact of random cross-validation — near-identical analogues land on both sides of the
fold, the forest interpolates between them, and the score describes a regime no production
model is ever in. §1 makes exactly this argument about reporting the random split as
generalisation. I had not noticed it applies just as hard to *choosing* a model.

With whole scaffolds held out of each fold, the order reverses and holds:

| endpoint | forest MAE | ridge MAE | picked | seeds agreeing |
|---|---|---|---|---|
| solubility | 0.871–0.894 | **0.831–0.842** | **ridge** | 5 / 5 |
| logD 7.4 | **0.660–0.676** | 0.815–0.819 | **forest** | 5 / 5 |
| EGFR potency | **0.650** | 0.977 | **forest** | — |

Two details that decide whether this works:

- **It compares MAE, not ρ.** Solubility was a *tie* on ranking and a loss on error, so a
  selector reading ρ keeps the forest and is wrong in precisely the case it exists for.
- **The in-training margin is thin.** 0.038 on solubility, against a held-out gap of 0.322.
  The direction is stable across all five seeds, but the selector is reading a much
  quieter signal than the one §18 measured, and on a smaller endpoint it could easily be
  a coin flip. It is a decision procedure, not an oracle.

### What shipping the ridge costs: its uncertainty is worse

The ridge is a linear model with no ensemble to disagree with itself, so σ comes from
refitting on bootstrap resamples. Measured the same way as everything else — Spearman of
σ against \|error\| — it is worse than the forest's on both endpoints:

| | forest | ridge |
|---|---|---|
| solubility, random | +0.273 | +0.188 |
| **solubility, scaffold** | +0.157 | **−0.026** |
| logD, scaffold | +0.324 | +0.168 |

On the endpoint we are shipping it for, the ridge's σ is **−0.026 on held-out scaffolds** —
no information, very slightly inverted. So the trade is a real one and it should be stated
in both directions: **better predictions, worse uncertainty.** A solubility campaign gets
the more accurate model and should run at β = 0 — greedy or expected improvement — because
every other acquisition function spends part of its batch on a σ that is not ranking
anything. The recommender already reads β from the worst split, which handles this, but it
is worth saying out loud rather than leaving it to be discovered.

### How it is wired

`CycleSettings.model` is `"auto"` by default, with `"forest"` and `"ridge"` to force it.
`"auto"` decides **once**, at 60 measured compounds, and refits that class each cycle
afterwards: the endpoint does not change mid-campaign, and `choose_model` is an order of
magnitude more work than one refit, so paying it every cycle would buy an answer that
cannot move. Below 60 compounds it stays on the forest and records no choice rather than
fixing the model for a whole run on a coin flip. The decision is in the cycle report.

### A guard one layer down, from another session's audit

Another session auditing a Cell Painting failure found its cause was a scaffold function
that swallowed an RDKit exception and returned `""` for every compound — so every compound
became its own group and a scaffold split silently became a split on row order, through a
full 38-protein experiment. `dmta.chembl.murcko` ends in the same shape:
`except Exception: return "unknown"`.

Everything in this section now rests on scaffold identity, so I checked rather than took
it on report: aspirin → `c1ccccc1`, triethylamine → `acyclic`, and 8,978 real EGFR
compounds give 3,294 distinct scaffolds with **zero** `"unknown"` and a largest group of
518. §1's and §18's scaffold numbers are sound as published. Two tests now pin it, and a
third pins the safe failure — if scaffolds ever do collapse to one group, `choose_model`
reports that it could not compare instead of silently returning whichever model it fitted
last.
