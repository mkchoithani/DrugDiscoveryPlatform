# PostEra's published methods, against what this platform actually does

Six papers, all with Alpha A. Lee as senior author. For each: what it does, whether we
have it, and whether it is worth having.

**Updated 2026-11-25.** Three things have changed since the first pass: paper 6 has been
implemented and measured (and rejected), paper 2 has been used to audit our own headline
result (which did not survive), and paper 3 has now been read in full — which overturned
what this file previously said about it. Corrections are marked where they sit.

| # | Paper | Venue | Read? | In the platform? |
|---|---|---|---|---|
| 1 | Learning probabilistic pharmacophore maps enables efficient hit finding | Digital Discovery 2026, 5, 5001 | full text (CC-BY) | **implemented, then rejected** on measurement |
| 2 | Understanding ML models trained on DNA-encoded libraries | Digital Discovery 2026, 5, 4976 | full text | **yes, as an auditor** — and it failed our own benchmark |
| 3 | Confusion Without Measure: units and best practices in reporting preclinical ADME/PK | J. Med. Chem. 2026, 69, 12723 | **full text** | **no** — and *not* applicable yet; see the correction |
| 4 | Structure-Guided Discovery of Selective PLK3 Inhibitors | ACS Med. Chem. Lett. 2026, 17, 1266 | metadata only (paywalled) | n/a — a campaign, not a method |
| 5 | SARS-CoV-2 Mpro inhibitors using a synthesis-directed de novo design model | Chem. Commun. 2021, 57, 5909 | full text (CC-BY) | **no** — closest in spirit to our loop |
| 6 | Ligand biological activity predicted by cleaning positive and negative chemical correlations | PNAS 2019, 116, 3373 | abstract only | **implemented, then rejected** on measurement |

---

## 1. Probabilistic pharmacophore maps (2026)

**What it does.** Learns, from known protein–ligand complexes, *where* interaction
features (donors, acceptors, hydrophobes, aromatics, ionisables) ought to sit in a
binding site, conditioned on the pocket residues. The predicted pharmacophore then
filters billions of purchasable molecules through a vector database — sublinear in
library size — and aligns candidates in the pocket. Virtual screening is reduced from
scoring every pose to identifying the right pharmacophore. Prospectively tested on
AKR1C1.

**Why it matters here.** It is a direct answer to the failure recorded in
[RESULTS.md](RESULTS.md) §4–5: docking gave us 0.0–0.37× enrichment across three
structures. Their framing is that docking is expensive and pose-dependent, and that
co-folding models carry training-set bias — so they moved to a coarser, cheaper
representation instead.

**It also corrected us.** [ROADMAP.md](ROADMAP.md) recommended Boltz-2 for affinity on
shortlists. This paper summarises evidence that Boltz-2's affinity predictions are
insensitive to pose quality and to mutations of key pocket residues, barely degrade when
pose information is removed, and failed to pair actives to the right target on a
LIT-PCBA-derived benchmark. That recommendation now carries a correction.

**Could we build it?** Partly. Inference code and model weights are released
(`github.com/postera-ai/phc4-field-generator`, weights on Zenodo); **training and
benchmarking code are proprietary and not released.** We already have the adjacent
pieces: structure preparation, binding-site definition from the crystal ligand, and
pocket-residue computation (the 3D viewer computes the 58 residues lining the 2RGP
pocket). The missing part — a model mapping pocket residues to a probabilistic feature
map — is a serious build, not a weekend.

**Implemented 2026-11-25, and not adopted.** Their inference ran on our 2RGP receptor in
four CPU-minutes. We wrote the missing screening half (`dmta/pharmacophore.py`: hot-spot
reduction, RDKit feature annotation, geometric clique matching). On 1,500 real EGFR
compounds it reached **1.24× enrichment** — above docking's 0.37×, but below the 1.5× bar
and **below ranking compounds by heavy-atom count** (1.37×), which costs nothing. The
time-split result looked like a win at 1.65× against the forest's 1.37×; that is 7 actives
against 6 out of 41, and a gap that size arises 65% of the time by chance.

**Where the test is unfair to them, and it matters:** the paper claims *hit finding* —
actives out of a large library, tested prospectively on AKR1C1. We tested *potency ranking
within a mature congeneric series*, where nearly every compound is already active. Those
are different jobs. The narrow, honest claim is that predicted pharmacophore maps do not
help rank the next batch in a programme. A fair test of their own claim would be actives
against property-matched decoys, which we do not run. Full detail, including four caveats
that belong to my implementation rather than their model, in [RESULTS.md](RESULTS.md) 10.

**Verdict: the most promising unimplemented idea turned out not to help here** — but this
is the one rejection of the four where I would not generalise beyond our own task.

---

## 2. DEL Simulator (2026)

**What it does.** An open-source, structure-aware simulator of DNA-encoded library
experiments — synthesis, selection, sequencing — built so that compounds have a *defined
structure–affinity mapping*. That gives computational ground truth against which ML
models trained on DEL data can be benchmarked honestly. They use it to show when
disynthon aggregation helps, and report the counterintuitive result that model
performance can *fall* as library size grows.

**Do we have it?** No — we use no DEL data. **But we built the same instrument for the
same reason.** `dmta/assay.py`'s `SimulatedAssay` is a hidden structure–activity
landscape whose only purpose is to benchmark selection strategies against known truth,
and it caught our UCB acquisition performing worse than random. Their argument — that
you cannot evaluate a model on noisy, biased real data without a ground-truth harness —
is the argument this entire repository is built on.

Their counterintuitive finding has a sibling in ours: exploration (β) *hurt*
monotonically at small budgets, and our 100-compound qualification rule evaporated on
held-out data.

**Implemented 2026-11-25, as an auditor rather than a feature.** `dmta/del_landscape.py`
loads their library and affinities as platform compounds; `del_benchmark.py` runs our own
acquisition comparison in their landscape and ours under one protocol. We use none of
their DEL-specific machinery (yield, selection rounds, PCR, NGS readout) — only the
ground truth, because the thing under test is our selection loop.

**It found what it was built to find, and the news was bad.** Our 2.75× enrichment
becomes **0.71×** in their landscape — worse than random screening. The mechanism is not
that UCB is a poor acquisition function: their landscape is essentially unlearnable from
Morgan fingerprints (ρ ≈ 0 with 2,000 training compounds), and a loop cannot beat random
selection with a blind model. What it really exposed is that **our own simulator is far
too easy** — ρ 0.89 from 100 compounds, against a median of 0.56 across 24 real
programmes. Full result, including an honest weakness in the audit's own control arm, in
[RESULTS.md](RESULTS.md) 9.

**Verdict: the single most valuable paper of the six for us**, and it never became a
feature. It falsified a number we were publishing.

---

## 3. Units and reporting in preclinical ADME/PK (2026)

**Read in full on 2026-11-25.** This section previously called it "the most immediately
applicable of the six" and said it was "the one whose subject we are already
reimplementing by trial and error". **Both claims were wrong, and they were made from the
title.** The correction is worth more than the original guess.

**What it actually is.** An *Editorial* on reporting conventions for preclinical **ADME
and PK** data: in vivo pharmacokinetics (dose, clearance, volume of distribution) and in
vitro ADME (microsomal and hepatocyte stability, protein binding), plus preclinical
efficacy. It ends with per-assay reporting checklists (Appendix 2), a table of common unit
conversions (Appendix 3), and a worked set of confusions — grains versus grams, pounds
entered as kilograms, bound versus unbound concentrations.

**What it is not about**, all of which this file previously implied:

- **Not about potency data.** IC50/Ki curation is barely touched; "potency" appears mainly
  as a unit list under preclinical efficacy.
- **Not about censored values.** It says nothing about `> X` or `< X` readings, assay
  limits, or how to use them. Our Open Source Malaria problem — 51 of 123 values censored
  — is simply not this paper's subject.
- **Not about replicate disagreement.** No guidance on aggregating replicates that differ.
- **Not about machine learning.** Standardised data sets for methodology development are
  the stated motivation, but no modelling content follows.

**So it does not apply to `dmta/ingest.py` as it stands**, because the platform models
potency and does not measure a single ADME or PK endpoint. It becomes relevant the day we
add one, and then it is a genuinely good specification to build against.

**Two things in it are useful now, both small.** First, its insistence that **salt and
solid form of the batch be recorded** supports keeping per-value provenance rather than
silently normalising, which our ingest already does when it strips salts. Second — and
this one lands squarely on real data we have — its advice to *avoid reporting only a
percentage at a single time point*, because such numbers "depend greatly on assay
conditions, and are hard to compare across assay setups". That is precisely the
`60% @ 120000` and `25-80% @ 80000` readings in the Open Source Malaria file. Our choice to
treat those as unusable rather than convert them is the right one, and now it has a
citation.

**The censored-value work stands on its own.** It was the one concrete change this paper
was supposed to justify; it does not, and it should be done anyway, because 51 of 123 OSM
compounds carry information we are currently discarding. See
[IMPLEMENTING_PAPERS.md](IMPLEMENTING_PAPERS.md).

---

## 4. Selective PLK3 inhibitors, with Pfizer (2026)

Paywalled; metadata only. Author list spans PostEra and Pfizer Research and Development
— an output of their partnership rather than a method paper. Potentially interesting as
*data* if the supporting information contains compounds with potencies, which would be
another real programme for `evaluate`.

---

## 5. Synthesis-directed de novo design (2021)

**What it does.** The closest published work to our loop. Starting from **42 compounds
with measurable IC50 and 515 inactives**, a model designed 5 new compounds predicted to
beat the incumbent, *together with predicted synthesis routes*. All 5 were made; **3 had
measurable activity**, and the best matched the strongest compound in the training set
while sitting on a different scaffold and showing activity in a live-virus assay.

Two design choices match ours closely:

- **Synthesis is part of the design**, not a filter afterwards — our library is
  synthesisable by construction for the same reason.
- They explicitly reject thresholding into active/inactive, because it destroys the
  ranking information needed to tell whether a compound beats the incumbent. Our loop
  ranks on a continuous potency scale for that reason, and the dataset-relative
  threshold fix in §4 came from the same realisation.

Where they go further: **de novo generation**. We enumerate a fixed combinatorial
library; they generate novel structures and predict routes for them. That is the natural
upgrade once enumeration stops covering the space of interest.

This is also the COVID Moonshot — the dataset we already replay in
[RESULTS.md](RESULTS.md) §3. Training set and generated molecules are open at
`github.com/wjm41/mpro-rank-gen`.

**Verdict: the clearest blueprint for our next capability**, and prospectively validated
(5 designed, 5 made, 3 active) rather than retrospectively.

---

## 6. Cleaning positive and negative chemical correlations (PNAS 2019)

**What it does.** Uses random matrix theory to separate real chemical signal from
undersampling noise, in the regime where known actives and inactives are vastly
outnumbered by possible chemical features. Prospectively found four confirmed M1
muscarinic agonists. The authors emphasise that the method is **interpretable** and
yielded chemical insight about binding modes.

**Relevance.** This is a direct alternative to `dmta/models.py`. We use a random forest
on 2048-bit Morgan fingerprints with 60–900 training compounds — exactly the
"features ≫ molecules" regime the paper addresses. And we have just shown
([RESULTS.md](RESULTS.md) §7) that a graph neural network does **not** beat that forest:
median ρ 0.510 vs 0.425 over 18 programmes. An RMT-denoised linear model is a different
bet from a bigger network, and it comes with interpretability the forest lacks.

**Implemented 2026-11-25 in `dmta/rmt_model.py`, and not adopted.** Over the same 24
programmes: forest median ρ **0.562**, best RMT variant **0.389**, worse by 0.098 (just
outside twice its standard error) and 27× slower to fit. It does win outright on 5 of 24.

**Verdict: the better-motivated next model was still worse.** Two things are worth keeping
from the exercise. The method is ambiguous in a way the paper does not resolve — inverting
the cleaned correlation matrix *whitens*, damping the signal directions RMT just found,
and solving inside the signal subspace instead is much better (ρ 0.389 vs 0.258). And the
first implementation was broken in a way that looked exactly like a negative result:
see [RESULTS.md](RESULTS.md) 8.

---

## What this review changed

1. **Boltz-2 recommendation corrected** in ROADMAP.md, citing paper 1.
2. **Pharmacophore screening replaces "try another docking structure"** as the route to
   structure-based signal, since docking has already failed three times here.
3. **Paper 6 was implemented and rejected.** It was the better-motivated model than
   Chemprop, and it still lost to the random forest. The model is not the bottleneck.
4. **Paper 2 was implemented as an auditor and falsified our own headline number.** The
   2.75× does not survive an independently authored landscape, and the reason is that our
   simulator is too easy. This is the most useful thing in the whole review.
5. **Paper 3 was read, and this file's claim about it was wrong.** It concerns ADME/PK
   reporting, a data domain the platform does not yet touch at all. It blocks nothing.
6. **Paper 5 is the template** for the de novo capability, if enumeration is ever
   outgrown — and it is now the most promising unimplemented item alongside paper 1.

Four of the six have now been acted on. **Paper 5 (de novo design) is the only one left
that would add a capability**; paper 4 remains data rather than method.

**Four attempts to beat a random forest on 2D fingerprints have now failed** — docking,
Chemprop, the RMT model and pharmacophore maps. Three of those are more sophisticated by
any reasonable measure. At this data scale the remaining value is not in the model or in
the 3D representation; it is in what the loop is allowed to propose, which is what paper 5
is about.

**The pattern across all three implementations is the same.** Each one was adopted or
rejected on a measurement, not on what its authors reported — and in two of the three the
first result was wrong for a reason that had nothing to do with the paper.
