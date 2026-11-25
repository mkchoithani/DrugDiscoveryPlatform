# PostEra's published methods, against what this platform actually does

Six papers, all with Alpha A. Lee as senior author. For each: what it does, whether we
have it, and whether it is worth having. The short answer is that **none of the six is
implemented here** — but two of them change decisions already recorded elsewhere in
this repository, and one of them we independently reinvented.

| # | Paper | Venue | Read? | In the platform? |
|---|---|---|---|---|
| 1 | Learning probabilistic pharmacophore maps enables efficient hit finding | Digital Discovery 2026, 5, 5001 | full text (CC-BY) | **no** |
| 2 | Understanding ML models trained on DNA-encoded libraries | Digital Discovery 2026, 5, 4976 | full text | **no** — but same philosophy, arrived at independently |
| 3 | Confusion Without Measure: units and best practices in reporting preclinical ADME/PK | J. Med. Chem. 2026, 69, 12723 | metadata only (paywalled) | **no** — most directly applicable of the six |
| 4 | Structure-Guided Discovery of Selective PLK3 Inhibitors | ACS Med. Chem. Lett. 2026, 17, 1266 | metadata only (paywalled) | n/a — a campaign, not a method |
| 5 | SARS-CoV-2 Mpro inhibitors using a synthesis-directed de novo design model | Chem. Commun. 2021, 57, 5909 | full text (CC-BY) | **no** — closest in spirit to our loop |
| 6 | Ligand biological activity predicted by cleaning positive and negative chemical correlations | PNAS 2019, 116, 3373 | abstract only | **no** — a direct alternative to our model |

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
(`github.com/postera-ai/phc4-eld-generator`, weights on Zenodo); **training and
benchmarking code are proprietary and not released.** We already have the adjacent
pieces: structure preparation, binding-site definition from the crystal ligand, and
pocket-residue computation (the 3D viewer computes the 58 residues lining the 2RGP
pocket). The missing part — a model mapping pocket residues to a probabilistic feature
map — is a serious build, not a weekend.

**Verdict: the most promising unimplemented idea for us**, because it targets the exact
capability we tested and rejected.

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

**Verdict: no code to adopt, but strong independent corroboration of the approach.**
Fully open at `github.com/postera-ai/del-simulator` if we ever touch DEL data.

---

## 3. Units and reporting in preclinical ADME/PK (2026)

Paywalled; I have title, author (Matthew C. Robinson, PostEra) and venue only, and will
not summarise content I have not read.

**Why it is the most immediately applicable of the six.** `dmta/ingest.py` already
fights precisely this problem, and the evidence is in our own runs:

- Units inferred from a column name, then from the value distribution (nM vs µM vs p-scale)
- Censored values — **51 of 123 Open Source Malaria compounds** were "> X", i.e. too weak to measure
- Readings that are not concentrations at all: `60% @ 120000`, `25-80% @ 80000`
- Replicate measurements disagreeing by more than 1 log unit
- A hardcoded "active ≥ pIC50 8" that silently reported 0% on a programme whose best compound was 7.71

A best-practice guide for exactly these decisions would be worth reading properly
before the ingest rules are finalised. **Worth obtaining the full text.**

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

**Verdict: a better-motivated next model than Chemprop was**, on this evidence.

---

## What this review changes

1. **Boltz-2 recommendation corrected** in ROADMAP.md, citing paper 1.
2. **Pharmacophore screening replaces "try another docking structure"** as the route to
   structure-based signal, since docking has already failed three times here.
3. **Paper 6 displaces Chemprop** as the model worth trying next — Chemprop has now been
   measured and did not win.
4. **Paper 3 should be read in full before the ingest rules are frozen.** It is the one
   whose subject we are already reimplementing by trial and error.
5. **Paper 5 is the template** for the de novo capability, if enumeration is ever
   outgrown.

Nothing here is implemented yet. This file exists so that none of it is quietly
forgotten, and so that the two corrections above are recorded next to the decisions
they overturn.
