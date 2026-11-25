# What it takes to actually do what the papers say

[PAPERS.md](PAPERS.md) says what the six papers do and that none of them was in the
platform. This file is the other half: **for each one, what we build, what we need from
you, what it costs, and the known-answer test it must pass before it ships.**

Every row obeys the rule the rest of the repository obeys: a capability is not adopted
because a paper reports it working. It is adopted when it beats a test whose answer we
already knew. That rule is what killed docking here (0.0–0.37× enrichment across three
structures) and Chemprop (median ρ 0.425 against the forest's 0.510), and it applies to
these six exactly as it applied to those.

## The short version

**Updated 2026-11-25.** Items 6, 2 and 3 are closed. Two were rejected on measurement and
one turned out not to apply; the remaining work is paper 1 and paper 5.

| # | Paper | Status | GPU cost | Needs from you |
|---|---|---|---|---|
| 6 | RMT linear model | ✅ built, **measured, rejected** — worse than the forest | $0 | — |
| 2 | DEL Simulator | ✅ built as an auditor — **falsified our own 2.75×** | $0 | — |
| 3 | ADME/PK reporting | ✅ read in full — **does not apply yet**; my earlier claim was wrong | $0 | — |
| 1 | Pharmacophore maps | 2 days → 2 weeks | ~$0 | nothing |
| 5 | De novo design | 2–3 weeks | $0 | a chemist to judge output |
| 4 | PLK3 campaign | hours | $0 | the supporting information |

Total GPU spend to do all of it: **approximately nothing.** Every one of these is CPU
work. That is the single most useful thing to know — the papers that look like they need
a GPU cluster do not, and the one capability that genuinely wanted a GPU (docking) is
the one we measured and rejected.

---

## ✅ DONE: DEL Simulator (paper 2) — it audited our results, and they failed

**Outcome, 2026-11-25.** Built in `dmta/del_landscape.py` and `del_benchmark.py`, one day
as estimated. **The test failed, which is the best thing that happened this week.** Our
2.75× became **0.71×** in their landscape — worse than random screening. The mechanism is
that their landscape is unlearnable from fingerprints (ρ ≈ 0 with 2,000 training
compounds), so the loop's model is blind there; what the exercise really proved is that
**our own simulator is far too easy** (ρ 0.89 from 100 compounds versus a real-programme
median of 0.56). The README's 2.75× is now labelled an upper bound, and the loop's
real-data evidence — which never used the simulator — stands unchanged. Full write-up with
the audit's own weaknesses in [RESULTS.md](RESULTS.md) 9.

**Why it was first.** It costs nothing, it is fully open
(`github.com/postera-ai/del-simulator`), and it tests something we currently have no way
to test: **whether our headline result is real or an artefact of our own simulator.**

We report that UCB selection beats random by 2.75× ([RESULTS.md](RESULTS.md) §2). That
was measured inside `dmta/assay.py`'s `SimulatedAssay` — a landscape we wrote. If the
2.75× is a property of that landscape rather than of the method, we would never know. The
DEL Simulator is a second, independently authored ground-truth landscape. Running our
acquisition functions in theirs is the cleanest possible check.

**Build**
1. `pip install` their simulator; generate a library with a known structure–affinity map.
2. Adapter so `dmta/select.py` can select from their compounds — a `Compound` shim, since
   selection only needs SMILES plus a measurable value.
3. Re-run the acquisition comparison: UCB, EI, Thompson, greedy, uncertainty, random.

**Known-answer test.** Does UCB still beat random, and by roughly how much? **Pass:** the
ordering of acquisition functions survives. **Fail:** it does not, and the 2.75× in
RESULTS.md gets a caveat saying it is simulator-specific. Either outcome is worth having;
the second is worth more.

**Cost:** $0, CPU. **Effort:** 1–2 days.

---

## ✅ READ, and the plan below was built on a wrong premise (paper 3)

**Correction, 2026-11-25.** I had not read this paper when I wrote the plan below, and I
guessed its contents from its title. Having now read it: it is an Editorial on reporting
conventions for **ADME and PK** data — clearance, volume of distribution, microsomal and
hepatocyte stability, protein binding, preclinical efficacy — with per-assay checklists
and a unit-conversion appendix. It says **nothing** about censored values, nothing about
replicate disagreement, and almost nothing about potency. The platform measures no ADME or
PK endpoint, so **it does not apply to `ingest.py` at all yet.** See
[PAPERS.md](PAPERS.md) 3 for what it does say, including the one piece of its advice that
does land on our data.

**The censored-value work below is still worth doing.** It was the change this paper was
supposed to justify. It does not justify it — but 51 of 123 Open Source Malaria compounds
still carry information we are throwing away, and that was true before any paper said so.

**Censored values are being thrown away.** Open Source Malaria: **51 of 123 compounds
report "> X"** — too weak to measure. We currently handle those by exclusion. That is
throwing away the most informative possible statement about a compound ("this one is
definitely bad") and it biases the training set towards actives.

**Build**
1. A censored likelihood (Tobit regression): a value recorded as `> 10 µM` contributes
   the probability that the true value exceeds 10, not a point. Straightforward for a
   linear model — and we now have one (the RMT model), which is why this is suddenly easy.
2. Make unit declaration explicit rather than inferred. `ingest.py` currently guesses from
   the column name, then from the value distribution. A guess that is right 95% of the
   time is a silent 5% error. Guess, then **show the user and make them confirm** — the
   upload page already has the inspect step for this.
3. A provenance field per value: file, column, original string, inferred unit, transform
   applied. Reports already take a `provenance` argument; this extends it per datum.

**Known-answer test.** On Open Source Malaria, scaffold split: does a censored-aware fit
beat the drop-them fit on ρ and on enrichment? **This is measurable today.** The paper
would tell us which other conventions to encode; it is not needed to start.

**From you:** nothing — you supplied the paper, and it turned out not to be needed for
this. Items 2 and 3 of the build list below are now optional polish rather than
paper-driven requirements; item 1, the censored likelihood, is the one that pays.

**Cost:** $0. **Effort:** 1–2 days for the censored likelihood alone.

---

## 3rd: Pharmacophore maps (paper 1) — the replacement for docking

Docking failed here three times. This is the paper's own answer to why, and it is the
most promising unimplemented idea we have. It splits into a cheap test and an expensive
build, and **we should only pay for the second if the first passes.**

### Phase A — inference only (about 2 days)

Their inference code and weights are released (`github.com/postera-ai/phc4-eld-generator`,
weights on Zenodo); **training and benchmarking code are not.** So the first version uses
their model as a black box.

We already own most of the surrounding machinery: `dmta/structure.py` fetches and cleans a
PDB entry and derives the binding box from the crystal ligand, and the 3D viewer already
computes the residues lining the pocket (58 of them for 2RGP). What is missing:

1. Pocket → their model's expected input format; out comes a predicted feature map.
2. Ligand side: RDKit ETKDG conformers, then RDKit's own feature factory
   (`BaseFeatures.fdef`) to annotate donors, acceptors, aromatics, hydrophobes and
   ionisables. All standard RDKit, all already a dependency.
3. Score each compound by how well its features can be aligned onto the predicted map.

**Known-answer test — the same one docking failed.** The 8,978-compound EGFR ChEMBL set,
scaffold split. Two bars, and it must clear both:
- enrichment in the top decile **> 1.5×** (docking managed 0.37×, 0.0× and 0.0×);
- some real signal on the **time split**, where the 2D forest collapses to ρ 0.11. That
  split is the whole reason to want a structure-based method at all. A method that only
  works where the 2D model already works has bought us nothing.

### Phase B — only if Phase A passes (about 2 weeks)

Vector-database search so this scales the way the paper describes: encode pharmacophores
as 3- or 4-point fingerprints, index with FAISS or hnswlib, and search sublinearly in
library size. This is the part that turns a scoring function into a screen over billions.

**Cost:** ~$0. Conformers for 9k compounds is about a CPU-hour; their model runs on one
pocket in seconds. **Effort:** 2 days for the test, 2 weeks for the scaled version.

---

## 4th: De novo design (paper 5) — the real capability upgrade

The closest published work to our loop: from 42 measured compounds and 515 inactives,
their model designed 5 compounds *with predicted routes*; all 5 were made and 3 were
active. Where we enumerate a fixed 924-compound library, they generate novel structures.

**Three ways to build it, and I recommend the cheapest.**

| Route | What it is | Honest assessment |
|---|---|---|
| **(a) Dynamic enumeration** | Keep the reaction-SMARTS approach, but expand the building-block set hugely and enumerate lazily on demand | **Recommended.** Routes stay free by construction — the thing the paper is actually about. Maybe 80% of the benefit for 10% of the work |
| (b) Generative model + route prediction | A SMILES or graph generator, then AiZynthFinder (open, MIT) for routes | The faithful reimplementation. Weeks, and route prediction is its own failure surface |
| (c) Their released code | `github.com/wjm41/mpro-rank-gen` | Training set and generated molecules are open. Worth reading before building either of the above |

Route (a) is a change to `dmta/library.py` rather than a new subsystem: more scaffolds,
more amines, more boronic acids, enumerated lazily and scored in batches. Every product
still carries a route, a yield estimate, a cost and a lead time, because that is how the
library is built.

**Known-answer test — retrospective, on data we already have.** On COVID Moonshot, hide
everything after a cut-off date. Let the generator propose. Then ask: among compounds that
the real programme *did* subsequently make, were the ones our loop proposed more potent
than the ones it did not? We already have the time-split machinery for this.

**From you:** this is the first item that needs a chemist, not a test. The published
practitioner guidance is that batch size should match **reviewer capacity** — someone has
to look at the proposals and say "no, not that one". Without that, generation produces
plausible nonsense at scale.

**Cost:** $0. **Effort:** 2–3 weeks for something credible.

---

## Opportunistic: PLK3 (paper 4)

Not a method — a campaign, with Pfizer. Useful only if the supporting information contains
a table of compounds with potencies, in which case it is another real programme for
`evaluate`, which takes a CSV and produces a report. **Hours, if the SI is accessible.**
Needs the same ACS access as paper 3.

---

## Done: RMT linear model (paper 6)

Implemented in `dmta/rmt_model.py` this session, with the method's central claim as a
unit test: random data must produce **zero** signal eigenvalues, and the cleaned model
must beat plain least squares on identical features. Both hold. The measured comparison
against the forest and Chemprop is in [RESULTS.md](RESULTS.md) §8.

---

## What I need from you, in one list

1. **Nothing, to continue.** Paper 1 (pharmacophore maps) needs no access and no budget;
   its Phase A is two days and faces the test docking already failed.
2. **A chemist's half-hour per batch**, before de novo design is worth starting. Not
   before — the loop as it stands proposes from a library a chemist already sanctioned.
3. **A decision on paper 5's route** (a), (b) or (c) when we get there. I recommend (a).
4. **The supporting information for paper 4**, if you want another real programme to
   evaluate. Optional.

Nothing on this list needs a GPU, and nothing is blocked on access any more.

## What the three closed items actually taught us

Not one of them became a feature, and that is the point:

- **Two models tried, neither better than a random forest** (Chemprop, then RMT). The
  model is not where the remaining value is.
- **Our own benchmark was flattering us**, and only an outside landscape could show it.
- **A paper I was confident about turned out not to apply**, and I only found that out by
  reading it rather than its title.

The honest summary of the literature review is therefore: it corrected three of our own
claims and added no capability. That is a better return than it sounds.
