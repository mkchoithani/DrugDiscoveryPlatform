# What the funded cohorts have that we do not, and what to do about it

A feature list derived from the companies in [COMPETITORS.md](COMPETITORS.md) that are
*not* our competitors — twelve US, then six Indian — plus a build plan.
Written 2026-11-26.

The second cohort is the one that matters commercially: it contains a customer we can
name and a gatekeeper we cannot route around. Its items are **F9–F12**, and two of them
outrank most of what the US list produced. A third group — Benchling, Latch Bio and
Applied Compute — sits at the *agent* layer rather than the molecule layer; its items are
**F13–F15**, and F13 is the cheapest thing on this page.

**The rule this document obeys:** nothing on this list ships because a funded company has
it. Each item carries the known-answer test it must pass first, and the condition under
which we delete it. That rule has already killed four capabilities here
([RESULTS.md](RESULTS.md)), and a feature list is exactly where it is most likely to be
quietly suspended.

**One caveat before the list.** Our top blocker is not a missing feature. It is that no
compound this loop chose has ever been made. Everything below improves a platform that
still has no prospective result. If only one thing happens next, it should be that — not
anything on this page.

---

## The whole list on one page

For comparison work. Plain-language version; the detail for each is further down.

| # | In plain words | Who has it | We have | Effort | Status |
|---|---|---|---|---|---|
| **F1** | Predict ADME — will the drug dissolve, survive the liver — not just potency | Every clinical company | partly | days | **DONE, split verdict.** logD passes (ρ 0.757 vs cLogP 0.219). Solubility fails — ship a ridge instead |
| **F2** | Check a route survives real manufacturing scale, not just that it exists | Resilience, WuXi | no | days | planned; expect it to lose to counting steps |
| **F3** | Show *what was sacrificed* in a compound, not one blended score | all | no | days | planned |
| **F4** | Judge a molecule by what it does to a cell, not by what it looks like | Eikon, Altos, NewLimit | no | ~4 wks | **planned — highest upside.** Only untried fix for our σ collapse |
| **F5** | Say out loud where our model stops working (natural-product chemistry) | Enveda, LifeMine | no | days | planned as a disclosure |
| **F6–F8** | Peptides, antibodies, conjugates | Chai, Earendil, Isomorphic, Parabilis | no | months | **recommend against.** A second platform, not a feature |
| **F9** | Stop assuming every drug is an oral pill — antibiotics break our scoring | exposed by Bugworks | **defect** | days | in progress; a bug, not a gap |
| **F10** | A tamper-proof record of who changed what, that an auditor accepts | Zifo | no | weeks | required before any pilot |
| **F11** | Plug into the customer's lab software | Zifo | no | — | **don't build.** Partnership surface |
| **F12** | Put the cost of each answer on the report, next to the accuracy | Immuneel, Cellogen (as strategy) | **true, unstated** | days | planned |
| **F13** | Prove the molecules never leave the customer's machine | Applied Compute / Palo Alto | **true, unstated** | days | **cheapest item here** |
| **F14** | Use our examiner to grade *other people's AI claims* | nobody | machinery exists | months | written down, not started |
| **F15** | The Benchling connector specifically | Benchling (1,300+ labs) | no | — | build when a customer asks |

**Three of these we already have and have never said:** F12 (our forest matches Chemprop
at 96× less compute), F13 (the loop makes no network calls), and the examiner itself.
**One is a defect:** F9. **One is nearly free and nobody else claims it:** F13.

And the thing no feature fixes: no compound this loop selected has ever been made.

---

## First, the honest cut: what cannot be copied

Four of the twelve have capabilities that are not features and cannot be built by us at
any budget. Listing them as a roadmap item would be dishonest.

| Company | Capability | Why it is not available to us |
|---|---|---|
| Resilience | Sterile biologics manufacturing, two sites | A factory. $750M of plant. Not software. |
| Kailera, Eikon, Rakuten | Clinical-stage pipelines, Phase 2/3 trials | Requires an asset, an IND and a clinical organisation |
| Enveda | Proprietary mass-spec metabolite library | Years of wet-lab data collection; the library *is* the moat |
| LifeMine | Fungal strain collection + Avatar-Rx genomics | Same — a physical collection, and target-side rather than compound-side |
| Altos, NewLimit, Retro | Primary-cell and humanised-mouse screening | Wet lab |

What is left after that cut is genuinely transferable, and there is more of it than I
expected.

---

## Tier 1 — build now. Public known answer exists, CPU only.

### F1. Measured ADMET endpoints, not just computed physchem

**Who has it:** every clinical company in the list — Kailera, Eikon, Parabilis, Rakuten.
It is why they have programmes and we have benchmarks.

**What we have:** `scoring.py` optimises measured potency plus *computed* physchem (MW,
cLogP, TPSA, QED, SA). **What we lack:** any endpoint that was measured in a lab other
than the potency one. RESULTS.md names this exactly — *"no selectivity, no ADMET, no tox —
the objectives a real programme dies on."*

**Known-answer test:** the [TDC ADMET Benchmark Group](https://tdcommons.ai/benchmark/admet_group/overview/)
— 22 standardised datasets, 80/20 scaffold splits, published leaderboards. Start with
aqueous solubility (AqSolDB, ~9,982 compounds) and microsomal clearance.

**Pass condition:** our model must (a) land within range of the published leaderboard on
the scaffold split, and (b) report σ reliability per split the way `chembl-validate`
already does for potency. If σ on an ADME endpoint is as uninformative as it is on our
time split, the endpoint is reported and *not* given to an acquisition function.

**Cost:** CPU, hours. **Verdict risk:** low — this is the best-characterised task in the
field.

**Result, 2026-11-26 (RESULTS.md §18): a split verdict.** logD passes — ρ 0.757 on
held-out scaffolds against cLogP's 0.219, and the forest still beats cLogP by +0.478,
CI [+0.407, +0.549], with every descriptor stripped out. Solubility does **not**: it ties
a twelve-descriptor ridge (Δρ −0.001, CI [−0.033, +0.032]) and loses on error, MAE 1.053
against 0.731. We ship the ridge for that endpoint. Neither public set carries dates, so
there is no time split — the column that matters most is still untested here. Note the 2026 critique that scaffold splits hide ADMET failures
([arXiv 2607.10729](https://arxiv.org/pdf/2607.10729)); we should run our time split too,
which the leaderboard does not.

### F2. Manufacturability as a constraint, not a cost number

**Who has it:** Resilience, as a business. Rakuten and Earendil live it through their
WuXi/LOTTE manufacturing agreements.

**What we have:** `library.py` gives every compound a route, a reagent cost and a lead
time, and `select.py` plans batches that share intermediates. **What we lack:** anything
about whether a route survives *scale*. Step count at scale, hazardous reagent classes,
chiral centre count, and chromatography-free feasibility are all invisible to us.

**Known-answer test:** hold out whole route classes and ask whether our makeability score
separates them. The external anchor is Enamine REAL's reported >80% synthesis success in
3–4 weeks — our score should rank the REAL-space subset above our exotic reactions.

**Pass condition:** the score must rank held-out route classes correctly better than
counting steps does. If counting steps wins, we ship counting steps and delete the score —
this is precisely how the pharmacophore work ended (beaten by counting heavy atoms).

**Cost:** CPU, days.

### F3. Show the trade-off, not the composite

**Who has it:** every drug company implicitly. A programme is a negotiation between
potency, exposure and safety; a single desirability number hides it.

**What we have:** desirability curves collapsing to one score. **What we lack:** any view
of where a compound was *sacrificed*, which is the first thing a chemist asks.

**Known-answer test:** none needed — this is a reporting change, not a claim. It therefore
ships without an examiner test, and this document should say so plainly rather than
inventing a benchmark to look rigorous.

**Cost:** CPU, days. `report.py` already renders self-contained HTML.
---

## Tier 1b — from the Indian cohort. These outrank most of Tier 2.

Added 2026-11-26, from the second cohort in [COMPETITORS.md](COMPETITORS.md). They sit in
Tier 1 because two of them are cheap, and one of them is a defect rather than a feature.

### F9. Property profiles that are not oral-drug profiles — **BUILT 2026-11-26, defect confirmed**

**Who exposed it:** Bugworks, developing Gram-negative antibacterials.

**The problem:** `scoring.py` applies desirability curves over MW, cLogP, TPSA and QED.
Those curves encode oral small-molecule drug-likeness. **Gram-negative antibacterials do
not live in that space**, and the rules that do govern them are different ones: the eNTRy
rules, derived from accumulation measurements on 180+ compounds in *E. coli*, call for a
non-sterically-hindered ionizable nitrogen (primary amines best), globularity ≤ 0.25, and
≤ 5 rotatable bonds. Their authors are explicit that these are empirical guidelines with
common exceptions, not filters.

So today we would hand an antibacterial programme a profile tuned for a different disease
class and score their best compounds badly. The architecture is already right — the
profile is data, not code, by an explicit design decision in LANDSCAPE.md — but **we ship
exactly one profile and present it as neutral.**

**Known-answer test, and it is a good one:** run the marketed antibacterials through our
default profile. **If we reject drugs that already exist and are on the market, our profile
is wrong** — and the number of them we reject is the size of the error. Then repeat with an
eNTRy-rules profile and check it rejects fewer.

**Pass condition:** the antibacterial profile must reject materially fewer marketed
antibacterials than the default does. If both reject the same number, the profile
machinery is not doing what we claim and the finding belongs in RESULTS.md.

**Cost:** CPU, days. Globularity and rotatable-bond counts are already available through
RDKit; the ionizable-nitrogen test is a SMARTS match.

**Why this is first:** it is the only item on either list that is a **defect we can
demonstrate against public data today**, it costs days, and it converts a prospect's first
objection into our own published finding.

**Outcome ([RESULTS.md](RESULTS.md) §19):** confirmed. The shipped profile rejects
marketed antibacterials **+19.6%** more often than other marketed drugs, CI
[+10.8%, +28.7%], at ≤600 Da where weight is not the explanation. Daptomycin scores
exactly zero; vancomycin, amikacin and streptomycin are all marked down. A repaired
profile ships as `PROFILES["antibacterial"]` and takes held-out rejection from 39.3% to
4.9%. Two failed repair attempts are recorded, and the second produced the project's
**sixth verdict-function error** — a saturation check measured on the property terms
alone, when potency is what ranks, which would have shipped a profile that marks down
ciprofloxacin. It also fits our standard exactly — the
pharmacophore work died because counting heavy atoms beat it, and this is the same kind of
test pointed at ourselves.

### F10. A real audit trail, validated — the entry condition, not a feature

**Who has it:** Zifo, as a business line (validation, documentation, audit readiness,
governance), across a reported 190+ organisations.

**What we have:** a SQLite store that records campaigns, batches, measurements and model
runs — the raw material. **What we lack:** immutability, per-user attribution, a
reason-for-change on every record, and any documentation a quality unit would accept.

**Known-answer test:** none — this is compliance, not science, and inventing a benchmark
for it would be exactly the dishonesty this document warns about. The test is an audit
checklist, passed or failed.

**Cost:** weeks. Unglamorous and non-negotiable. Already in Phase 1; this cohort gives it
a name and shows there is a services industry built on it.

### F11. Instrument, ELN and LIMS integration

**Who has it:** Zifo, as its core competence.

**The honest position:** "the lab is a CSV" remains the right first contract, and we
should not build connectors speculatively. But a CSV is the thing a customer's informatics
team has to wrap, and that wrapping is Zifo's business model. **Recommendation: do not
build this; treat it as the partnership surface.** An examiner is a natural product for a
validation-and-governance firm to resell, and they have the distribution we never will.

### F12. Cost per answer, reported as a first-class number

**Who made it a strategy:** Immuneel and Cellogen, whose pitch is a CAR-T course at
$60–70k against a global $500–700k. In this market the funded differentiator is cost per
outcome.

**What we have:** the fact, buried as a footnote — the forest is level with Chemprop at
**96× less compute**, and docking, pharmacophore and RMT were each rejected for a few hours
of CPU.

**What to add:** compute cost, wall-clock and cost-per-compound-scored on every evaluation
report, next to the accuracy numbers. `report.py` already renders the artefact.

**Known-answer test:** none needed — it is measurement and reporting, not a claim.

**Cost:** days.

---

## Tier 1c — the agent layer: Applied Compute, Benchling, Latch

Added 2026-11-26. **This is a different layer from everything above, and it is the first
group that is neither competitor nor customer but a *platform we would sit on top of*.**

### What is actually being done there

| | what it is | layer |
|---|---|---|
| **Benchling** | Cloud R&D platform — ELN, LIMS, registry — across a reported 1,300+ biotechs including Merck, Moderna and Sanofi. Benchling AI adds agents grounded in the customer's own structured R&D data, with MCP-based connectors. **Science Foundries** (Oct 2026) lets agents submit designs and order their own experiments by web or API. | the system of record |
| **Latch Bio** | No-code cloud bioinformatics: versioned containerised pipelines, RNA-seq, CRISPR, structure prediction. Positions as "the AI agent for biology data analysis", and publishes benchmark work across antibodies, single-cell, spatial and biosecurity. | the analysis layer |
| **Applied Compute** | Post-training and inference of **open-weight** models on a customer's proprietary data, in the customer's own or a sovereign environment. Acquired by Palo Alto Networks on 18 Nov 2026 for **$4.5B in stock**, closing FY3Q 2027. | the model layer |

None of the three decides which molecule to make. **They automate the work around the
decision; we make the decision.** So this is not a competitive threat — it is the first
credible distribution surface either cohort has produced.

### The claim in the source, and what is wrong with how it was measured

The pasted note reports that Applied Compute's custom model beat frontier models by
**over 11 percentage points on Benchling's protocol benchmarks**, on the reasoning that
scientists rarely publish how they troubleshoot failed protocols, so public models have
never seen the data.

**The reasoning is exactly ours.** It is the known-answer argument moved one layer up: find
a task where the answer is known internally and absent publicly, and measure against it.
We should take it seriously for that reason.

**But look at who owns what.** The benchmark belongs to Benchling. The model belongs to
Applied Compute. The partnership belongs to both. The number was published on the
vendor's own feed, and one search did not surface it anywhere independent. There is no
split discipline stated, no confidence interval, no baseline description, and no account
of what a tie would have looked like.

That is not an accusation of dishonesty — it is the ordinary structure of a vendor
benchmark, and we have published our own version of the same mistake five times and
corrected it each time ([RESULTS.md](RESULTS.md) §16, §17). **It is also, precisely, the
market we said we were in.** A $4.5B acquisition has just made agent-performance claims in
biotech a category with real money behind it and no independent referee.

### F13. Sovereignty as a measured property, not a marketing line

**Who made it valuable:** Palo Alto paid $4.5B substantially for the ability to train and
serve models in a customer's own environment. The pasted note names the reason: molecular
data, IP and patient records cannot go to a public model provider.

**What we already have, by accident:** a random forest on CPU, a local SQLite store, and a
loop that calls nothing. I checked — **only two modules in the codebase open a network
connection at all** (`chembl.py` fetching activity data, `structure.py` fetching a PDB),
and both are explicit, separate fetch steps. The Design–Make–Test cycle itself does not
touch the network.

**For a customer whose compound structures *are* the asset, that is a feature we have never
once claimed.**

**Known-answer test — and this one can be genuinely proven rather than asserted:** run a
full campaign cycle with the socket layer blocked, and assert it completes. A test that
fails the moment any module tries to open a connection turns "your molecules never leave
your machine" from a promise into something with a test name attached to it.

**Pass condition:** a complete cycle — enumerate, score, select, plan, ingest results,
retrain — with sockets disabled, and the ChEMBL and PDB fetchers clearly marked as the
only opt-in network steps.

**Cost:** days. **This is the single cheapest credibility gain on any of these lists**, and
unlike every other item here it needs no new science — only a test and a sentence we are
already entitled to write.

### F14. Point the examiner at agent claims

**The observation:** our verdict machinery is not molecule-specific. The paired bootstrap,
the three-valued verdict, the scaffold/time split discipline, the pre-specification
reported *alongside* the result — none of it knows what is being compared. We built it to
examine property models; §16 and §17 show it examining our own reasoning.

An "11 points better on our own benchmark" claim is the same shape as "cleaning beats no
cleaning" and would go through the same machine: what is the control, what is the split,
what is the interval, and what would a tie have looked like?

**Honest limits.** Agent benchmarks are not AUC over compounds. Pass rates on protocol
tasks are the natural unit; our bootstrap resamples compounds, so the resampling unit
would have to change. That is real work, not a rename. And we have no agent benchmark of
our own and no access to Benchling's.

**Recommendation: do not build this yet — write it down as the direction.** It is the
clearest commercial expression of the thing [COMPETITORS.md](COMPETITORS.md) concluded we
should sell, and it now has a $4.5B category forming around it. But it is a product bet,
and Phase 0 outranks it.

### F15. The Benchling surface, which changes F11's answer

F11 said instrument and ELN integration was a partnership surface and not a build, with
Zifo as the example. Benchling sharpens that: **a specific platform, 1,300+ organisations,
MCP connectors, and agents that can already order experiments through an API.**

"The lab is a CSV" was the right first contract and should stay. But the API we said we
would not build now exists on someone else's side of the line, and our worklist-out /
results-in loop is the exact shape it accepts.

**Recommendation: still do not build it speculatively** — build it when a customer on
Benchling asks. The change is that we now know what to build *against*, which we did not
before.


---

## Revised order

F9 and F12 move ahead of most of Tier 2, because both are days of work and both change what
we can say to a customer we can actually name.

| | item | why it moved |
|---|---|---|
| **Phase 1** | F1 ADME · **F9 non-oral profiles** · **F13 sovereignty, proven by test** · F3 trade-off view · **F12 cost reporting** · F10 audit trail | F9 is a demonstrable defect; F13 is days of work for a property we already have and never claimed; F12 is a strategy this market already rewards |
| **Phase 2** | F4 Cell Painting · F2 manufacturability | unchanged — F4 remains the highest-value experiment |
| **Phase 3** | F5 natural-product domain | unchanged |
| **Not scheduled** | F6–F8 modalities · **F11 / F15 integrations** · **F14 examiner for agents** | F11 and F15 are a partnership surface, not a build; F14 is the product bet, written down rather than started |

The Phase 0 caveat is unchanged and now has more force, not less. A discovery organisation
with $43M has fewer compounds to waste than Lilly does, so **"what has this loop
prospectively found?" arrives sooner from this customer, not later.**

---

## Tier 2 — build next. Needs new data plumbing; public data exists.

### F4. A phenotypic readout — the one with real upside

**Who has it:** Eikon, whose whole discovery engine is live-cell microscopy. Altos,
NewLimit and Retro all read cell state rather than binding.

**Why this one matters more than the rest:** our central unsolved finding is that σ
collapses on the time split (0.44 random → 0.36 scaffold → **0.12 time**, inverting to
−0.13 on μ-opioid), because a fingerprint model only knows resemblance to what it has
seen. Docking was the proposed fix and failed four times; Chemprop tied. **A morphological
profile is not a structural resemblance** — it is the only candidate fix left that is
different in kind rather than degree.

**Known-answer test:** the [JUMP Cell Painting dataset](https://www.biorxiv.org/content/10.1101/2023.03.23.534023.full.pdf)
(cpg0016; >116,750 compounds in U2OS, public on the AWS Cell Painting Gallery), with
JUMP-MOA as the positive-control plate. The test is the one we already run: does adding a
morphological profile beat the fingerprint forest **on the time split**, with a paired
interval over targets?

**Pass condition:** a paired CI that excludes zero on the time split. A tie is a tie —
Chemprop got one and we kept the forest on cost. Beware the known batch-effect problem;
any gain must survive batch correction or it is a measurement of the plates.

**Cost:** modest GPU for the embedding — from [ROADMAP.md](ROADMAP.md)'s catalogue read, an
RTX 4090 at ~$0.34/hr community. Call it tens of dollars, plus real engineering time for
the data plumbing.

**Honest risk:** we would be joining compound structure to a phenotype measured on *other*
compounds, in a different assay, in one cell line. That is a weaker link than it sounds,
and the test above is designed to expose it rather than hide it.

### F5. Natural-product chemistry support — the library, not the model

**Who has it:** Enveda (PRISM, MS2Mol), LifeMine.

**The honest version:** building an MS/MS-to-structure model is Enveda's entire company
and we would build a worse one. What we *can* fix cheaply is that natural-product-like
chemistry — macrocycles, high sp³ fraction, many chiral centres — is handled badly by
Morgan fingerprints and absent from our reaction-SMARTS library.

**Known-answer test:** run `chembl-validate` on a natural-product-derived target set and
compare scaffold-split ρ against our synthetic targets. If it collapses, we have found a
real applicability-domain limit and should say so in the product rather than quietly
scoring those compounds.

**Cost:** CPU, days. **This is a limitation disclosure as much as a feature**, which makes
it exactly our kind of work.

---

## Tier 3 — modality expansion. Large, and mostly a different product.

### F6. Peptides · F7. Antibodies and bispecifics · F8. Conjugates

**Who has it:** Chai (Chai-2/3, de novo antibodies), Earendil (bispecifics, Genentech),
Isomorphic (small molecules, antibodies, peptides, glues in one engine), Parabilis
(helical peptides), Rakuten (antibody-dye conjugates).

**The honest assessment:** this is not a feature, it is a second platform. Every layer we
own is small-molecule-specific — Morgan fingerprints, reaction SMARTS enumeration, RDKit
descriptors, SA score. A sequence has none of those. The *only* parts that transfer are
the ones that are genuinely ours: the campaign store, the acquisition functions, the
split discipline, the calibration check and the verdict function.

**That transfer is the interesting observation.** Our examiner is modality-agnostic. The
loop is not.

**Known-answer test, if we ever do it:** [SKEMPI 2.0](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6361233/)
— 7,085 manually curated mutations with measured binding changes, filterable to
antibody/antigen — or AB-Bind. Watch the inequality entries ("non-binding", "weaker
than"): we already have machinery for exactly that problem in `censored.py`, which is a
genuine and unexpected fit.

**Cost:** months, not days. **Recommendation: do not start this.** It competes directly
with three companies holding $7.7B, $787M and $630M, on their chosen ground.

---

## The plan

Each phase gates the next. No phase starts until the previous one's verdict is written
into RESULTS.md, including when the verdict is "delete it".

### Phase 0 — before any of this (the thing that actually blocks us)
One target, one chemist, one batch of 24 made and tested. Not a build. Unchanged as the
top item in ROADMAP.md, and nothing below substitutes for it.

### Phase 1 — credibility, ~2 weeks of work, $0 of compute
- **F1** ADMET on the TDC group, with our own time split added to their scaffold split.
- **F3** trade-off reporting in `report.py`.
- Auth, multi-user and an audit trail — not a feature from this list, but the floor for
  any pilot, and it blocks a customer conversation that F1 would otherwise unlock.

**Gate:** F1 passes on at least two endpoints, with σ reliability reported per split.
If σ is uninformative on both, we report the endpoints and keep them out of acquisition —
and that negative is itself publishable under our own standard.

### Phase 2 — the real experiment, ~4 weeks, tens of dollars of GPU
- **F4** JUMP Cell Painting against the time split.
- **F2** manufacturability scoring, in parallel since it shares nothing with F4.

**Gate:** F4's paired CI on the time split. This is the highest-value experiment available
to us, because a win would fix the finding that currently limits the whole platform, and a
loss closes the last cheap hypothesis for it — which is worth knowing and worth publishing
either way.

### Phase 3 — scope honesty, ~1 week
- **F5** natural-product applicability domain, reported as a limit.

### Phase 4 — not scheduled
- **F6–F8** modality expansion. Revisit only if a customer pays for it, and only after
  Phase 0 exists. Until then the right answer is that our examiner is modality-agnostic
  and our loop is not.

---

## What this list is worth, stated plainly

Five of these fifteen items are strong: **F13** because the property is already true of
the codebase and needs only a test to become a claim; **F9** because it is a defect we can
demonstrate against public data in days; **F1** because programmes die on ADME and we do
not model it; **F4** because it is the only untried fix for our central measured weakness;
and **F3** because it is cheap and chemists ask for it immediately. **F12** is nearly free
and worth doing beside them.

**F5** is mostly a disclosure. **F10** is compliance, not science, and non-negotiable.
**F2** may well be beaten by counting steps, and we should expect that. **F6–F8** should
not be built, and **F11** is a partnership surface rather than a build.

And none of the fifteen changes the fact that no compound this loop selected has ever been
made. A longer feature list is the most comfortable way to avoid that sentence, which is
the main reason this document leads with it.

## Sources

- [TDC ADMET Benchmark Group](https://tdcommons.ai/benchmark/admet_group/overview/) · [TDC paper, NeurIPS 2021](https://zitniklab.hms.harvard.edu/publications/papers/TDC-neurips21-main.pdf) · [Beyond Scaffold Splits (2026)](https://arxiv.org/pdf/2607.10729)
- [JUMP Cell Painting dataset](https://www.biorxiv.org/content/10.1101/2023.03.23.534023.full.pdf) · [consortium GitHub](https://github.com/jump-cellpainting) · [batch-correction work](https://arxiv.org/pdf/2601.22331)
- [SKEMPI 2.0](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6361233/)
- Agent layer, for F13–F15: [Palo Alto / Applied Compute, $4.5B](https://www.constellationr.com/insights/news/palo-alto-networks-acquires-applied-compute-45-billion) · [8-K press release](https://www.sec.gov/Archives/edgar/data/0001327567/000132756726000030/ex992pressrelease_november.htm) · [Benchling AI](https://www.benchling.com/ai) · [Benchling AI connectors](https://www.benchling.com/news/benchling-launches-ai-connectors-to-power-data-ecosystem) · [LatchBio](https://www.verdict.co.uk/latchbio-unveils-cloud-bioinformatics-platform/) · [Latch on verifiable environments](https://ai.engineer/orgs/latchbio)
  *The 11-percentage-point figure is from the vendor's own feed and did not surface independently in one search.*
- eNTRy rules, for F9: [converting Gram-positive-only compounds](https://pubmed.ncbi.nlm.nih.gov/29446459/) · [eNTRyway and the three criteria](https://pmc.ncbi.nlm.nih.gov/articles/PMC8715754) · [worked examples](https://pmc.ncbi.nlm.nih.gov/articles/PMC8612397) · [exceptions are common](https://revive.gardp.org/?p=14643)
- Company capabilities as cited in [COMPETITORS.md](COMPETITORS.md)
