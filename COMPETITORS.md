# The funded cohort, against us

Fourteen companies from a 2026 funding table, checked one at a time against what this
platform actually is. Written 2026-11-26.

Funding figures are as given in the source table and were **not** independently verified;
where my own search contradicted them it is noted. Product descriptions are from the
searches cited at the end, and vendor claims remain the vendors' own.

---

## The first finding: twelve of the fourteen are not competitors

The table mixes two completely different kinds of company, and the distinction decides
everything about what we should do.

| | what it is | how it makes money | competes with us? |
|---|---|---|---|
| **Isomorphic Labs** | AI drug-design engine (IsoDDE) + own pipeline | Pharma partnerships: upfront + milestones, plus wholly-owned programmes | **Adjacent** — same problem, no software sold |
| **Chai Discovery** | Chai-2/Chai-3, de novo antibody design models | Model access to selected partners (Novartis, argenx reported) | **Adjacent** — biologics, not small molecule |
| **Enveda** | PRISM: mass-spec/metabolomics discovery from natural chemistry | Own pipeline; two Phase readouts reported 2026 | **No** — different discovery input entirely |
| **Earendil Labs** | AI biologics; bispecifics | Genentech $55M upfront / >$1.5B milestones; Sanofi licence | **No** — biologics asset company |
| **Eikon Therapeutics** | Live-cell microscopy platform, internally applied | Late-stage oncology pipeline (EIK1001/1003/1005) | **No** — a clinical oncology company |
| **LifeMine** | Avatar-Rx, fungal genomics | GSK alliance, 50/50 to IND | **No** — different modality source |
| **Altos Labs** | Cellular rejuvenation research | Nothing yet — no disclosed candidate or IND as of Sep 2026 | **No** |
| **NewLimit** | Epigenetic reprogramming; AI on screening data | Nothing yet — first human study not before 2027 | **No** |
| **Retro Biosciences** | Autophagy, microglia, HSC reprogramming | Lead oral autophagy programme in Phase I | **No** |
| **Kailera** | Four clinical-stage GLP-1 obesity candidates | Drug asset company; Phase 3 from 1H 2027 | **No** |
| **Parabilis Medicines** | Helical peptide oncology (could not verify independently) | Drug asset company | **No** |
| **Rakuten Medical** | Alluminox photoimmunotherapy | Drug asset company | **No** |
| **Resilience** | Biologics CDMO, sterile drug product | Contract manufacturing | **No** — a supplier |
| **Isomorphic / Chai** | *see above* | | |

**So the table is a map of where capital went, not a map of our market.** Eleven of the
fourteen are drug companies. One is a factory. The two nearest neighbours — Isomorphic
and Chai — do not sell software either: Isomorphic runs partnerships and its own
pipeline, Chai gives model access to selected partners. Nobody in this list is going to
show up in a procurement process against us, because nobody in this list sells a thing a
med-chem team can buy.

The competitors we actually have are in [LANDSCAPE.md](LANDSCAPE.md): Schrödinger,
Iktos, PostEra, XtalPi, Optibrium. **None of them is in this funding table.** That is
itself the signal — see the next section.

## The second finding: the money is going to assets, not to tools

Every large, verifiable deal in this cohort is structured the same way:

| deal | upfront | milestones |
|---|---|---|
| Earendil / Genentech (Sep 2026) | $55M | >$1.5B |
| Isomorphic / Lilly (third-party estimate) | ~$45M | ~$1.7B |
| Isomorphic / Novartis (third-party estimate) | ~$37.5M | ~$1.2B+ |
| LifeMine / GSK (2022) | $70M incl. equity | shared cost to IND, then royalties |

The upfront is 2–4% of the headline. **The capital is being paid for a share of a future
drug, not for a platform.** A software licence cannot be structured this way, which is
why the software vendors are absent from a table built on round size. It does not mean
the software market is bad; it means round size is the wrong axis for finding it.

Two consequences for us:

1. If the goal is a round that looks like these, the product has to be an *asset*, and we
   have no asset, no target, no chemist and no lab. That is a different company.
2. If the goal is revenue from med-chem teams, this table tells us nothing about our
   competition and we should stop reading it as if it did.

## Where each near-neighbour is genuinely ahead of us

Only three comparisons are worth making in detail.

**Isomorphic Labs.** Structure prediction, binding-affinity estimation, pocket finding,
generative design and multi-objective optimisation in one engine, across small molecules,
antibodies, peptides and glues. We have multi-objective optimisation and the loop. We do
*not* have structure prediction or affinity prediction — and §4–6 of
[RESULTS.md](RESULTS.md) is the record of us trying the cheap version of that (docking)
and measuring it to be worthless on our target. They are not ahead of us because they
bought better software; they are ahead because they can afford the expensive version.
Their programmes are still preclinical as of November 2026, so neither of us has a drug.

**Chai Discovery.** Reported roughly half of ~50 targets yielding validated hits from
fewer than 20 designs each. That is a *prospective* number — designs made and tested. Our
2.37× is retrospective on compounds someone else already chose to make. This is the single
most important gap in this document and it is treated on its own below.

**Enveda.** Two positive early clinical readouts in 2026 from a platform-discovered
molecule. A platform with a readout is a different asset class from a platform with a
benchmark, however good the benchmark.

## Where we are ahead, and it is not nothing

No company in this table publishes what we publish. Across the three AI-platform
companies I could research, I found **no published negative result, no failed internal
capability, and no corrected prior claim.** Chai's figure is reported in a press release;
Isomorphic's engine is described in company articles; Enveda's platform in trade coverage.

We have:

- seven capabilities examined across eight claims, with three rejections, two ties, and
  the test that produced each ([RESULTS.md](RESULTS.md));
- five of our own published claims corrected, and five errors in our own verdict logic
  found and pinned with tests;
- docking rejected four times on our own target rather than shipped because the field
  uses it;
- a paired confidence interval attached to every comparison, and a tie reported as a tie.

In a market where the buyer cannot evaluate the claims — which is this market — that is a
real asset. It is also, today, an asset with no price tag attached to it.

---

# What is actually blocking us from selling

In priority order. The first one dominates the rest.

### 1. Nothing prospective. This is the gap.

Every number we have comes from revealing measurements that already existed. The loop has
never chosen a compound that was then made and tested. Chai's "half of 50 targets" and
Enveda's Phase readouts are different in kind, not degree, and a sophisticated buyer will
ask within the first ten minutes.

We have said this honestly in RESULTS.md — *"every compound in these libraries was made by
someone who thought it worth making, so even the random baseline sits on an enriched deck"*
— which protects our integrity and does nothing for our sales.

**What closes it:** one target, one chemist, one batch of 24 made and tested. This is
already the top blocker in [ROADMAP.md](ROADMAP.md). It does not need more compute.

### 2. We optimise the endpoint that programmes do not die on

The loop optimises measured potency plus computed physchem. No selectivity, no ADMET, no
tox. RESULTS.md states it plainly: *"the objectives a real programme dies on."* A buyer with
a live programme has a potency model already; what they lack is the trade-off against
everything else. Our multi-objective profile is the right architecture for this and is
currently pointed at the wrong objectives.

**What closes it:** one public ADME endpoint (solubility, microsomal clearance) run
through the same examiner, with the same scaffold/time split discipline. Cheap. CPU-only.

### 3. Our own best finding is unsolved and we are selling into it

σ ranks error at 0.44 random, 0.36 scaffold, **0.12 time** — and inverts to −0.13 on
μ-opioid. A production model lives in that third column. We currently handle this by
*telling the user* (`chembl-validate` recommends β from the worst split). That is the right
thing to do and it is not a fix.

Docking was the proposed fix and it failed four times. Chemprop tied. So we have measured
two of the three obvious answers and neither works.

**The honest positioning:** this is not a defect to hide before launch; it is the most
valuable thing we know, and nobody else in this cohort has published anything like it. But
it has to be positioned deliberately rather than discovered by a prospect.

### 4. No chemist, no lab, no catalogue

`make` is a CSV out and results in. That is correct for a first engagement and it caps us
at advisory. Iktos owns a robotic lab running 100 reactions/day; Resilience exists because
making things is hard. An Enamine-class catalogue licence is still unbought, which keeps
our "synthesisable by construction" space small.

### 5. Nothing a pharma IT department will accept

Single-process uvicorn, SQLite, no authentication, no multi-user, no audit trail, no
tenancy. The science is further along than the software. This blocks a pilot, not a demo.

---

# Where I think the product goes

**Recommendation: sell the examiner, not the loop.**

The loop is a commodity — LANDSCAPE.md already found that everyone converged on it and
differentiates on which arm they own. We cannot out-spend Isomorphic on prediction or
out-build Iktos on robotics.

But every company in that funding table has the same problem, and none of them has solved
it in public: **they cannot independently tell whether their own model works.** Isomorphic
publishes engine descriptions. Chai publishes a press release. Enveda publishes trade
coverage. The buyers of all three — Lilly, Novartis, J&J, Merck, Genentech, GSK, Sanofi —
are writing nine-figure milestone cheques against claims with no paired interval anywhere
near them.

What we have actually built and proven is a machine for converting a plausible claim into a
measured verdict with a confidence interval, including when the verdict goes against us.
Seven capabilities, eight claims, five self-corrections. That is the product. The DMTA
platform is the thing it was built to examine.

Three ways to take that to market, in the order I would try them:

1. **Validation as a service, aimed at the pharma side of these deals.** Someone with a
   $1.5B milestone exposure has an obvious reason to want an independent scaffold- and
   time-split verdict on the platform they are buying into. This needs no lab, no
   catalogue, no chemist — only what we have.
2. **The benchmark harness as the licensed artefact.** `chembl-validate`, the split
   discipline, the calibration check and the three-valued verdict function, sold to
   internal teams who want to hold their own models to it. Small, defensible, and it
   plays to the one thing we are measurably better at.
3. **The full DMTA platform, after a prospective result.** Not before. Gap 1 has to close
   first or we are selling a retrospective number into a market that has prospective ones.

What I would **not** do is position against the companies in that table on their terms.
They raised that money for assets. We do not have an asset, and the fastest way to look
small is to be compared on one.

## The next three things, concretely

1. One ADME endpoint through the examiner — cheapest credibility gain available, CPU-only.
2. Auth, multi-user and an audit trail — the floor for any pilot, independent of science.
3. One chemist and one target, for a single prospective batch of 24.

Items 1 and 2 are things we can do now. Item 3 is a conversation, not a build, and it has
been the top blocker in ROADMAP.md for some time.

---

## Sources

- [Isomorphic Labs — partnerships](https://www.isomorphiclabs.com/partnerships) · [drug design engine](https://isomorphiclabs.com/articles/the-isomorphic-labs-drug-design-engine-unlocks-a-new-frontier) · [J&J collaboration](https://www.isomorphiclabs.com/articles/isomorphic-labs-enters-into-a-research-collaboration-with-johnson-johnson) · [Merck, after the raise](https://www.fiercebiotech.com/biotech/isomorphic-labs-adds-merck-ai-drug-discovery-roster-days-after-5b-raise) · [Sacra profile](https://sacra.com/c/isomorphic-labs)
  *Note: sources disagree on the 2026 round — a May 2026 [press release](https://www.isomorphiclabs.com/press/isomorphic-labs-funding) says $2.1B Series B, November 2026 [trade coverage](https://www.drugdiscoverytrends.com/isomorphic-labs-adds-merck-as-an-oncology-partner-after-5-billion-raise/) says a $5B Series C. The source table lists both.*
- [Chai Discovery — Chai-2 announcement](https://www.businesswire.com/news/home/20250630307418/en/Chai-Discovery-Unveils-Chai-2-Breakthrough-Achieving-Fully-De-Novo-Antibody-Design-With-AI) · [Contrary Research](https://research.contrary.com/company/chai-discovery)
- [Enveda — Series E](https://secure.businesswire.com/news/home/20260923027166/en/Enveda-Raises-%24311-Million-From-Leading-AI-and-Biotech-Investors-to-Bring-Pharma-Into-the-21st-Century) · [pipeline coverage](https://www.pharmaceutical-technology.com/newsletters/enveda-nabs-311m-to-progress-nature-derived-drug-pipeline)
- [Earendil / Genentech](https://www.biopharminternational.com/view/earendil-labs-1-5-billion-bispecific-antibody-pact-genentech) · [funding](https://www.biopharminternational.com/view/earendil-labs-to-scale-ai-driven-biologics-platform-with-787-million-funding)
- [Eikon — Q1 2026 results](https://www.barchart.com/story/news/1857266/eikon-therapeutics-reports-first-quarter-2026-financial-results-and-provides-clinical-updates) · [ESMO abstracts](https://eikontx.gcs-web.com/news-releases/news-release-details/eikon-therapeutics-announces-seven-abstracts-accepted)
- [LifeMine / GSK alliance](https://insights.citeline.com/SC146115/LifeMine-Gets-Funding-Inks-Discovery-Alliance-With-GSK)
- [Altos, NewLimit, Retro — comparison](https://sacra.com/c/altos-labs/) · [NewLimit](https://research.contrary.com/company/newlimit)
- [Kailera — Q3 2026](https://www.streetinsider.com/dr/news.php?id=27289696) · [company site](https://www.kailera.com/)
- [Resilience — manufacturing](https://resilience.com/manufacturing) · [CDMO coverage](https://www.pharmamanufacturing.com/sector/contract-manufacturing/article/55389494/resilience-continues-to-find-its-footing-in-cdmo-space-amid-headquarters-relocation)
- [Rakuten Medical / LOTTE Biologics](https://www.aap.com.au/aapreleases/cision20260113ae59192)
- [XtalPi — drug discovery, active learning and DMTA](https://en.xtalpi.com/drug-discovery/) · [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger,_Inc.)
- [Query Matters: How Selection Strategies Influence Active Learning in Drug Discovery, *JCIM*](https://pubs.acs.org/jcisd8/article/66/6/3288/5080393/Query-Matters-How-Selection-Strategies-Influence)

**Parabilis Medicines could not be verified.** Two searches returned nothing on the
company or its helical-peptide platform; the row above reflects the source table only.
