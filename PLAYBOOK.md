# Playbook: turning this into something that works

Written after measuring what this platform can and cannot do. Every claim below is
backed by a number in [RESULTS.md](RESULTS.md), not by ambition.

---

## 1. First, be honest about what the evidence supports

Three things are now measured facts about this system:

| What we tested | Result | What it means |
|---|---|---|
| Ranking compounds in known chemotypes | **ρ 0.815, 3.51× enrichment** | Genuinely strong |
| Choosing what to make next (the loop) | **2.37× hit rate vs random** | Genuinely useful |
| Predicting compounds from *unseen* chemotypes | **ρ 0.11, R² −0.66** | Does not work |
| 3D docking to rescue that | **0.0–0.37×, three attempts** | Does not work |

So the honest product claim is **not** "our AI discovers new drugs."

It is:

> **"Give us an active chemistry programme. We will get you to your best compound in
> fewer experiments than you would otherwise need — and we will prove it on your own
> historical data before you pay us anything."**

That claim is narrower, defensible, falsifiable, and supported by every number above.
The broader claim is what gets a company killed eighteen months in, when the compounds
it promised do not work.

**Optimisation within a series, not invention of new series.** Say it plainly. Many
buyers want exactly that, because it is where their money actually burns.

---

## 2. The one move that matters: get real data

Everything built so far runs on public ChEMBL. Public data is a demo, not a product.
It has mixed assays, mixed labs, mixed protein variants, and every competitor has it.

**A pharma or biotech partner's own assay data is both the moat and the barrier.**

You do not need much. One active programme with 300–1,000 compounds and their measured
results is enough to run the entire retrospective test.

How to get it without asking anyone to trust you first:

- Offer the retrospective test **free**, under an NDA, on data they already have.
- They send a CSV: compound SMILES, measured value, date.
- You return a report: *"Of your last 400 compounds, our loop would have found your
  top 10 after 120 experiments instead of 400."*
- They can check every claim against what actually happened. Nothing to take on faith.

This is the whole sales motion, and **you already have the machinery** —
`chembl_benchmark.py` does exactly this; it just needs to accept their CSV instead of
ChEMBL, and produce a report rather than terminal output.

---

## 3. What to build next, in order

### Priority 1 — Turn the retrospective test into the product (1–2 weeks)

Today it is a script that prints to a terminal. It needs to be:

- **Upload a CSV** → compound id, SMILES, result, date
- **Automatic data hygiene** → the cleaning in `chembl.py` already handles salts,
  replicate disagreement, censored values. Report what was dropped and why.
- **The three splits** → random, scaffold, time. Show all three. **Especially show the
  time split**, even when it fails. A vendor who shows you their own failure mode is
  the one you trust.
- **The retrospective loop** → "we would have found your best N compounds in X
  experiments instead of Y"
- **A PDF or shareable page** they can take to their boss

This is the smallest thing that is genuinely saleable.

### Priority 2 — ADMET, not potency (2–3 weeks)

Real programmes rarely die from lack of potency. They die from solubility, metabolism,
toxicity, hERG. The multi-objective scoring framework in `scoring.py` is already built
for this — it currently scores MW, cLogP, TPSA, QED, SA.

Add real ADMET predictions as further objectives. Two routes:

- **Train your own** on public TDC / ChEMBL ADMET endpoints (free, CPU)
- **Buy them in** — Inductive Bio exposes ADMET models over MCP, SandboxAQ exposes
  potency screening. Both connect in minutes.

**And grade the bought-in ones against known answers before trusting them**, exactly as
we graded docking. That discipline is now a reusable asset, not a one-off.

This is almost certainly worth more than docking ever was, and it is a far smaller lift.

### Priority 3 — Chemprop (3–4 days)

A graph neural network to replace the random forest. Free, CPU, no GPU. It has one job:
beat ρ 0.815 and 3.51× on the scaffold split. If it does not, keep the random forest and
stop paying attention to model architecture.

### Priority 4 — Only now, scale (weeks, with budget)

Bigger libraries, Enamine REAL, GPU inference, foundation models. **None of this is
worth starting until a real customer's data has gone through priorities 1–3.** Scale
multiplies whatever you have; multiplying something unproven just produces more of it.

---

## 4. What to stop doing

- **Docking.** Three structures, no signal, decision recorded. Revisit only if a
  customer's own target gives a different answer on the same test.
- **GPU spend.** Nothing in the measured results currently needs a GPU.
- **Foundation models.** Interesting; not the bottleneck. The bottleneck is data.
- **Chasing "new chemotype" discovery.** The time split says plainly that neither the
  model nor docking can do this. Selling it would be selling the wrong town.

---

## 5. The rule that makes the whole thing work

> **No feature ships until it has beaten a test whose answers you already knew.**

This is not a process nicety. It is the product's only real differentiator, and it has
already earned its keep once:

- Docking passed its smoke test (lapatinib −11.03, ethanol −2.85)
- It would have demoed beautifully — rotatable 3D protein, molecules in the pocket
- It scored **0.37× — worse than random**
- Cost of finding out: **70 minutes of CPU and $0**
- Cost of *not* finding out: a year of lab work steered towards weaker compounds

Write the rule on the wall. Every vendor model, every new structure, every new
architecture goes through it. When something fails, say so publicly — in a field full
of overclaiming, being the team that publishes its own negative results is a position
nobody else wants and everybody trusts.

---

## 6. 30 / 60 / 90 days

**Days 1–30 — one real dataset.**
Sign an NDA with one partner. Get one programme's compounds and results. Run the three
splits and the retrospective loop. Decide, on the number, whether the loop helps on
*their* chemistry. Build the CSV-in, report-out flow while you wait for the data.

**Days 31–60 — make the report a product.**
Shareable report. ADMET objectives added and graded. Chemprop benchmarked against the
random forest. Second partner dataset, so you know whether result one was luck.

**Days 61–90 — one prospective cycle.**
The real test: your loop picks a batch, the partner's chemists make it, the assay comes
back. Retrospective success is encouraging; prospective success is the business.
Expect worse numbers than retrospective — every honest group reports this.

---

## 7. What can still go wrong, said out loud

- **Retrospective is easier than prospective.** Every compound in a historical set was
  made by someone who thought it worth making. Your real-world numbers will be lower.
- **Data access may simply not happen.** No partner, no product. This is the single
  biggest risk, and it is commercial, not technical.
- **One target is not a platform.** EGFR worked. Two or three more programmes are needed
  before "it works" is a claim rather than an anecdote.
- **Chemists must stay in the loop.** Batch size should match what a chemist can review
  in a week, not what the algorithm can generate. A system that overrules its users
  gets switched off, however good it is.

---

## 8. The one-line version

**Get one real dataset. Prove the loop on it against answers that already exist. Sell
that proof, not the promise. Let the examiner kill every feature that cannot pass it —
including the impressive ones.**
