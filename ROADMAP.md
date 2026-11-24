# Roadmap: foundation models, docking, and scale

What it takes to add the three capabilities, what each one costs, and the order to do
them in. GPU prices are live from Runpod's `list-gpu-types` catalogue read (Secure and
Community cloud, 2026-11-24).

**The rule that governs all of it:** each stage must prove it adds signal *on the
scaffold split* before you pay for the next one. The platform already measures that
([RESULTS.md](RESULTS.md)), so there is no excuse for buying GPU time on faith.

---

## A. Chemistry foundation model

### What it actually is

A large neural network pre-trained on hundreds of millions of molecules, which has
learned general "chemistry sense" and is then fine-tuned on your few thousand measured
compounds. The pitch: it brings knowledge your 9,000 data points cannot contain.

### What to use

| Option | What it is | Cost to try |
|---|---|---|
| **Chemprop (D-MPNN)** | Not a foundation model — a graph neural net trained from scratch on your data. The standard strong baseline. | **CPU, free.** ~20 min for 9k compounds |
| **ChemFM-3B** | 3-billion-parameter chemical foundation model; the authors report ~7% average improvement over prior state of the art across 34 property datasets | A LoRA fine-tune: **1–3 GPU-hours** |
| **MoLFormer / ChemBERTa** | Older, smaller transformer models pre-trained on SMILES strings. Cheap to fine-tune. | **<1 GPU-hour** |
| **TabPFN + molecular features** | A tabular foundation model used in-context; a 2026 preprint reports 50 wins across 58 tasks | **<1 GPU-hour** |

### Hardware and cost

Fine-tuning a 3B model with LoRA needs ~24 GB of VRAM; 80 GB is comfortable. From the
catalogue read, with current stock:

| GPU | VRAM | Community $/hr | Secure $/hr | Availability now |
|---|---|---|---|---|
| `NVIDIA GeForce RTX 4090` | 24 GB | **0.34** | 0.89 | LOW — EU-RO-1, EUR-IS-1, EUR-IS-2, US-IL-1 |
| `NVIDIA GeForce RTX 5090` | 32 GB | 0.69 | 1.39 | LOW — EUR-IS-1 |
| `NVIDIA L40S` | 48 GB | 0.79 | 1.09 | LOW — OC-AU-1, US-TX-4 |
| `NVIDIA A100 80GB PCIe` | 80 GB | **1.19** | 2.09 | LOW — CA-MTL-3 |
| `NVIDIA H100 80GB HBM3` | 80 GB | 2.69 | 4.39 | LOW — 7 data centres |

**A full fine-tune experiment costs about $2–$6.** This is not a budget decision; it is
an evidence decision.

### The honest caveat

The April 2026 benchmark survey found **four different models each winning on different
ADME endpoints, and rankings changing completely under time-based splitting** — the same
collapse we measured on EGFR. No foundation model fixes the extrapolation problem. Expect
a better scaffold-split ρ, not a different kind of model.

### What to do

1. Run Chemprop on CPU first. If it does not beat the random forest on your scaffold
   split, a foundation model probably will not either, and you have spent nothing.
2. If it does, rent one `RTX 4090` by the hour and LoRA fine-tune ChemFM on the same
   split. Keep whichever wins, by the measured number.

---

## B. Docking and physics

### What it actually is

Everything so far judges a molecule by its 2D structure. Docking puts a 3D model of the
molecule into a 3D model of the protein and scores the physical fit. It is the one thing
that can rank a chemotype the model has never seen — exactly where our time split failed.

### The three tiers

| Tier | Tool | Speed | What it buys |
|---|---|---|---|
| **Fast docking** | **Uni-Dock** (GPU) | Developers report **>1M compounds/day on one GPU**, and >1000× a single Vina CPU core | Triage: throw away the obvious non-binders |
| **Accurate docking** | **GNINA** (CNN rescoring) | Much slower — one benchmark: ~58 min vs ~8.6 min for AutoDock-GPU on a T4 | Better poses on the shortlist |
| **Affinity prediction** | **Boltz-2**, or **Boltzina** | Boltz-2: **~20 s per compound per GPU**. Boltzina skips structure prediction and uses Vina poses, far faster | A real binding-affinity estimate |

### Hardware and cost, from the same catalogue read

- **Screening 1 million compounds with Uni-Dock:** roughly one GPU-day. On
  `NVIDIA GeForce RTX 4090` at **$0.34/hr community**, that is about **$8**.
- **Boltz-2 on a 1,000-compound shortlist:** ~5.5 GPU-hours at 20 s each. On
  `NVIDIA L40S` (48 GB) at **$0.79/hr community**, about **$4.50**.
- Our current 924-compound library would dock in **minutes**.

Docking is far cheaper than people assume. The expensive part is not the GPU.

### What you must supply

A 3D structure of the protein. For EGFR there are many experimental structures in the
PDB; for a novel target you predict one (Boltz-2 or AlphaFold) and accept the extra
uncertainty. **Choosing the right structure — which conformation, which mutation, where
the pocket is — is a medicinal-chemistry judgement call, not a compute problem.** This is
the step where the project needs a structural chemist, not a bigger GPU.

### What to do

1. Pick and prepare an EGFR structure from the PDB.
2. Run Uni-Dock over the existing 8,978-compound ChEMBL set — a few minutes of GPU time.
3. **Test whether the docking score actually ranks the real potencies.** If it does, add
   it as a ninth objective in the scoring profile. If it does not, you have learned that
   for $10 instead of finding out after building a pipeline around it.

---

## C. Scale to millions (and billions)

### The thing nobody tells you

**You do not dock a billion molecules.** The published billion-scale methods all avoid it:

| Approach | Scale reached | How |
|---|---|---|
| **V-SYNTHES2** | 36 billion | Docks fragments ("synthons") first, builds up only the winners — needed only ~3.8M actual dockings |
| **V-SYNTHES / xREAL** | 173 billion | 16× bigger space for only **10–15% more compute** |
| **AdaptiveFlow (ATG-VS)** | 69 billion ready-to-dock | ML-guided adaptive screening; authors report **up to 1,000-fold cost reduction** versus exhaustive search |
| **HIDDEN GEM** | 37 billion | Benchmarked on 16 targets; authors claim the highest enrichment among accelerated methods |

All of them are the same idea the platform already runs: **a cheap model triages, an
expensive method confirms, and results feed back.** Scaling is not a new architecture —
it is the existing loop with a bigger candidate pool and a cheaper first filter.

### What has to change in our code

| Now | At a million | At a billion |
|---|---|---|
| Python list of 924 compounds in memory | Parquet on disk + batched scoring | Synthon-based: never enumerate the full space |
| Fingerprint every candidate each cycle | Precompute once, memory-map | Fingerprint the synthons, not the products |
| Score the full library every cycle | Score a sampled candidate pool | Fragment-first search, build up winners |
| SQLite | SQLite still fine (it only stores what you *made*) | Unchanged — you still only make hundreds |

Fingerprinting a billion molecules with RDKit is roughly **280 CPU-core-hours** — a few
hours on a 64-core box, and a one-off. The model inference is the part to move to GPU.

### The commercial dependency

Enamine REAL is a commercial catalogue (~78 billion in their navigator; BioSolveIT report
**>80% synthesis success in 3–4 weeks**). You will need a licence or a supplier
relationship. The open alternative is to enumerate your own space from purchasable
building blocks — which is exactly what `dmta/library.py` already does, just smaller.

---

## The order to do it in

| Stage | Spend | What it proves |
|---|---|---|
| **0. Your own data** | $0 | Whether a loop is worth starting at all. The scaffold-split ρ decides it |
| **1. Chemprop on CPU** | $0 | Whether a neural model beats the random forest here |
| **2. Foundation model fine-tune** | ~$5 | Whether pre-trained chemistry knowledge adds anything on your target |
| **3. Docking on the existing library** | ~$10 | Whether 3D adds signal 2D does not — the only known route past the time-split failure |
| **4. Scale the candidate pool** | $50–200/campaign | Only worth it once 1–3 have proven out |
| **5. Boltz-2 on shortlists** | ~$5/1,000 | Final ranking before committing lab money |

Stages 0–3 together cost **under $20 of GPU time** and answer every important question.
Stage 4 is the only one that needs real budget, and by then you will know whether to
spend it.

### What each stage needs from you, not from compute

- **Stage 0:** your assay data, or a decision on which public target to stand in for it.
- **Stage 3:** a structural chemist's call on the protein structure and binding site.
- **Stage 4:** an Enamine (or equivalent) catalogue relationship.
- **Throughout:** a chemist who will look at the proposed batch and say "no, not that
  one" — the published practitioner guidance is that batch size should match **reviewer
  capacity**, not what the algorithm can generate.

---

## Sources

- [Benchmark survey of molecular property prediction in the foundation-model era (ChemRxiv, Apr 2026)](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.15002208)
- [ChemFM, *Communications Chemistry*](https://www.nature.com/articles/s42004-025-01793-8)
- [Tabular foundation models for molecular properties (arXiv 2026)](https://arxiv.org/pdf/2604.16123)
- [Uni-Dock: GPU-accelerated ultralarge virtual screening](https://www.researchgate.net/publication/370306163_Uni-Dock_GPU-Accelerated_Docking_Enables_Ultralarge_Virtual_Screening)
- [Boltzina: docking-guided binding prediction with Boltz-2 (arXiv)](https://arxiv.org/abs/2508.17555)
- [Benchmarking DiffDock, AutoDock-GPU, GNINA on LIT-PCBA (arXiv)](https://arxiv.org/pdf/2605.01681)
- [V-SYNTHES2, *Nature*](https://www.nature.com/articles/s44386-026-00053-6) · [V-SYNTHES at 173B](https://www.nature.com/articles/s44386-026-00067-0)
- [AI-enhanced adaptive virtual screening, *Nature Biotechnology*](https://www.nature.com/articles/s41587-026-03217-x)
- [Enamine REAL Space Navigator](https://enamine.net/library-synthesis/real-compounds/real-space-navigator) · [BioSolveIT REAL Space](https://www.biosolveit.de/chemical-spaces/real-space/)
- GPU prices and availability: Runpod `list-gpu-types` catalogue read, 2026-11-24
