# GPU docking worker

A Runpod Serverless worker that docks a batch of molecules against a receptor.

The platform runs on CPU and stays up. This worker is rented per job: it wakes, docks
a batch, returns scores, and shuts down. **Nothing bills while idle**, which is the
whole reason docking a million compounds costs about $8 rather than a monthly GPU bill.

## Why this is build-your-own

I checked the Runpod Hub: **121 listings, none for chemistry.** Searching "boltz",
"protein" and "molecul" all return zero. The Hub is LLM, image, audio and video. So
there is nothing to fork here — but the worker itself is small, and this is all of it.

## Build and deploy

```bash
docker build --platform linux/amd64 -t <your-registry>/dmta-docking:1 worker/
```

```bash
docker push <your-registry>/dmta-docking:1
```

Then create a Serverless endpoint from that image. Two settings people get wrong:

- **Set a `timeout`.** Runpod's v2 API returns `"timeout": 0` when it is omitted, and a
  zero timeout will bite you on a long batch.
- **Pin the CUDA floor** (`gpu.minCudaVersion`, e.g. `"12.8"`) so the image cannot land
  on a host too old to run it.

Keep max workers low until you have watched one real batch end to end.

## Calling it

```json
{
  "input": {
    "receptor_pdb": "<contents of data/structures/2RGP_receptor.pdb>",
    "reference_pdb": "<contents of data/structures/2RGP_ref_ligand.pdb>",
    "compounds": [{"id": "CHEMBL1234", "smiles": "CCOc1cc2ncnc(Nc3ccccc3)c2cc1"}],
    "exhaustiveness": 8
  }
}
```

It returns `scores`, plus `failed_embed` and `failed_dock` as **separate lists**. That
separation matters more than it looks: a compound that could not be processed must never
be confused with a compound that docked badly, or the loop will quietly learn that
"awkward to build in 3D" means "poor binder".

It also returns which `engine` actually ran. If the GPU engine fails and it falls back
to CPU, you will see `"smina (unidock failed)"` rather than a mysteriously large bill.

## Sizing it

From the Runpod catalogue read on 2026-11-24, GPUs that suit this job:

| GPU | VRAM | Community $/hr | Secure $/hr |
|---|---|---|---|
| `NVIDIA GeForce RTX 4090` | 24 GB | 0.34 | 0.89 |
| `NVIDIA GeForce RTX 5090` | 32 GB | 0.69 | 1.39 |
| `NVIDIA L40S` | 48 GB | 0.79 | 1.09 |

Docking is not memory-hungry — a 4090 is plenty. Save the 48 GB cards for Boltz-2.

## Before you trust a single score from this

Run `dock_validate.py` on the structure first. The worker will happily return confident
numbers for a receptor prepared wrongly, and nothing in the response will look unusual.
The validation against compounds of known potency is the only thing standing between you
and a pipeline that is precisely, expensively wrong.
