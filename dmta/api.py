"""HTTP service: the platform as an API, with the UI served alongside it.

    uvicorn dmta.api:app --reload --port 8000
"""
from __future__ import annotations

import io
from functools import lru_cache
from pathlib import Path

from fastapi import Body, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .assay import SimulatedAssay, write_batch_csv
from .campaign import Campaign, CycleSettings, fresh_library
from .library import CATALOGUE, REACTIONS, SCAFFOLDS
from .scoring import CHEMBL_PROFILE, profile_as_dicts
from .select import ACQUISITIONS, diversity_report
from .store import Store

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "campaigns.db"

app = FastAPI(title="DMTA platform", version="1.0",
              description="Design–Make–Test loop for medicinal chemistry")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

store = Store(DB)
_live: dict[str, Campaign] = {}


def _build(name: str, target: str, settings: CycleSettings) -> Campaign:
    """A campaign over a ChEMBL target replays real measurements; otherwise it
    runs over the enumerated make-on-demand library against a simulated assay."""
    if str(target).startswith("CHEMBL"):
        from .chembl import load_compounds
        library = load_compounds(str(target).split()[0])
        return Campaign(name=name, target=target, settings=settings, store=store,
                        library=library, profile=list(CHEMBL_PROFILE))
    return Campaign(name=name, target=target, settings=settings, store=store)


def get_campaign(name: str) -> Campaign:
    if name not in _live:
        row = store.get_campaign(name)
        if row is None:
            raise HTTPException(404, f"no campaign named {name!r}")
        _live[name] = _build(name, row["target"], CycleSettings(**row["settings"]))
    return _live[name]


# ------------------------------------------------------------------ models

class NewCampaign(BaseModel):
    name: str = Field(..., min_length=1, max_length=64)
    target: str = "simulated kinase"
    batch_size: int = 16
    acquisition: str = "ucb"
    beta: float = 0.25
    diversity: float = 0.15
    budget: float | None = None
    max_per_scaffold: int | None = None
    assay_noise: float = 0.2
    oracle_seed: int = 42
    hit_threshold: float = 7.5


class Results(BaseModel):
    """Measurements coming back from a real lab."""
    endpoint: str = "pIC50"
    values: dict[str, float]


# ------------------------------------------------------------------- routes

@app.get("/api/health")
def health():
    return {"ok": True, "library_size": len(fresh_library()), "acquisitions": list(ACQUISITIONS)}


@app.get("/api/chemistry")
def chemistry():
    """What the platform can actually make."""
    return {
        "reactions": [{"key": k, "name": r.name, "smarts": r.smarts,
                       "typical_yield": r.typical_yield, "reliability": r.reliability}
                      for k, r in REACTIONS.items()],
        "scaffolds": [s.__dict__ for s in SCAFFOLDS],
        "reagents": {role: [r.__dict__ for r in rs] for role, rs in CATALOGUE.items()},
        "profile": profile_as_dicts(),
    }


@app.get("/api/targets")
def targets():
    from .chembl import SUGGESTED_TARGETS, cache_path
    return [{"id": k, "name": v, "cached": cache_path(k).exists()}
            for k, v in SUGGESTED_TARGETS.items()]


# ------------------------------------------------- structures & 3D viewing

@app.get("/api/structures/{pdb_id}")
def structure_detail(pdb_id: str):
    """Everything the 3D viewer needs: the manifest plus the atoms themselves."""
    from .structure import prepare
    try:
        s = prepare(pdb_id)
    except Exception as exc:
        raise HTTPException(404, f"could not prepare {pdb_id}: {exc}")
    d = s.as_dict()
    d["receptor_data"] = Path(s.receptor_pdb).read_text()
    d["reference_ligand_data"] = Path(s.reference_ligand_pdb).read_text()
    return d


def _pose_facts(pdb_id: str):
    """Scores and measurements read from what was actually computed on disk —
    not from whatever a live campaign happens to hold in memory."""
    from .docking import _load_cache
    cache = _load_cache()
    scores: dict[str, float] = {}
    prefix = f"{pdb_id.upper()}|"
    for key, value in cache.items():
        if key.startswith(prefix):
            scores[key.rsplit("|", 1)[-1]] = value          # inchikey -> score
    return scores


@lru_cache(maxsize=8)
def _chembl_index(target: str = "CHEMBL203") -> dict[str, dict]:
    """inchikey -> {id, smiles, pchembl} straight from the cached CSV.

    Deliberately does NOT build Compound objects: that recomputes descriptors for
    every row and turned a viewer page load into a minute of CPU.
    """
    import csv as _csv
    from .chembl import cache_path
    path = cache_path(target)
    if not path.exists():
        return {}
    out = {}
    with open(path, newline="") as fh:
        for row in _csv.DictReader(fh):
            out[row["inchikey"]] = {"id": row["chembl_id"], "smiles": row["smiles"],
                                    "pchembl": float(row["pchembl"])}
    return out


@app.get("/api/structures/{pdb_id}/pose/{compound_id}")
def compound_pose(pdb_id: str, compound_id: str, campaign: str | None = None):
    """The docked 3D pose of one compound, if it has been docked."""
    from .docking import pose_path
    comp = None
    if campaign:
        comp = get_campaign(campaign).by_id.get(compound_id)
    if comp is None:
        for c in _live.values():
            comp = c.by_id.get(compound_id)
            if comp is not None:
                break
    record = None
    if comp is None:                        # fall back to the reference dataset
        record = next((r for r in _chembl_index().values() if r["id"] == compound_id), None)
        if record is None:
            raise HTTPException(404, f"unknown compound {compound_id}")
    inchikey = comp.inchikey if comp is not None else next(
        k for k, r in _chembl_index().items() if r["id"] == compound_id)
    path = pose_path(pdb_id, inchikey)
    if path is None:
        raise HTTPException(404, f"{compound_id} has not been docked against {pdb_id}")
    from .chem import descriptors
    facts = _pose_facts(pdb_id)
    smiles = comp.smiles if comp is not None else record["smiles"]
    measured = (comp.measured if comp is not None and comp.measured is not None
                else (getattr(comp, "truth", None) if comp is not None else record["pchembl"]))
    return {
        "compound_id": compound_id,
        "smiles": smiles,
        "docking_score": (comp.docking_score if comp is not None and comp.docking_score is not None
                          else facts.get(inchikey)),
        "measured": measured,
        "predicted": comp.predicted if comp is not None else None,
        "desc": comp.desc if comp is not None else descriptors(smiles),
        "pose_sdf": path.read_text(),
    }


@app.get("/api/structures/{pdb_id}/docked")
def docked_compounds(pdb_id: str, campaign: str | None = None, limit: int = 200):
    """Which compounds have a pose available to view."""
    from .docking import POSES
    folder = POSES / pdb_id.upper()
    if not folder.exists():
        return []
    keys = {f.stem for f in folder.glob("*.sdf")}
    facts = _pose_facts(pdb_id)
    index = _chembl_index()
    pool = []
    for key in keys:
        record = index.get(key)
        if record is None:
            continue
        pool.append({"id": record["id"], "smiles": record["smiles"],
                     "docking_score": facts.get(key), "measured": record["pchembl"]})
    pool.sort(key=lambda r: (r["docking_score"] is None, r["docking_score"] or 0))
    return pool[:limit]


@app.get("/api/campaigns")
def list_campaigns():
    return store.list_campaigns()


@app.post("/api/campaigns", status_code=201)
def create_campaign(body: NewCampaign):
    if store.get_campaign(body.name):
        raise HTTPException(409, f"campaign {body.name!r} already exists")
    settings = CycleSettings(
        batch_size=body.batch_size, acquisition=body.acquisition, beta=body.beta,
        diversity=body.diversity, budget=body.budget, max_per_scaffold=body.max_per_scaffold,
        assay_noise=body.assay_noise, oracle_seed=body.oracle_seed, hit_threshold=body.hit_threshold)
    if body.acquisition not in ACQUISITIONS:
        raise HTTPException(400, f"acquisition must be one of {ACQUISITIONS}")
    c = _build(body.name, body.target, settings)
    _live[body.name] = c
    return c.state()


@app.delete("/api/campaigns/{name}", status_code=204)
def delete_campaign(name: str):
    row = store.get_campaign(name)
    if not row:
        raise HTTPException(404, "no such campaign")
    store.delete_campaign(row["id"])
    _live.pop(name, None)


@app.get("/api/campaigns/{name}/state")
def state(name: str):
    c = get_campaign(name)
    s = c.state()
    if hasattr(c.assay, "ceiling"):
        s["oracle"] = {
            "ceiling": round(c.assay.ceiling(), 3),
            "top5pct_threshold": c.assay.threshold_at(0.05),
            "note": ("real ChEMBL measurements, revealed only when a compound is made"
                     if c.assay.name == "retrospective"
                     else "simulated landscape — for benchmarking only"),
        }
    s["data_source"] = c.assay.name
    return s


@app.post("/api/campaigns/{name}/propose")
def propose(name: str, size: int | None = None):
    """DESIGN + MAKE: the next batch and its synthesis plan, without testing it."""
    c = get_campaign(name)
    batch = c.propose(size)
    return {
        "cycle": c.cycle,
        "plan": c.plan(batch),
        "diversity": diversity_report(batch),
        "compounds": [x.as_dict() for x in batch],
    }


@app.post("/api/campaigns/{name}/cycle")
def run_cycle(name: str, size: int | None = None):
    """Full loop against the configured assay provider."""
    c = get_campaign(name)
    try:
        return c.run_cycle(size)
    except RuntimeError as e:
        raise HTTPException(409, str(e))


@app.post("/api/campaigns/{name}/measurements")
def submit(name: str, body: Results = Body(...)):
    """TEST with real data: hand back results for compounds the lab actually made."""
    c = get_campaign(name)
    unknown = [k for k in body.values if k not in c.by_id]
    if unknown:
        raise HTTPException(400, f"unknown compound ids: {unknown[:5]}")
    rows = []
    for cid, value in body.values.items():
        comp = c.by_id[cid]
        comp.measured = float(value)
        comp.cycle_made = c.cycle
        rows.append({"compound_id": cid, "smiles": comp.smiles, "value": float(value),
                     "endpoint": body.endpoint, "predicted": comp.predicted,
                     "uncertainty": comp.uncertainty})
    if c.store:
        c.store.add_measurements(c.campaign_id, None, rows, source="external")
    report = c._retrain()
    return {"accepted": len(rows), "tested_total": len(c.measured), "model": report}


@app.get("/api/campaigns/{name}/compounds")
def compounds(name: str, status: str = Query("tested", pattern="^(tested|untested|all)$"),
              limit: int = 50):
    c = get_campaign(name)
    pool = {"tested": c.measured, "untested": c.untested, "all": c.library}[status]
    key = (lambda x: x.measured) if status == "tested" else (lambda x: x.score or 0)
    return [x.as_dict() for x in sorted(pool, key=key, reverse=True)[:limit]]


@app.get("/api/campaigns/{name}/worklist.csv")
def worklist(name: str, size: int | None = None):
    """The CSV a chemistry team receives: what to make, how, and where to type results."""
    c = get_campaign(name)
    batch = c.propose(size)
    path = ROOT / f"worklist-{name}-cycle{c.cycle}.csv"
    write_batch_csv(batch, path)
    return FileResponse(path, media_type="text/csv", filename=path.name)


@app.get("/api/campaigns/{name}/export.csv")
def export(name: str):
    c = get_campaign(name)
    buf = io.StringIO()
    buf.write("compound_id,smiles,cycle,measured,predicted,uncertainty,score,mw,clogp,tpsa,qed,cost\n")
    for x in sorted(c.measured, key=lambda y: y.measured, reverse=True):
        buf.write(f"{x.id},{x.smiles},{x.cycle_made},{x.measured},{x.predicted},{x.uncertainty},"
                  f"{x.score},{x.desc['mw']},{x.desc['clogp']},{x.desc['tpsa']},{x.desc['qed']},{x.cost}\n")
    buf.seek(0)
    return StreamingResponse(buf, media_type="text/csv",
                             headers={"Content-Disposition": f'attachment; filename="{name}-results.csv"'})


# the UI is static files served from the same origin
UI = ROOT / "ui"
if UI.exists():
    app.mount("/", StaticFiles(directory=UI, html=True), name="ui")
