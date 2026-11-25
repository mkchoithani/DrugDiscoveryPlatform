"""The checks that matter: valid chemistry, honest models, working persistence."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from dmta.assay import CsvAssay, SimulatedAssay, write_batch_csv
from dmta.campaign import Campaign, CycleSettings, fresh_library
from dmta.chem import alerts, descriptors, feature_vector, fingerprint, tanimoto
from dmta.library import REACTIONS, enumerate_library, group_by_intermediate
from dmta.models import PropertyModel, regression_metrics
from dmta.select import select_batch
from dmta.store import Store


# ----------------------------------------------------------------- chemistry

def test_descriptors_match_known_values():
    caffeine = descriptors("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    assert caffeine["mw"] == pytest.approx(194.19, abs=0.05)
    assert caffeine["clogp"] == pytest.approx(-1.03, abs=0.05)
    aspirin = descriptors("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin["mw"] == pytest.approx(180.16, abs=0.05)
    assert aspirin["hbd"] == 1


def test_every_enumerated_product_is_a_valid_molecule():
    lib = enumerate_library()
    assert len(lib) > 500
    for c in lib[:200]:
        assert Chem.MolFromSmiles(c.smiles) is not None


def test_products_are_unique():
    lib = enumerate_library()
    assert len({c.inchikey for c in lib}) == len(lib)


def test_reactions_actually_transform():
    acid = Chem.MolFromSmiles("OC(=O)c1ccccc1")
    amine = Chem.MolFromSmiles("C1CCNCC1")
    product = REACTIONS["amide_coupling"].run(acid, amine)
    assert product is not None
    assert Chem.MolToSmiles(product) == Chem.MolToSmiles(Chem.MolFromSmiles("O=C(c1ccccc1)N1CCCCC1"))


def test_amide_coupling_refuses_an_amide_nitrogen():
    acid = Chem.MolFromSmiles("OC(=O)c1ccccc1")
    amide = Chem.MolFromSmiles("CC(=O)NC")        # already acylated — not a nucleophile here
    assert REACTIONS["amide_coupling"].run(acid, amide) is None


def test_structural_alerts_fire():
    assert "epoxide" in alerts("C1OC1")
    assert alerts("CCO") == []


def test_library_carries_a_route_for_every_compound():
    for c in enumerate_library()[:100]:
        assert len(c.route) == 2
        assert c.route[0].product_smiles == c.intermediate_smiles
        assert 0 < c.overall_yield < 1
        assert c.cost > 0


# ---------------------------------------------------------------- similarity

def test_tanimoto_bounds():
    a, b = fingerprint("c1ccccc1"), fingerprint("c1ccccc1")
    assert tanimoto(a, b) == pytest.approx(1.0)
    assert 0 <= tanimoto(fingerprint("CCO"), fingerprint("c1ccncc1C(=O)N")) < 0.3


# -------------------------------------------------------------------- models

def test_model_learns_a_known_signal():
    """Train on cLogP — a property fully determined by structure. A model that
    cannot recover this is broken, regardless of what it does on assay data."""
    lib = fresh_library()[:300]
    for c in lib:
        c.measured = c.desc["clogp"]
    m = PropertyModel(n_estimators=120).fit(lib)
    assert m.cv_r2 is not None and m.cv_r2 > 0.5
    m.apply(lib[:50])
    assert all(c.predicted is not None and c.uncertainty > 0 for c in lib[:50])


def test_uncertainty_is_higher_off_the_training_set():
    lib = fresh_library()
    train = [c for c in lib if c.scaffold_id == "S1"][:80]
    for c in train:
        c.measured = c.desc["clogp"]
    m = PropertyModel(n_estimators=120).fit(train)
    far = [c for c in lib if c.scaffold_id == "S4"][:80]
    _, sigma_near = m.predict(train)
    _, sigma_far = m.predict(far)
    assert sigma_far.mean() > sigma_near.mean()


def test_regression_metrics():
    r = regression_metrics([1, 2, 3, 4], [1, 2, 3, 4])
    assert r["r2"] == 1.0 and r["mae"] == 0.0 and r["spearman"] == 1.0


# ----------------------------------------------------------------- selection

def test_batch_respects_size_and_uniqueness():
    lib = fresh_library()
    for c in lib:
        c.score = float(np.random.default_rng(0).random())
    batch = select_batch(lib, size=12, acquisition="greedy")
    assert len(batch) == 12
    assert len({c.id for c in batch}) == 12


def test_budget_constraint_is_enforced():
    lib = fresh_library()
    for c in lib:
        c.score = 0.5
    batch = select_batch(lib, size=20, acquisition="greedy", budget=150.0)
    assert sum(c.cost for c in batch) <= 150.0


def test_diversity_pressure_reduces_similarity():
    lib = fresh_library()
    for c in lib:
        c.score = 1.0 - c.desc["mw"] / 1000     # a smooth, clusterable objective
    tight = select_batch(lib, size=10, acquisition="greedy", diversity=0.0)
    spread = select_batch(lib, size=10, acquisition="greedy", diversity=0.5)
    from dmta.select import diversity_report
    assert diversity_report(spread)["mean_pairwise_tanimoto"] < diversity_report(tight)["mean_pairwise_tanimoto"]


def test_shared_intermediates_save_steps():
    lib = fresh_library()
    batch = [c for c in lib if c.intermediate_key == lib[0].intermediate_key][:6]
    plan = group_by_intermediate(batch)
    assert plan["shared_steps"] < plan["naive_steps"]
    assert plan["steps_saved"] == plan["naive_steps"] - plan["shared_steps"]


# ------------------------------------------------------------------ campaign

def test_cycle_makes_progress_and_records_history():
    c = Campaign("t-cycle", settings=CycleSettings(batch_size=8, oracle_seed=42))
    s0 = c.run_cycle()
    s1 = c.run_cycle()
    assert s0["cycle"] == 0 and s1["cycle"] == 1
    assert s1["tested_total"] == 16
    assert s1["prospective"] is not None       # second batch was model-chosen
    assert s0["prospective"] is None           # first batch could not be
    assert len({x.id for x in c.measured}) == 16


def test_no_compound_is_made_twice():
    c = Campaign("t-dupe", settings=CycleSettings(batch_size=10, oracle_seed=7))
    for _ in range(4):
        c.run_cycle()
    ids = [x.id for x in c.measured]
    assert len(ids) == len(set(ids)) == 40


def test_active_learning_beats_random_on_the_same_landscape():
    def found(strategy, diversity):
        c = Campaign(f"t-{strategy}", settings=CycleSettings(
            batch_size=16, acquisition=strategy, diversity=diversity, oracle_seed=42))
        target = c.assay.top_fraction(0.05)
        for _ in range(5):
            c.run_cycle()
        return sum(1 for x in c.measured if x.id in target)

    assert found("ucb", 0.15) > found("random", 0.0)


def test_persistence_round_trip():
    with tempfile.TemporaryDirectory() as tmp:
        db = Path(tmp) / "t.db"
        store = Store(db)
        a = Campaign("persist", settings=CycleSettings(batch_size=8), store=store)
        a.run_cycle()
        a.run_cycle()
        before = {x.id: x.measured for x in a.measured}

        b = Campaign("persist", settings=CycleSettings(batch_size=8), store=Store(db))
        after = {x.id: x.measured for x in b.measured}
        assert before == after
        assert b.cycle == 2
        assert b.model is not None


# ---------------------------------------------------------------- lab handoff

def test_worklist_csv_round_trip():
    with tempfile.TemporaryDirectory() as tmp:
        c = Campaign("t-csv", settings=CycleSettings(batch_size=6))
        batch = c.propose()
        path = write_batch_csv(batch, Path(tmp) / "work.csv")
        text = path.read_text().strip().splitlines()
        assert len(text) == 7                      # header + 6 compounds
        # lab fills the empty pIC50 column and hands it back
        filled = Path(tmp) / "done.csv"
        filled.write_text(text[0] + "\n" + "\n".join(
            row[:row.rfind(",") + 1] + "7.5" for row in text[1:]) + "\n")
        results = CsvAssay(filled).measure(batch)
        assert set(results) == {x.id for x in batch}
        assert all(v == 7.5 for v in results.values())


def test_csv_assay_rejects_missing_results():
    with tempfile.TemporaryDirectory() as tmp:
        c = Campaign("t-csv2", settings=CycleSettings(batch_size=4))
        batch = c.propose()
        path = Path(tmp) / "partial.csv"
        path.write_text("compound_id,pIC50\n" + f"{batch[0].id},7.0\n")
        with pytest.raises(ValueError):
            CsvAssay(path).measure(batch)


# --------------------------------------------------------------- ChEMBL data

from dmta.chembl import cache_path, dataset_summary, load_compounds, murcko, standardize
from dmta.validate import evaluate, random_split, scaffold_split, time_split

HAVE_CHEMBL = cache_path("CHEMBL203").exists()
needs_data = pytest.mark.skipif(not HAVE_CHEMBL, reason="run `dmta.cli chembl-fetch` first")


def test_standardize_strips_salts_and_charges():
    # sodium benzoate → benzoic acid
    assert standardize("[Na+].[O-]C(=O)c1ccccc1") == "O=C(O)c1ccccc1"
    assert standardize("CC") is None            # too small to be a compound
    assert standardize("not a smiles") is None


def test_murcko_strips_side_chains():
    assert murcko("CC(=O)Oc1ccccc1C(=O)O") == "c1ccccc1"
    assert murcko("CCO") == "acyclic"


@needs_data
def test_chembl_dataset_is_clean():
    lib = load_compounds("CHEMBL203")
    assert len(lib) > 1000
    assert len({c.inchikey for c in lib}) == len(lib)        # deduplicated
    for c in lib[:200]:
        assert Chem.MolFromSmiles(c.smiles) is not None
        assert 3.0 <= c.truth <= 12.0                        # plausible pChEMBL
        assert "." not in c.smiles                           # salts stripped


@needs_data
def test_splits_are_disjoint_and_scaffold_split_is_honest():
    lib = load_compounds("CHEMBL203")[:2000]
    tr, te = random_split(lib)
    assert not ({c.id for c in tr} & {c.id for c in te})
    tr, te = scaffold_split(lib)
    assert not ({c.scaffold_id for c in tr} & {c.scaffold_id for c in te})
    tr, te, cutoff = time_split(lib)
    assert max(c.year for c in tr) < cutoff <= min(c.year for c in te)


@needs_data
def test_model_has_real_signal_on_scaffold_split():
    """The claim the platform rests on: it can rank unseen chemotypes."""
    lib = load_compounds("CHEMBL203")[:3000]
    tr, te = scaffold_split(lib)
    m = evaluate(tr, te, n_estimators=100)
    assert m["spearman"] > 0.5
    assert m["enrichment"] > 1.5


@needs_data
def test_retrospective_loop_never_peeks():
    """A compound's measurement must not exist until the loop selects it."""
    from dmta.campaign import Campaign
    from dmta.scoring import CHEMBL_PROFILE
    lib = load_compounds("CHEMBL203")[:1500]
    c = Campaign("t-retro", library=lib, profile=list(CHEMBL_PROFILE),
                 settings=CycleSettings(batch_size=12, hit_threshold=8.0))
    assert all(x.measured is None for x in c.library)
    s = c.run_cycle()
    assert sum(1 for x in c.library if x.measured is not None) == 12
    assert s["prospective"] is None
    for x in c.measured:
        assert x.measured == x.truth        # real values, replayed faithfully


def test_replayed_history_matches_a_live_one():
    """A rehydrated campaign must expose the same history shape as a live one —
    missing keys here silently break every chart and report downstream."""
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as tmp:
        db = Path(tmp) / "h.db"
        live = Campaign("replay-shape", settings=CycleSettings(batch_size=8), store=Store(db))
        live.run_cycle(); live.run_cycle()
        reloaded = Campaign("replay-shape", settings=CycleSettings(batch_size=8), store=Store(db))
        assert len(reloaded.history) == len(live.history)
        for a, b in zip(live.history, reloaded.history):
            assert set(a) == set(b), "history keys differ between live and replayed"
            for key in ("cycle", "made", "tested_total", "best_measured", "hits_in_batch"):
                assert a[key] == b[key], f"{key}: {a[key]} != {b[key]}"
            assert all(v is not None for v in
                       (b["tested_total"], b["best_measured"], b["made"]))


CACHED_TARGETS = [t for t in ("CHEMBL203", "CHEMBL279") if cache_path(t).exists()]


@pytest.mark.parametrize("target", CACHED_TARGETS or ["none"])
def test_pipeline_holds_for_every_cached_target(target):
    """Whatever target is loaded, the dataset must satisfy the same invariants —
    a pipeline that only works on the target it was written against is useless."""
    if target == "none":
        pytest.skip("no cached ChEMBL datasets")
    lib = load_compounds(target)
    assert len(lib) > 500
    assert len({c.inchikey for c in lib}) == len(lib)
    assert len({c.scaffold_id for c in lib}) > 100       # Murcko actually varied
    assert all(4.0 <= c.truth <= 12.0 for c in lib)
    assert all(c.year is None or 1980 <= c.year <= 2030 for c in lib)


# ------------------------------------------------- uncertainty trustworthiness

from dmta.validate import recommend_beta, sigma_verdict


def test_sigma_verdict_boundaries():
    assert sigma_verdict(0.45)[0] == "usable"
    assert sigma_verdict(0.30)[0] == "usable"
    assert sigma_verdict(0.20)[0] == "weak"
    assert sigma_verdict(0.05)[0] == "unreliable"
    assert sigma_verdict(-0.13)[0] == "unreliable"      # measured on mu opioid
    assert sigma_verdict(None)[0] == "unknown"


def test_beta_recommendation_uses_the_worst_split_not_the_best():
    """mu opioid: scaffold σ is fine (0.36) but time σ is inverted (−0.13).
    Recommending on the scaffold split alone would tell the user to explore on an
    uncertainty estimate that is worse than useless in production."""
    report = {"scaffold": {"sigma_error_spearman": 0.361},
              "time": {"sigma_error_spearman": -0.129}}
    beta, why = recommend_beta(report)
    assert beta == 0.0
    assert "greedy" in why or "buys nothing" in why


def test_beta_recommendation_never_exceeds_what_was_measured():
    """The β sweep only ever supported values up to 0.25; the recommender must not
    invent confidence beyond the evidence."""
    generous = {"scaffold": {"sigma_error_spearman": 0.9},
                "time": {"sigma_error_spearman": 0.9}}
    beta, _ = recommend_beta(generous)
    assert beta <= 0.25
    assert recommend_beta({})[0] == 0.0                 # no evidence → no exploration


# ------------------------------------------------------- structure & docking

from dmta import docking as dk
from dmta.structure import DATA as STRUCT_DATA, list_ligands, pick_reference_ligand, prepare

HAVE_2RGP = (STRUCT_DATA / "2RGP.pdb").exists()
needs_structure = pytest.mark.skipif(not HAVE_2RGP, reason="2RGP not downloaded")
needs_engine = pytest.mark.skipif(not dk.engine_available(), reason="no docking engine installed")


@needs_structure
def test_reference_ligand_is_the_inhibitor_not_a_buffer_component():
    ref = pick_reference_ligand("2RGP")
    assert ref["resname"] == "HYZ"          # the co-crystallised inhibitor
    assert ref["atoms"] == 33
    names = {l["resname"] for l in list_ligands("2RGP")}
    assert {"PO4", "HOH"} & names           # junk is present in the file…
    assert ref["resname"] not in ("PO4", "HOH")   # …but never chosen as the site


@needs_structure
def test_preparation_strips_junk_and_records_it():
    s = prepare("2RGP")
    assert s.resolution == 2.0
    assert s.n_protein_atoms > 2000
    assert s.removed.get("HOH", 0) > 50     # waters gone, and counted
    assert s.removed.get("PO4", 0) == 10
    text = Path(s.receptor_pdb).read_text()
    assert "HOH" not in text and "HYZ" not in text
    assert (STRUCT_DATA / "2RGP_manifest.json").exists()


@needs_structure
def test_binding_box_surrounds_the_crystal_ligand():
    s = prepare("2RGP", padding=5.0)
    coords = []
    for line in Path(s.reference_ligand_pdb).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    for axis in range(3):
        values = [c[axis] for c in coords]
        lo, hi = s.site.center[axis] - s.site.size[axis] / 2, s.site.center[axis] + s.site.size[axis] / 2
        assert lo <= min(values) and max(values) <= hi, "ligand falls outside its own box"


def test_embed_3d_makes_a_real_conformer():
    mol = dk.embed_3d("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
    assert mol is not None and mol.GetNumConformers() == 1
    pos = mol.GetConformer().GetPositions()
    assert pos.shape[1] == 3
    assert np.ptp(pos[:, 2]) > 0.1            # genuinely 3D, not flat
    assert any(a.GetSymbol() == "H" for a in mol.GetAtoms())


def test_embed_3d_returns_none_for_rubbish():
    assert dk.embed_3d("this is not a molecule", "junk") is None


def test_missing_engine_raises_rather_than_silently_skipping():
    """A missing docking engine must be loud — a silent zero would poison the loop."""
    original = dk._find_engine
    dk._find_engine = lambda: None
    try:
        with pytest.raises(dk.DockingUnavailable):
            dk.require_engine()
    finally:
        dk._find_engine = original


@needs_engine
@needs_structure
def test_docking_separates_a_real_drug_from_a_non_binder():
    from dmta.chem import descriptors, inchikey
    from dmta.library import Compound

    def mk(cid, smi):
        return Compound(id=cid, smiles=smi, inchikey=inchikey(smi), scaffold_id="x",
                        intermediate_key=cid, intermediate_smiles=smi, route=[],
                        desc=descriptors(smi), alerts=[], overall_yield=1.0,
                        route_confidence=None, cost=0.0, lead_days=0)

    drug = mk("gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1")
    junk = mk("ethanol", "CCO")
    result = dk.dock([drug, junk], prepare("2RGP"), exhaustiveness=4)
    assert set(result.scores) == {"gefitinib", "ethanol"}
    assert result.scores["gefitinib"] < -6.0          # plausible binder
    assert result.scores["gefitinib"] < result.scores["ethanol"] - 2.0
    assert dk.apply_scores([drug, junk], result) == 2
    assert drug.docking_score == result.scores["gefitinib"]


# ---------------------------------------------------------- CSV ingest & API

from dmta.ingest import detect_units, load_csv


def test_ingest_handles_a_messy_real_export(tmp_path):
    """Wrong column names, nM units, salts, censored values, junk rows."""
    csv_path = tmp_path / "export.csv"
    csv_path.write_text(
        "Corporate ID,Structure,IC50 (nM),Date Tested\n"
        "A-1,CC(=O)Oc1ccccc1C(=O)O,120,2021-03-01\n"
        "A-2,[Na+].[O-]C(=O)c1ccccc1,45,2021-04-01\n"       # a salt
        "A-3,COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1,8,2021-05-01\n"
        "A-4,CCO,> 10000,2021-06-01\n"                       # censored
        "A-5,NOTASMILES,15,2021-07-01\n"                     # unparsable
        "A-6,CC(=O)Oc1ccccc1C(=O)O,,2021-08-01\n")           # no result
    compounds, report = load_csv(csv_path)

    assert report.columns["smiles"] == "Structure"
    assert report.columns["value"] == "IC50 (nM)"
    assert report.columns["date"] == "Date Tested"
    assert report.units == "nM"
    assert report.dropped.get("censored value (>, <, ~)") == 1
    assert report.dropped.get("no result recorded") == 1
    assert any("unreadable" in k or "too small" in k for k in report.dropped)
    assert all("." not in c.smiles for c in compounds)       # salts stripped
    assert all(c.year is not None for c in compounds)
    # 120 nM -> pIC50 6.92
    aspirin = next(c for c in compounds if "OC(=O)C" in c.smiles or "C(=O)O" in c.smiles)
    assert 6.0 < aspirin.truth < 8.0


def test_units_are_detected_from_name_and_from_values():
    assert detect_units("pIC50", [7.2, 8.1]) == "log"
    assert detect_units("IC50 (uM)", [1.2]) == "uM"
    assert detect_units("IC50 (nM)", [120.0]) == "nM"
    assert detect_units("activity", [7.1, 6.8, 8.2]) == "log"     # looks like a p-scale
    assert detect_units("activity", [1200.0, 45.0]) == "nM"       # looks like nanomolar


def test_ingest_refuses_a_file_it_cannot_use(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("name,colour\nthing,red\n")
    with pytest.raises(ValueError, match="no SMILES column"):
        load_csv(bad)


def test_api_routes_are_not_shadowed_by_the_static_mount():
    """The UI is mounted at '/', which swallows any route declared after it.
    This caught a real 405 on the upload endpoint."""
    from fastapi.testclient import TestClient
    import dmta.api as api
    client = TestClient(api.app)
    res = client.post("/api/evaluate/inspect", files={"file": ("x.csv", b"a,b\n1,2\n", "text/csv")})
    assert res.status_code != 405, "static mount is shadowing the API routes again"
    # inspect answers even when it recognises nothing: the page then asks the user
    # to pick the columns, rather than rejecting their file outright.
    assert res.status_code == 200
    assert set(res.json()["missing"]) == {"smiles", "value"}


def test_api_inspect_detects_columns():
    from fastapi.testclient import TestClient
    import dmta.api as api
    client = TestClient(api.app)
    body = (b"Corporate ID,Structure,IC50 (nM),Date Tested\n"
            b"A-1,CC(=O)Oc1ccccc1C(=O)O,120,2021-03-01\n")
    res = client.post("/api/evaluate/inspect", files={"file": ("export.csv", body, "text/csv")})
    assert res.status_code == 200
    data = res.json()
    assert data["detected"]["smiles"] == "Structure"
    assert data["detected"]["value"] == "IC50 (nM)"
    assert data["units"] == "nM"
    assert data["missing"] == []


# --- RMT-denoised linear model (PNAS 2019) ------------------------------------

def test_mp_edge_rejects_pure_noise():
    """The whole method rests on this: random data must yield no signal at all."""
    import numpy as np
    from dmta.rmt_model import clip_spectrum
    rng = np.random.default_rng(0)
    n, p = 100, 300
    Z = rng.standard_normal((n, p))
    Z = (Z - Z.mean(axis=0)) / Z.std(axis=0)
    eig = np.linalg.svd(Z, full_matrices=False)[1] ** 2 / n
    signal, noise, edge = clip_spectrum(eig, n, p)
    assert len(signal) == 0, f"{len(signal)} noise eigenvalues mistaken for signal"
    assert 0.5 < noise < 1.5
    assert edge > 1.0


def test_mp_edge_finds_a_planted_factor():
    import numpy as np
    from dmta.rmt_model import clip_spectrum
    rng = np.random.default_rng(1)
    n, p = 120, 400
    Z = rng.standard_normal((n, 1)) @ (rng.standard_normal((1, p)) * 0.8)
    Z = Z + rng.standard_normal((n, p))
    Z = (Z - Z.mean(axis=0)) / Z.std(axis=0)
    eig = np.linalg.svd(Z, full_matrices=False)[1] ** 2 / n
    signal, _, edge = clip_spectrum(eig, n, p)
    assert len(signal) >= 1
    assert signal[0] > edge


def test_clipped_spectrum_preserves_the_trace():
    """Clipping redistributes eigenvalues; it must not create or destroy variance."""
    import numpy as np
    from dmta.rmt_model import clip_spectrum
    rng = np.random.default_rng(2)
    n, p = 80, 250
    Z = rng.standard_normal((n, p)) + rng.standard_normal((n, 1)) * 0.7
    Z = (Z - Z.mean(axis=0)) / Z.std(axis=0)
    eig = np.linalg.svd(Z, full_matrices=False)[1] ** 2 / n
    signal, noise, _ = clip_spectrum(eig, n, p)
    assert abs(signal.sum() + noise * (p - len(signal)) - p) < 1e-6


def test_rmt_model_fits_and_predicts_every_compound():
    from dmta.library import enumerate_library
    from dmta.assay import SimulatedAssay
    from dmta.rmt_model import RmtModel
    lib = enumerate_library()[:160]
    assay = SimulatedAssay(lib, seed=3)
    train, test = lib[:120], lib[120:]
    for c, v in zip(train, assay.measure(train).values()):
        c.measured = v
    model = RmtModel().fit(train)
    mu, sigma = model.predict(test)
    assert len(mu) == len(test) == len(sigma)
    assert model.n_signal >= 1, "no signal found in real chemistry — suspicious"
    assert all(s > 0 for s in sigma)
    r = model.report()
    assert r["model"] == "rmt-linear-project" and r["n_features"] > 0
    assert r["mp_edge"] > 1.0, "the bulk edge collapsed — the variance fit diverged again"
    assert 0.0 < r["noise_level"] < 2.0


def test_rmt_beats_plain_least_squares_on_the_same_features():
    """The paper's central claim, as a test: inverting the raw covariance is worse."""
    import numpy as np
    from dmta.library import enumerate_library
    from dmta.assay import SimulatedAssay
    from dmta.chem import feature_vector
    from dmta.rmt_model import RmtModel, ols_baseline_metrics
    from dmta.models import regression_metrics
    lib = enumerate_library()[:200]
    assay = SimulatedAssay(lib, seed=4)
    train, test = lib[:150], lib[150:]
    measured = assay.measure(train)
    for c in train:
        c.measured = measured[c.id]
    held = assay.measure(test)
    truth = [held[c.id] for c in test]

    mu, _ = RmtModel().fit(train).predict(test)
    rmt = regression_metrics([float(v) for v in mu], truth)
    X = np.array([feature_vector(c.smiles) for c in train], dtype=float)
    y = np.array([c.measured for c in train], dtype=float)
    Xt = np.array([feature_vector(c.smiles) for c in test], dtype=float)
    ols = ols_baseline_metrics(X, y, Xt, np.array(truth))
    assert rmt["mae"] < ols["mae"], f"rmt {rmt['mae']} vs ols {ols['mae']}"


def test_rmt_handles_a_constant_feature_block():
    """Identical molecules make every fingerprint bit constant; must not divide by zero."""
    from dmta.library import enumerate_library
    from dmta.rmt_model import RmtModel
    lib = enumerate_library()[:12]
    same = [lib[0].__class__(**{**lib[0].__dict__}) for _ in range(10)]
    for i, c in enumerate(same):
        c.measured = 6.0 + i * 0.01
    model = RmtModel().fit(same)
    mu, sigma = model.predict(same[:3])
    assert len(mu) == 3 and all(s > 0 for s in sigma)


# --- DEL Simulator landscape (Digital Discovery 2026) -------------------------

def _write_del_pair(tmp_path, rows, affinities):
    """rows: (nsynthon_id, smiles, relative_fraction). affinities: one per row, in order."""
    lib = tmp_path / "library.csv"
    aff = tmp_path / "affinity_seh_512.csv"
    lib.write_text("nsynthon_id,smiles,relative_fraction\n"
                   + "".join(f"{a},{b},{c}\n" for a, b, c in rows))
    aff.write_text("pKd\n" + "".join(f"{v}\n" for v in affinities))
    return lib, aff


def test_del_join_is_positional_and_correct(tmp_path):
    """The affinity file has no identifier, so row i must map to library row i."""
    from dmta.del_landscape import load_del_landscape
    rows = [("1_1_1", "CCO", 1.0), ("1_1_2", "CCCO", 1.0), ("1_2_1", "c1ccccc1", 1.0)]
    lib, aff = _write_del_pair(tmp_path, rows, [5.0, 6.0, 7.0])
    got = {c.smiles: c.truth for c in load_del_landscape(lib, aff, n=None)}
    assert got == {"CCO": 5.0, "CCCO": 6.0, "c1ccccc1": 7.0}


def test_del_refuses_a_misaligned_pair(tmp_path):
    """A length mismatch means the two files came from different runs — never guess."""
    import pytest
    from dmta.del_landscape import load_del_landscape
    rows = [("1_1_1", "CCO", 1.0), ("1_1_2", "CCCO", 1.0)]
    lib, aff = _write_del_pair(tmp_path, rows, [5.0, 6.0, 7.0])
    with pytest.raises(ValueError, match="positional"):
        load_del_landscape(lib, aff, n=None)


def test_del_keeps_the_major_product_of_a_duplicated_molecule(tmp_path):
    """Yield noise emits the same molecule several times; keep the intended product."""
    from dmta.del_landscape import load_del_landscape
    rows = [("1_1_1", "CCO", 0.1), ("1_1_1", "CCO", 0.9), ("1_1_2", "CCCO", 1.0)]
    lib, aff = _write_del_pair(tmp_path, rows, [5.0, 8.0, 6.0])
    got = {c.smiles: c.truth for c in load_del_landscape(lib, aff, n=None)}
    assert got["CCO"] == 8.0, "kept a trace side product over the major one"
    assert len(got) == 2


def test_del_dense_sampling_is_denser_than_uniform(tmp_path):
    """The two sample shapes must actually differ, or the three-arm test is vacuous."""
    from dmta.del_landscape import load_del_landscape
    rows, affinities = [], []
    for c1 in range(8):
        for c3 in range(8):
            rows.append((f"{c1}_0_{c3}", f"{'C' * (c1 + 1)}{'O' * (c3 + 1)}", 1.0))
            affinities.append(4.0 + 0.1 * c1 + 0.01 * c3)
    lib, aff = _write_del_pair(tmp_path, rows, affinities)
    uniform = load_del_landscape(lib, aff, n=16, seed=0, sampling="uniform")
    dense = load_del_landscape(lib, aff, n=16, seed=0, sampling="dense")
    assert len({c.intermediate_key for c in dense}) < len({c.intermediate_key for c in uniform})


def test_noisy_replay_assay_keeps_truth_clean():
    from dmta.del_landscape import NoisyReplayAssay
    from dmta.library import enumerate_library
    lib = enumerate_library()[:20]
    for i, c in enumerate(lib):
        c.truth = 5.0 + i * 0.1
    assay = NoisyReplayAssay(lib, noise=0.5, seed=1)
    readings = assay.measure(lib)
    assert all(assay.truth(c.id) == c.truth for c in lib), "truth was perturbed"
    assert any(abs(readings[c.id] - c.truth) > 1e-6 for c in lib), "no noise was added"
    assert len(assay.top_fraction(0.25)) == 5


# --- predicted pharmacophore maps (Digital Discovery 2026) --------------------

def _synthetic_field(tmp_path, shape=(6, 6, 6), spacing=2.0, hot=None):
    """A field we control, so the reduction can be checked against known answers."""
    import numpy as np
    from dmta.pharmacophore import FEATURE_CHANNELS
    n = shape[0] * shape[1] * shape[2]
    probs = np.full((n, len(FEATURE_CHANNELS)), 0.01, dtype=np.float32)
    for flat_index, channel, value in (hot or []):
        probs[flat_index, channel] = value
    path = tmp_path / "phc4_field.npz"
    np.savez(path, probs=probs, origin=np.zeros(3), shape=np.array(shape),
             spacing=np.array([spacing]), features=np.array(list(FEATURE_CHANNELS)))
    return path


def test_field_keeps_only_confident_points(tmp_path):
    from dmta.pharmacophore import load_field
    path = _synthetic_field(tmp_path, hot=[(0, 0, 0.95), (10, 0, 0.3)])
    field = load_field(path, threshold=0.9)
    assert len(field.spots) == 1
    assert field.spots[0].feature == "Donor"
    assert field.spots[0].probability > 0.9


def test_field_suppresses_neighbouring_maxima(tmp_path):
    """Two confident points 2 A apart are one pharmacophore feature, not two."""
    from dmta.pharmacophore import load_field
    path = _synthetic_field(tmp_path, spacing=2.0, hot=[(0, 0, 0.98), (1, 0, 0.97)])
    assert len(load_field(path, threshold=0.9, min_separation=3.0).spots) == 1
    assert len(load_field(path, threshold=0.9, min_separation=1.0).spots) == 2


def test_field_rejects_an_inconsistent_file(tmp_path):
    import numpy as np
    import pytest
    from dmta.pharmacophore import FEATURE_CHANNELS, load_field
    path = tmp_path / "bad.npz"
    np.savez(path, probs=np.zeros((10, len(FEATURE_CHANNELS)), dtype=np.float32),
             origin=np.zeros(3), shape=np.array([6, 6, 6]), spacing=np.array([1.0]),
             features=np.array(list(FEATURE_CHANNELS)))
    with pytest.raises(ValueError, match="rows but shape"):
        load_field(path)

    path2 = tmp_path / "bad2.npz"
    np.savez(path2, probs=np.zeros((8, 2), dtype=np.float32), origin=np.zeros(3),
             shape=np.array([2, 2, 2]), spacing=np.array([1.0]),
             features=np.array(["Donor", "Wingding"]))
    with pytest.raises(ValueError, match="unexpected feature channels"):
        load_field(path2)


def test_feature_vocabulary_matches_rdkit():
    """The model's channel names must be real RDKit families, or nothing aligns."""
    from dmta.pharmacophore import FEATURE_CHANNELS, _feature_factory
    families = set(_feature_factory().GetFeatureFamilies())
    missing = set(FEATURE_CHANNELS) - families
    assert not missing, f"not RDKit feature families: {sorted(missing)}"


def test_matching_is_geometric_not_a_feature_count(tmp_path):
    """The point of the clique constraint: right features, wrong distances, no match."""
    import numpy as np
    from dmta.pharmacophore import load_field, match_score
    # three Donor spots on a 6 A grid, far apart
    path = _synthetic_field(tmp_path, spacing=6.0,
                            hot=[(0, 0, 0.95), (1, 0, 0.95), (7, 0, 0.95)])
    field = load_field(path, threshold=0.9, min_separation=3.0)
    assert len(field.spots) == 3

    spread = field.positions                     # exactly on the spots
    tight = np.array([[0, 0, 0], [0.4, 0, 0], [0, 0.4, 0]], dtype=float)
    types = ["Donor"] * 3
    matched_spread = match_score(field, types, spread, tolerance=0.5)[1]
    matched_tight = match_score(field, types, tight, tolerance=0.5)[1]
    assert matched_spread == 3, "features sitting on the spots should match them all"
    assert matched_tight <= 1, "a cluster matched spots 6 A apart — geometry ignored"


def test_unembeddable_molecule_scores_none(tmp_path):
    """A compound with no 3D structure must not silently score zero."""
    from dmta.pharmacophore import conformer_features, load_field, score_compound
    assert conformer_features("not a molecule at all") == []
    field = load_field(_synthetic_field(tmp_path, hot=[(0, 0, 0.95)]), threshold=0.9)
    assert score_compound(field, "not a molecule at all")["score"] is None


def test_bigger_molecule_can_match_more(tmp_path):
    """Sanity check on direction, and a reminder that size is the confound here."""
    from dmta.pharmacophore import load_field, score_compound
    path = _synthetic_field(tmp_path, spacing=3.0, hot=[
        (0, 4, 0.95), (5, 4, 0.95), (10, 5, 0.95), (40, 1, 0.95), (80, 0, 0.95)])
    field = load_field(path, threshold=0.9, min_separation=2.0)
    small = score_compound(field, "c1ccccc1", n_conformers=2)
    big = score_compound(field, "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
                         n_conformers=2)
    assert big["matched"] >= small["matched"]
