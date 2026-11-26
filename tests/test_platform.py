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
    assert plan["n_compounds"] == len(batch)
    assert plan["routed"] is True
    assert plan["n_shared_intermediates"] == 1          # all six share one


def test_the_plan_reports_the_compound_count_itself():
    """The UI used to infer it as shared_steps - len(groups).

    For a batch of measured compounds that is exactly zero -- each one keys its own
    group, so the two counts are equal -- and the panel read "reactions for 0
    compounds". The count has to come from the batch, and the plan has to say when
    nothing is shared, so a display cannot dress fifteen separate syntheses up as a
    parallel plan that saved nothing.
    """
    lib = fresh_library()
    batch = lib[:15]
    for c in batch:                  # as ChEMBL and CSV ingest build them
        c.route = []
        c.intermediate_key = c.id
    plan = group_by_intermediate(batch)
    assert plan["n_compounds"] == 15
    assert plan["routed"] is False
    assert plan["n_shared_intermediates"] == 0
    assert len(plan["groups"]) == 15
    assert plan["shared_steps"] - len(plan["groups"]) == 0   # the old expression
    assert plan["steps_saved"] == 0


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
        "A-4,Cc1ccccc1C(=O)NCc1ccccc1,> 10000,2021-06-01\n"   # censored, kept
        "A-7,CCO,> 10000,2021-09-01\n"                       # censored, too small
        "A-5,NOTASMILES,15,2021-07-01\n"                     # unparsable
        "A-6,CC(=O)Oc1ccccc1C(=O)O,,2021-08-01\n")           # no result
    compounds, report = load_csv(csv_path)

    assert report.columns["smiles"] == "Structure"
    assert report.columns["value"] == "IC50 (nM)"
    assert report.columns["date"] == "Date Tested"
    assert report.units == "nM"
    # Censored results are KEPT now, with their bound and direction, rather than
    # dropped — ">10000 nM" says the compound is weak, which is a measurement.
    assert report.censored == 1
    assert report.censored_left == 1 and report.censored_right == 0
    # ...but a censored row whose structure is unusable is still dropped, and the
    # reason must name the structure rather than the number.
    assert report.dropped.get("censored structure unreadable or too small") == 1
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


def test_the_uncleaned_control_is_really_the_raw_ridged_inverse():
    """The RMT screening test is only meaningful if its control is honest.

    `rmt_screen_test.py` attributes every difference between clean=True and
    clean=False to the Marchenko-Pastur clipping. That attribution holds only if
    clean=False solves against (C_raw + lambda I), the raw sample correlation
    matrix with a ridge -- not against some partly-cleaned matrix wearing the
    control's name. One wrong noise level would make the control a second
    cleaned model and the whole comparison would measure nothing.

    So this computes the ridged inverse densely, the slow and obvious way, and
    requires the model's weights to match it.
    """
    import numpy as np
    from dmta.rmt_model import RmtModel

    rng = np.random.default_rng(7)
    n, p, lam = 30, 120, 0.4                      # p > n, so C_raw is singular
    X = rng.normal(size=(n, p))
    y = X @ rng.normal(size=p) * 0.1 + rng.normal(size=n)

    model = RmtModel(mode="clip", clean=False, shrinkage=lam).fit_matrix_only(X, y)

    Z = (X - X.mean(axis=0)) / X.std(axis=0)
    C = Z.T @ Z / n
    b = Z.T @ (y - y.mean()) / n
    expected = np.linalg.solve(C + lam * np.eye(p), b)

    assert np.allclose(model.weights_, expected, atol=1e-8), (
        "clean=False is not solving against the raw ridged correlation matrix, so "
        "it is not a control for the cleaning")


def test_the_uncleaned_control_refuses_to_run_without_a_ridge():
    """At shrinkage 0 the control is singular; it must say so, not return noise.

    With p > n the raw correlation matrix has exact zero eigenvalues. Inverting
    it without a ridge gives weights of order 1e16 and predictions that look like
    a catastrophic negative result for the *uncleaned* arm -- which would hand
    the comparison to RMT for a reason that has nothing to do with chemistry.
    """
    import numpy as np
    import pytest
    from dmta.rmt_model import RmtModel

    rng = np.random.default_rng(1)
    X = rng.normal(size=(20, 90))
    y = rng.normal(size=20)
    with pytest.raises(ValueError, match="singular"):
        RmtModel(mode="clip", clean=False).fit_matrix_only(X, y)


def test_cleaning_still_happens_when_it_is_asked_for():
    """The control must not be the default. clean=True clips; clean=False does not."""
    import numpy as np
    from dmta.rmt_model import RmtModel

    rng = np.random.default_rng(3)
    n, p = 40, 200
    direction = rng.normal(size=p)
    scores = rng.normal(size=n)
    X = np.outer(scores, direction) + rng.normal(size=(n, p))   # one real factor
    y = scores + rng.normal(size=n) * 0.1

    cleaned = RmtModel(mode="clip", clean=True).fit_matrix_only(X, y)
    raw = RmtModel(mode="clip", clean=False, shrinkage=1.0).fit_matrix_only(X, y)

    assert cleaned.edge > 1.0, "a Marchenko-Pastur edge of zero means no clipping"
    assert 0 < cleaned.n_signal < min(n, p), (
        "clipping should keep some directions and discard most of them")
    assert raw.edge == 0.0, "the control must not compute a clipping edge at all"
    assert raw.n_signal == min(n, p), (
        f"the control must keep every one of the {min(n, p)} directions the data "
        f"spans, but kept {raw.n_signal}")


def test_the_screening_verdict_will_not_call_a_tie_a_rejection():
    """RESULTS 17's real numbers, pinned against the function that reads them.

    The tie and the loss in that section differ by the modelling mode alone, and
    the whole point of the section is that reporting only the pre-specified mode
    would have turned a tie into a published rejection. So both are pinned here:
    the arm that loses must say LOSES, the arm that ties must say TIE and must
    volunteer how large an undetected effect would still fit.
    """
    import numpy as np
    from rmt_screen_test import read_verdict

    # project minus uncleaned, the nine targets of section 17
    tie = np.array([-0.025, +0.036, -0.025, +0.028, -0.007,
                    -0.045, -0.020, -0.002, -0.042])
    lines = " ".join(read_verdict(tie, -0.028, +0.006, 0.300, 0.180, "cleaning"))
    assert "TIE, NOT A REJECTION" in lines
    assert "0.028" in lines, "a tie must state the effect it could not rule out"
    assert "LOSES" not in lines and "WINS" not in lines

    # clip minus uncleaned, same targets: the interval clears zero
    loss = np.array([-0.029, +0.001, -0.085, -0.142, -0.082,
                     -0.049, -0.086, -0.088, -0.084])
    lines = " ".join(read_verdict(loss, -0.096, -0.046, 0.008, 0.039, "cleaning"))
    assert "LOSES" in lines and "TIE" not in lines

    # and the mirror image must not be read as a loss
    lines = " ".join(read_verdict(-loss, +0.046, +0.096, 0.008, 0.039, "cleaning"))
    assert "WINS" in lines and "LOSES" not in lines


def test_a_tie_interval_that_barely_contains_zero_is_still_a_tie():
    """No fudge factor. An interval touching zero does not get rounded to a win.

    The temptation in section 17 was real: [-0.028, +0.006] is nearly all on one
    side, and "basically a loss" was an available sentence. It is not available to
    the function.
    """
    import numpy as np
    from rmt_screen_test import read_verdict

    lines = " ".join(read_verdict(np.array([-0.05, -0.04, -0.03, -0.02, +0.001]),
                                  -0.049, +0.0005, 0.06, 0.20, "cleaning"))
    assert "TIE, NOT A REJECTION" in lines
    assert "0.049" in lines


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


# --- censored values: compounds too weak to measure ---------------------------

def test_censoring_direction_flips_on_a_concentration_scale():
    """The trap: "IC50 > 10 uM" means WEAK, which is a LOW pIC50.

    pIC50 = -log10(molar) decreases as concentration rises, so a ">" on a
    concentration becomes a "<" on the p-scale. On a column that is already log, the
    inequality stands. Getting this backwards inverts the dataset.
    """
    from dmta.censored import censoring_on_p_scale
    weak = censoring_on_p_scale(">", 5.0, reported_in_log_units=False)
    assert weak.direction == "left", "IC50 > 10uM was read as potent"
    potent = censoring_on_p_scale("<", 9.0, reported_in_log_units=False)
    assert potent.direction == "right"
    # already on a log scale: no flip
    assert censoring_on_p_scale(">", 8.0, reported_in_log_units=True).direction == "right"
    assert censoring_on_p_scale("<", 5.0, reported_in_log_units=True).direction == "left"


def test_truncated_expectation_respects_the_bound():
    """A compound known to be below a bound must never be imputed above it."""
    import numpy as np
    from dmta.censored import _truncated_expectation
    mu = np.array([9.0, 3.0, 4.0, 9.0])          # first is wildly too potent
    sigma = np.array([0.5, 0.5, 0.5, 0.5])
    bound = np.array([5.0, 5.0, 7.0, 7.0])
    left = np.array([True, True, False, False])
    out = _truncated_expectation(mu, sigma, bound, left)
    assert out[0] <= 5.0, "left-censored compound imputed above its bound"
    assert out[1] <= 5.0
    assert out[2] >= 7.0, "right-censored compound imputed below its bound"
    assert out[3] >= 7.0
    # a prediction already inside the region should barely move
    assert abs(out[1] - 3.0) < 1.0


def test_ingest_keeps_censored_rows_with_their_bound(tmp_path):
    from dmta.ingest import censored_compounds, load_csv
    body = ("smiles,IC50 (nM)\n"
            "CC(=O)Oc1ccccc1C(=O)O,120\n"
            "c1ccc(cc1)S(=O)(=O)Nc1ccccc1,45\n"
            "COc1ccc(CCN)cc1,>10000\n"
            "CCCCc1ccc(cc1)C(=O)O,> 50000\n"
            "Cc1ccccc1C(=O)NCc1ccccc1,<0.5\n")
    path = tmp_path / "assay.csv"
    path.write_text(body)
    compounds, report = load_csv(path)
    assert report.compounds == 2
    assert report.censored == 3, f"censored rows were dropped: {report.dropped}"
    assert report.censored_left == 2 and report.censored_right == 1
    cens = censored_compounds(report)
    assert all(c.truth is None for c in cens), "a bound was recorded as a measurement"
    assert all(getattr(c, "censoring", None) is not None for c in cens)
    # none of the censored molecules may leak into the observed set
    assert not ({c.inchikey for c in cens} & {c.inchikey for c in compounds})
    assert any("censored" in w for w in report.warnings)


def test_a_real_measurement_beats_a_bound_for_the_same_molecule(tmp_path):
    from dmta.ingest import censored_compounds, load_csv
    path = tmp_path / "dup.csv"
    path.write_text("smiles,IC50 (nM)\n"
                    "CC(=O)Oc1ccccc1C(=O)O,120\n"
                    "CC(=O)Oc1ccccc1C(=O)O,>10000\n"
                    "c1ccc(cc1)S(=O)(=O)Nc1ccccc1,45\n")
    compounds, report = load_csv(path)
    assert report.compounds == 2
    assert report.censored == 0, "kept a bound for a compound that was measured"
    assert censored_compounds(report) == []


def test_censored_model_pulls_weak_compounds_down():
    """EM must place censored compounds below their bound, not at the training mean."""
    from dmta.assay import SimulatedAssay
    from dmta.censored import CensoredModel, Censoring
    from dmta.library import enumerate_library
    lib = enumerate_library()[:220]
    assay = SimulatedAssay(lib, seed=11)
    values = assay.measure(lib)
    bound = sorted(values.values())[int(len(lib) * 0.4)]
    observed = [c for c in lib if values[c.id] >= bound]
    censored = [c for c in lib if values[c.id] < bound]
    for c in observed:
        c.measured = values[c.id]
    assert len(observed) > 20 and len(censored) > 10

    model = CensoredModel(iterations=3).fit(
        observed, censored, [Censoring(bound=bound, direction="left")] * len(censored))
    assert model.imputations_ is not None
    assert (model.imputations_ <= bound + 1e-6).all(), "imputed above the bound"
    r = model.report()
    assert r["n_censored"] == len(censored) and r["model"] == "censored-em"


# --- hit-finding metrics (actives vs decoys) ----------------------------------

def _labelled(n=2000, n_actives=100):
    import numpy as np
    labels = np.zeros(n, bool)
    labels[:n_actives] = True
    return labels


def test_hit_finding_metrics_on_known_rankings():
    """Perfect, random and inverted rankings have known answers. Both of this
    module's original bugs showed up here: BEDROC returned 12.6 instead of 1.0, and
    an enrichment of exactly 0.0 was reported as None."""
    import numpy as np
    from dmta.hitfinding import hit_finding_report
    labels = _labelled()
    rng = np.random.default_rng(0)
    perfect = np.where(labels, 1.0, 0.0) + rng.normal(0, 1e-9, len(labels))

    good = hit_finding_report(perfect, labels)
    assert good["roc_auc"] == 1.0
    assert good["ef_1pct"] == 20.0                    # the ceiling at 5% actives
    assert 0.99 <= good["bedroc_20"] <= 1.0, good["bedroc_20"]

    chance = hit_finding_report(rng.normal(0, 1, len(labels)), labels)
    assert 0.45 < chance["roc_auc"] < 0.55
    assert chance["ef_1pct"] is not None, "an enrichment of 0.0 was reported as None"
    assert chance["ef_1pct"] < 3.0

    worst = hit_finding_report(-perfect, labels)
    assert worst["roc_auc"] == 0.0
    assert worst["bedroc_20"] <= 0.01


def test_bedroc_is_bounded():
    """BEDROC must live in [0, 1]; an unbounded value means the rescaling is wrong."""
    import numpy as np
    from dmta.hitfinding import bedroc
    labels = _labelled(1000, 50)
    rng = np.random.default_rng(1)
    for scores in (np.where(labels, 1.0, 0.0), rng.normal(0, 1, 1000),
                   -np.where(labels, 1.0, 0.0), np.zeros(1000)):
        value = bedroc(scores, labels)
        assert -0.01 <= value <= 1.01, value


def test_tied_scores_do_not_earn_credit():
    """A method that gives every compound the same score is exactly at chance."""
    import numpy as np
    from dmta.hitfinding import roc_auc
    labels = _labelled(500, 50)
    assert abs(roc_auc(np.zeros(500), labels) - 0.5) < 1e-9


def test_hit_finding_needs_both_classes():
    import numpy as np
    import pytest
    from dmta.hitfinding import roc_auc
    with pytest.raises(ValueError, match="both actives and decoys"):
        roc_auc(np.arange(10.0), np.ones(10, bool))


def test_load_ism_reads_dude_format(tmp_path):
    from dmta.hitfinding import load_ism
    path = tmp_path / "actives_final.ism"
    path.write_text("CC(=O)Oc1ccccc1C(=O)O 101343 CHEMBL63786\n"
                    "c1ccccc1 C04955856\n")
    rows = load_ism(path)
    assert rows == [("CC(=O)Oc1ccccc1C(=O)O", "CHEMBL63786"),
                    ("c1ccccc1", "C04955856")]


def test_every_cli_script_parses_and_has_an_entrypoint():
    """A script with no __main__ guard runs, prints nothing, and exits 0.

    That happened: an edit to hitfinding_benchmark.py replaced the tail of the file
    and took the entrypoint with it, so the benchmark "succeeded" three times in a
    row while producing an empty log. Silent success is the worst failure mode there
    is, and it costs one test to rule out.
    """
    import ast
    import pathlib
    root = pathlib.Path(__file__).resolve().parent.parent
    scripts = sorted(p for p in root.glob("*.py") if not p.name.startswith("_"))
    assert len(scripts) > 5, f"expected several CLI scripts, found {len(scripts)}"
    for path in scripts:
        source = path.read_text()
        ast.parse(source)                       # a syntax error must fail here
        tree = ast.parse(source)
        has_main = any(isinstance(n, ast.FunctionDef) and n.name == "main"
                       for n in tree.body)
        if not has_main:
            continue                            # not a CLI script; nothing to check
        assert '__name__ == "__main__"' in source, (
            f"{path.name} defines main() but never calls it — it will exit 0 "
            "having done nothing")


def test_saturated_enrichment_is_flagged():
    """EF is capped at 1/base_rate; past that the column cannot compare anything.

    On a 26.5%-active test set the cap is 3.77 and the top 1% is ten molecules, so
    three methods came back at 3.77, 3.77 and 3.74 — ten, ten and 9.9 actives out of
    ten — and a verdict function compared them.
    """
    import numpy as np
    from dmta.hitfinding import enrichment_ceiling, hit_finding_report
    labels = np.zeros(1021, bool)
    labels[:271] = True                               # 26.5% active
    assert abs(enrichment_ceiling(labels) - 3.77) < 0.01

    perfect = hit_finding_report(np.where(labels, 1.0, 0.0), labels)
    assert perfect["ef_1pct_saturated"] is True
    assert perfect["ef_1pct"] <= perfect["ef_ceiling"] + 1e-9

    rare = np.zeros(15327, bool)
    rare[:77] = True                                  # 0.5% active, lots of headroom
    assert enrichment_ceiling(rare) > 100
    modest = hit_finding_report(
        np.where(rare, 1.0, 0.0) + np.random.default_rng(0).normal(0, 2, len(rare)), rare)
    assert modest["ef_1pct_saturated"] is False


def test_sign_test_matches_the_binomial_by_hand():
    """The crude test, pinned, because it is the one that gets quoted.

    "The forest won 11 of 18" sounds decisive and is not: p = 0.24. These values are
    computed from the binomial directly, so a regression in the helper shows up as a
    wrong p-value rather than a plausible one.
    """
    from model_paired_test import sign_test
    assert sign_test(5, 10) == 1.0                      # a dead heat
    assert round(sign_test(11, 18), 3) == 0.481         # two-sided
    assert round(sign_test(10, 10), 4) == 0.002         # a clean sweep
    # symmetric: losing 11 of 18 is as unlikely as winning 11 of 18
    assert sign_test(11, 18) == sign_test(7, 18)


def test_signed_rank_finds_a_real_shift_and_ignores_noise():
    import numpy as np
    from model_paired_test import signed_rank_p
    rng = np.random.default_rng(0)
    # a consistent shift of +0.3 on 20 pairs is unmistakable
    shifted = rng.normal(0.3, 0.05, 20)
    assert signed_rank_p(shifted, rng, resamples=4000) < 0.01
    # the same spread with no shift must not be called significant
    noise = rng.normal(0.0, 0.05, 20)
    assert signed_rank_p(noise, rng, resamples=4000) > 0.05
    # all-zero differences cannot support any claim
    assert signed_rank_p(np.zeros(8), rng, resamples=100) == 1.0


def test_the_paired_test_refuses_to_call_a_tie_a_rejection():
    """The bug this whole script exists to prevent.

    Chemprop was summarised as "rejected" on a median gap of -0.026 against a
    between-programme spread of 0.135. Fed differences of that shape, the verdict
    must come back as no evidence of a difference, and must say so in words that
    cannot be quoted as a rejection.
    """
    import json
    import subprocess
    import sys
    import pathlib
    import numpy as np
    root = pathlib.Path(__file__).resolve().parent.parent
    rng = np.random.default_rng(1)
    diffs = rng.normal(-0.026, 0.135, 18)
    payload = {"models": ["forest", "contender"], "programmes": [
        {"document": f"P{i}", "n": 200,
         "scores": {"forest": {"spearman": 0.5, "mae": 0.4},
                    "contender": {"spearman": 0.5 + float(d), "mae": 0.4}}}
        for i, d in enumerate(diffs)]}
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(payload, fh)
        path = fh.name
    out = subprocess.run([sys.executable, "model_paired_test.py", path],
                         cwd=root, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "NO EVIDENCE OF A DIFFERENCE" in out.stdout, out.stdout
    assert "a tie, not a rejection" in out.stdout
    assert "cost decision" in out.stdout
    # and it must not claim the gap is zero
    assert "still inside the interval" in out.stdout


def test_a_plateau_is_not_read_as_a_trend():
    """The fourth verdict-function error, pinned with the real curve.

    chemprop_scale.py first fitted ONE log-linear slope through every training size.
    On the measured EGFR curve that gives +0.125 of rho per decade, which extrapolates
    to a crossover at n = 7,111 — while the measured difference at n = 7,173 was
    -0.033, nowhere near zero. The slope was carried entirely by the n = 250 point.

    Same class of error as comparing a saturated enrichment factor: a summary
    statistic used where it is not valid.
    """
    from chemprop_scale import read_trend
    ns = [250, 500, 1000, 2000, 4000, 7173]
    real = [-0.233, -0.114, -0.082, -0.026, -0.042, -0.033]
    out = "\n".join(read_trend(ns, real))
    assert "TWO REGIMES" in out and "PLATEAUED" in out, out
    assert "must NOT be extrapolated" in out
    assert "-0.034" in out                      # the level, not the slope
    assert "CROSSED" not in out
    assert "STILL CLOSING" not in out           # the wrong verdict it used to print

    # a genuinely closing curve must still be called closing
    closing = [-0.233, -0.180, -0.130, -0.090, -0.050, -0.012]
    out = "\n".join(read_trend(ns, closing))
    assert "STILL CLOSING" in out, out

    # and one that crosses must be called crossed, whatever the slopes do
    crossed = [-0.233, -0.114, -0.040, +0.010, +0.030, +0.044]
    assert "HAS CROSSED" in "\n".join(read_trend(ns, crossed))

    # a flat curve from the start is a different finding from a plateau
    flat = [-0.030, -0.034, -0.029, -0.033, -0.031, -0.034]
    out = "\n".join(read_trend(ns, flat))
    assert "FLAT" in out and "never the explanation" in out, out


# ------------------------------------------------------- the ADME control arm

from adme_test import CLOGP_AT, HEAVY_AT, N_PHYS, check_feature_layout, read_verdict


def test_the_fingerprint_only_arm_really_excludes_clogp():
    """Section 18 turns on one claim: that X[:, :-12] contains no cLogP.

    The whole ADME comparison is "does the model know chemistry, or is it reading
    the trivial baseline off its own input". If stripping the physchem block left
    cLogP behind -- because the block moved, or changed length -- the control arm
    would silently contain the thing it exists to exclude, and the verdict would
    be meaningless while still looking fine.
    """
    smiles = "CCOc1ccccc1C(=O)O"
    vec = feature_vector(smiles)
    d = descriptors(smiles)
    phys = vec[-N_PHYS:]
    assert np.isclose(phys[CLOGP_AT], d["clogp"] / 5.0, atol=1e-5)
    assert np.isclose(vec[HEAVY_AT], d["heavy_atoms"] / 40.0, atol=1e-5)

    fp_only = vec[:-N_PHYS]
    assert len(fp_only) == len(vec) - N_PHYS
    assert np.isclose(fp_only, d["clogp"] / 5.0, atol=1e-5).sum() == 0, (
        "a column equal to cLogP/5 survived the strip"
    )
    # a fingerprint is bits; a physchem column is not
    assert set(np.unique(fp_only)) <= {0.0, 1.0}


def test_the_feature_layout_guard_fires_when_the_layout_moves(monkeypatch):
    """The guard has to fail loudly, not pass quietly, or it is decoration."""
    check_feature_layout()                      # the real layout passes
    import adme_test
    monkeypatch.setattr(adme_test, "CLOGP_AT", 0)
    with pytest.raises(SystemExit, match="feature layout changed"):
        adme_test.check_feature_layout()


def test_the_adme_verdict_will_not_call_a_tie_a_win():
    """Same three-valued rule as sections 16 and 17, pinned again here."""
    assert "WINS" in read_verdict(0.04, 0.01, 0.07, "the forest")
    assert "LOSES" in read_verdict(-0.04, -0.07, -0.01, "the forest")
    tie = read_verdict(0.004, -0.028, 0.036, "the forest")
    assert "TIE" in tie and "WINS" not in tie
    assert "0.036" in tie              # reports the effect still inside the interval


def test_spearman_uses_mid_ranks_and_matches_scipy():
    """_rank had no tie handling: it ordered equal values arbitrarily.

    Section 18 found it because a constant predictor then has real rank variance,
    the 'is this constant' guard never fires, and the correlation comes back as a
    small random number. The visible symptom was a bootstrap point estimate
    outside its own confidence interval.
    """
    scipy_stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(0)
    for n, n_tied in ((60, 0), (60, 20), (200, 80)):
        obs = rng.normal(size=n)
        pred = obs + rng.normal(size=n) * 0.7
        if n_tied:
            pred[:n_tied] = pred[0]
        ours = regression_metrics(list(pred), list(obs))["spearman"]
        ref = round(float(scipy_stats.spearmanr(pred, obs).statistic), 3)
        assert ours == ref, f"n={n} ties={n_tied}: {ours} != {ref}"


def test_a_constant_prediction_has_no_rank_correlation():
    """Undefined, not zero — and emphatically not a small random number."""
    rng = np.random.default_rng(1)
    obs = list(rng.normal(size=50))
    assert regression_metrics([3.0] * 50, obs)["spearman"] is None
    # the floor still has a defined error, which is the point of keeping it
    assert regression_metrics([3.0] * 50, obs)["mae"] is not None


def test_ties_do_not_move_the_potency_numbers():
    """The fix must not silently restate section 1.

    Forest predictions on distinct compounds are almost never exactly equal, so
    mid-ranks and ordinal ranks agree. If this ever fails, the published split
    numbers need re-checking, not the test.
    """
    rng = np.random.default_rng(2)
    obs = rng.normal(size=300)
    pred = obs + rng.normal(size=300) * 0.5          # no exact ties
    assert len(set(pred)) == len(pred)
    mid = regression_metrics(list(pred), list(obs))["spearman"]
    order = pred.argsort()
    ordinal = np.empty(len(pred));  ordinal[order] = np.arange(len(pred))
    ro = _rank_free_spearman(ordinal, obs)
    assert abs(mid - ro) < 0.002, (mid, ro)


def _rank_free_spearman(rank_p: np.ndarray, o: np.ndarray) -> float:
    from dmta.models import _rank
    rank_o = _rank(o)
    sd = rank_p.std() * rank_o.std()
    return float(((rank_p - rank_p.mean()) * (rank_o - rank_o.mean())).mean() / sd)


# ------------------------------------------------- F13: the molecules stay here

import socket as _socket
from pathlib import Path as _Path

from dmta.egress import EgressAttempt, is_local_address, no_egress
from egress_test import FETCHERS, scan_for_network_use


def test_the_guard_fires_on_a_remote_address():
    """Without this, every other test in this block passes against a no-op guard."""
    with no_egress():
        with pytest.raises(EgressAttempt):
            _socket.create_connection(("93.184.216.34", 80), timeout=1)
        with pytest.raises(EgressAttempt):
            _socket.getaddrinfo("example.com", 80)


def test_the_guard_leaves_local_traffic_alone():
    """The claim is "nothing leaves the machine", not "no sockets exist".

    A guard that blocks loopback would stop the customer running the UI, and one
    that blocks AF_UNIX would stop scikit-learn's own worker processes -- so it
    would fail for reasons that have nothing to do with data leaving, and whoever
    hit that would switch the guard off rather than debug it.
    """
    with no_egress():
        listener = _socket.socket(_socket.AF_INET, _socket.SOCK_STREAM)
        listener.bind(("127.0.0.1", 0))
        listener.listen(1)
        client = _socket.create_connection(listener.getsockname(), timeout=2)
        client.close()
        listener.close()
        a, b = _socket.socketpair(_socket.AF_UNIX, _socket.SOCK_STREAM)
        a.close()
        b.close()


def test_an_unparseable_destination_is_treated_as_remote():
    assert is_local_address("127.0.0.1") and is_local_address("::1")
    assert is_local_address("localhost") and is_local_address("")
    assert not is_local_address("example.com")
    assert not is_local_address("8.8.8.8")
    assert not is_local_address(None)          # the safe default for "I cannot tell"
    assert not is_local_address(("tuple",))


def test_a_full_cycle_runs_with_every_remote_address_blocked():
    """The claim itself: Design -> Make -> Test contacts nothing."""
    with no_egress() as contacted:
        c = Campaign("t-egress", settings=CycleSettings(batch_size=8, oracle_seed=11))
        s0 = c.run_cycle()
        s1 = c.run_cycle()
    assert contacted == []
    assert s1["tested_total"] == 16 and s0["cycle"] == 0


def test_only_the_two_opt_in_fetchers_touch_the_network():
    """This is the test that keeps the claim true after today.

    The cycle not contacting anything is a fact about the code as it stands; one
    `import requests` in select.py would quietly end it, and the run-time test above
    would still pass as long as that import was not reached on the happy path. So
    the scan is static and over-sensitive on purpose: a new network-capable module
    has to be added to FETCHERS deliberately, by someone who reads why.
    """
    found = set(scan_for_network_use(_Path("dmta")))
    assert found == FETCHERS, (
        f"network use changed: unexpected {sorted(found - FETCHERS)}, "
        f"gone {sorted(FETCHERS - found)}. If this is deliberate and the module is an "
        f"opt-in fetch the user invokes by name, add it to FETCHERS in egress_test.py "
        f"and say so in RESULTS.md section 20."
    )


def test_the_scaffold_split_is_deterministic_and_never_shares_a_scaffold():
    """Section 21's split, and the precedence bug it was written with.

    The first version of this function read

        if n + len(grp) <= target or rng.random() < 0.5 and n < target:

    which Python groups as `A or (B and C)` -- a coin flip inside a split that the
    whole comparison is run on. Two runs of the same experiment would then be two
    different experiments, and the paired intervals across proteins would silently
    mix them. This pins both properties the split has to have: the same input gives
    the same split, and no Murcko scaffold appears on both sides, because a shared
    scaffold leaks the training chemistry into the held-out set and turns a
    generalisation test back into a random split.
    """
    from jump_annotated_test import scaffold_key, scaffold_split

    # Six ring families of descending size plus two acyclic compounds. The split
    # fills training with the largest families first, so a fixture needs enough
    # distinct ring systems for the held-out side to contain a real scaffold --
    # with too few, every test-side compound is acyclic, the overlap check compares
    # two empty sets and passes without testing anything.
    smiles = [
        "CC(=O)Oc1ccccc1C(=O)O", "CCOc1ccccc1C(=O)O", "COc1ccccc1C(=O)O",
        "c1ccc2ccccc2c1", "Cc1ccc2ccccc2c1", "CCc1ccc2ccccc2c1",
        "C1CCNCC1", "CC1CCNCC1", "CCC1CCNCC1",
        "c1ccncc1", "Cc1ccncc1",
        "c1ccc(-c2ccccc2)cc1", "Cc1ccc(-c2ccccc2)cc1",
        "c1ccc2[nH]ccc2c1", "Cc1ccc2[nH]ccc2c1",
        "CCN(CC)CC", "CCCN(CC)CC",
    ]
    a = scaffold_split(smiles)
    b = scaffold_split(smiles)
    assert (a == b).all(), "the split must not change between two identical calls"

    # These two assertions are why this test exists. scaffold_key once returned ""
    # for EVERY molecule -- rdkit.Chem.Scaffolds.MurckoScaffold raises "bad any cast"
    # on all inputs in RDKit 2026.09.2, and a bare `except: return ""` turned that
    # total failure into a silent one. Every compound then became its own group and
    # the "scaffold split" was a split on row order, which no assertion about
    # disjointness can catch: distinct keys are trivially disjoint. So check first
    # that real scaffolds exist and that compounds actually share them.
    keys = [scaffold_key(s) for s in smiles]
    assert all(k != "unknown" for k in keys), "every fixture SMILES must parse"
    real = [k for k in keys if k != "acyclic"]
    assert len(set(real)) < len(real), \
        "no two fixture compounds share a scaffold -- scaffold_key is not grouping"

    train_scaffolds = {k for k, t in zip(keys, a) if t and k != "acyclic"}
    test_scaffolds = {k for k, t in zip(keys, a) if not t and k != "acyclic"}
    assert train_scaffolds and test_scaffolds, \
        "both sides must hold real scaffolds or the overlap check is vacuous"
    shared = train_scaffolds & test_scaffolds
    assert not shared, f"scaffolds on both sides of the split: {shared}"


# ------------------------------------------- choosing the model by measurement

from dmta.models import DescriptorRidge, _grouped_oof_mae, choose_model


def _label_library(fn):
    lib = fresh_library()[:260]
    for c in lib:
        c.measured = float(fn(c))
    return lib


def test_the_selector_keeps_the_forest_when_the_signal_is_not_in_the_descriptors():
    """Known answer: a target built from fingerprint structure, not physchem.

    The measured value keys off a substructure the descriptors cannot see, so a
    twelve-number linear model has nothing to fit and the forest must win.
    """
    lib = _label_library(lambda c: 7.0 + 1.5 * ("c1ccccc1" in c.smiles)
                         + 0.9 * c.smiles.count("N") % 3)
    _, pick = choose_model(lib)
    assert pick["chosen"] == "forest", pick
    assert pick["scoring"] == "scaffold-grouped out-of-fold MAE"


def test_the_selector_ships_the_ridge_when_the_endpoint_is_linear_in_descriptors():
    """Known answer: section 18's solubility shape, made exact.

    A clean linear function of cLogP and MW over a wide range is what a ridge
    extrapolates and a forest cannot -- the forest's leaves stop at the edge of
    its training data. If the selector cannot find the ridge here it will not
    find it on a real endpoint either.
    """
    lib = _label_library(lambda c: -0.9 * c.desc["clogp"] - 0.011 * c.desc["mw"] + 2.0)
    _, pick = choose_model(lib)
    assert pick["chosen"] == "ridge", pick
    assert pick["ridge_mae"] < pick["forest_mae"]


def test_the_selector_folds_never_share_a_scaffold():
    """The fix that made the selector work, pinned.

    The first version cross-validated on random folds and picked the forest 5
    times out of 5 on solubility -- the honest answer to the wrong question,
    because random folds put near-identical analogues on both sides. If grouping
    silently regresses, the selector goes back to being wrong in exactly the
    case it exists for.
    """
    seen: list[set[str]] = []

    class Spy:
        def fit(self, train, target="measured"):
            seen.append({c.scaffold_id for c in train})
            self.train = train
            return self

        def predict(self, test):
            self.held = {c.scaffold_id for c in test}
            assert not (seen[-1] & self.held), "a scaffold was in train and test at once"
            return np.zeros(len(test)), np.ones(len(test))

    lib = _label_library(lambda c: c.desc["mw"] / 100.0)
    out = _grouped_oof_mae(Spy, lib, "measured")
    assert out is not None and len(seen) >= 2


def test_the_ridge_reports_a_usable_uncertainty():
    """A linear model has no ensemble to disagree with itself.

    Without the bootstrap there is no sigma at all, and every acquisition
    function except greedy spends part of its batch on sigma.
    """
    lib = _label_library(lambda c: -0.9 * c.desc["clogp"] + 2.0)
    m = DescriptorRidge(n_boot=20, seed=0).fit(lib)
    mu, sigma = m.predict(lib[:40])
    assert len(mu) == 40 and (sigma > 0).all()
    r = m.report()
    assert r["kind"] == "descriptor_ridge" and r["n_bootstrap"] > 0
    m.apply(lib[:5])
    assert all(c.predicted is not None and c.uncertainty is not None for c in lib[:5])


def test_a_campaign_honours_a_forced_model_and_will_not_guess_on_too_little_data():
    forced = Campaign("t-ridge", settings=CycleSettings(batch_size=8, model="ridge"))
    forced.run_cycle()
    assert forced.model_choice["chosen"] == "ridge"
    assert isinstance(forced.model, DescriptorRidge)

    auto = Campaign("t-auto", settings=CycleSettings(batch_size=8, model="auto"))
    auto.run_cycle()
    assert len(auto.measured) < Campaign.MIN_FOR_CHOICE
    assert auto.model_choice is None, "decided the model on a handful of compounds"
    assert isinstance(auto.model, PropertyModel)


# ------------------------------------- the Cell Painting controls (section 21)

def test_the_formula_prefilter_cannot_discard_a_true_match():
    """jump_match.py skips ~95% of molecules without computing their InChIKey.

    It is allowed to do that only because two molecules with the same InChIKey
    necessarily have the same molecular formula, so a formula miss is a guaranteed
    InChIKey miss. If that ever stopped holding, the matcher would silently find
    fewer compounds, and the size of the overlap is the number section 21's whole
    verdict rests on.
    """
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    pairs = [
        ("CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1OC(C)=O"),      # aspirin, two writings
        ("C[C@H](N)C(=O)O", "C[C@H](N)C(=O)O"),                   # with stereo
        ("[Na+].[Cl-]", "[Cl-].[Na+]"),                           # a salt, reordered
        ("c1ccccc1", "C1=CC=CC=C1"),                              # kekulised
    ]
    for a, b in pairs:
        ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
        assert ma is not None and mb is not None, (a, b)
        if Chem.MolToInchiKey(ma) == Chem.MolToInchiKey(mb):
            assert CalcMolFormula(ma) == CalcMolFormula(mb), (
                f"{a} and {b} share an InChIKey but not a formula -- the pre-filter "
                f"in jump_match.py would drop this true match")


def test_the_scrambled_arm_keeps_the_columns_and_destroys_the_pairing():
    """The control's whole job is to be identical except for the information.

    Section 21 turns on the scrambled arm costing the forest the same 0.12 of AUC
    that the real profile did. That is only evidence of dilution if the scrambling
    is real: a permutation that left most rows in place would still carry the true
    compound-to-profile link, would score like the real arm, and the agreement
    would look like proof while being an artefact of a broken control.
    """
    rng = np.random.default_rng(0)
    n, p = 400, 20
    P = rng.normal(size=(n, p))
    perm = np.random.default_rng(0).permutation(n)
    S = P[perm]

    assert np.allclose(np.sort(P, axis=0), np.sort(S, axis=0)), \
        "a permutation of rows must leave every column's values untouched"
    assert (perm == np.arange(n)).mean() < 0.05, \
        "the permutation left too many rows in place to be a control"
    moved = (np.linalg.norm(P - S, axis=1) > 1e-9).mean()
    assert moved > 0.95, f"only {moved:.0%} of rows actually moved"


def test_rank_fusion_is_immune_to_one_arm_gaining_columns():
    """Late fusion exists because concatenation is not immune.

    Concatenating the profile lost 0.12 of AUC to pure noise, so it could not
    answer whether the profile carries information. Rank-averaging two separately
    fitted models cannot dilute either one: it depends only on the ORDER of each
    input, so a change of scale in one arm -- which is what a differently sized
    model produces -- leaves the fused ranking untouched.
    """
    def ranks(v):
        order = np.empty(len(v), float)
        order[np.argsort(v, kind="stable")] = np.arange(len(v))
        return order / max(len(v) - 1, 1)

    rng = np.random.default_rng(1)
    a, b = rng.random(200), rng.random(200)
    fused = (ranks(a) + ranks(b)) / 2
    fused_rescaled = (ranks(0.01 * a + 5.0) + ranks(b)) / 2
    assert np.allclose(fused, fused_rescaled), \
        "rank fusion must depend only on the order of each arm"


def test_the_cell_painting_verdict_will_not_call_a_tie_a_rejection():
    """Section 21's real numbers, through the verdict function.

    Every comparison that matters there is a tie, and the one arm that "loses"
    loses to its own scrambled control by the same margin -- which is the finding,
    not a defeat for Cell Painting. If this function ever starts calling those
    intervals rejections, the section's conclusion inverts: we would publish "a
    morphological profile does not predict activity" on an interval containing zero.
    """
    from jump_cellpaint_test import verdict

    # profile against its own permutation: +0.015, CI [-0.082, +0.114]
    out = verdict(0.015, -0.082, 0.114, "profile")
    assert "TIE" in out and "not a measurement of no effect" in out, out

    # real profile against scrambled, same column count: -0.001, CI [-0.092, +0.095]
    assert "TIE" in verdict(-0.001, -0.092, 0.095, "both")

    # fused against fingerprint: -0.043, CI [-0.088, +0.003] -- barely a tie, still one
    assert "TIE" in verdict(-0.043, -0.088, 0.003, "fused")

    # and the one real loss: concatenation against the fingerprint alone
    out = verdict(-0.123, -0.198, -0.047, "both")
    assert "LOSES" in out and "TIE" not in out, out


def test_murcko_returns_real_scaffolds_and_not_a_silent_fallback():
    """murcko() ends in `except Exception: return "unknown"`.

    If the RDKit call it wraps ever breaks the way MurckoScaffold did in
    2026.09.2 -- raising on every input rather than on odd ones -- every
    compound would come back "unknown", collapse into a single group, and both
    the scaffold split and the model selector's scaffold-grouped folds would
    silently become splits on row order. Nothing would raise. This is the
    known-answer guard for that, one layer under section 18's control arm.
    """
    assert murcko("CC(=O)Oc1ccccc1C(=O)O") == "c1ccccc1"      # aspirin
    assert murcko("CCN(CC)CC") == "acyclic"                   # no ring system
    assert murcko("not a smiles") == "unknown"                # the fallback works

    lib = fresh_library()[:200]
    keys = [c.scaffold_id for c in lib]
    assert "unknown" not in keys, "the scaffold fallback fired on a clean library"
    assert len(set(keys)) > 1, "every compound landed in one scaffold group"


def test_the_selector_says_so_rather_than_guessing_when_scaffolds_collapse():
    """The safe failure, if the guard above ever does fire.

    One group means no held-out fold is possible. The selector must report that
    it could not compare, not quietly return whichever model it fitted last.
    """
    lib = fresh_library()[:80]
    for c in lib:
        c.measured = 7.0
        c.scaffold_id = "unknown"
    _, pick = choose_model(lib)
    assert pick["chosen"] == "forest"
    assert "too few scaffolds" in pick["reason"]
