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
