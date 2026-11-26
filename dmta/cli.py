"""Command line for the platform — the path a computational chemist actually uses.

    python3 -m dmta.cli list
    python3 -m dmta.cli propose  --campaign kinaseX --size 16 --out worklist.csv
    python3 -m dmta.cli submit   --campaign kinaseX --results worklist.csv
    python3 -m dmta.cli simulate --campaign kinaseX --cycles 8 --batch 16
    python3 -m dmta.cli export   --campaign kinaseX --out results.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .assay import CsvAssay, SimulatedAssay, write_batch_csv
from .campaign import Campaign, CycleSettings
from .scoring import CHEMBL_PROFILE
from .select import ACQUISITIONS
from .cost import RATES as COST_RATES
from .store import Store

DB = Path(__file__).resolve().parent.parent / "campaigns.db"


def load(name: str, store: Store, chembl_target: str | None = None, **overrides) -> Campaign:
    row = store.get_campaign(name)
    settings = CycleSettings(**{**(row["settings"] if row else {}),
                                **{k: v for k, v in overrides.items() if v is not None}})
    target = chembl_target or (row or {}).get("target", "simulated kinase")
    library = profile = None
    if chembl_target or str(target).startswith("CHEMBL"):
        from .chembl import load_compounds
        code = chembl_target or str(target).split()[0]
        library = load_compounds(code)
        profile = list(CHEMBL_PROFILE)
    return Campaign(name=name, target=target, settings=settings, store=store,
                    library=library, **({"profile": profile} if profile else {}))


def cmd_list(args, store):
    rows = store.list_campaigns()
    if not rows:
        print("no campaigns yet — run `simulate` or `propose` to start one")
        return
    print(f"{'campaign':<24}{'target':<24}{'cycles':>8}{'tested':>8}  created")
    for r in rows:
        print(f"{r['name']:<24}{(r['target'] or ''):<24}{r['cycles']:>8}{r['measured']:>8}  {r['created_at']}")


def cmd_propose(args, store):
    c = load(args.campaign, store, batch_size=args.size, acquisition=args.acquisition)
    batch = c.propose(args.size)
    plan = c.plan(batch)
    out = Path(args.out or f"worklist-{args.campaign}-cycle{c.cycle}.csv")
    write_batch_csv(batch, out)
    print(f"cycle {c.cycle}: {len(batch)} compounds proposed")
    print(f"  {plan['shared_steps']} reactions ({plan['naive_steps']} one-at-a-time, "
          f"{plan['steps_saved']} saved via {len(plan['groups'])} shared intermediates)")
    print(f"  reagent spend ${plan['reagent_cost']:.0f}, longest lead {plan['longest_lead_days']} days")
    print(f"  worklist written to {out} — fill the pIC50 column and run `submit`")


def cmd_submit(args, store):
    c = load(args.campaign, store)
    pending = {x.id for x in c.untested}
    with open(args.results, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r.get(args.column, "").strip()]
    accepted = []
    for r in rows:
        cid = r["compound_id"].strip()
        if cid in pending:
            comp = c.by_id[cid]
            comp.measured = float(r[args.column])
            comp.cycle_made = c.cycle
            accepted.append({"compound_id": cid, "smiles": comp.smiles, "value": comp.measured,
                             "predicted": comp.predicted, "uncertainty": comp.uncertainty})
    if not accepted:
        sys.exit("no new results found in that file")
    batch_id = store.add_batch(c.campaign_id, c.cycle, "external",
                               c.plan([c.by_id[a["compound_id"]] for a in accepted]), {})
    store.add_measurements(c.campaign_id, batch_id, accepted, source="lab")
    report = c._retrain()
    print(f"accepted {len(accepted)} results — {len(c.measured)} compounds now measured")
    print(f"  model: CV R² {report['cv_r2']}, MAE {report['cv_mae']}, σ calibration {report['sigma_calibration']}")
    best = max(c.measured, key=lambda x: x.measured)
    print(f"  best so far: {best.id} at {best.measured} ({best.smiles})")


def cmd_simulate(args, store):
    c = load(args.campaign, store, chembl_target=getattr(args, "chembl", None),
             batch_size=args.batch, acquisition=args.acquisition,
             oracle_seed=args.seed, assay_noise=args.noise)
    oracle = c.assay if hasattr(c.assay, "top_fraction") else None
    target = oracle.top_fraction(0.05) if oracle else set()
    print(f"campaign {args.campaign}: {len(c.library)} compounds, "
          f"{len(c.measured)} already tested, acquisition {c.settings.acquisition}")
    for _ in range(args.cycles):
        s = c.run_cycle()
        p = s["prospective"]
        print(f"  cycle {s['cycle']:>2}: made {s['made']:>3} | tested {s['tested_total']:>4}"
              f" | best {s['best_measured']:>5.2f} | hits {s['hits_in_batch']:>2}"
              f" | {s['plan']['shared_steps']:>3} rxns (vs {s['plan']['naive_steps']})"
              + (f" | prospective ρ {p['spearman']:>5.2f} MAE {p['mae']:.2f}" if p else " | seed batch"))
    if oracle:
        found = sum(1 for x in c.measured if x.id in target)
        print(f"\nfound {found} of the library's best {len(target)} compounds "
              f"by making {len(c.measured)} of {len(c.library)} "
              f"({found / max(1, len(c.measured) * 0.05):.2f}x random expectation)")
        print(f"best measured {max(x.measured for x in c.measured):.2f} "
              f"against a library ceiling of {oracle.ceiling():.2f}")


def cmd_export(args, store):
    c = load(args.campaign, store)
    out = Path(args.out or f"{args.campaign}-results.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["compound_id", "smiles", "cycle", "measured", "predicted", "uncertainty",
                    "score", "mw", "clogp", "tpsa", "qed", "sa_score", "cost_usd"])
        for x in sorted(c.measured, key=lambda y: y.measured, reverse=True):
            w.writerow([x.id, x.smiles, x.cycle_made, x.measured, x.predicted, x.uncertainty,
                        x.score, x.desc["mw"], x.desc["clogp"], x.desc["tpsa"],
                        x.desc["qed"], x.desc["sa_score"], x.cost])
    print(f"wrote {len(c.measured)} measured compounds to {out}")


def cmd_chembl_fetch(args, store):
    from .chembl import SUGGESTED_TARGETS, build_dataset, dataset_summary, load_compounds
    if args.list:
        print("well-populated targets:")
        for k, v in SUGGESTED_TARGETS.items():
            print(f"  {k:<12} {v}")
        return
    build_dataset(args.target, args.type, refresh=args.refresh)
    info = dataset_summary(load_compounds(args.target, args.type))
    print()
    for k, v in info.items():
        print(f"  {k:<16} {v}")


def cmd_evaluate(args, store):
    """The customer-facing flow: their CSV in, an honest report out."""
    from .evaluate import evaluate_dataset
    from .ingest import load_csv
    from .report import write_report

    print(f"reading {args.csv}")
    compounds, ingest = load_csv(
        args.csv, smiles_column=args.smiles_column, value_column=args.value_column,
        id_column=args.id_column, date_column=args.date_column, units=args.units)
    print(f"  columns used: {ingest.columns}")
    print(f"  units: {ingest.units}")
    print(f"  {ingest.rows_read} rows -> {ingest.compounds} usable compounds, "
          f"{ingest.scaffolds} chemotypes")
    for reason, n in ingest.dropped.items():
        print(f"    dropped {n:>5}  {reason}")
    for w in ingest.warnings:
        print(f"  ! {w}")

    ev = evaluate_dataset(compounds, ingest, batch_size=args.batch,
                          cycles=args.cycles, top_n=args.top, verbose=True,
                          cost_rate=args.cost_rate)

    print(f"\nVERDICT: {ev.verdict['headline']}")
    print(f"  {ev.verdict['detail']}")
    if ev.verdict.get("caveat"):
        print(f"  caveat: {ev.verdict['caveat']}")
    print(f"\nCOST: {ev.cost.summary(ev.cost_rate)}")

    out = write_report(ev, args.out or "evaluation-report.html",
                       title=args.title, customer=args.customer,
                       provenance=args.provenance)
    print(f"\nreport written to {out}")
    return ev


def cmd_structure(args, store):
    """Prepare a PDB entry for docking and show exactly what was done to it."""
    from .structure import list_ligands, prepare
    if args.ligands:
        for l in list_ligands(args.pdb):
            print(f"  {l['resname']:<6} chain {l['chain']:<2} res {l['resseq']:<6} {l['atoms']:>4} atoms")
        return
    s = prepare(args.pdb, chain=args.chain, ligand_resname=args.ligand,
                padding=args.padding, refresh=args.refresh)
    print(f"{s.pdb_id}: {s.title}")
    print(f"  resolution   {s.resolution} A, chain {s.chain}, {s.n_protein_atoms} protein atoms")
    print(f"  removed      {s.removed}")
    print(f"  site ligand  {s.site.reference_ligand} ({s.site.n_ligand_atoms} atoms)")
    print(f"  box centre   {s.site.center}  size {s.site.size}")
    print(f"  receptor     {s.receptor_pdb}")


def cmd_chembl_validate(args, store):
    from .chembl import load_compounds
    from .validate import full_report
    lib = load_compounds(args.target, args.type)
    print(f"{args.target}: {len(lib)} compounds, "
          f"{len({c.scaffold_id for c in lib})} Murcko scaffolds\n")
    full_report(lib)
    print("\nrandom is optimistic; scaffold is the honest generalisation test;"
          "\ntime is what the model faces in production.")


def main(argv=None):
    p = argparse.ArgumentParser(prog="dmta", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", default=str(DB))
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("list", help="show campaigns").set_defaults(fn=cmd_list)

    q = sub.add_parser("propose", help="design a batch and write a worklist")
    q.add_argument("--campaign", required=True)
    q.add_argument("--size", type=int, default=None)
    q.add_argument("--acquisition", choices=ACQUISITIONS, default=None)
    q.add_argument("--out")
    q.set_defaults(fn=cmd_propose)

    q = sub.add_parser("submit", help="load assay results from a CSV")
    q.add_argument("--campaign", required=True)
    q.add_argument("--results", required=True)
    q.add_argument("--column", default="pIC50")
    q.set_defaults(fn=cmd_submit)

    q = sub.add_parser("simulate", help="run cycles against the simulated assay")
    q.add_argument("--campaign", required=True)
    q.add_argument("--cycles", type=int, default=5)
    q.add_argument("--batch", type=int, default=None)
    q.add_argument("--acquisition", choices=ACQUISITIONS, default=None)
    q.add_argument("--seed", type=int, default=None)
    q.add_argument("--noise", type=float, default=None)
    q.add_argument("--chembl", help="run retrospectively over a real ChEMBL target, e.g. CHEMBL203")
    q.set_defaults(fn=cmd_simulate)

    q = sub.add_parser("chembl-fetch", help="download and clean a ChEMBL target dataset")
    q.add_argument("--target", default="CHEMBL203")
    q.add_argument("--type", default="IC50")
    q.add_argument("--refresh", action="store_true")
    q.add_argument("--list", action="store_true", help="show suggested targets")
    q.set_defaults(fn=cmd_chembl_fetch)

    q = sub.add_parser("chembl-validate", help="random / scaffold / time split report")
    q.add_argument("--target", default="CHEMBL203")
    q.add_argument("--type", default="IC50")
    q.set_defaults(fn=cmd_chembl_validate)

    q = sub.add_parser("evaluate", help="evaluate a customer CSV and write a report")
    q.add_argument("--csv", required=True, help="their export: structures and results")
    q.add_argument("--out", help="where to write the HTML report")
    q.add_argument("--title", default="Retrospective evaluation")
    q.add_argument("--customer", help="name shown on the report")
    q.add_argument("--provenance", help="where the data came from and who ran the analysis; "
                                        "always set this on a report you will share")
    q.add_argument("--smiles-column"); q.add_argument("--value-column")
    q.add_argument("--id-column"); q.add_argument("--date-column")
    q.add_argument("--units", choices=["log", "nM", "uM"],
                   help="override the detected potency units")
    q.add_argument("--batch", type=int, help="compounds per cycle (default ~4% of the set)")
    q.add_argument("--cycles", type=int)
    q.add_argument("--top", type=int, default=10, help="size of the 'best compounds' target set")
    q.add_argument("--cost-rate", default="local", choices=sorted(COST_RATES),
                   help="hourly rate to price this run at; 'local' means hardware you "
                        "already own, whose marginal cash cost is zero")
    q.set_defaults(fn=cmd_evaluate)

    q = sub.add_parser("structure", help="prepare a PDB entry for docking")
    q.add_argument("--pdb", default="2RGP")
    q.add_argument("--chain", default=None)
    q.add_argument("--ligand", default=None, help="residue name to define the site")
    q.add_argument("--padding", type=float, default=5.0)
    q.add_argument("--ligands", action="store_true", help="just list what is in the file")
    q.add_argument("--refresh", action="store_true")
    q.set_defaults(fn=cmd_structure)

    q = sub.add_parser("export", help="write measured compounds to CSV")
    q.add_argument("--campaign", required=True)
    q.add_argument("--out")
    q.set_defaults(fn=cmd_export)

    args = p.parse_args(argv)
    args.fn(args, Store(args.db))


if __name__ == "__main__":
    main()
