"""Read a customer's own CSV of compounds and results.

Nobody's export looks like ours. Columns are named differently, potency is in
different units, duplicates and salts are everywhere, and some rows are censored
(">10000 nM" means "we stopped measuring", not "10000").

Everything dropped is counted and explained. A cleaning step that quietly discards
rows is how an evaluation ends up flattering itself.
"""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

from rdkit import Chem

from .chem import descriptors, inchikey
from .chembl import murcko, standardize
from .library import Compound

# Column names seen in real exports, in rough order of preference.
SMILES_COLUMNS = ("smiles", "canonical_smiles", "structure", "smiles_string", "mol", "molecule")
VALUE_COLUMNS = ("pic50", "pchembl", "pchembl_value", "pki", "pkd", "pec50",
                 "activity", "value", "result", "ic50", "ic50_nm", "ic50_um",
                 "ki", "ki_nm", "kd", "kd_nm", "ec50", "ec50_nm", "potency")
ID_COLUMNS = ("compound_id", "id", "name", "compound", "molecule_id", "corporate_id", "reg_no")
DATE_COLUMNS = ("date", "year", "created", "synthesised", "synthesized", "registered",
                "document_year", "test_date", "assay_date")

CENSORED = re.compile(r"^\s*[<>~≤≥]")
LOG_UNITS = ("pic50", "pchembl", "pchembl_value", "pki", "pkd", "pec50")
NANOMOLAR = ("ic50", "ic50_nm", "ki", "ki_nm", "kd", "kd_nm", "ec50", "ec50_nm", "potency",
             "activity", "value", "result")
MICROMOLAR = ("ic50_um",)


@dataclass
class IngestReport:
    source: str
    rows_read: int = 0
    compounds: int = 0
    columns: dict = field(default_factory=dict)
    units: str = ""
    dropped: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    value_range: tuple | None = None
    date_range: tuple | None = None
    scaffolds: int = 0

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _pick(fieldnames: list[str], candidates: tuple[str, ...]) -> str | None:
    lowered = {f.lower().strip(): f for f in fieldnames}
    for want in candidates:
        if want in lowered:
            return lowered[want]
    for want in candidates:                     # then substring, e.g. "IC50 (nM)"
        for low, original in lowered.items():
            if want in low:
                return original
    return None


def _to_plog(raw: float, units: str) -> float | None:
    """Convert a potency to the p-scale (pIC50 etc.), where bigger = more potent."""
    if units == "log":
        return raw
    if raw <= 0:
        return None
    molar = raw * 1e-9 if units == "nM" else raw * 1e-6
    return -math.log10(molar)


def detect_units(column: str, values: list[float]) -> str:
    """Log units or a concentration? Decided from the column name, then the data."""
    low = column.lower().strip()
    for name in LOG_UNITS:
        if name in low:
            return "log"
    if "um" in low or "µm" in low or "micro" in low:
        return "uM"
    if "nm" in low or "nano" in low:
        return "nM"
    if values:
        mid = median(values)
        if 0 < mid <= 20:        # pIC50 values live here; nanomolar rarely does
            return "log"
    return "nM"


def load_csv(path: str | Path, smiles_column: str | None = None,
             value_column: str | None = None, id_column: str | None = None,
             date_column: str | None = None, units: str | None = None,
             max_replicate_spread: float = 1.0,
             mw_range: tuple[float, float] = (120.0, 900.0)) -> tuple[list[Compound], IngestReport]:
    path = Path(path)
    report = IngestReport(source=str(path))

    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError("the file has no header row")
        fields = list(reader.fieldnames)
        smi_col = smiles_column or _pick(fields, SMILES_COLUMNS)
        val_col = value_column or _pick(fields, VALUE_COLUMNS)
        id_col = id_column or _pick(fields, ID_COLUMNS)
        date_col = date_column or _pick(fields, DATE_COLUMNS)
        if smi_col is None:
            raise ValueError(f"no SMILES column found. Columns present: {fields}")
        if val_col is None:
            raise ValueError(f"no activity column found. Columns present: {fields}")
        rows = list(reader)

    report.rows_read = len(rows)
    report.columns = {"smiles": smi_col, "value": val_col, "id": id_col, "date": date_col}

    raw_numbers = []
    for r in rows:
        text = (r.get(val_col) or "").strip()
        if text and not CENSORED.match(text):
            try:
                raw_numbers.append(float(re.sub(r"[^\d.eE+-]", "", text)))
            except ValueError:
                pass
    report.units = units or detect_units(val_col, raw_numbers)

    dropped: dict[str, int] = defaultdict(int)
    measurements: dict[str, list[float]] = defaultdict(list)
    dates: dict[str, list[int]] = defaultdict(list)
    kept: dict[str, dict] = {}

    for i, r in enumerate(rows):
        raw_value = (r.get(val_col) or "").strip()
        if not raw_value:
            dropped["no result recorded"] += 1
            continue
        if CENSORED.match(raw_value):
            dropped["censored value (>, <, ~)"] += 1
            continue
        try:
            number = float(re.sub(r"[^\d.eE+-]", "", raw_value))
        except ValueError:
            dropped["result not a number"] += 1
            continue
        value = _to_plog(number, report.units)
        if value is None or not math.isfinite(value):
            dropped["result could not be converted"] += 1
            continue

        smiles = standardize((r.get(smi_col) or "").strip())
        if smiles is None:
            dropped["structure unreadable or too small"] += 1
            continue
        mol = Chem.MolFromSmiles(smiles)
        from rdkit.Chem import Descriptors
        mw = Descriptors.MolWt(mol)
        if not (mw_range[0] <= mw <= mw_range[1]):
            dropped[f"MW outside {mw_range[0]:.0f}–{mw_range[1]:.0f}"] += 1
            continue

        key = inchikey(smiles)
        measurements[key].append(value)
        if date_col:
            raw_date = str(r.get(date_col) or "")
            # Prefer YYYY-MM so a programme that ran over two years still splits
            # sensibly in time; fall back to the bare year.
            m = re.search(r"((?:19|20)\d{2})[-/](\d{1,2})", raw_date)
            if m:
                dates[key].append(int(m.group(1)) + (int(m.group(2)) - 1) / 12)
            else:
                m = re.search(r"(19|20)\d{2}", raw_date)
                if m:
                    dates[key].append(float(m.group(0)))
        kept.setdefault(key, {
            "id": (r.get(id_col) or "").strip() if id_col else f"ROW{i + 1}",
            "smiles": smiles,
        })

    compounds: list[Compound] = []
    for key, rec in kept.items():
        values = measurements[key]
        if len(values) > 1 and (max(values) - min(values)) > max_replicate_spread:
            dropped[f"replicates disagree by >{max_replicate_spread} log"] += 1
            continue
        smiles = rec["smiles"]
        c = Compound(
            id=rec["id"] or key[:14], smiles=smiles, inchikey=key,
            scaffold_id=murcko(smiles), intermediate_key=rec["id"] or key,
            intermediate_smiles=smiles, route=[], desc=descriptors(smiles), alerts=[],
            overall_yield=1.0, route_confidence=None, cost=0.0, lead_days=0,
        )
        c.truth = round(median(values), 3)
        c.year = min(dates[key]) if dates[key] else None
        compounds.append(c)

    if not compounds:
        raise ValueError("no usable compounds survived cleaning — check the columns and units")

    seen_ids: dict[str, int] = defaultdict(int)
    for c in compounds:                       # ids must be unique for the report to make sense
        seen_ids[c.id] += 1
        if seen_ids[c.id] > 1:
            c.id = f"{c.id}#{seen_ids[c.id]}"

    report.compounds = len(compounds)
    report.dropped = dict(sorted(dropped.items(), key=lambda kv: -kv[1]))
    values = [c.truth for c in compounds]
    report.value_range = (round(min(values), 2), round(max(values), 2))
    years = [c.year for c in compounds if c.year]
    report.date_range = ((round(min(years), 2), round(max(years), 2)) if years else None)
    report.scaffolds = len({c.scaffold_id for c in compounds})

    if report.date_range is None:
        report.warnings.append(
            "No usable date column, so the time split — the one that matters most for "
            "predicting compounds nobody has made yet — cannot be run.")
    if len(compounds) < 150:
        report.warnings.append(
            f"Only {len(compounds)} compounds. Results below are indicative; "
            "300+ gives a much firmer answer.")
    if report.units != "log":
        report.warnings.append(
            f"Values read as {report.units} and converted to a p-scale "
            f"(p = -log10[molar]). Pass --units if that is wrong.")
    spread = report.value_range[1] - report.value_range[0]
    if spread < 1.5:
        report.warnings.append(
            f"Potency spans only {spread:.1f} log units. With so little variation there "
            "may be nothing for a model to learn.")
    return compounds, report
