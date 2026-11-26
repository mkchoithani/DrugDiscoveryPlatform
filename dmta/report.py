"""Render an evaluation as one self-contained HTML file.

No external assets, so it can be emailed, archived, or printed to PDF and still look
the same in five years. It prints light, because this gets forwarded to people who
read it on paper.

The report shows failures as prominently as successes. A vendor report that only
contains good news is the one a customer should not trust.
"""
from __future__ import annotations

import datetime as _dt
import html
from pathlib import Path

from .evaluate import Evaluation

CSS = """
*{box-sizing:border-box}
body{margin:0;background:#f6f7f9;color:#14181f;font:15px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.page{max-width:860px;margin:0 auto;padding:44px 36px 80px;background:#fff;min-height:100vh}
h1{font-size:25px;margin:0 0 6px;letter-spacing:-.3px}
h2{font-size:17px;margin:38px 0 12px;padding-bottom:7px;border-bottom:2px solid #14181f}
h3{font-size:14px;margin:22px 0 8px;text-transform:uppercase;letter-spacing:.06em;color:#5b6572}
p{margin:10px 0}
.sub{color:#5b6572;margin:0 0 4px}
.verdict{border-left:4px solid #1f6feb;background:#f0f5ff;padding:18px 20px;margin:22px 0;border-radius:0 7px 7px 0}
.verdict.good{border-color:#1a7f37;background:#f0f9f2}
.verdict.warn{border-color:#9a6700;background:#fff8e8}
.verdict.bad{border-color:#b42318;background:#fdf3f2}
.verdict strong{font-size:17px;display:block;margin-bottom:6px}
.headline{display:flex;gap:26px;flex-wrap:wrap;margin:22px 0}
.stat{flex:1;min-width:140px;border:1px solid #e3e6ea;border-radius:8px;padding:14px 16px}
.stat b{display:block;font-size:26px;line-height:1.15;font-variant-numeric:tabular-nums}
.stat span{color:#5b6572;font-size:12px}
table{width:100%;border-collapse:collapse;margin:12px 0;font-size:14px}
th,td{text-align:left;padding:8px 10px;border-bottom:1px solid #e9ecef}
th{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:#5b6572}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
tr.good td{background:#f4fbf6}
tr.bad td{background:#fdf5f4}
.note{background:#f6f7f9;border-radius:7px;padding:13px 16px;font-size:13.5px;color:#39414d;margin:14px 0}
.warn{background:#fff8e8;border-left:3px solid #9a6700}
.muted{color:#5b6572}
svg{display:block;width:100%;height:auto;margin:14px 0}
.axis{stroke:#c9ced6;stroke-width:1}
.grid{stroke:#eceff3;stroke-width:1}
.lbl{fill:#5b6572;font:10px ui-monospace,Menlo,monospace}
footer{margin-top:46px;padding-top:16px;border-top:1px solid #e9ecef;color:#5b6572;font-size:12.5px}
code{background:#f0f2f5;padding:1px 5px;border-radius:4px;font-size:13px}
@media print{body{background:#fff}.page{padding:0}h2{break-after:avoid}}
"""


def _year_label(value) -> str:
    """2021.25 -> 2021-04, 2021 -> 2021."""
    try:
        year = int(value)
        month = round((float(value) - year) * 12) + 1
        return f"{year}-{month:02d}" if month > 1 or float(value) != year else str(year)
    except (TypeError, ValueError):
        return str(value)


def _esc(v) -> str:
    return html.escape(str(v))


def _fmt(v, dp=3) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{dp}f}"
    return str(v)


def _line_chart(series: list[tuple[str, list[float], str]], x_values: list[float],
                x_label: str, y_label: str, height: int = 230) -> str:
    """Minimal inline SVG — no chart library, nothing to break later."""
    W, H, L, R, T, B = 820, height, 56, 18, 16, 40
    all_y = [v for _, ys, _ in series for v in ys if v is not None]
    if not all_y or not x_values:
        return "<p class='muted'>not enough data to plot</p>"
    y_lo, y_hi = min(all_y), max(all_y)
    pad = (y_hi - y_lo) * 0.12 or 0.5
    y_lo, y_hi = y_lo - pad, y_hi + pad
    x_lo, x_hi = min(x_values), max(x_values)
    sx = lambda v: L + (0 if x_hi == x_lo else (v - x_lo) / (x_hi - x_lo)) * (W - L - R)
    sy = lambda v: H - B - (0 if y_hi == y_lo else (v - y_lo) / (y_hi - y_lo)) * (H - T - B)

    parts = [f'<svg viewBox="0 0 {W} {H}" role="img">']
    for i in range(5):
        y = y_lo + i * (y_hi - y_lo) / 4
        parts.append(f'<line class="grid" x1="{L}" y1="{sy(y):.1f}" x2="{W-R}" y2="{sy(y):.1f}"/>')
        parts.append(f'<text class="lbl" x="{L-8}" y="{sy(y)+3:.1f}" text-anchor="end">{y:.2f}</text>')
    parts.append(f'<line class="axis" x1="{L}" y1="{H-B}" x2="{W-R}" y2="{H-B}"/>')
    parts.append(f'<line class="axis" x1="{L}" y1="{T}" x2="{L}" y2="{H-B}"/>')
    for v in x_values:
        parts.append(f'<text class="lbl" x="{sx(v):.1f}" y="{H-B+16}" text-anchor="middle">{v:g}</text>')
    parts.append(f'<text class="lbl" x="{(W+L)/2:.0f}" y="{H-6}" text-anchor="middle">{_esc(x_label)}</text>')
    parts.append(f'<text class="lbl" x="14" y="{(H-B+T)/2:.0f}" text-anchor="middle" '
                 f'transform="rotate(-90 14 {(H-B+T)/2:.0f})">{_esc(y_label)}</text>')
    for name, ys, colour in series:
        pts = [(sx(x), sy(y)) for x, y in zip(x_values, ys) if y is not None]
        if not pts:
            continue
        d = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(pts))
        parts.append(f'<path d="{d}" fill="none" stroke="{colour}" stroke-width="2.5"/>')
        for x, y in pts:
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{colour}"/>')
        parts.append(f'<text class="lbl" x="{pts[-1][0]-6:.1f}" y="{pts[-1][1]-9:.1f}" '
                     f'text-anchor="end" fill="{colour}">{_esc(name)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _splits_table(splits: dict) -> str:
    rows = []
    meaning = {
        "random": "Optimistic. Close analogues of each test compound sit in the training set.",
        "scaffold": "Whole chemotypes held out. The honest test of generalisation.",
        "time": "Trained on earlier compounds, tested on later ones. What production actually faces.",
    }
    for name in ("random", "scaffold", "time"):
        m = splits.get(name, {})
        if "error" in m or not m:
            rows.append(f"<tr><td><b>{name}</b></td><td colspan='5' class='muted'>"
                        f"{_esc(m.get('error', 'not run'))}</td>"
                        f"<td class='muted'>{_esc(meaning[name])}</td></tr>")
            continue
        rho = m.get("spearman")
        cls = "good" if (rho or 0) >= 0.5 else ("bad" if (rho or 0) < 0.3 else "")
        hit = f"{m['top_decile_hit_rate']:.0%} vs {m['base_hit_rate']:.0%}"
        if m.get("active_threshold") is not None:
            hit += f"<br><span class='muted' style='font-size:11px'>active ≥ {m['active_threshold']}</span>"
        rows.append(
            f"<tr class='{cls}'><td><b>{name}</b>"
            + (f" <span class='muted'>(from {_year_label(m['cutoff_year'])})</span>"
               if m.get("cutoff_year") else "")
            + f"</td><td class='n'>{m['n_train']}</td><td class='n'>{m['n_test']}</td>"
            f"<td class='n'>{_fmt(m.get('spearman'))}</td><td class='n'>{_fmt(m.get('r2'))}</td>"
            f"<td class='n'>{hit}</td><td class='muted'>{_esc(meaning[name])}</td></tr>")
    return ("<table><thead><tr><th>Split</th><th class='n'>Train</th><th class='n'>Test</th>"
            "<th class='n'>ρ</th><th class='n'>R²</th><th class='n'>Top 10% hit rate</th>"
            "<th>What it means</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _cost_block(ev) -> str:
    """What the analysis cost, beside what it found.

    Measured seconds and estimated money are rendered differently on purpose: the
    estimate carries its rate and that rate's source inline, so the page cannot be
    read as having observed a price it only assumed.
    """
    from .cost import COMPOUND_COST_SOURCES, COMPOUND_COST_USD

    d = ev.cost.as_dict(ev.cost_rate)
    m, est, ass = d["measured"], d["estimate"], d["assumed"]
    if not m["steps"]:
        return ""

    rows = "".join(
        f"<tr><td>{_esc(st['name'])}</td><td class='n'>{st['wall_s']:.2f}</td>"
        f"<td class='n'>{st['cpu_s']:.2f}</td></tr>" for st in m["steps"])
    rows += (f"<tr><td><b>total</b></td><td class='n'><b>{m['wall_s']:.2f}</b></td>"
             f"<td class='n'><b>{m['cpu_s']:.2f}</b></td></tr>")

    if est is None:
        money = ("<p style='margin-bottom:0'>No hourly rate was supplied, so no cost is "
                 "shown. Time is measured; money is an assumption, and this page does "
                 "not invent one.</p>")
    else:
        lo, hi = d["compounds_not_made"]
        lo_c, hi_c = min(COMPOUND_COST_USD), max(COMPOUND_COST_USD)
        per = (est["per_unit"] or {}).get("verdicts")
        money = (
            f"<p><b>Estimated cost: {_esc(_money_html(est['usd']))}</b>"
            + (f" for the answer above" if per is not None else "")
            + f", assuming <b>${ass['usd_per_hour']:.2f}/hour</b> "
              f"({_esc(ass['label'])}). That rate is an assumption, not a measurement: "
              f"{_esc(ass['source'])}.</p>"
            f"<p style='margin-bottom:0'>Making and testing <em>one</em> compound costs "
            f"roughly ${lo_c:,.0f}–${hi_c:,.0f}"
            + (f", so this entire analysis cost under a thousandth of a single compound."
               if hi < 0.001 and hi > 0 else
               f", so this entire analysis cost about {lo:.3g}–{hi:.3g} of one compound."
               if hi > 0 else
               ". This ran on hardware already owned, so its marginal cash cost was nil.")
            + "</p>")

    if d.get("warning"):
        money += f"<p style='margin-bottom:0'><b>Note.</b> {_esc(d['warning'])}</p>"
    sources = "".join(f"<li>{_esc(x)}</li>" for x in COMPOUND_COST_SOURCES)
    scored = d["units"].get("compounds_scored")
    aside = ("" if not scored else
             f"<p style='margin-bottom:0'><small>Spread across the {scored:,} compounds "
             f"scored this is a fraction of a cent each, which is true and not the "
             f"number that decides anything. One run answers one question; that is what "
             f"the figure above prices.</small></p>")

    return f"""
<h2>5 · What this answer cost</h2>
<table><thead><tr><th>Step</th><th class="n">Wall seconds</th>
<th class="n">CPU seconds</th></tr></thead><tbody>{rows}</tbody></table>
<div class="note">{money}{aside}
<p style="margin-bottom:0"><small>Per-compound cost sources (both are single-project
self-reports, not industry benchmarks):<ul>{sources}</ul></small></p>
</div>
"""


def _money_html(usd: float) -> str:
    if usd == 0:
        return "$0.00"
    if usd < 0.01:
        return f"${usd:.5f}"
    return f"${usd:,.2f}"


def render(ev: Evaluation, title: str = "Retrospective evaluation",
           customer: str | None = None, provenance: str | None = None) -> str:
    """provenance: where the data came from and who made this analysis. Always set it
    when a report may be shared, so the page cannot be mistaken for the data owner's
    own document."""
    ing, loop, splits, verd = ev.ingest, ev.loop, ev.splits, ev.verdict
    today = _dt.date.today().isoformat()

    kind = ("good" if verd["headline"].startswith("the loop should")
            else "bad" if verd["headline"].startswith("not yet")
            else "warn" if verd["headline"].startswith("promising") else "")

    cycles = list(range(1, len(loop["best_curve"]) + 1))
    experiments = [c * loop["batch_size"] for c in cycles]
    chart = _line_chart(
        [("this loop", loop["best_curve"], "#1a7f37"),
         ("random order", loop["random_curve"], "#8a929d")],
        experiments, "compounds made", "best potency found so far")
    found_chart = _line_chart(
        [(f"top {loop['top_n']} found", [float(x) for x in loop["found_curve"]], "#1f6feb")],
        experiments, "compounds made", f"how many of the best {loop['top_n']} found", height=200)

    dropped = "".join(
        f"<tr><td>{_esc(k)}</td><td class='n'>{v}</td></tr>" for k, v in ing.dropped.items())
    warnings = "".join(f"<div class='note warn'>{_esc(w)}</div>" for w in ing.warnings)

    saved = ""
    if loop["experiments_for_best"]:
        pct = loop["experiments_for_best"] / loop["experiments_available"]
        saved = (f"The single best compound in this dataset was reached after "
                 f"<b>{loop['experiments_for_best']} experiments</b> — "
                 f"{pct:.0%} of the {loop['experiments_available']} that were actually run.")
    else:
        saved = (f"The best compound was not reached within {loop['tested_total']} experiments "
                 f"({loop['fraction_tested']:.0%} of the programme).")

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_esc(title)}</title><style>{CSS}</style></head>
<body><div class="page">

<h1>{_esc(title)}</h1>
<p class="sub">{_esc(customer or 'Dataset')} · {ing.compounds} compounds · generated {today}</p>

<div class="verdict {kind}">
  <strong>{_esc(verd['headline'].capitalize())}</strong>
  {_esc(verd['detail'])}
  {f"<p class='muted' style='margin-bottom:0'>{_esc(verd['caveat'])}</p>" if verd.get('caveat') else ""}
</div>

<div class="headline">
  <div class="stat"><b>{loop['found_of_top']}/{loop['top_n']}</b>
    <span>of your best compounds found</span></div>
  <div class="stat"><b>{loop['tested_total']}</b>
    <span>experiments used ({loop['fraction_tested']:.0%} of the programme)</span></div>
  <div class="stat"><b>{_fmt(verd['lift_vs_random'], 2)}×</b>
    <span>better than testing in random order</span></div>
  <div class="stat"><b>{_fmt(splits.get('scaffold', {}).get('spearman'))}</b>
    <span>ρ on unseen chemotypes</span></div>
</div>

{f'<div class="note"><b>About this analysis.</b> {_esc(provenance)}</div>' if provenance else ""}

<h2>1 · What we read from your file</h2>
<table><tbody>
<tr><td>Rows in the file</td><td class="n">{ing.rows_read}</td></tr>
<tr><td>Usable compounds after cleaning</td><td class="n"><b>{ing.compounds}</b></td></tr>
<tr><td>Distinct chemotypes (Murcko scaffolds)</td><td class="n">{ing.scaffolds}</td></tr>
<tr><td>Potency range (p-scale, higher = more potent)</td><td class="n">{ing.value_range[0]} – {ing.value_range[1]}</td></tr>
<tr><td>Dates covered</td><td class="n">{_esc(f"{_year_label(ing.date_range[0])} – {_year_label(ing.date_range[1])}" if ing.date_range else "no date column")}</td></tr>
<tr><td>Columns used</td><td class="n">{_esc(", ".join(f"{k}={v}" for k, v in ing.columns.items() if v))}</td></tr>
<tr><td>Units detected</td><td class="n">{_esc(ing.units)}</td></tr>
</tbody></table>

<h3>What we removed, and why</h3>
{f"<table><thead><tr><th>Reason</th><th class='n'>Rows</th></tr></thead><tbody>{dropped}</tbody></table>"
 if dropped else "<p class='muted'>Nothing was removed.</p>"}
{warnings}

<h2>2 · Can a model rank your compounds?</h2>
<p>Three tests, from easiest to hardest. Only the last two matter for deciding whether
to use a model on compounds that do not exist yet.</p>
{_splits_table(splits)}
<div class="note">ρ (Spearman) measures whether the model puts compounds in the right
<em>order</em> — the thing a design loop actually consumes. R² measures absolute accuracy
and is routinely poor even when ranking is good, because the model is asked about
compounds unlike anything it was trained on.</div>

<h2>3 · Replaying your programme</h2>
<p>Starting from nothing, the system picked {loop['batch_size']} compounds at a time.
It was shown a result only after choosing that compound — exactly the position a live
project is in. Everything below is compared against what your programme actually found.</p>
{chart}
<p>{saved}</p>
{found_chart}
<table><tbody>
<tr><td>Best potency in the dataset</td><td class="n">{loop['best_possible']}</td></tr>
<tr><td>Best the loop reached</td><td class="n"><b>{loop['best_found']}</b></td></tr>
<tr><td>Best expected from random order, same budget</td><td class="n">{loop['random_expected_best']}</td></tr>
<tr><td>Top {loop['top_n']} compounds found</td><td class="n"><b>{loop['found_of_top']}</b></td></tr>
<tr><td>Expected by chance, same budget</td><td class="n">{loop['random_found_expected']}</td></tr>
<tr><td>Ranking agreement in later cycles (ρ)</td><td class="n">{_fmt(loop['late_rho'])}</td></tr>
</tbody></table>

<h2>4 · How to read this honestly</h2>
<div class="note">
<p style="margin-top:0"><b>A replay flatters any method.</b> Every compound in this file was made
by a chemist who thought it worth making. A live project also faces compounds nobody
chose, so expect real numbers below these.</p>
<p><b>Ranking, not invention.</b> Where the time split is weak, the model optimises within
chemistry it has seen; it does not propose new chemotypes.</p>
<p style="margin-bottom:0"><b>One dataset is one dataset.</b> A second programme is what turns
this from an anecdote into a claim.</p>
</div>

{_cost_block(ev)}
<footer>
Generated by the DMTA platform from {_esc(Path(ing.source).name)} on {today}.
Model: random forest on Morgan fingerprints and physicochemical descriptors,
{ev.settings['n_estimators']} trees, uncertainty calibrated against out-of-fold residuals.
Every figure here is reproducible from the same file with
<code>python3 -m dmta.cli evaluate --csv &lt;file&gt;</code>.
</footer>
</div></body></html>"""


def write_report(ev: Evaluation, path: str | Path, title: str = "Retrospective evaluation",
                 customer: str | None = None, provenance: str | None = None) -> Path:
    path = Path(path)
    path.write_text(render(ev, title=title, customer=customer, provenance=provenance),
                    encoding="utf-8")
    return path
