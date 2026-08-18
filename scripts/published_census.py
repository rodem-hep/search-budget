#!/usr/bin/env python3
import os, csv, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def _p(*a): return os.path.join(ROOT, *a)

rows = list(csv.DictReader(open(_p("data", "published_spectra.csv"))))
papers = {a for r in rows for a in r["arxiv"].split()}
by_fam = collections.Counter(r["family"] for r in rows)
pap_fam = collections.Counter()
for r in rows:
    pap_fam[r["family"]] += int(r["n_papers"])

STATUS = ("current", "ageing", "stale")
by_status = collections.Counter(r["status"] for r in rows)
run3 = [r for r in rows if r["run3"] == "yes"]
stale = sorted((r for r in rows if r["status"] == "stale"), key=lambda r: r["last_year"])
on_axis = [r for r in stale if r["budget_axis"] != "-"]
off_axis = [r for r in stale if r["budget_axis"] == "-"]
observables = {r["observable"] for r in rows}

print(f"catalogued spectra      : {len(rows)}")
print(f"distinct bump observables: {len(observables)}")
print(f"search papers           : {len(papers)}  ({min(r['first_year'] for r in rows)}"
      f"-{max(r['last_year'] for r in rows)})")
print(f"recency                 : {by_status['current']} current (2024+), "
      f"{by_status['ageing']} ageing (2019-2023), {by_status['stale']} stale (pre-2019)")
print(f"with a Run-3 result     : {len(run3)}")
print("\nby family:")
for f, n in by_fam.most_common():
    print(f"  {f:42s} {n:3d} spectra   {pap_fam[f]:3d} papers")
print(f"\nnot updated since before 2019 ({len(stale)}): "
      f"{len(on_axis)} on an axis already counted, {len(off_axis)} would add one")
for r in stale:
    ax = r["budget_axis"]
    print(f"  {r['last_year']}  {r['spectrum']:52s} {r['observable']:18s} "
          f"{'-> ' + ax if ax != '-' else 'adds an axis'}")

with open(_p("results", "tables", "published_census.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["family", "spectra", "papers", "current", "ageing", "stale", "run3"])
    for fam in sorted(by_fam):
        sel = [r for r in rows if r["family"] == fam]
        w.writerow([fam, len(sel), pap_fam[fam]] +
                   [sum(1 for r in sel if r["status"] == s) for s in STATUS] +
                   [sum(1 for r in sel if r["run3"] == "yes")])
    w.writerow(["TOTAL", len(rows), len(papers)] +
               [by_status[s] for s in STATUS] + [len(run3)])

md = f"""# The publication-side census

Every ATLAS resonance search, grouped by the mass spectrum it scans. Assembled from the
collaboration's publication record on INSPIRE-HEP and curated by hand; the per-spectrum rows, with
their arXiv references, are in `data/published_spectra.csv`. Scope matches the budget's: bump hunts
for new states, so hadron-spectroscopy measurements are out even where they are bump hunts.

This is the complement to `SEARCH_BUDGET.md`. That counts the spectra public BSM models motivate
(46 canonical mass axes); this counts the searches that have actually been published
({len(rows)} entries over {len(papers)} papers). **The two use different bases and must not be
summed**: the publication record separates entries that share a mass axis when they are different
analyses ({len(observables)} distinct bump observables appear across the {len(rows)} entries), while the
budget merges them onto one axis and counts resolution elements along it.

| | |
|---|---|
| catalogued spectra | **{len(rows)}** |
| search papers | **{len(papers)}** ({min(r['first_year'] for r in rows)}-{max(r['last_year'] for r in rows)}) |
| current (latest paper 2024 or later) | {by_status['current']} |
| ageing (2019-2023) | {by_status['ageing']} |
| **stale (nothing since before 2019)** | **{by_status['stale']}** |
| with a published Run-3 (13.6 TeV) result | **{len(run3)}** |

## By final state

| family | spectra | papers | current | ageing | stale | Run-3 |
|---|--:|--:|--:|--:|--:|--:|
"""
for fam in sorted(by_fam):
    sel = [r for r in rows if r["family"] == fam]
    md += (f"| {fam} | {len(sel)} | {pap_fam[fam]} | "
           + " | ".join(str(sum(1 for r in sel if r["status"] == s)) for s in STATUS)
           + f" | {sum(1 for r in sel if r['run3'] == 'yes')} |\n")
md += (f"| **total** | **{len(rows)}** | **{len(papers)}** | "
       + " | ".join(f"**{by_status[s]}**" for s in STATUS) + f" | **{len(run3)}** |\n")

md += f"""
## Not revisited since before 2019

These carry Run-1 or early-Run-2 sensitivity. What the budget says about them splits them in two.

**{len(on_axis)} of the {len(stale)} sit on a mass axis that is already counted in `N`, so revisiting one costs
nothing in trials** -- the discovery bar for the re-run is the bar the program already pays, and the
whole cost is analysis effort. The remaining {len(off_axis)} fall on no axis in the budget's 46, so re-running
one extends the axis count rather than reusing it, and it is priced like any other new spectrum.

| last published | spectrum | observable | counted axis |
|--:|---|---|---|
"""
for r in stale:
    ax = f"`{r['budget_axis']}`" if r["budget_axis"] != "-" else "**adds an axis**"
    md += f"| {r['last_year']} | {r['spectrum']} | `{r['observable']}` | {ax} |\n"

md += f"""
## Already re-run at 13.6 TeV

| spectrum | observable |
|---|---|
"""
for r in sorted(run3, key=lambda r: r["spectrum"]):
    md += f"| {r['spectrum']} | `{r['observable']}` |\n"

md += """
Source: `scripts/published_census.py` from `data/published_spectra.csv`. What this record costs in
trials, entry by entry: `results/overviews/CENSUS_BUDGET.md` (`scripts/census_budget.py`).
"""
open(_p("results", "overviews", "PUBLISHED_CENSUS.md"), "w").write(md)
print("\nwrote results/tables/published_census.csv, results/overviews/PUBLISHED_CENSUS.md")
