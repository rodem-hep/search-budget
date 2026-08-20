import collections

from .. import io, paths
from ..core.bump_observables import (RES_DEFAULT, SCAN, canon, n_s, res, scan_segments,
                                     z_local_for_global5 as z5)
from ..core.lee import band, merge_segments
from ..registry import stage

FIXED = "fixed"
OFF = "-"


@stage(
    name="census-budget",
    group="census",
    summary="the publication record priced in trials, per search and per axis",
    outputs=["tables/census_budget.csv", "overviews/CENSUS_BUDGET.md"],
    inputs=["data/published_spectra.csv"],
    needs=["tables/search_budget.csv"],
)
def main(options=None):
    rows = io.read_rows(paths.data("published_spectra.csv"))

    def axes_of(row):
        return [] if row["budget_axis"] == OFF else [canon(a.strip())
                                                     for a in row["budget_axis"].split(";")]

    for r_ in rows:
        for a in axes_of(r_):
            assert a in SCAN, f"{r_['spectrum']}: unknown axis {a!r}"

    def segments(row, axis):
        s = row["scan_GeV"]
        if s == FIXED:
            return []
        if s:
            return [tuple(float(x) for x in seg.split("-")) for seg in s.split("+")]
        return scan_segments(axis) if axis else []

    entries = []
    for row in rows:
        axes = axes_of(row)
        for axis in (axes or [None]):
            r = res(axis) if axis else RES_DEFAULT
            segs = segments(row, axis)
            if row["scan_GeV"] == FIXED:
                ns = 1.0
            else:
                ns = sum(n_s(lo, hi, r) for lo, hi in segs)
            entries.append((row, axis, r, segs, ns))

    priced = [e for e in entries if e[4] > 0]
    unpriced = [e for e in entries if e[4] == 0]
    on_axis = [e for e in priced if e[1]]
    off_axis = [e for e in priced if not e[1]]
    n_fixed = sum(1 for e in priced if e[0]["scan_GeV"] == FIXED)

    N_entry = sum(e[4] for e in priced)

    by_axis = collections.defaultdict(list)
    for row, axis, _r, segs, _ns in priced:
        if axis:
            by_axis[axis] += segs
    union = {a: merge_segments(s) for a, s in by_axis.items()}
    N_axes = sum(sum(n_s(lo, hi, res(a)) for lo, hi in segs) for a, segs in union.items())
    N_axes += sum(e[4] for e in off_axis)
    covered, uncovered = sorted(union), sorted(set(SCAN) - set(union))
    no_axis_declared = [e for e in unpriced
                        if e[0]["observable"].strip().lower().startswith("many")]
    no_range = [e for e in unpriced if e not in no_axis_declared]
    observables = len(covered) + len(off_axis) + len(unpriced)
    N_off = sum(e[4] for e in off_axis)

    model = {}
    try:
        for r_ in io.read_rows(paths.table("search_budget.csv")):
            model[r_["observable"]] = float(r_["ns_scan"])
    except FileNotFoundError:
        pass
    N_model = sum(model.values()) or None

    N_uncov = sum(model.get(a, 0.0) for a in uncovered)
    why_union = "" if not N_model else f"""**Why the union basis lands on the model space.** Counting each axis once gives {N_axes:,.0f}
against the model space's {N_model:,.0f}, a gap of {abs(N_axes-N_model)/N_model:.1%} -- closer than
it has any right to be, because two omissions cancel. The union misses the {len(uncovered)} model
axes with no published search, worth {N_uncov:,.0f} looks on the
model side, and it adds {len(off_axis)} scanned axes outside the {len(SCAN)}, worth {N_off:,.0f}. Charging
neither leaves {N_axes-N_off:,.0f} against {N_model - N_uncov:,.0f}
on the same {len(covered)} axes. The robust statement is the ten-percent one, and it is not a
tautology: the model space takes its windows from published searches but its axes from the models,
so the agreement says the program scans almost exactly the axes the models motivate."""

    print(f"census entries: {len(rows)}   priced pairs: {len(priced)} "
          f"({len(on_axis)} on a budget axis, {len(off_axis)} off-axis, {n_fixed} fixed-mass)")
    print(f"unpriced      : {len(unpriced)} ({len(no_range)} publish no range, "
          f"{len(no_axis_declared)} declare no single axis)")
    print(f"observables   : {observables} = {len(covered)} budget axes + {len(off_axis)} off-axis "
          f"with a range + {len(unpriced)} unchargeable")
    print("published searches :", band(N_entry))
    print("axes scanned (union):", band(N_axes), f"[{len(covered)} of {len(SCAN)} axes covered]")
    if N_model:
        print(f"model side (reference): N = {N_model:,.0f} over {len(model)} axes, "
              f"Z_local = {z5(N_model):.2f}")
    print("\ntop contributors (published searches):")
    for row, axis, r, segs, ns in sorted(priced, key=lambda e: -e[4])[:8]:
        print(f"  {row['spectrum'][:46]:46s} {str(axis):16s} r={r:5.3f}  n_s={ns:5.0f}")
    print(f"\naxes with no published search ({len(uncovered)}): {', '.join(uncovered)}")

    def segstr(segs):
        return "+".join(f"{lo:g}-{hi:g}" for lo, hi in segs) or FIXED

    table = []
    for row, axis, r, segs, ns in entries:
        src = (FIXED if row["scan_GeV"] == FIXED else
               "census" if row["scan_GeV"] else
               "axis window" if axis else "unpriced")
        table.append([row["family"], row["spectrum"], row["status"], axis or OFF,
                      segstr(segs) if segs or row["scan_GeV"] == FIXED else "",
                      src, f"{r:g}", f"{ns:.1f}"])
    io.write_rows(paths.table("census_budget.csv"),
                  ["family", "spectrum", "status", "budget_axis", "window_GeV", "window_from",
                   "r", "n_s"], table)

    def md_entries():
        lines = ["| family | published search | axis | window [GeV] | from | r | n_s |",
                 "|---|---|---|---|---|--:|--:|"]
        for row, axis, r, segs, ns in entries:
            src = (FIXED if row["scan_GeV"] == FIXED else
                   "census" if row["scan_GeV"] else "axis" if axis else "**unpriced**")
            win = segstr(segs) if segs else ("single mass" if row["scan_GeV"] == FIXED else "-")
            lines.append(f"| {row['family']} | {row['spectrum']} | "
                         f"{'`'+axis+'`' if axis else '-'} | {win} | {src} | {r:g} | "
                         f"{ns:.0f} |")
        lines.append(f"| | **total ({len(rows)} entries, {len(priced)} priced pairs)** | | | | | "
                     f"**{N_entry:,.0f}** |")
        return "\n".join(lines)

    def ns_union(a):
        return sum(n_s(lo, hi, res(a)) for lo, hi in union[a])

    def md_axes():
        lines = ["| axis | scanned range(s) [GeV] | r | n_s | model-side n_s |",
                 "|---|---|--:|--:|--:|"]
        for a in sorted(union, key=lambda a: -ns_union(a)):
            ref = f"{model[a]:.0f}" if a in model else "-"
            lines.append(f"| `{a}` | {segstr(union[a])} | {res(a):g} | {ns_union(a):.0f} | "
                         f"{ref} |")
        return "\n".join(lines)

    def md_levels():
        lines = ["| basis | units | N_trials | band (r x0.5..x2) | Z_local for 5s global | band |",
                 "|---|--:|--:|---|--:|---|"]
        rows_ = [(f"**published searches** ({len(rows)} census entries)", len(priced), N_entry),
                 ("axes scanned (union of the ranges)", len(covered) + len(off_axis), N_axes)]
        if N_model:
            rows_.append(("model-motivated axes (reference, `SEARCH_BUDGET.md`)",
                          len(model), N_model))
        for name, u, N in rows_:
            lines.append(f"| {name} | {u} | **{N:,.0f}** | {N*0.5:,.0f}-{N*2:,.0f} | "
                         f"**{z5(N):.2f}** | {z5(N*0.5):.2f}-{z5(N*2):.2f} |")
        return "\n".join(lines)

    md = f"""# The publication record, priced in trials

`SEARCH_BUDGET.md` counts the spectra public BSM models motivate. This prices the other side: the
searches ATLAS has actually **published**, taken from the census (`data/published_spectra.csv`,
{len(rows)} entries over {len({a for r_ in rows for a in r_['arxiv'].split()})} papers) and run through
the same rule, `n_s = (1/r) ln(M_hi/M_lo)`, `Z_local = sqrt(25 + 2 ln N)`.

## Method
Every census entry carries the canonical budget axis it scans (`budget_axis`; several when the
entry scans several; `-` when it falls on none of the {len(SCAN)}) and the range it scanned (`scan_GeV`,
transcribed from the published range the census records). The resolution `r` is that axis' -- the
budget's one physics input, unchanged here; an entry on no axis is priced at the default
`r = {RES_DEFAULT:g}`. **Where the census does not record a range, the axis' own published window is the
fallback**, which is why the two bases below differ.

| basis | what one look is |
|---|---|
| published searches | one census entry on one axis, over the range that entry scanned |
| axes scanned | one axis, over the **union** of every published range on it, counted once |

A fixed-mass search (LFV `Z` and `tau` decays, exclusive `H`/`Z -> quarkonium + gamma`,
`H -> Z gamma`) scans nothing and contributes exactly one look.

## Summary
{md_levels()}

**Reading it.** Pricing every published ATLAS resonance search over the range it actually scanned
gives **N = {N_entry:,.0f}**, so a 5 sigma global discovery in the published program needs a local
**{z5(N_entry):.2f} sigma**. Counting each axis once instead, over the union of everything published on
it, gives **N = {N_axes:,.0f}** and **{z5(N_axes):.2f} sigma**. The two bases differ by a factor
{N_entry/N_axes:.1f} in N and {z5(N_entry)-z5(N_axes):.2f} sigma in the bar -- the same lesson as the
model side: the answer is nearly independent of how finely the program is sliced.

## Coverage
The census reaches {len(covered)} of the {len(SCAN)} canonical axes. {len(uncovered)} carry no published
search at all: {', '.join(f'`{a}`' for a in uncovered)}.

{len(off_axis)} published entries fall on no canonical axis and are priced on their own published range at
`r = {RES_DEFAULT:g}`; {len(unpriced)} more fall on no axis **and** carry no chargeable range, so
they are listed below as unpriced and are missing from `N` -- `N` is a lower bound by that much. Of
those {len(unpriced)}, {len(no_range)} published no scanned range and {len(no_axis_declared)} declare
no single axis to scan (the anomaly-detection and generic multi-body entries, whose observable is
"many").

The {observables} distinct bump observables the census scans are therefore
{len(covered)} of the {len(SCAN)} budget axes, {len(off_axis)} scanned axes outside them, and
{len(unpriced)} carrying nothing chargeable.

{why_union}

## Per axis (union of every published range)
{md_axes()}

## Per published search
{md_entries()}

## Assumptions & caveats
- **The axis assignment is a curated judgement**, recorded per row in `data/published_spectra.csv`
  (`budget_axis`) so it can be checked and changed. Entries whose final state spans several axes
  (LFV dilepton, leptoquark pair, single VLQ, doubly charged Higgs) contribute one look per axis.
- **The window fallback over-counts.** {sum(1 for r_ in rows if not r_['scan_GeV'])} entries carry no published range in
  the census and inherit their axis' full window, so on the published-searches basis several
  analyses on one axis are each charged the whole axis. The union basis is free of this and is the
  conservative reading; the truth sits between them.
- **Off-axis entries are priced at the default resolution** `r = {RES_DEFAULT:g}`, not at a measured one.
- **Scope**: the same as the budget's -- invariant/transverse-mass bump hunts. Census entries that
  are themselves multi-spectrum scans (the anomaly-detection and generic multi-body papers) and the
  displaced programs fall on no axis and stay unpriced: their trials belong in the combinatorial
  count of `scaled_scan.py`, not here.
- Everything the model-side budget assumes about resolution elements, correlations and the
  narrow-resonance approximation applies unchanged; see `SEARCH_BUDGET.md`.

Source: `searchbudget/stages/census_budget.py` from `data/published_spectra.csv` ->
`results/tables/census_budget.csv`. The census itself: `results/overviews/PUBLISHED_CENSUS.md`.
The model-side budget: `results/overviews/SEARCH_BUDGET.md`.
"""
    io.write_text(paths.overview("CENSUS_BUDGET.md"), md)
    print("\nwrote results/tables/census_budget.csv, results/overviews/CENSUS_BUDGET.md")
