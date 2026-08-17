#!/usr/bin/env python3
"""The search budget: how many independent bump hunts the public BSM resonance space implies,
and what local significance a 5 sigma GLOBAL discovery in that program therefore costs.

Reads nothing but the two data-free modules (bump_observables.py, public_obs_map.py) -- no
sample catalogue, no dataset identifiers. Pure standard library, so it runs anywhere.

Writes results/tables/search_budget.csv, results/tables/search_budget_selections.csv and
results/overviews/SEARCH_BUDGET.md (the headline). Figures: budget_plots.py, budget_waterfall.py.
"""
import os, csv, math, collections, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def _p(*a): return os.path.join(ROOT, *a)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import (floor, res, canon, CANON_ORDER, scan_segments, scan_source,
                              ns_scan, ns_achievable, z_local_for_global5 as z5, SQRTS)
from public_obs_map import PUBLIC_OBS, NSEL, NSEL_DEFAULT, nsel, WIDTH, nonpeak_only

# ---------------------------------------------------------------- the public channel set
# One canonical observable per scanned mass axis; a model contributes to every axis it peaks in.
pub_models = collections.defaultdict(set)
for model, obss in PUBLIC_OBS.items():
    for o in obss:
        pub_models[canon(o)].add(model)

order = [o for o in CANON_ORDER if o == canon(o) and o in pub_models]
order += [o for o in sorted(pub_models) if o not in order]

# ---------------------------------------------------------------- the three granularity levels
incl = {o: ns_scan(o) for o in order}                       # published window, 1 look / spectrum
sel  = {o: nsel(o) * incl[o] for o in order}                # x published event selections
env  = {o: ns_achievable(o) for o in order}                 # reference: floor -> sqrt(s)

N_incl, n_incl = sum(incl.values()), len(order)
N_sel,  n_sel  = sum(sel.values()),  sum(nsel(o) for o in order)
N_env          = sum(env.values())
N_full, N_full_lo, N_full_hi = 5e4, 1e4, 1e5                # full ATLAS BSM program (literature)

# ---------------------------------------------------------------- the narrow-resonance assumption
wclass = collections.Counter(c for c, _ in WIDTH.values())
nonpeak_axes = [o for o in order if nonpeak_only(o, pub_models[o])]
N_peak = N_incl - sum(incl[o] for o in nonpeak_axes)

# ---------------------------------------------------------------- console
def band(N):
    return (f"N = {N:,.0f}  (r x0.5..x2 -> {N*0.5:,.0f}-{N*2:,.0f});  "
            f"Z_local(5s global) = {z5(N):.2f}  ({z5(N*0.5):.2f}-{z5(N*2):.2f})")

print(f"public bump spectra: {n_incl}   public model classes: {len(PUBLIC_OBS)}")
print("inclusive (1/spectrum) :", band(N_incl))
print("event selections       :", band(N_sel), f"[{n_sel} channels]")
print("kinematic envelope     :", band(N_env), "(reference bound)")
print("\ntop contributors (inclusive n_s):")
for o in sorted(order, key=lambda o: -incl[o])[:8]:
    print(f"  {o:16s} r={res(o):5.3f}  n_s={incl[o]:5.0f}  x{nsel(o)} selections")

# ---------------------------------------------------------------- CSV: per spectrum
def segs(o): return "+".join(f"{lo:g}-{hi:g}" for lo, hi in scan_segments(o))

with open(_p("results", "tables", "search_budget.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["observable", "r", "floor_GeV", "scan_window_GeV", "ns_scan", "ns_envelope",
                "n_models_public", "n_event_selections", "ns_with_selections", "scan_source"])
    for o in order:
        w.writerow([o, res(o), f"{floor(o):g}", segs(o), f"{incl[o]:.1f}", f"{env[o]:.1f}",
                    len(pub_models[o]), nsel(o), f"{sel[o]:.1f}", scan_source(o)])

# ---------------------------------------------------------------- CSV: event selections
with open(_p("results", "tables", "search_budget_selections.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["observable", "n_event_selections", "ns_inclusive", "ns_with_selections", "channels"])
    for o in order:
        w.writerow([o, nsel(o), f"{incl[o]:.1f}", f"{sel[o]:.1f}",
                    NSEL.get(o, (NSEL_DEFAULT, "default"))[1]])

# ---------------------------------------------------------------- markdown
def md_spectra():
    lines = ["| # | spectrum (bump observable) | r | scan window [GeV] | n_s | n_s envelope | "
             "#models | #selections | n_s x sel |",
             "|--:|---|--:|---|--:|--:|--:|--:|--:|"]
    for n, o in enumerate(order, 1):
        lines.append(f"| {n} | `{o}` | {res(o):g} | {' + '.join(f'{lo:g}-{hi:g}' for lo, hi in scan_segments(o))} | "
                     f"{incl[o]:.0f} | {env[o]:.0f} | {len(pub_models[o])} | {nsel(o)} | {sel[o]:.0f} |")
    lines.append(f"| | **total ({n_incl} spectra)** | | | **{N_incl:,.0f}** | {N_env:,.0f} | | "
                 f"**{n_sel}** | **{N_sel:,.0f}** |")
    return "\n".join(lines)

def md_levels():
    lines = ["| granularity | # spectra | N_trials | band (r x0.5..x2) | Z_local for 5s global | band |",
             "|---|--:|--:|---|--:|---|"]
    for name, nsp, N in [("inclusive (1 spectrum / observable)", n_incl, N_incl),
                         ("**published event selections**", n_sel, N_sel),
                         ("kinematic envelope (reference bound)", n_incl, N_env),
                         ("full ATLAS BSM program (literature)", "-", N_full)]:
        lines.append(f"| {name} | {nsp} | **{N:,.0f}** | {N*0.5:,.0f}-{N*2:,.0f} | "
                     f"**{z5(N):.2f}** | {z5(N*0.5):.2f}-{z5(N*2):.2f} |")
    return "\n".join(lines)

md = f"""# Search budget: model coverage x clustering -> the LEE price of the program

**One bump hunt = one invariant-mass spectrum.** Every model that peaks in the same spectrum is
tested by the same search, so the *number of searches* = the number of distinct bump observables,
and the global-significance cost follows the Look-Elsewhere relation
`Z_global^2 = Z_local^2 - 2 ln N`.

Public information only: the channel set comes from `public_obs_map.PUBLIC_OBS` (public BSM model
classes -> the spectra they populate) and the windows from `bump_observables.SCAN` (published ATLAS
search families, one source note per channel). No sample catalogue enters.

Scope: **invariant/transverse-mass bump hunts only**. The dedicated non-bump programs (displaced
HNL/Zd/RPV/ALP, MET+jet, dE/dx monopole, EFT off-shell N) carry their own look-elsewhere effect and
are **not** summed here.

## Method
Constant fractional resolution `sigma_M = r*M` -> effective independent looks across a scanned
window `[M_lo, M_hi]` is `n_s = (1/r) ln(M_hi/M_lo)` (summed over disjoint scan segments);
`N_trials = sum_spectra n_s`; the local significance for a **5 sigma global** discovery is
`Z_local = sqrt(25 + 2 ln N_trials)`. The fractional resolution `r` is the one real physics input,
so every headline is quoted as a **band**: `r` x2 (-> N x0.5) and x0.5 (-> N x2).

**Window convention.** Every spectrum is counted over its **published-search scan window**
(`bump_observables.SCAN`, with a per-channel source note), never over a range inferred from which
signal samples happen to exist -- that would understate sparsely sampled channels and make any
comparison against a sample catalogue circular.

**What counts as one spectrum.** Lepton flavour and b-jet content are part of the mass axis, so
`m(ee)` and `m(mumu)`, and likewise `m(eb)` and `m(ej)`, are separate spectra with their own
resolution and window. Because flavour is an axis, it is deliberately **not** also one of the
event-selection multipliers below; counting both would double-count it.

## Trials per spectrum
`n_s` = published-search window (the headline); `n_s envelope` = kinematic reference
(analyzable floor -> sqrt(s) = {SQRTS/1000:g} TeV); `#selections` = published event selections that
scan this axis separately (`public_obs_map.NSEL`).

{md_spectra()}

## Summary
{md_levels()}

**Reading it.** Covering every bump channel that public BSM models motivate costs
**N_trials = {N_incl:,.0f}** over {n_incl} spectra: a local 5 sigma degrades to
~{math.sqrt(max(0.0, 25 - 2*math.log(N_incl))):.1f} sigma global, and a 5 sigma-global discovery
needs local **~{z5(N_incl):.2f} sigma**. Slicing at the granularity real searches actually use
({n_sel} event selections) raises N to {N_sel:,.0f} and the bar only to {z5(N_sel):.2f}. Because N
enters through `ln N`, every level sits within ~0.1 sigma of 6.5 -- **breadth is cheap**, and the
budget is extremely robust to counting choices.

## Assumptions & caveats
Each of these is varied and priced in `results/overviews/BUDGET_UNCERTAINTY.md`
(`scripts/budget_uncertainty.py`), which carries the band on every number above and is where the
`r` x0.5..x2 band quoted here is only one line: the largest term there is not a physics input but
the convention that turns a resolution element into an independent look.

- **Resolution dict `r` is coarse** (per-channel central values; headline carries the x0.5..x2
  band). A factor-2 error in `r` moves `Z_local` by ~+-0.1 sigma only (it enters via `ln N`).
- **Scan windows are hand-curated** from published ATLAS search families (source column in the
  CSV), extended where the search family's signal grid prepares a wider scan.
- **Same-axis merges**: `m(HH) 4b` counts as an extra *event selection* of the `m(HH)` axis, not as
  an independent spectrum. The `(Zd)` dark-photon axes (0.3-400 GeV) and the high-mass dilepton
  axes (150 GeV-8 TeV) remain separate spectra (different selections; their 150-400 GeV overlap
  double-counts ~50 looks, a <2% excess).
- **Lepton flavour is a spectrum, not a selection**: `m(ee)` and `m(mumu)` are counted separately
  with their OWN resolution (EM-calorimeter ~1.5% vs a sagitta-limited muon measurement averaging
  ~5% over this window), because they are separate analyses with separate triggers.
- `NSEL` is hand-curated from published ATLAS channel counts. Channels within one search are partly
  correlated and often statistically combined, so treating each as an independent scan is a mild
  over-count: the true N sits between the inclusive and the selections level.
- `n_s` ignores cross-channel correlations (conservative: slight over-count) and uses the
  fixed-resolution-element approximation of Gross-Vitells (the up-crossing refinement adds a mild
  Z-dependence).
- **Narrow-resonance assumption.** `n_s` counts *detector* resolution elements, which is the right
  step size only while the natural width stays below `r`. `public_obs_map.WIDTH` places all
  {len(WIDTH)} public model classes: **{wclass['narrow']} narrow** on every axis they populate,
  **{wclass['benchmark']}** narrow only at the benchmark ATLAS publishes (Z'_SSM at 3% against
  `r`=0.015 on `m(ee)`, the DM-mediator coupling grid, single-VLQ, U1 at the flavour-anomaly point),
  **{wclass['broad']}** already broader than `r` there (KK gluon 15-30% vs `r`=0.08, coloron/axigluon,
  composite/NJL), and **{wclass['nonpeak']}** with no Breit-Wigner peak at all (QBH thresholds,
  ADD/HEIDI continua, pair-produced Type-III/VLL, toponium at 2 m_t). A signal wider than `r`
  correlates neighbouring mass points, so counting resolution elements **over**-counts independent
  looks: the bias is conservative, `Z_local` too strict rather than too loose.
  Only {len(nonpeak_axes)} axes are motivated *exclusively* by non-peaking models
  ({', '.join(f'`{o}`' for o in nonpeak_axes)}); dropping both takes N from {N_incl:,.0f} to
  {N_peak:,.0f} and `Z_local` from {z5(N_incl):.2f} to {z5(N_peak):.2f}, so the whole question is
  worth {z5(N_incl)-z5(N_peak):.3f} sigma. `NONPEAK_ON` records the (model, axis) pairs where an
  otherwise-narrow class does not peak, including the H/A interference with SM ttbar.

Source: `scripts/search_budget.py` -> `results/tables/search_budget.csv`,
`results/tables/search_budget_selections.csv`. Figures: `scripts/budget_plots.py`,
`scripts/budget_waterfall.py`. Windows and their sources: `scripts/bump_observables.py` (SCAN).
Excess bookkeeping: `results/overviews/EXCESS_COUNTING.md`.
Uncertainty budget: `results/overviews/BUDGET_UNCERTAINTY.md`.
"""
open(_p("results", "overviews", "SEARCH_BUDGET.md"), "w").write(md)
print("\nwrote results/tables/search_budget.csv, results/tables/search_budget_selections.csv, "
      "results/overviews/SEARCH_BUDGET.md")
