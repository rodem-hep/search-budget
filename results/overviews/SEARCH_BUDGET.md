# Search budget: model coverage x clustering -> the LEE price of the program

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
(analyzable floor -> sqrt(s) = 13.6 TeV); `#selections` = published event selections that
scan this axis separately (`public_obs_map.NSEL`).

| # | spectrum (bump observable) | r | scan window [GeV] | n_s | n_s envelope | #models | #selections | n_s x sel |
|--:|---|--:|---|--:|--:|--:|--:|--:|
| 1 | `m(gammagamma)` | 0.01 | 66-110 + 150-5000 | 402 | 583 | 5 | 2 | 803 |
| 2 | `m(mumu) (Zd)` | 0.02 | 0.3-400 | 360 | 536 | 2 | 1 | 360 |
| 3 | `m(ee) SS` | 0.015 | 200-1300 | 125 | 454 | 2 | 1 | 125 |
| 4 | `m(ee) (Zd)` | 0.025 | 1-400 | 240 | 381 | 1 | 1 | 240 |
| 5 | `m(egamma)` | 0.015 | 200-5000 | 215 | 328 | 1 | 1 | 215 |
| 6 | `m(ee)` | 0.015 | 150-8000 | 265 | 321 | 7 | 1 | 265 |
| 7 | `m(emu) SS` | 0.025 | 200-1300 | 75 | 272 | 2 | 1 | 75 |
| 8 | `m(emu) LFV` | 0.03 | 100-8000 | 146 | 234 | 4 | 1 | 146 |
| 9 | `m(ej)` | 0.03 | 200-5000 | 107 | 234 | 3 | 4 | 429 |
| 10 | `m(muj)` | 0.04 | 200-5000 | 80 | 176 | 3 | 4 | 322 |
| 11 | `m(mumu) SS` | 0.04 | 200-1300 | 47 | 170 | 2 | 1 | 47 |
| 12 | `m(eb)` | 0.04 | 200-2000 | 58 | 170 | 3 | 2 | 115 |
| 13 | `m(eZ)` | 0.03 | 100-1100 | 80 | 164 | 1 | 1 | 80 |
| 14 | `m(mugamma)` | 0.03 | 200-5000 | 107 | 164 | 1 | 1 | 107 |
| 15 | `m(Vgamma)` | 0.03 | 220-6000 | 110 | 141 | 3 | 2 | 220 |
| 16 | `m(cb) dijet` | 0.05 | 20-1000 | 78 | 136 | 1 | 2 | 156 |
| 17 | `m(mub)` | 0.05 | 200-2000 | 46 | 136 | 3 | 2 | 92 |
| 18 | `m(muZ)` | 0.04 | 100-1100 | 60 | 123 | 1 | 1 | 60 |
| 19 | `m(jgamma)` | 0.04 | 500-7000 | 66 | 113 | 2 | 1 | 66 |
| 20 | `m(mumu)` | 0.05 | 150-8000 | 80 | 96 | 7 | 1 | 80 |
| 21 | `m(jj)` | 0.05 | 200-8000 | 74 | 84 | 8 | 2 | 148 |
| 22 | `m(3j)` | 0.05 | 200-1800 | 44 | 84 | 2 | 1 | 44 |
| 23 | `m(eejj)` | 0.05 | 400-7000 | 57 | 76 | 2 | 1 | 57 |
| 24 | `m(VV)` | 0.06 | 200-6000 | 57 | 74 | 12 | 9 | 510 |
| 25 | `m(Vh)` | 0.06 | 200-5000 | 54 | 70 | 2 | 6 | 322 |
| 26 | `m(bj)` | 0.06 | 200-6000 | 57 | 70 | 1 | 2 | 113 |
| 27 | `m(bb)` | 0.1 | 15-62 + 450-6000 | 40 | 68 | 4 | 3 | 120 |
| 28 | `m(multi)` | 0.06 | 3000-10500 | 21 | 64 | 2 | 1 | 21 |
| 29 | `multilepton` | 0.12 | 50-10000 | 44 | 62 | 3 | 4 | 177 |
| 30 | `m(Wb)` | 0.08 | 800-3000 | 17 | 61 | 1 | 1 | 17 |
| 31 | `mT(ev)` | 0.1 | 150-7000 | 38 | 60 | 2 | 1 | 38 |
| 32 | `m(tauj)` | 0.1 | 200-5000 | 32 | 56 | 3 | 2 | 64 |
| 33 | `m(mumujj)` | 0.07 | 400-7000 | 41 | 54 | 2 | 1 | 41 |
| 34 | `m(tb)` | 0.08 | 180-6000 | 44 | 54 | 3 | 3 | 131 |
| 35 | `m(tt)/m(jj)` | 0.08 | 400-2000 | 20 | 53 | 2 | 2 | 40 |
| 36 | `m(HH)` | 0.08 | 250-6000 | 40 | 50 | 4 | 6 | 238 |
| 37 | `m(tautau)` | 0.12 | 35-6000 | 43 | 50 | 4 | 3 | 129 |
| 38 | `m(etau) LFV` | 0.12 | 100-3000 | 28 | 50 | 4 | 1 | 28 |
| 39 | `m(mutau) LFV` | 0.12 | 100-3000 | 28 | 50 | 4 | 1 | 28 |
| 40 | `m(tW)` | 0.08 | 500-3000 | 22 | 49 | 1 | 2 | 45 |
| 41 | `m(Ht)` | 0.08 | 1000-3000 | 14 | 48 | 1 | 1 | 14 |
| 42 | `m(taub)` | 0.12 | 200-2000 | 19 | 47 | 3 | 2 | 38 |
| 43 | `mT(taunu)` | 0.12 | 200-5000 | 27 | 47 | 1 | 1 | 27 |
| 44 | `m(tt)` | 0.08 | 350-6000 | 36 | 46 | 9 | 4 | 142 |
| 45 | `m(ttZ)/m(Zt)` | 0.08 | 1000-4000 | 17 | 43 | 1 | 2 | 35 |
| 46 | `mT(muv)` | 0.15 | 150-7000 | 26 | 40 | 2 | 1 | 26 |
| | **total (46 spectra)** | | | **3,685** | 6,442 | | **94** | **6,597** |

## Summary
| granularity | # spectra | N_trials | band (r x0.5..x2) | Z_local for 5s global | band |
|---|--:|--:|---|--:|---|
| inclusive (1 spectrum / observable) | 46 | **3,685** | 1,842-7,370 | **6.44** | 6.33-6.54 |
| **published event selections** | 94 | **6,597** | 3,298-13,193 | **6.53** | 6.42-6.63 |
| kinematic envelope (reference bound) | 46 | **6,442** | 3,221-12,884 | **6.52** | 6.42-6.63 |
| full ATLAS BSM program (literature) | - | **50,000** | 25,000-100,000 | **6.83** | 6.73-6.93 |

**Reading it.** Covering every bump channel that public BSM models motivate costs
**N_trials = 3,685** over 46 spectra: a local 5 sigma degrades to
~2.9 sigma global, and a 5 sigma-global discovery
needs local **~6.44 sigma**. Slicing at the granularity real searches actually use
(94 event selections) raises N to 6,597 and the bar only to 6.53. Because N
enters through `ln N`, every level sits within ~0.1 sigma of 6.5 -- **breadth is cheap**, and the
budget is extremely robust to counting choices.

## Assumptions & caveats
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

Source: `scripts/search_budget.py` -> `results/tables/search_budget.csv`,
`results/tables/search_budget_selections.csv`. Figures: `scripts/budget_plots.py`,
`scripts/budget_waterfall.py`. Windows and their sources: `scripts/bump_observables.py` (SCAN).
Excess bookkeeping: `results/overviews/EXCESS_COUNTING.md`.
