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
| 1 | `m(gammagamma)` | 0.01 | 66-110 + 150-5000 | 402 | 583 | 7 | 2 | 803 |
| 2 | `m(mumu) (Zd)` | 0.02 | 0.3-400 | 360 | 536 | 3 | 1 | 360 |
| 3 | `m(ee) SS` | 0.015 | 200-1300 | 125 | 454 | 4 | 1 | 125 |
| 4 | `m(ee) (Zd)` | 0.025 | 1-400 | 240 | 381 | 2 | 1 | 240 |
| 5 | `m(egamma)` | 0.015 | 200-5000 | 215 | 328 | 1 | 1 | 215 |
| 6 | `m(ee)` | 0.015 | 150-8000 | 265 | 321 | 7 | 1 | 265 |
| 7 | `m(emu) SS` | 0.025 | 200-1300 | 75 | 272 | 4 | 1 | 75 |
| 8 | `m(emu) LFV` | 0.03 | 100-8000 | 146 | 234 | 5 | 1 | 146 |
| 9 | `m(ej)` | 0.03 | 200-5000 | 107 | 234 | 4 | 4 | 429 |
| 10 | `m(muj)` | 0.04 | 200-5000 | 80 | 176 | 4 | 4 | 322 |
| 11 | `m(mumu) SS` | 0.04 | 200-1300 | 47 | 170 | 4 | 1 | 47 |
| 12 | `m(eb)` | 0.04 | 200-2000 | 58 | 170 | 3 | 2 | 115 |
| 13 | `m(eZ)` | 0.03 | 100-1100 | 80 | 164 | 4 | 1 | 80 |
| 14 | `m(mugamma)` | 0.03 | 200-5000 | 107 | 164 | 1 | 1 | 107 |
| 15 | `m(Vgamma)` | 0.03 | 220-6000 | 110 | 141 | 3 | 2 | 220 |
| 16 | `m(cb) dijet` | 0.05 | 20-1000 | 78 | 136 | 1 | 2 | 156 |
| 17 | `m(mub)` | 0.05 | 200-2000 | 46 | 136 | 3 | 2 | 92 |
| 18 | `m(muZ)` | 0.04 | 100-1100 | 60 | 123 | 4 | 1 | 60 |
| 19 | `m(jgamma)` | 0.04 | 500-7000 | 66 | 113 | 3 | 1 | 66 |
| 20 | `m(bZ)` | 0.05 | 200-2000 | 46 | 98 | 1 | 1 | 46 |
| 21 | `m(mumu)` | 0.05 | 150-8000 | 80 | 96 | 8 | 1 | 80 |
| 22 | `m(bgamma)` | 0.05 | 200-5000 | 64 | 90 | 1 | 1 | 64 |
| 23 | `m(jj)` | 0.05 | 200-8000 | 74 | 84 | 11 | 2 | 148 |
| 24 | `m(3j)` | 0.05 | 200-1800 | 44 | 84 | 3 | 1 | 44 |
| 25 | `m(gammajj)` | 0.05 | 200-4000 | 60 | 84 | 2 | 1 | 60 |
| 26 | `m(ejj)` | 0.05 | 250-6000 | 64 | 80 | 5 | 1 | 64 |
| 27 | `m(eejj)` | 0.05 | 400-7000 | 57 | 76 | 2 | 1 | 57 |
| 28 | `m(VV)` | 0.06 | 200-6000 | 57 | 74 | 13 | 9 | 510 |
| 29 | `m(Vh)` | 0.06 | 200-5000 | 54 | 70 | 2 | 6 | 322 |
| 30 | `m(bj)` | 0.06 | 200-6000 | 57 | 70 | 1 | 2 | 113 |
| 31 | `m(jV)` | 0.06 | 500-7000 | 44 | 70 | 2 | 1 | 44 |
| 32 | `m(bb)` | 0.1 | 15-62 + 450-6000 | 40 | 68 | 7 | 3 | 120 |
| 33 | `m(tgamma)` | 0.06 | 300-3000 | 38 | 67 | 1 | 1 | 38 |
| 34 | `m(multi)` | 0.06 | 3000-10500 | 21 | 64 | 2 | 1 | 21 |
| 35 | `multilepton` | 0.12 | 50-10000 | 44 | 62 | 4 | 4 | 177 |
| 36 | `m(Wb)` | 0.08 | 800-3000 | 17 | 61 | 1 | 1 | 17 |
| 37 | `m(taugamma)` | 0.08 | 200-5000 | 40 | 61 | 1 | 1 | 40 |
| 38 | `mT(ev)` | 0.1 | 150-7000 | 38 | 60 | 2 | 1 | 38 |
| 39 | `m(mujj)` | 0.07 | 250-6000 | 45 | 57 | 6 | 1 | 45 |
| 40 | `m(tauj)` | 0.1 | 200-5000 | 32 | 56 | 3 | 2 | 64 |
| 41 | `m(mumujj)` | 0.07 | 400-7000 | 41 | 54 | 3 | 1 | 41 |
| 42 | `m(tb)` | 0.08 | 180-6000 | 44 | 54 | 3 | 3 | 131 |
| 43 | `m(tt)/m(jj)` | 0.08 | 400-2000 | 20 | 53 | 2 | 2 | 40 |
| 44 | `m(tj)` | 0.08 | 350-3000 | 27 | 53 | 2 | 1 | 27 |
| 45 | `m(eH)` | 0.08 | 200-1500 | 25 | 53 | 2 | 1 | 25 |
| 46 | `m(muH)` | 0.08 | 200-1500 | 25 | 53 | 2 | 1 | 25 |
| 47 | `m(et)` | 0.08 | 200-2000 | 29 | 53 | 2 | 1 | 29 |
| 48 | `m(mut)` | 0.08 | 200-2000 | 29 | 53 | 2 | 1 | 29 |
| 49 | `m(HH)` | 0.08 | 250-6000 | 40 | 50 | 7 | 6 | 238 |
| 50 | `m(jH)` | 0.08 | 400-3000 | 25 | 50 | 1 | 1 | 25 |
| 51 | `m(tautau)` | 0.12 | 35-6000 | 43 | 50 | 5 | 3 | 129 |
| 52 | `m(etau) LFV` | 0.12 | 100-3000 | 28 | 50 | 5 | 1 | 28 |
| 53 | `m(mutau) LFV` | 0.12 | 100-3000 | 28 | 50 | 5 | 1 | 28 |
| 54 | `m(tW)` | 0.08 | 500-3000 | 22 | 49 | 2 | 2 | 45 |
| 55 | `m(tauV)` | 0.1 | 100-1100 | 24 | 49 | 3 | 1 | 24 |
| 56 | `m(Ht)` | 0.08 | 1000-3000 | 14 | 48 | 1 | 1 | 14 |
| 57 | `m(taub)` | 0.12 | 200-2000 | 19 | 47 | 3 | 2 | 38 |
| 58 | `mT(taunu)` | 0.12 | 200-5000 | 27 | 47 | 2 | 1 | 27 |
| 59 | `m(tt)` | 0.08 | 350-6000 | 36 | 46 | 9 | 4 | 142 |
| 60 | `m(ttZ)/m(Zt)` | 0.08 | 1000-4000 | 17 | 43 | 1 | 2 | 35 |
| 61 | `m(tbj)` | 0.1 | 400-3000 | 20 | 40 | 1 | 1 | 20 |
| 62 | `mT(muv)` | 0.15 | 150-7000 | 26 | 40 | 2 | 1 | 26 |
| 63 | `m(taujj)` | 0.1 | 300-5000 | 28 | 38 | 1 | 1 | 28 |
| | **total (63 spectra)** | | | **4,319** | 7,491 | | **111** | **7,231** |

## Summary
| granularity | # spectra | N_trials | band (r x0.5..x2) | Z_local for 5s global | band |
|---|--:|--:|---|--:|---|
| inclusive (1 spectrum / observable) | 63 | **4,319** | 2,160-8,638 | **6.46** | 6.35-6.57 |
| **published event selections** | 111 | **7,231** | 3,615-14,462 | **6.54** | 6.43-6.65 |
| kinematic envelope (reference bound) | 63 | **7,491** | 3,746-14,982 | **6.55** | 6.44-6.65 |
| full ATLAS BSM program (literature) | - | **50,000** | 25,000-100,000 | **6.83** | 6.73-6.93 |

**Reading it.** Covering every bump channel that public BSM models motivate costs
**N_trials = 4,319** over 63 spectra: a local 5 sigma degrades to
~2.9 sigma global, and a 5 sigma-global discovery
needs local **~6.46 sigma**. Slicing at the granularity real searches actually use
(111 event selections) raises N to 7,231 and the bar only to 6.54. Because N
enters through `ln N`, every level sits within ~0.1 sigma of 6.5 -- **breadth is cheap**, and the
budget is extremely robust to counting choices.

## Assumptions & caveats
Each of these is varied and priced in `results/overviews/BUDGET_UNCERTAINTY.md`
(`searchbudget/stages/budget_uncertainty.py`), which carries the band on every number above and is
where the `r` x0.5..x2 band quoted here is only one line: the largest term there is not a physics
input but the convention that turns a resolution element into an independent look.

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
  58 public model classes: **36 narrow** on every axis they populate,
  **16** narrow only at the benchmark ATLAS publishes (Z'_SSM at 3% against
  `r`=0.015 on `m(ee)`, the DM-mediator coupling grid, single-VLQ, U1 at the flavour-anomaly point),
  **3** already broader than `r` there (KK gluon 15-30% vs `r`=0.08, coloron/axigluon,
  composite/NJL), and **3** with no Breit-Wigner peak at all (QBH thresholds,
  ADD/HEIDI continua, pair-produced Type-III/VLL, toponium at 2 m_t). A signal wider than `r`
  correlates neighbouring mass points, so counting resolution elements **over**-counts independent
  looks: the bias is conservative, `Z_local` too strict rather than too loose.
  Only 1 axis is motivated
  *exclusively* by non-peaking models
  (`m(multi)`); dropping it
  takes N from 4,319 to
  4,298 and `Z_local` from 6.46 to 6.46, so the whole question is
  worth 0.001 sigma. `NONPEAK_ON` records the (model, axis) pairs where an
  otherwise-narrow class does not peak, including the H/A interference with SM ttbar.

Source: `searchbudget/stages/search_budget.py` -> `results/tables/search_budget.csv`,
`results/tables/search_budget_selections.csv`. Figures: `searchbudget/stages/budget_plots.py`,
`searchbudget/stages/budget_waterfall.py`. Windows and their sources:
`searchbudget/core/bump_observables.py` (SCAN).
Excess bookkeeping: `results/overviews/EXCESS_COUNTING.md`.
Uncertainty budget: `results/overviews/BUDGET_UNCERTAINTY.md`.
