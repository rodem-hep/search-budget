# The publication-side census

Every ATLAS resonance search, grouped by the mass spectrum it scans. Assembled from the
collaboration's publication record on INSPIRE-HEP and curated by hand; the per-spectrum rows, with
their arXiv references, are in `data/published_spectra.csv`. Scope matches the budget's: bump hunts
for new states, so hadron-spectroscopy measurements are out even where they are bump hunts.

This is the complement to `SEARCH_BUDGET.md`. That counts the spectra public BSM models motivate
(63 canonical mass axes); this counts the searches that have actually been published
(86 entries over 290 papers). **The two use different bases and must not be
summed**: the publication record separates entries that share a mass axis when they are different
analyses (62 distinct bump observables appear across the 86 entries), while the
budget merges them onto one axis and counts resolution elements along it.

| | |
|---|---|
| catalogued spectra | **86** |
| search papers | **290** (2010-2026) |
| current (latest paper 2024 or later) | 29 |
| ageing (2019-2023) | 38 |
| **stale (nothing since before 2019)** | **19** |
| with a published Run-3 (13.6 TeV) result | **6** |

## By final state

| family | spectra | papers | current | ageing | stale | Run-3 |
|---|--:|--:|--:|--:|--:|--:|
| Diboson | 8 | 24 | 0 | 4 | 4 | 0 |
| Generic & other | 10 | 24 | 3 | 5 | 2 | 0 |
| Higgs pairs & extended scalars | 14 | 46 | 7 | 6 | 1 | 2 |
| Jets & hadronic | 11 | 31 | 4 | 4 | 3 | 0 |
| Leptons | 14 | 59 | 7 | 5 | 2 | 1 |
| Leptoquarks, excited & heavy fermions | 7 | 31 | 3 | 3 | 1 | 2 |
| Photons | 7 | 15 | 1 | 3 | 3 | 0 |
| Top & heavy quarks | 8 | 44 | 3 | 2 | 3 | 0 |
| Vector boson + X | 7 | 16 | 1 | 6 | 0 | 1 |
| **total** | **86** | **290** | **29** | **38** | **19** | **6** |

## Not revisited since before 2019

These carry Run-1 or early-Run-2 sensitivity. What the budget says about them splits them in two.

**14 of the 19 sit on a mass axis that is already counted in `N`, so revisiting one costs
nothing in trials** -- the discovery bar for the re-run is the bar the program already pays, and the
whole cost is analysis effort. The remaining 5 fall on no axis in the budget's 63, so re-running
one extends the axis count rather than reusing it, and it is priced like any other new spectrum.

| last published | spectrum | observable | counted axis |
|--:|---|---|---|
| 2012 | Diphoton, extra dimensions | `m_gamgam` | `m(gammagamma)` |
| 2012 | Resonant top + jet | `m_tq` | `m(tj)` |
| 2013 | Same-sign dimuon (strong gravity) | `m_mumu (SS)` | `m(mumu) SS` |
| 2014 | Xb -> Upsilon pi pi | `m_Upsilon pipi` | **adds an axis** |
| 2015 | Three photons | `m_3gamma` | **adds an axis** |
| 2015 | Heavy Higgs -> ZZ / WW | `m_ZZ, m_WW` | `m(VV)` |
| 2015 | Fourth-generation / heavy quark -> Wq, Zq | `m_Wq, m_Zq` | `m(bZ)` |
| 2015 | Wt (excited/single VLQ) | `m_Wt` | `m(tW)` |
| 2015 | Heavy lepton -> Z + lepton | `m_Zl` | `m(eZ); m(muZ)` |
| 2016 | Diboson combination | `m_VV` | `m(VV)` |
| 2016 | Triboson WWW | `m_WWW` | **adds an axis** |
| 2017 | Photon + jet | `m_(gamma j)` | `m(jgamma)` |
| 2017 | Pair-produced dijet resonances (4 jets) | `m_jj (paired)` | `m(jj)` |
| 2017 | WW -> e-mu | `m_WW` | `m(VV)` |
| 2018 | Three-jet / multijet (3-quark res.) | `m_jjj` | `m(3j)` |
| 2018 | Photon-jets (collimated gamma pairs) | `m_gamgam-jets` | **adds an axis** |
| 2018 | SFOS dilepton + MET (edge) | `m_ll` | `m(ee); m(mumu)` |
| 2018 | HH -> bbWW / WWWW / gamgamWW | `m_HH` | `m(HH)` |
| 2018 | B_s invariant mass structure (X(5568)) | `m_Bs pi` | **adds an axis** |

## Already re-run at 13.6 TeV

| spectrum | observable |
|---|---|
| Dilepton (ee/mumu, high mass) | `m_ll` |
| H -> Z gamma | `m_llgam` |
| HH + top quarks (ttHH) | `m_HH` |
| Quantum black holes -> lepton + jet | `m_lj` |
| Resonant single leptoquark | `m_lj` |
| X -> S H (scalar + Higgs) | `m_SH` |

Source: `searchbudget/stages/published_census.py` from `data/published_spectra.csv`. What this
record costs in trials, entry by entry: `results/overviews/CENSUS_BUDGET.md`
(`searchbudget/stages/census_budget.py`).
