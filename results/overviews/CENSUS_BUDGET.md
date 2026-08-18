# The publication record, priced in trials

`SEARCH_BUDGET.md` counts the spectra public BSM models motivate. This prices the other side: the
searches ATLAS has actually **published**, taken from the census (`data/published_spectra.csv`,
86 entries over 290 papers) and run through
the same rule, `n_s = (1/r) ln(M_hi/M_lo)`, `Z_local = sqrt(25 + 2 ln N)`.

## Method
Every census entry carries the canonical budget axis it scans (`budget_axis`; several when the
entry scans several; `-` when it falls on none of the 46) and the range it scanned (`scan_GeV`,
transcribed from the published range the census records). The resolution `r` is that axis' -- the
budget's one physics input, unchanged here; an entry on no axis is priced at the default
`r = 0.06`. **Where the census does not record a range, the axis' own published window is the
fallback**, which is why the two bases below differ.

| basis | what one look is |
|---|---|
| published searches | one census entry on one axis, over the range that entry scanned |
| axes scanned | one axis, over the **union** of every published range on it, counted once |

A fixed-mass search (LFV `Z` and `tau` decays, exclusive `H`/`Z -> quarkonium + gamma`,
`H -> Z gamma`) scans nothing and contributes exactly one look.

## Summary
| basis | units | N_trials | band (r x0.5..x2) | Z_local for 5s global | band |
|---|--:|--:|---|--:|---|
| **published searches** (86 census entries) | 100 | **7,710** | 3,855-15,420 | **6.55** | 6.44-6.65 |
| axes scanned (union of the ranges) | 50 | **3,672** | 1,836-7,344 | **6.44** | 6.33-6.54 |
| model-motivated axes (reference, `SEARCH_BUDGET.md`) | 46 | **3,685** | 1,842-7,369 | **6.44** | 6.33-6.54 |

**Reading it.** Pricing every published ATLAS resonance search over the range it actually scanned
gives **N = 7,710**, so a 5 sigma global discovery in the published program needs a local
**6.55 sigma**. Counting each axis once instead, over the union of everything published on
it, gives **N = 3,672** and **6.44 sigma**. The two bases differ by a factor
2.1 in N and 0.11 sigma in the bar -- the same lesson as the
model side: the answer is nearly independent of how finely the program is sliced.

## Coverage
The census reaches 40 of the 46 canonical axes. 6 carry no published
search at all: `m(cb) dijet`, `m(eb)`, `m(mub)`, `m(multi)`, `m(tt)/m(jj)`, `multilepton`.

10 published entries fall on no canonical axis and are priced on their own published range at
`r = 0.06`; 12 more fall on no axis **and** carry no chargeable range, so
they are listed below as unpriced and are missing from `N` -- `N` is a lower bound by that much. Of
those 12, 9 published no scanned range and 3 declare
no single axis to scan (the anomaly-detection and generic multi-body entries, whose observable is
"many").

The 62 distinct bump observables the census scans are therefore
40 of the 46 budget axes, 10 scanned axes outside the 46, and
12 carrying nothing chargeable.

**Why the union basis lands on the model space.** Counting each axis once gives 3,672
against the model space's 3,685, a gap of 0.3% -- closer than
it has any right to be, because two omissions cancel. The union misses the 6 model
axes with no published search, worth 267 looks on the
model side, and it adds 10 scanned axes outside the 46, worth 251. Charging
neither leaves 3,422 against 3,418
on the same 40 axes. The robust statement is the ten-percent one, and it is not a
tautology: the model space takes its windows from published searches but its axes from the models,
so the agreement says the program scans almost exactly the axes the models motivate.

## Per axis (union of every published range)
| axis | scanned range(s) [GeV] | r | n_s | model-side n_s |
|---|---|--:|--:|--:|
| `m(gammagamma)` | 66-110+150-5000 | 0.01 | 402 | 402 |
| `m(mumu) (Zd)` | 0.3-400 | 0.02 | 360 | 360 |
| `m(ee)` | 150-8000 | 0.015 | 265 | 265 |
| `m(ee) (Zd)` | 1-400 | 0.025 | 240 | 240 |
| `m(egamma)` | 200-5000 | 0.015 | 215 | 215 |
| `m(emu) LFV` | 100-8000 | 0.03 | 146 | 146 |
| `m(ee) SS` | 200-1300 | 0.015 | 125 | 125 |
| `m(Vgamma)` | 220-6800 | 0.03 | 114 | 110 |
| `m(ej)` | 200-5000 | 0.03 | 107 | 107 |
| `m(mugamma)` | 200-5000 | 0.03 | 107 | 107 |
| `m(mumu)` | 35-75+150-8000 | 0.05 | 95 | 80 |
| `m(jj)` | 100-8200 | 0.05 | 88 | 74 |
| `m(muj)` | 200-5000 | 0.04 | 80 | 80 |
| `m(eZ)` | 100-1100 | 0.03 | 80 | 80 |
| `m(emu) SS` | 200-1300 | 0.025 | 75 | 75 |
| `m(jgamma)` | 500-7000 | 0.04 | 66 | 66 |
| `m(muZ)` | 100-1100 | 0.04 | 60 | 60 |
| `m(eejj)` | 400-7000 | 0.05 | 57 | 57 |
| `m(bj)` | 200-6000 | 0.06 | 57 | 57 |
| `m(VV)` | 200-6000 | 0.06 | 57 | 57 |
| `m(Vh)` | 200-5000 | 0.06 | 54 | 54 |
| `m(bb)` | 12-160+450-6000 | 0.1 | 52 | 40 |
| `m(mumu) SS` | 200-1300 | 0.04 | 47 | 47 |
| `m(3j)` | 200-1800 | 0.05 | 44 | 44 |
| `m(tb)` | 180-6000 | 0.08 | 44 | 44 |
| `m(mumujj)` | 400-7000 | 0.07 | 41 | 41 |
| `m(HH)` | 250-6000 | 0.08 | 40 | 40 |
| `mT(ev)` | 150-7000 | 0.1 | 38 | 38 |
| `m(tt)` | 350-6000 | 0.08 | 36 | 36 |
| `m(tauj)` | 200-5000 | 0.1 | 32 | 32 |
| `m(etau) LFV` | 100-3000 | 0.12 | 28 | 28 |
| `m(mutau) LFV` | 100-3000 | 0.12 | 28 | 28 |
| `mT(taunu)` | 200-5000 | 0.12 | 27 | 27 |
| `mT(muv)` | 150-7000 | 0.15 | 26 | 26 |
| `m(tW)` | 500-3000 | 0.08 | 22 | 22 |
| `m(ttZ)/m(Zt)` | 1000-4000 | 0.08 | 17 | 17 |
| `m(Wb)` | 800-3000 | 0.08 | 17 | 16 |
| `m(tautau)` | 15-85 | 0.12 | 14 | 43 |
| `m(Ht)` | 1000-3000 | 0.08 | 14 | 14 |
| `m(taub)` | 1500-3000 | 0.12 | 6 | 19 |

## Per published search
| family | published search | axis | window [GeV] | from | r | n_s |
|---|---|---|---|---|--:|--:|
| Jets & hadronic | Dijet (inclusive, high mass) | `m(jj)` | 1800-8200 | census | 0.05 | 30 |
| Jets & hadronic | Photon + jet | `m(jgamma)` | 500-7000 | axis | 0.04 | 66 |
| Jets & hadronic | Dijet (low mass, trigger-level) | `m(jj)` | 100-250 | census | 0.05 | 18 |
| Jets & hadronic | Dijet + ISR photon/jet | `m(jj)` | 200-650 | census | 0.05 | 24 |
| Jets & hadronic | Dijet, b-tagged | `m(bb)` | 15-62+450-6000 | axis | 0.1 | 40 |
| Jets & hadronic | Dijet, b-tagged | `m(bj)` | 200-6000 | axis | 0.06 | 57 |
| Jets & hadronic | Pair-produced dijet resonances (4 jets) | `m(jj)` | 200-8000 | axis | 0.05 | 74 |
| Jets & hadronic | Three-jet / multijet (3-quark res.) | `m(3j)` | 200-1800 | axis | 0.05 | 44 |
| Jets & hadronic | Z + light hadronic resonance (H decay) | - | 0.5-3.5 | census | 0.06 | 32 |
| Jets & hadronic | Dijet + isolated lepton | `m(jj)` | 250-6000 | census | 0.05 | 64 |
| Jets & hadronic | Dijet from dark quarks | `m(jj)` | 200-8000 | axis | 0.05 | 74 |
| Jets & hadronic | Boosted light hadronic res. + photon | - | 20-100 | census | 0.06 | 27 |
| Photons | Diphoton (high mass) | `m(gammagamma)` | 160-2800 | census | 0.01 | 286 |
| Photons | Diphoton (low mass / boosted) | `m(gammagamma)` | 66-110+150-5000 | axis | 0.01 | 402 |
| Photons | Diphoton, extra dimensions | `m(gammagamma)` | 2520-3920 | census | 0.01 | 44 |
| Photons | Photon-jets (collimated gamma pairs) | - | 2-10 | census | 0.06 | 27 |
| Photons | Diphoton + forward protons (ALP) | `m(gammagamma)` | 150-1600 | census | 0.01 | 237 |
| Photons | Three photons | - | - | **unpriced** | 0.06 | 0 |
| Photons | H -> Za -> ll gamgam | - | 0.1-33 | census | 0.06 | 97 |
| Leptons | Ditau | `m(tautau)` | 20-85 | census | 0.12 | 12 |
| Leptons | Dilepton (ee/mumu, high mass) | `m(ee)` | 150-8000 | axis | 0.015 | 265 |
| Leptons | Dilepton (ee/mumu, high mass) | `m(mumu)` | 150-8000 | axis | 0.05 | 80 |
| Leptons | Lepton-flavour-violating dilepton (e-mu/e-tau/mu-tau) | `m(emu) LFV` | 100-8000 | axis | 0.03 | 146 |
| Leptons | Lepton-flavour-violating dilepton (e-mu/e-tau/mu-tau) | `m(etau) LFV` | 100-3000 | axis | 0.12 | 28 |
| Leptons | Lepton-flavour-violating dilepton (e-mu/e-tau/mu-tau) | `m(mutau) LFV` | 100-3000 | axis | 0.12 | 28 |
| Leptons | Four leptons / 4mu | `m(ee) (Zd)` | 5-81 | census | 0.025 | 111 |
| Leptons | Four leptons / 4mu | `m(mumu) (Zd)` | 5-81 | census | 0.02 | 139 |
| Leptons | Lepton + MET (W-prime) | `mT(ev)` | 150-7000 | census | 0.1 | 38 |
| Leptons | Lepton + MET (W-prime) | `mT(muv)` | 150-7000 | census | 0.15 | 26 |
| Leptons | LFV Z decays | - | single mass | fixed | 0.06 | 1 |
| Leptons | Dimuon (low / intermediate mass) | `m(mumu)` | 35-75 | census | 0.05 | 15 |
| Leptons | SFOS dilepton + MET (edge) | `m(ee)` | 150-8000 | axis | 0.015 | 265 |
| Leptons | SFOS dilepton + MET (edge) | `m(mumu)` | 150-8000 | axis | 0.05 | 80 |
| Leptons | Same-sign dimuon (strong gravity) | `m(mumu) SS` | 200-1300 | axis | 0.04 | 47 |
| Leptons | Collimated / displaced lepton pairs | `m(ee) (Zd)` | 1-400 | axis | 0.025 | 240 |
| Leptons | Collimated / displaced lepton pairs | `m(mumu) (Zd)` | 0.3-400 | axis | 0.02 | 360 |
| Leptons | Tau + MET | `mT(taunu)` | 200-5000 | axis | 0.12 | 27 |
| Leptons | Trilepton resonance | `m(eZ)` | 100-1100 | census | 0.03 | 80 |
| Leptons | Trilepton resonance | `m(muZ)` | 100-1100 | census | 0.04 | 60 |
| Leptons | Periodic signals in dielectron/diphoton | `m(ee)` | 150-8000 | axis | 0.015 | 265 |
| Leptons | Periodic signals in dielectron/diphoton | `m(gammagamma)` | 66-110+150-5000 | axis | 0.01 | 402 |
| Leptons | LFV tau -> 3 mu | - | single mass | fixed | 0.06 | 1 |
| Diboson | ZZ/ZW -> llqq, vvqq, 4l | `m(VV)` | 200-6000 | axis | 0.06 | 57 |
| Diboson | WZ/WW fully leptonic | `m(VV)` | 200-6000 | axis | 0.06 | 57 |
| Diboson | WW/WZ semileptonic (lvqq) | `m(VV)` | 300-5000 | census | 0.06 | 47 |
| Diboson | Heavy Higgs -> ZZ / WW | `m(VV)` | 200-6000 | axis | 0.06 | 57 |
| Diboson | Diboson fully hadronic (boson-tagged jets) | `m(VV)` | 1300-6000 | census | 0.06 | 25 |
| Diboson | Diboson combination | `m(VV)` | 500-3000 | census | 0.06 | 30 |
| Diboson | WW -> e-mu | `m(VV)` | 200-6000 | axis | 0.06 | 57 |
| Diboson | Triboson WWW | - | - | **unpriced** | 0.06 | 0 |
| Vector boson + X | VH (llbb, lvbb, vvbb) | `m(Vh)` | 220-5000 | census | 0.06 | 52 |
| Vector boson + X | Z-gamma | `m(Vgamma)` | 220-3400 | census | 0.03 | 91 |
| Vector boson + X | V-gamma, hadronic V | `m(Vgamma)` | 1000-6800 | census | 0.03 | 64 |
| Vector boson + X | H -> Z gamma | - | single mass | fixed | 0.06 | 1 |
| Vector boson + X | VH fully hadronic | `m(Vh)` | 1500-5000 | census | 0.06 | 20 |
| Vector boson + X | Higgs + photon | - | 700-4000 | census | 0.06 | 29 |
| Vector boson + X | High-pT Z + X resonances | - | - | **unpriced** | 0.06 | 0 |
| Top & heavy quarks | ttbar resonances | `m(tt)` | 400-5000 | census | 0.08 | 32 |
| Top & heavy quarks | Single VLQ: Wb / Ht / Zt / Hb | `m(Wb)` | 800-3000 | axis | 0.08 | 17 |
| Top & heavy quarks | Single VLQ: Wb / Ht / Zt / Hb | `m(Ht)` | 1000-3000 | axis | 0.08 | 14 |
| Top & heavy quarks | Single VLQ: Wb / Ht / Zt / Hb | `m(ttZ)/m(Zt)` | 1000-4000 | axis | 0.08 | 17 |
| Top & heavy quarks | tb resonances (W-prime, charged Higgs) | `m(tb)` | 500-6000 | census | 0.08 | 31 |
| Top & heavy quarks | Fourth-generation / heavy quark -> Wq, Zq | - | - | **unpriced** | 0.06 | 0 |
| Top & heavy quarks | Wt (excited/single VLQ) | `m(tW)` | 500-3000 | axis | 0.08 | 22 |
| Top & heavy quarks | Resonant top + jet | - | - | **unpriced** | 0.06 | 0 |
| Top & heavy quarks | Top-philic resonances | `m(tt)` | 350-6000 | axis | 0.08 | 36 |
| Top & heavy quarks | Dark mesons -> tb | `m(tb)` | 180-6000 | axis | 0.08 | 44 |
| Higgs pairs & extended scalars | H -> aa (exotic Higgs decays) | `m(bb)` | 15-60 | census | 0.1 | 14 |
| Higgs pairs & extended scalars | H -> aa (exotic Higgs decays) | `m(tautau)` | 15-60 | census | 0.12 | 12 |
| Higgs pairs & extended scalars | HH -> bbbb | `m(HH)` | 1000-5000 | census | 0.08 | 20 |
| Higgs pairs & extended scalars | X -> S H (scalar + Higgs) | `m(HH)` | 300-6000 | census | 0.08 | 37 |
| Higgs pairs & extended scalars | A -> ZH | `m(Vh)` | 200-5000 | axis | 0.06 | 54 |
| Higgs pairs & extended scalars | HH -> bb gamgam | `m(HH)` | 250-6000 | axis | 0.08 | 40 |
| Higgs pairs & extended scalars | HH -> bbWW / WWWW / gamgamWW | `m(HH)` | 250-6000 | axis | 0.08 | 40 |
| Higgs pairs & extended scalars | HH -> bb tautau | `m(HH)` | 251-1600 | census | 0.08 | 23 |
| Higgs pairs & extended scalars | X -> H Y (anomaly detection / XH) | `m(HH)` | 1500-6000 | census | 0.08 | 17 |
| Higgs pairs & extended scalars | Charged Higgs -> WZ / WH | `m(VV)` | 250-3000 | census | 0.06 | 41 |
| Higgs pairs & extended scalars | Charged Higgs -> WZ / WH | `m(Vh)` | 250-3000 | census | 0.06 | 41 |
| Higgs pairs & extended scalars | bbA -> bbbb (b-associated) | `m(bb)` | 12-100 | census | 0.1 | 21 |
| Higgs pairs & extended scalars | HH + vector boson | `m(HH)` | 260-1000 | census | 0.08 | 17 |
| Higgs pairs & extended scalars | Triple Higgs (6b) | - | - | **unpriced** | 0.06 | 0 |
| Higgs pairs & extended scalars | t -> qX, X -> bb | `m(bb)` | 20-160 | census | 0.1 | 21 |
| Higgs pairs & extended scalars | HH + top quarks (ttHH) | `m(HH)` | 250-6000 | axis | 0.08 | 40 |
| Leptoquarks, excited & heavy fermions | Leptoquark pair -> lepton + jet | `m(ej)` | 200-5000 | axis | 0.03 | 107 |
| Leptoquarks, excited & heavy fermions | Leptoquark pair -> lepton + jet | `m(muj)` | 200-5000 | axis | 0.04 | 80 |
| Leptoquarks, excited & heavy fermions | Leptoquark pair -> lepton + jet | `m(tauj)` | 200-5000 | axis | 0.1 | 32 |
| Leptoquarks, excited & heavy fermions | Third-generation LQ -> b tau / t tau | `m(taub)` | 1500-3000 | census | 0.12 | 6 |
| Leptoquarks, excited & heavy fermions | Excited leptons (e*, mu*, tau*) | `m(egamma)` | 200-5000 | axis | 0.015 | 215 |
| Leptoquarks, excited & heavy fermions | Excited leptons (e*, mu*, tau*) | `m(mugamma)` | 200-5000 | axis | 0.03 | 107 |
| Leptoquarks, excited & heavy fermions | Excited leptons (e*, mu*, tau*) | `m(ej)` | 200-5000 | axis | 0.03 | 107 |
| Leptoquarks, excited & heavy fermions | Excited leptons (e*, mu*, tau*) | `m(muj)` | 200-5000 | axis | 0.04 | 80 |
| Leptoquarks, excited & heavy fermions | Heavy neutrino / right-handed W (lljj) | `m(eejj)` | 400-7000 | axis | 0.05 | 57 |
| Leptoquarks, excited & heavy fermions | Heavy neutrino / right-handed W (lljj) | `m(mumujj)` | 400-7000 | axis | 0.07 | 41 |
| Leptoquarks, excited & heavy fermions | Quantum black holes -> lepton + jet | `m(ej)` | 200-5000 | axis | 0.03 | 107 |
| Leptoquarks, excited & heavy fermions | Quantum black holes -> lepton + jet | `m(muj)` | 200-5000 | axis | 0.04 | 80 |
| Leptoquarks, excited & heavy fermions | Resonant single leptoquark | `m(ej)` | 200-5000 | axis | 0.03 | 107 |
| Leptoquarks, excited & heavy fermions | Resonant single leptoquark | `m(muj)` | 200-5000 | axis | 0.04 | 80 |
| Leptoquarks, excited & heavy fermions | Heavy lepton -> Z + lepton | `m(eZ)` | 114-176 | census | 0.03 | 14 |
| Leptoquarks, excited & heavy fermions | Heavy lepton -> Z + lepton | `m(muZ)` | 114-176 | census | 0.04 | 11 |
| Generic & other | Quarkonium + photon (exclusive H/Z decays) | - | single mass | fixed | 0.06 | 1 |
| Generic & other | Doubly charged Higgs (same-sign ll / WW) | `m(ee) SS` | 200-1300 | axis | 0.015 | 125 |
| Generic & other | Doubly charged Higgs (same-sign ll / WW) | `m(emu) SS` | 200-1300 | axis | 0.025 | 75 |
| Generic & other | Doubly charged Higgs (same-sign ll / WW) | `m(mumu) SS` | 200-1300 | axis | 0.04 | 47 |
| Generic & other | Heavy neutral lepton (displaced/prompt mass) | - | 8-65 | census | 0.06 | 35 |
| Generic & other | Anomaly detection (multilepton) | - | - | **unpriced** | 0.06 | 0 |
| Generic & other | Generic multi-body invariant masses | - | - | **unpriced** | 0.06 | 0 |
| Generic & other | Two-body masses, unsupervised anomaly detection | - | - | **unpriced** | 0.06 | 0 |
| Generic & other | Dark photon in rare Z decays | `m(ee) (Zd)` | 1-400 | axis | 0.025 | 240 |
| Generic & other | Dark photon in rare Z decays | `m(mumu) (Zd)` | 0.3-400 | axis | 0.02 | 360 |
| Generic & other | Displaced diphoton / dielectron | - | - | **unpriced** | 0.06 | 0 |
| Generic & other | B_s invariant mass structure (X(5568)) | - | - | **unpriced** | 0.06 | 0 |
| Generic & other | Xb -> Upsilon pi pi | - | - | **unpriced** | 0.06 | 0 |
| | **total (86 entries, 100 priced looks)** | | | | | **7,710** |

## Assumptions & caveats
- **The axis assignment is a curated judgement**, recorded per row in `data/published_spectra.csv`
  (`budget_axis`) so it can be checked and changed. Entries whose final state spans several axes
  (LFV dilepton, leptoquark pair, single VLQ, doubly charged Higgs) contribute one look per axis.
- **The window fallback over-counts.** 44 entries carry no published range in
  the census and inherit their axis' full window, so on the published-searches basis several
  analyses on one axis are each charged the whole axis. The union basis is free of this and is the
  conservative reading; the truth sits between them.
- **Off-axis entries are priced at the default resolution** `r = 0.06`, not at a measured one.
- **Scope**: the same as the budget's -- invariant/transverse-mass bump hunts. Census entries that
  are themselves multi-spectrum scans (the anomaly-detection and generic multi-body papers) and the
  displaced programs fall on no axis and stay unpriced: their trials belong in the combinatorial
  count of `scaled_scan.py`, not here.
- Everything the model-side budget assumes about resolution elements, correlations and the
  narrow-resonance approximation applies unchanged; see `SEARCH_BUDGET.md`.

Source: `scripts/census_budget.py` from `data/published_spectra.csv` ->
`results/tables/census_budget.csv`. The census itself: `results/overviews/PUBLISHED_CENSUS.md`.
The model-side budget: `results/overviews/SEARCH_BUDGET.md`.
