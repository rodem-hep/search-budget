# Models the catalogue would catch that no ATLAS search has scanned

Sections A-C are the first sweep, for models on *unscanned* axes. Section D is the second sweep, for
models on axes the map already carries: same sources, but asking whether a spectrum already in the
map is motivated by more classes than the map lists. Section D buys no trials by construction --
adding a class to a counted spectrum leaves `N` untouched -- so it changes only the model column of
the map. Section E is the third sweep, over the two-body pairs of `two_body_matrix` that carried no
axis at all; unlike D it does move `N`, and it is what raised the model-motivated-but-unscanned
count from 14 to 19.

A literature sweep (arXiv, INSPIRE, FeynRules model database; 2026-08-19) for published models
predicting a narrow resonance on an axis the combinatorial catalogue can form but no published
ATLAS bump hunt has scanned. Coverage is judged against `data/published_spectra.csv` (the census)
and means an ATLAS mass scan of the *resonant axis* — a counting analysis on the same final state
does not count, and CMS coverage does not count (noted in brackets where it exists, since a CMS
scan is at least evidence the axis is analysable). Hand-written: nothing regenerates this file.
Every arXiv id below was verified against its abstract page during the sweep.

Candidates from the model families already in `public_obs_map.py` are included only when the
resonant axis is new; they are marked "(new signature of a mapped family)".

All of the candidates below have since been folded into the model catalogue
(`public_obs_map.PUBLIC_OBS`/`WIDTH`, with the new axes in `bump_observables`); this file remains
as the per-model evidence record. The defining references the paper cites for each
literature-sourced class are declared in `public_obs_map.LITERATURE`, their metadata fetched by
`fetch-model-meta`, and the per-class provenance written to `results/tables/model_classes.csv`.

## A. New axes — no ATLAS scan at any energy

**The radiative excited-fermion trio: b* → bγ, t* → tγ, τ* → τγ.** ATLAS has scanned the
radiative decay of excited light quarks (m(jγ)) and of excited e/μ (m(ℓγ)), but never the
b-tagged, top, or tau siblings, though all three follow from the same Baur–Spira–Zerwas
Lagrangian (PRD 42 (1990) 815; NLO ℓ*→ℓV in 1210.8307). For f = −f′ the ℓ*→ℓγ coupling vanishes
and the covered ℓγ axes say nothing. Compositions bγ, tγ, τγ; narrow (Γ/m a few % at f = 1).
[CMS has scanned all three: bγ 2305.07998, tγ 2602.20477 (spin-3/2 benchmark 1208.5811),
τγ 2410.21137 — ATLAS has none.]

**Excited neutrino ν* → ℓjj** (hep-ph/0401066; 2606.24486 notes the absence of any dedicated
search). Mandatory partner of the excited charged leptons ATLAS does scan; dominant visible decay
ℓW(→jj). Composition ejj/mjj (±MET category). No scan at any LHC experiment.

**Excited muon μ* → μjj** (hep-ph/0212006 framework). ATLAS scanned the electron channel only
(1906.03204, 36 fb⁻¹); the muon channel has never been scanned at ATLAS. Composition mjj inside
mmjj. [CMS: 2001.04521.]

**Composite Majorana neutrino N_ℓ → ℓjj via contact interactions** (1510.07988). Resonant
production at compositeness-scale cross sections, same-sign Majorana signature — a different
family from the mapped Keung–Senjanovic W_R, whose census row charges only the 4-body m(ℓℓjj).
The 3-body m(ℓjj) sub-axis is unscanned at ATLAS. [CMS: 1706.08578, 2210.03082.]

**Leptogluons ℓ₈ → ℓg** (1211.6394; resonant lepton-gluon production 2212.06178). Coloured
lepton partners, QCD-sized cross sections. The 2-body m(ℓj) axis is charged in the census
(resonant single LQ, 2507.03650), so the pair-production double-(ℓj) peak and the eej/mmj
single-production category carry the new content; borderline between class A and C.

**Single vector-like B → bZ(→ℓℓ) at 13 TeV** (framework 1305.4172, 1306.0572; (B,Y) doublet has
BR(B→bZ) ≈ 50%). ATLAS scanned m(bZ) only at 8 TeV (1409.5500); no 13 TeV bZ(ℓℓ) resonance scan
exists at either experiment beyond CMS's 2.3 fb⁻¹ (1701.07409) — the axis is ~60× behind in
luminosity. Compositions bZ, jbZ, bbZ. (New signature of the mapped VLQ family.)

**Same-sign m(tt): colour-sextet diquarks** (vector sextet 1009.5379, scalar "diquark Higgs"
0709.1486, FCNC Z′ 0907.4112, 1011.4960). uu-valence s-channel production, LHC-unique. ATLAS and
CMS constrain same-sign tops only by counting; no mass scan exists anywhere. The catalogue's
charge-blind tt category covers it automatically.

**m(tt̄) in three-top and four-top categories** (triple-top: 1001.0221, 1901.04643, 1909.03998;
ttH-associated heavy scalar: 2507.05334, 1910.09581). Associated production evades the
gg→H→tt̄ interference problem of the inclusive m(tt) scan. ATLAS has no 3-top search of any kind
and only BDT counting for ttH/A→tttt (2211.01136). (New category of the mapped 2HDM family.)
[CMS 4-top resonance: 2604.14058, 2608.10148.]

**Paired-(bb) resonances in 4b: Manohar–Wise colour-octet scalars and symmetric TRSM
h₃ → h₂h₂ → (bb)(bb)** (hep-ph/0606172 + 0710.3133; TRSM 1908.08554). ATLAS's paired-dijet scan
is untagged and stale, and its X→SH→bbbb scans (2405.20926, 2607.18484) fix one leg to H(125) —
the equal-mass paired-bb axis with both masses free is unscanned. Composition bbbb, mass group
bb. (TRSM entry is a new axis of a mapped family.)

**b-philic Z′ → bb̄ in a 3b category** (1707.07016 and the hadrophobic corners of the
flavour-anomaly Z′ family). The census's b-associated m(bb) scan (bbA→4b) stops at 100 GeV;
0.2–2 TeV in bbb/jbb is open at ATLAS.

**Leptophobic Z′ → dark states → 4ℓ** (1111.0633). No ℓℓ decay by construction, so dilepton
scans are blind; the m(4ℓ) axis with free intermediate masses is unscanned at ATLAS (heavy
ZZ→4ℓ fixes both pairs to m_Z). Composition eemm/eeee/mmmm. [CMS: 8 TeV 1701.01345; 2604.14236
covers only sub-15 GeV intermediates.]

**Dark-Z sub-peak: heavy scalar → Z Z_d → (ℓℓ)(ℓℓ)** (1203.2947, 1412.0018). The m(ℓℓ) = m(Z_d)
sub-axis is scanned only inside H(125)→ZX (1802.03388, 15–55 GeV); any heavier parent falls in
the hole between that scan and heavy ZZ→4ℓ. Compositions eeZ/mmZ, mass group ℓℓ — the catalogue
forms exactly this.

**Warped KK cascades: γ+(jj), (γγ)+j, Z_lep+(jj), WWW** (1612.00047, 1711.09920). Double peaks
m(total) = KK resonance, m(jj/γγ/VV) = radion. ATLAS has scanned none of these cascade axes
(the census WWW row is stale and charges no axis). [CMS: trijet 2201.02140, W(WW) 2112.13090 /
2201.08476, g+(WW) 2410.17303 — the photon and leptonic-Z cascades are unscanned by anyone.]

**Low-scale technicolor ω_T/a_T → γ π_T(→bb̄/jj)** (0706.2339, hep-ph/9903369). Sub-GeV widths,
200–700 GeV; ATLAS's boosted-X+γ scans charge only the m(jj) sub-axis, never the total m(γjj).
Never carried from the Tevatron to the LHC.

**MFV RPV gluino → tbs** (1111.1239). Fully hadronic 3-body m(tbj) peak, pair-produced with 50%
same-sign tops; only counting constraints exist. Reconstructed width O(10%) — marginal but
huntable.

**Resonant smuon (RPV λ′₂₁₁): μ̃ → μχ̃⁰, χ̃⁰ → μjj** (hep-ph/0001224, 1201.5014). Double peak;
the census charges the 4-body m(μμjj) (W_R row) but not the 3-body m(μjj) sub-axis.
[CMS: 1811.09760, 2016 data only.]

## B. Existing axis, new flavour/association category

These resonate on a census-charged axis but in a category no ATLAS scan has looked at — the
catalogue's category split covers them with no new axis needed.

- **Third-family-hypercharge / "Plan B" Z′ → μμ with b-jets** (1809.01158, 1905.10327,
  2111.06691, 2607.18411): the b-associated category tags the flavour-non-universal coupling that
  inclusive dimuon scans dilute. [CMS: 2307.08708, m ≥ 350 GeV only.]
- **L_μ−L_τ "dressed in color"** (1403.1269): same category, μμ+b below 350 GeV — open at both
  experiments.
- **RPV ν̃_τ → eμ from bb̄ fusion** (hep-ph/0311327): eμ axis scanned inclusively; the
  b-associated category isolates λ′₃₃₃.
- **LFV leptoquark, single production → eμj** (2301.04119, 2103.12724): ℓj peaks in the mixed
  eμ+jet category; no LFV-LQ scan at any collider since HERA.
- **Flavoured muoquark, single production μ + LQ(→μb)** (2103.12504, 2103.12724): the m(μb) axis
  is now census-current (2507.03650); what remains is the μμb category and the anomaly-era
  motivation.

## C. Stale or ageing census axes these models re-motivate

- **m(tj)** ("Resonant top + jet", stale): flavored W′/Z′ with t–q couplings (1102.0018,
  0907.4112) and t* → tg (1208.5811) — the AFB-legacy mediators survive precisely where only
  counting analyses constrain them.
- **m(ℓZ)** ("Trilepton resonance", ageing; VLL row stale at 8 TeV): vector-like lepton doublets
  (1510.03456), type-III triplets as a mass scan rather than counting (2111.07949, 2006.04123),
  ℓ* → ℓZ (1210.8307) — ATLAS's 8 TeV VLL scan stops at 176 GeV.
- **ℓ±ℓ± below ~300 GeV** (doubly-charged row is ageing): Zee-Babu k±± (1402.4491, 2206.14833).
- **X→aa→4γ with a heavy parent**: ATLAS scanned it only at 8 TeV (1509.05051);
  hep-ph/0005308-type hidden sectors.

## Compositions chased without finding a resonant model

Pure mixed-flavour trileptons eem/emm and ℓ+3b (ebbb/mbbb): the natural candidates (bileptons,
doubly charged scalars) peak in the ℓ±ℓ± pair, an axis ATLAS already scans. m(tγ) had no
narrow-resonance model paper beyond the spin-3/2 t* benchmark CMS adopted. Resonant monotop was
rejected (peak involves the invisible leg — no invariant-mass bump).

## D. More classes on axes the map already lists

A second sweep (arXiv, INSPIRE, FeynRules model database, and the benchmark lists of the published
ATLAS and CMS searches on each axis; 2026-08-19), asking the complementary question: for each of the
56 spectra already in the map, is a public model class predicting a peak there missing from its row?
It adds two classes and eighteen class-spectrum pairs, `N` unchanged. Two are new classes with
defining references, declared in `public_obs_map.LITERATURE`; the rest attach a class the map already
carries to a further axis, and need no reference beyond the one that class already stands on.

**New classes.**

- **Bilepton, 331 model** (1806.04536, 1812.02723). The doubly charged gauge boson `Y^{++}` of
  SU(3)_c x SU(3)_L x U(1)_X decays to a like-sign lepton pair, so it peaks in exactly the
  `m(ll) SS` axes ATLAS already scans for `H^{++}`. A vector, not a scalar: distinct from the
  Georgi-Machacek, type-II seesaw and Zee-Babu doubly charged scalars the map already carries, and
  distinct in width (gauge coupling, a few %). The first sweep noted bileptons and dropped them for
  pointing at a scanned axis; that is precisely what this sweep keeps. Axes: `m(ee) SS`,
  `m(emu) SS`, `m(mumu) SS`.
- **String / Regge resonance** (0808.0497 dijet, 0804.2013 direct photon + jet). The first Regge
  excitations of the quark and the gluon in TeV-scale open-string constructions, with amplitudes
  independent of the compactification. A standing benchmark of the CMS dijet scan (excluded below
  7.4 TeV) with no UFO implementation in the database. Axes: `m(jj)`, `m(jgamma)`.

**Further axes for classes already in the map.** Each is a decay channel of a mapped class that a
published search interprets on an axis the map did not connect it to.

- `2HDM (general/typeII/CPV)` -> `mT(taunu)` (`H+ -> tau nu`, the type-II interpretation of the
  ATLAS charged-Higgs transverse-mass scan) and `m(HH)` (the heavy CP-even `H -> hh`, the spin-0
  benchmark of the di-Higgs searches).
- `MSSM/NMSSM/RPV SUSY` -> `m(HH)` (hMSSM `H -> hh`, excluded by the ATLAS `bbtautau` scan).
- `KK graviton (Gstar)` -> `m(HH)` (the bulk RS spin-2 benchmark of the di-Higgs searches) and
  `m(jj)` (RS gravitons, a CMS dijet benchmark).
- `Singlet scalar / SM+Scalars` -> `m(gammagamma)` (the generic narrow spin-0 `X -> gammagamma`,
  the archetypal diphoton benchmark; the class already carried `HH`, `VV` and `Vgamma`).
- `Minimal Z' / U(1)` -> `m(tautau)` (`Z'_SSM` and `Z'_NU`, scanned alongside the MSSM Higgs bosons
  in the ATLAS ditau search).
- `Vector-like quark (VLQ)` -> `m(tW)` (the charge-5/3 partner `X_{5/3} -> tW`, mass-reconstructed
  in the same-sign-lepton searches).
- `Color-octet scalar (MW)` -> `m(jj)` (the Manohar-Wise octet decaying to `gg`, a CMS dijet
  benchmark; the class already carried `m(bb)`).
- `RPV resonant slepton` -> `m(emu) LFV`, `m(etau) LFV`, `m(mutau) LFV` (the resonant tau sneutrino
  with R-parity-violating couplings, one of the three benchmarks of the ATLAS LFV dilepton search
  alongside the Z' and the quantum black hole already mapped there).
- `Heavy neutrino / HNL (prompt)` and `LRSM / Alt-LRSM` -> `m(ejj)`, `m(mujj)` (the heavy neutrino
  itself is a peak in `m(l jj)` inside the `lljj` final state the map charges on the 4-body axis --
  the same double-peak bookkeeping the map already applies to the RPV resonant slepton).

**Checked and not added.** Technicolor on `m(jj)` (the technirho dijet branching fraction is
suppressed relative to the `VV` and `gammajj` modes the map already carries); a heavy spin-2 on
`m(tt)` (no published ttbar-resonance search benchmarks it); 2HDM on `m(gammagamma)` (the branching
fraction is negligible away from the alignment limit).

## E. The two-body pairs that carried no axis

The two-body grid (`two_body_matrix.py`) classifies every pair of grid objects as *scanned* (a
published ATLAS search covers it), *axes unscanned* (the catalogue holds an axis there but nothing
published scans it), or *gap* (no axis at all). The eleven gap pairs had never been swept: sweep A
only looked at compositions that already carried a model. Sweeping them (arXiv, INSPIRE, FeynRules
model database, and the benchmark lists of the published searches on the neighbouring axes;
2026-08-19) turns up a public model for every one of them, and in each case a class the catalogue
already carries -- these are structural omissions, the flavour or decay-mode sibling of an axis the
map already holds.

Seven of the eleven have no published ATLAS scan of any kind *and* can be fitted in the
combinatorial scan, and became axes: the catalogue goes from 56 to 63 spectra, `N` from 4,118 to
4,319, `Z_local` from 6.45 to 6.46, and the model-motivated-but-unscanned count
(`unscanned_spectra.py`) from 14 to 19 fittable, 21 counting the two the scan cannot fit.
Resolutions follow the worse-leg rule of the resolution appendix; windows are stated per axis in
`bump_observables.SCAN` with their source. An eighth, `m(tauH)`, is motivated but **not** fittable
and was left out -- see below.

Counted instead in the scan's own units, where a spectrum is one axis in one category and counts as
motivated only when the category is a final state some model predicts with the mass on its resonant
sub-system (`scaled_scan.FINAL_STATES`), **164 of the 4,438 fittable spectra (3.7%)** are motivated,
over 59 axes and 76 distinct (final state, mass) pairs; **51 of them (1.15% of the scan)**, over 18
axes, sit on an axis no published ATLAS search scans (`unscanned_scan_units.csv`). Matching on the
mass composition alone -- asking only whether some model likes `m(jj)`, not whether it predicts the
category the `m(jj)` sits in -- would put the motivated share at 70% instead, and is wrong: it counts
`m(jj)` in a two-jet-plus-two-electron category as motivated.

| new axis | class(es) that predict a peak | why the axis, and its window |
|---|---|---|
| `m(jV)` | Excited quark `q*/b*`; VLQ | `q* -> qW/qZ` is in the same Baur-Spira-Zerwas Lagrangian as the mapped `q*->q gamma` and `q*->qg`; BR(Z) is 3-5% at `f_s = 1` and above 20% at `f_s = 0`. Window: the published `q*->q gamma` grid, 500 GeV-7 TeV. |
| `m(jH)` | VLQ | single light-flavour vector-like quark `Q -> qH`, the `H` sibling of the `qW`/`qZ` modes. Window 0.4-3 TeV, from the light-flavour VLQ searches. |
| `m(eH)`, `m(muH)` | VLL; Type-III seesaw | `L -> l h` and `Sigma -> l h` are the mandatory partners of the mapped `l Z` decays; 1112.3080 reconstructs exactly this mass. Window 0.2-1.5 TeV. |
| `m(tauV)` | Excited lepton `l*`; VLL; Type-III seesaw | the tau sibling of the mapped `m(eZ)`/`m(muZ)`; ATLAS's trilepton-resonance scan is e/mu only (1506.01291). Window 100-1100 GeV, the e/mu grid. Mapped to the **boosted** `V` composition, not the leptonic `Z` the e/mu axes use: with a hadronic tau the leptonic-Z composition carries a weight of 3.5e-8 and the scan cannot fit it, while `TV` fits at 30 looks. |
| `m(et)`, `m(mut)` | Scalar LQ `S1`; Vector LQ `U1` | leptoquarks with a third-generation quark leg decay to `t l`; the `R_K` triplet's `t,mu` component is the standard flavour-anomaly benchmark. ATLAS's LQ3 search covers `b tau` and `t tau` only. Window 0.2-2 TeV, the lepton + heavy-quark LQ grid. |

**The four left out, and why.** `taut` (`m(t tau)`), `bH` (`m(Hb)`) and `gaH` (`m(H gamma)`) are
motivated -- LQ3 `-> t tau`, single VLQ `B -> bH`, and the `H gamma` resonance -- but ATLAS has
*published a scan of each*: they appear in `published_spectra.csv` as the second observable of the
third-generation-LQ row, the fourth observable of the single-VLQ row, and the `Higgs + photon` row
(2008.05928, 0.7-4 TeV). They stay outside the axis set as three of the ten published observables the
budget does not carry as axes.

`tauH` is the fourth, and it is left out for the opposite reason: nothing scans it, VLL `L -> tau h`
and Type-III `Sigma -> tau h` motivate it, but it fails the statistics requirement. Its kinematic
floor sits at 165 GeV (the Higgs mass) and its one-event mass at 2.0 TeV, which leaves 24.6
resolution elements at `r = 0.102` against the 25 the yield model asks for -- a knife-edge miss that
the fittability systematic (15-50 elements) spans. The grid's own `gap_price` prices it at 29.6
looks only because it starts every gap pair at the generic 100 GeV rather than at the pair's
kinematic floor.

All four remain the `gap` cells of the grid and the whole content of the "axis set" systematic.
Adding them would take `N` to 4,470 and `Z_local` to 6.47.
