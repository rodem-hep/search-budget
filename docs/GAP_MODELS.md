# Models the catalogue would catch that no ATLAS search has scanned

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
