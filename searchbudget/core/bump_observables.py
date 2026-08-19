BUMP_DCLASS = {"D2v", "D2i", "DNb"}
SQRTS = 13600.0
FLOOR_DEFAULT = 150.0

OBS_MERGE = {"m(HH) 4b": "m(HH)", "m(SH)": "m(HH)"}
def canon(obs): return OBS_MERGE.get(obs, obs)

FLAV_SPLIT = {
    "m(ll)":         ["m(ee)", "m(mumu)"],
    "m(ll) (Zd)":    ["m(ee) (Zd)", "m(mumu) (Zd)"],
    "m(ll') LFV":    ["m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
    "m(l+-l+-) SS":  ["m(ee) SS", "m(emu) SS", "m(mumu) SS"],
    "m(lj)":         ["m(ej)", "m(muj)", "m(tauj)"],
    "mT(lv)":        ["mT(ev)", "mT(muv)"],
    "m(lljj)":       ["m(eejj)", "m(mumujj)"],
    "m(lgamma)":     ["m(egamma)", "m(mugamma)"],
    "m(lb)":         ["m(eb)", "m(mub)", "m(taub)"],
    "m(lZ)":         ["m(eZ)", "m(muZ)"],
}
def flav_channels(obs):
    return FLAV_SPLIT.get(obs, [obs])

def is_inclusive(obs): return obs in FLAV_SPLIT

FLOOR = {
    "m(ee)":         (110, "Z-peak"),
    "m(mumu)":       (110, "Z-peak"),
    "m(ee) (Zd)":    (1.0, "e-LJ (e+e- merges into one cluster below ~1 GeV)"),
    "m(mumu) (Zd)":  (0.3, "prompt dark-photon muon LJ / 4mu"),
    "m(tautau)":     ( 35, "bb-tagged a->tautau (off Z)"),
    "m(emu) LFV":    ( 12, "trig (no SM bkg)"),
    "m(etau) LFV":   ( 35, "hadronic-tau reco"),
    "m(mutau) LFV":  ( 35, "hadronic-tau reco"),
    "m(ee) SS":      ( 15, "di-electron (same-sign)"),
    "m(emu) SS":     ( 15, "e-mu (same-sign)"),
    "m(mumu) SS":    ( 15, "di-muon (same-sign)"),
    "multilepton":   (  8, "di-muon"),
    "m(jj)":         (200, "TLA+ISR"),
    "m(cb) dijet":   ( 15, "tt-tagged (Cc)"),
    "m(3j)":         (200, "multijet / TLA"),
    "m(bb)":         ( 15, "h->aa->4b (h-tagged)"),
    "m(gammagamma)": ( 40, "low-mass diphoton"),
    "m(jgamma)":     (150, "photon+jet"),
    "m(tt)":         (350, "2*m_t kin"),
    "m(tt)/m(jj)":   (200, "pair, jj/tt"),
    "m(tb)":         (175, "m_t+m_b kin"),
    "m(VV)":         (160, "2*m_V kin"),
    "m(Vh)":         (200, "A->Zh (off-shell V)"),
    "m(HH)":         (250, "2*m_h kin"),
    "m(HH) 4b":      (250, "2*m_h kin"),
    "m(SH)":         (170, "m_h + m_S kin (S down to ~20 GeV)"),
    "m(ttZ)/m(Zt)":  (450, "m_Vlt+m_t kin"),
    "m(ej)":         ( 12, "electron trig"),
    "m(muj)":        ( 12, "muon trig"),
    "m(tauj)":       ( 50, "hadronic-tau + jet trig"),
    "m(eb)":         ( 15, "electron + b-tagged jet"),
    "m(mub)":        ( 15, "muon + b-tagged jet"),
    "m(taub)":       ( 50, "hadronic-tau + b-tagged jet"),
    "m(bj)":         (200, "b-tagged + light jet (dijet trig)"),
    "m(eZ)":         (100, "Z-peak + electron (trilepton)"),
    "m(muZ)":        (100, "Z-peak + muon (trilepton)"),
    "mT(ev)":        ( 35, "electron/MET"),
    "mT(muv)":       ( 35, "muon/MET"),
    "mT(taunu)":     ( 50, "tau+MET"),
    "m(multi)":      (300, "multijet threshold"),
    "m(eejj)":       (300, "2e+2j kin (W_R)"),
    "m(mumujj)":     (300, "2mu+2j kin (W_R)"),
    "m(Vgamma)":     (200, "m_V + photon pT"),
    "m(egamma)":     (100, "electron+photon trig"),
    "m(mugamma)":    (100, "muon+photon trig"),
    "m(tW)":         (260, "m_t+m_W kin"),
    "m(Wb)":         (100, "m_W+b kin"),
    "m(Ht)":         (300, "m_h+m_t kin"),
    "m(bgamma)":     (150, "photon + b-tagged jet"),
    "m(tgamma)":     (250, "top tag + photon pT"),
    "m(taugamma)":   (100, "hadronic-tau + photon trig"),
    "m(ejj)":        (250, "electron + two jets kin"),
    "m(mujj)":       (250, "muon + two jets kin"),
    "m(taujj)":      (300, "hadronic-tau + two jets"),
    "m(bZ)":         (100, "Z-peak + b-jet"),
    "m(tj)":         (200, "m_t + jet kin"),
    "m(gammajj)":    (200, "photon + dijet (ISR-style trigger)"),
    "m(tbj)":        (250, "m_t + b + jet kin"),
    "m(jV)":         (200, "jet + boosted V (large-R V mass)"),
    "m(jH)":         (250, "jet + boosted h (large-R h mass)"),
    "m(eH)":         (200, "electron + h(bb) kin"),
    "m(muH)":        (200, "muon + h(bb) kin"),
    "m(tauV)":       (100, "hadronic tau + W/Z"),
    "m(et)":         (200, "electron + m_t kin"),
    "m(mut)":        (200, "muon + m_t kin"),
}

def floor(obs, default=FLOOR_DEFAULT):
    return float(FLOOR[obs][0]) if obs in FLOOR else float(default)

RESOLUTION = {
    "m(ee)": 0.015, "m(mumu)": 0.05,
    "m(ee) (Zd)": 0.025, "m(mumu) (Zd)": 0.02,
    "m(emu) LFV": 0.03, "m(etau) LFV": 0.12, "m(mutau) LFV": 0.12,
    "m(ee) SS": 0.015, "m(emu) SS": 0.025, "m(mumu) SS": 0.04,
    "m(ej)": 0.03, "m(muj)": 0.04, "m(tauj)": 0.10,
    "m(eb)": 0.04, "m(mub)": 0.05, "m(taub)": 0.12, "m(bj)": 0.06,
    "m(eZ)": 0.03, "m(muZ)": 0.04,
    "mT(ev)": 0.10, "mT(muv)": 0.15,
    "m(eejj)": 0.05, "m(mumujj)": 0.07,
    "m(egamma)": 0.015, "m(mugamma)": 0.03,
    "m(gammagamma)": 0.01, "m(jgamma)": 0.04,
    "m(jj)": 0.05, "m(3j)": 0.05, "m(cb) dijet": 0.05,
    "m(VV)": 0.06, "m(Vh)": 0.06, "m(multi)": 0.06,
    "m(tt)": 0.08, "m(tb)": 0.08, "m(tt)/m(jj)": 0.08,
    "m(HH)": 0.08, "m(HH) 4b": 0.08, "m(ttZ)/m(Zt)": 0.08,
    "m(bb)": 0.10,
    "m(tautau)": 0.12, "mT(taunu)": 0.12, "multilepton": 0.12,
    "m(Vgamma)": 0.03,
    "m(tW)": 0.08, "m(Wb)": 0.08, "m(Ht)": 0.08,
    "m(bgamma)": 0.05, "m(tgamma)": 0.06, "m(taugamma)": 0.08,
    "m(ejj)": 0.05, "m(mujj)": 0.07, "m(taujj)": 0.10,
    "m(bZ)": 0.05, "m(tj)": 0.08, "m(gammajj)": 0.05, "m(tbj)": 0.10,
    "m(jV)": 0.06, "m(jH)": 0.08,
    "m(eH)": 0.08, "m(muH)": 0.08, "m(tauV)": 0.10,
    "m(et)": 0.08, "m(mut)": 0.08,
}

RESOLUTION_SOURCE = {
    "e/gamma": ("EM calorimeter energy resolution", "arXiv:2309.05471"),
    "muon": ("track sagitta, rises with pT", "arXiv:2012.00578"),
    "light jet": ("jet energy scale and resolution", "arXiv:2007.02645"),
    "b-jet": ("as a jet, plus semileptonic neutrinos", "arXiv:1907.05120"),
    "tau_had": ("visible decay products only", "arXiv:1512.05955"),
    "large-R V/h/t/H": ("large-radius jet mass", "arXiv:2311.08885"),
    "MET": ("transverse masses, coarsest axes", "arXiv:2402.05858"),
}
RES_DEFAULT = 0.06

def res(obs): return RESOLUTION.get(obs, RES_DEFAULT)

SCAN = {
    "m(gammagamma)": ([(66, 110), (150, 5000)],
                      "low-mass diphoton 66-110 (JHEP 01 (2025) 053); spin-0/2 high-mass 150 GeV-5 TeV (incl. X->aa->4gamma)"),
    "m(ee)":         (150, 8000, "high-mass Drell-Yan dielectron scan (grid prepared to 8 TeV)"),
    "m(mumu)":       (150, 8000, "high-mass Drell-Yan dimuon scan (grid prepared to 8 TeV)"),
    "m(ee) (Zd)":    (1.0, 400,  "e-LJ 1-10 + H->ZdZd->4e 1-60 + heavy ZdZd grid to 400"),
    "m(mumu) (Zd)":  (0.3, 400,  "prompt dark-photon muon LJ 0.3-10 + H->ZdZd->4mu 1-60 + grid to 400"),
    "m(emu) LFV":    (100, 8000, "Z'->emu LFV dilepton (grid 100 GeV-8 TeV)"),
    "m(etau) LFV":   (100, 3000, "Z'->etau LFV (hadronic-tau channel, grid 100 GeV-3 TeV)"),
    "m(mutau) LFV":  (100, 3000, "Z'->mutau LFV (hadronic-tau channel, grid 100 GeV-3 TeV)"),
    "m(ee) SS":      (200, 1300, "H++ same-sign ee pair scan"),
    "m(emu) SS":     (200, 1300, "H++ same-sign e-mu pair scan"),
    "m(mumu) SS":    (200, 1300, "H++ same-sign mumu pair scan"),
    "multilepton":   (50, 10000, "type-III seesaw / VLL / prompt heavy-N mass hypotheses (grid to 10 TeV)"),
    "m(ej)":         (200, 5000, "leptoquark pair -> eq (2006.05872) + QBH jet-electron; grid to 5 TeV"),
    "m(muj)":        (200, 5000, "leptoquark pair -> muq (2006.05872) + QBH jet-muon; grid to 5 TeV"),
    "m(tauj)":       (200, 5000, "leptoquark pair -> tau q mass scan (LQ3; grid to 5 TeV)"),
    "m(eb)":         (200, 2000, "RPV stop -> b e (2406.18367) + b/c-tagged LQ pair (2006.05872)"),
    "m(mub)":        (200, 2000, "RPV stop -> b mu (2406.18367) + b/c-tagged LQ pair (2006.05872)"),
    "m(taub)":       (200, 2000, "3rd-generation LQ -> tau b"),
    "m(bj)":         (200, 6000, "b-tagged + light jet: b* -> b g, b-tagged dijet resonance "
                                 "(two-body ML anomaly scan 2307.01612)"),
    "m(eZ)":         (100, 1100, "RPV chargino/neutralino trilepton resonance e+Z (2011.10543)"),
    "m(muZ)":        (100, 1100, "RPV chargino/neutralino trilepton resonance mu+Z (2011.10543)"),
    "m(cb) dijet":   (20, 1000,  "t->H+(cb)b 60-160 + tb-associated heavy H+->cb grid to 1 TeV"),
    "m(jj)":         (200, 8000, "dijet: TLA+ISR 200-450, TLA 450-1800, high-mass to 8 TeV"),
    "m(3j)":         (200, 1800, "RPV gluino -> 3-jet paired resonance scan"),
    "m(bb)":         ([(15, 62), (450, 6000)],
                      "h->aa->4b low-mass 15-62; b-tagged dijet resonance 450 GeV-6 TeV"),
    "m(jgamma)":     (500, 7000, "q* -> q gamma (grid 500 GeV-7 TeV)"),
    "m(tt)":         (350, 6000, "ttbar resonance (resolved+boosted), from the 2*m_t threshold"),
    "m(tt)/m(jj)":   (400, 2000, "pair-produced coloron/sgluon, per-particle tt/jj legs"),
    "m(tb)":         (180, 6000, "W'->tb / heavy H+->tb (grid from 180 GeV)"),
    "m(VV)":         (200, 6000, "WW/WZ/ZZ diboson family (ggF+VBF; grid 200 GeV-6 TeV)"),
    "m(Vh)":         (200, 5000, "W/Z + h(bb) resonance"),
    "m(HH)":         (250, 6000, "X->HH bbbb/bbtautau/bbyy incl. boosted 4b"),
    "m(ttZ)/m(Zt)":  (1000, 4000,"VLQ T->Zt single/pair + W'->Vlt cascade"),
    "m(multi)":      (3000, 10500,"QBH / sphaleron multi-object threshold scan"),
    "mT(ev)":        (150, 7000, "W' -> e nu transverse mass"),
    "mT(muv)":       (150, 7000, "W' -> mu nu transverse mass"),
    "mT(taunu)":     (200, 5000, "W' -> tau nu"),
    "m(tautau)":     (35, 6000,  "a/H/A/Z' -> tautau (low-mass a/H grid from 35 GeV)"),
    "m(eejj)":       (400, 7000, "W_R -> eejj (Keung-Senjanovic), published scan 0.4-7 TeV"),
    "m(mumujj)":     (400, 7000, "W_R -> mumujj (Keung-Senjanovic), published scan 0.4-7 TeV"),
    "m(Vgamma)":     (220, 6000, "X -> Z gamma / W gamma (ll/qq + photon)"),
    "m(egamma)":     (200, 5000, "excited electron e* -> e gamma (7/8 TeV only, 1201.3293/1308.1364; no 13 TeV search)"),
    "m(mugamma)":    (200, 5000, "excited muon mu* -> mu gamma (7/8 TeV only, 1201.3293/1308.1364; no 13 TeV search)"),
    "m(tW)":         (500, 3000, "excited b* -> tW"),
    "m(Wb)":         (800, 3000, "single VLQ T/Y -> Wb"),
    "m(Ht)":         (1000, 3000,"single VLQ T -> Ht"),
    "m(bgamma)":     (200, 5000, "b* -> b gamma: no ATLAS search at any energy; the q*->q gamma "
                                 "grid with a b-tag (CMS 2305.07998 scans to ~6 TeV)"),
    "m(tgamma)":     (300, 3000, "excited top t* -> t gamma: no ATLAS search; spin-3/2 t* "
                                 "pair-production grid (CMS 2602.20477)"),
    "m(taugamma)":   (200, 5000, "tau* -> tau gamma: no ATLAS search at any energy; the e*/mu* "
                                 "grid (CMS 2410.21137 scans 0.175-5 TeV)"),
    "m(ejj)":        (250, 6000, "e*/nu* -> e qq: the e* system of the eejj search (1906.03204); "
                                 "composite-N grid to 6 TeV"),
    "m(mujj)":       (250, 6000, "mu*/nu* -> mu qq: no ATLAS search at any energy "
                                 "(CMS 2001.04521 scans the mumujj system); grid to 6 TeV"),
    "m(taujj)":      (300, 5000, "tau* -> tau qq via contact interaction (2303.09444)"),
    "m(bZ)":         (200, 2000, "VLQ b'/B -> bZ(ll): 7 TeV pair scan (1204.1265), single+pair "
                                 "at 8 TeV (1409.5500); no 13 TeV scan, grid to 2 TeV"),
    "m(tj)":         (350, 3000, "W'/Z' -> t + jet (1209.6593, 7 TeV only); flavoured-mediator "
                                 "grids (1102.0018)"),
    "m(gammajj)":    (200, 4000, "warped KK cascade gamma + dijet (1612.00047) and LSTC "
                                 "omega_T/a_T -> gamma pi_T (0706.2339): no search at any collider"),
    "m(tbj)":        (400, 3000, "MFV RPV gluino -> tbs (1111.1239): counting limits only, "
                                 "no mass scan anywhere"),
    "m(jV)":         (500, 7000, "q* -> qW/qZ, the gauge decays of the same excited-quark "
                                 "Lagrangian as the published q*->q gamma grid; no ATLAS search "
                                 "(CMS light-flavour VLQ partners 1708.02510 cover 0.4-2 TeV)"),
    "m(jH)":         (400, 3000, "single light-flavour vector-like quark Q -> qH; no search at "
                                 "any collider (CMS 1708.02510 scans the qW/qZ modes 0.4-2 TeV)"),
    "m(eH)":         (200, 1500, "vector-like lepton L -> e h and type-III seesaw Sigma -> e h; "
                                 "no ATLAS search, the lZ partner mode scanned over 100-1100 "
                                 "(1506.01291)"),
    "m(muH)":        (200, 1500, "vector-like lepton L -> mu h and type-III seesaw Sigma -> mu h; "
                                 "no ATLAS search, the lZ partner mode scanned over 100-1100 "
                                 "(1506.01291)"),
    "m(tauV)":       (100, 1100, "tau* -> tau Z/W, vector-like lepton and type-III seesaw in "
                                 "the tau channel: the e/mu trilepton grid (1506.01291, "
                                 "2011.10543); no ATLAS scan of the tau channel"),
    "m(et)":         (200, 2000, "leptoquark -> t e (S1/S3/U1 with a third-generation quark leg), "
                                 "on the lepton + heavy-quark LQ grid (2006.05872); no ATLAS "
                                 "scan, the LQ3 search covering b tau / t tau only"),
    "m(mut)":        (200, 2000, "leptoquark -> t mu, the R_K-motivated triplet component, on the "
                                 "lepton + heavy-quark LQ grid (2006.05872); no ATLAS scan "
                                 "(CMS covers mu b only)"),
}

def scan_segments(obs):
    v = SCAN[canon(obs)]
    if isinstance(v[0], (list, tuple)):
        return [(float(lo), float(hi)) for lo, hi in v[0]]
    return [(float(v[0]), float(v[1]))]

def scan(obs):
    seg = scan_segments(obs)
    return seg[0][0], seg[-1][1]

def scan_source(obs): return SCAN[canon(obs)][-1]

import math as _math
def n_s(m_lo, m_hi, r):
    if m_hi <= m_lo: return 0.0
    return _math.log(m_hi / m_lo) / r

def ns_scan(obs):
    return sum(n_s(lo, hi, res(obs)) for lo, hi in scan_segments(obs))

def ns_achievable(obs):
    return n_s(floor(obs), SQRTS, res(obs))

def z_local_for_global5(N):
    return _math.sqrt(25.0 + 2.0 * _math.log(N)) if N > 1 else 5.0

CANON_ORDER = sorted(FLOOR, key=lambda o: -ns_achievable(o))

for _parent, _leaves in FLAV_SPLIT.items():
    assert _parent not in SCAN, f"inclusive {_parent!r} must not keep a SCAN entry (double count)"
    assert _parent not in FLOOR, f"inclusive {_parent!r} must not keep a FLOOR entry"
    for _lf in _leaves:
        assert _lf in SCAN and _lf in FLOOR, f"leaf {_lf!r} missing a SCAN/FLOOR entry"
assert set(SCAN) == set(FLOOR) - set(OBS_MERGE), \
    f"SCAN/FLOOR key mismatch: {set(SCAN) ^ (set(FLOOR) - set(OBS_MERGE))}"
del _parent, _leaves, _lf

PCAP = {"P1": SQRTS, "P2": SQRTS/2, "P3": SQRTS-400, "P4": SQRTS/2, "P5": SQRTS, "P6": SQRTS}

def prod_ceiling(pmodes):
    return max((PCAP.get(p, SQRTS) for p in pmodes), default=SQRTS)

def fmt_mass(m):
    m = float(m)
    if m < 1000:
        return (f"{m:.1f} GeV" if m < 10 else f"{m:.0f} GeV")
    t = m / 1000.0
    return (f"{t:.0f} TeV" if abs(t - round(t)) < 0.05 else f"{t:.1f} TeV")

def fmt_range(lo, hi):
    return f"{fmt_mass(lo)}–{fmt_mass(hi)}"
