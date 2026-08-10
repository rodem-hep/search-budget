#!/usr/bin/env python3
"""Single source of truth for the BUMP observables. Import this; never redefine it locally.

  BUMP_DCLASS  decay-class codes that mark a genuine mass bump (D2v/D2i/DNb)
  SQRTS        kinematic ceiling [GeV]
  FLOOR        observable -> (analyzable-floor GeV, reason)
  SCAN         observable -> (lo, hi, source): the published-search scan window
  OBS_MERGE    raw decay-observable label -> canonical budget observable (same mass axis)
  FLAV_SPLIT   flavour-inclusive label -> the LEAF flavour channels it spans (m(ll) -> m(ee),
               m(mumu)); FLOOR/RESOLUTION/SCAN are keyed on leaves. Aggregate as
               flav_channels(canon(Dobs)).

Public information only: published scan windows, resolutions and kinematic floors. No side
effects on import beyond the self-consistency asserts at the bottom. Why the published window
and not a sample-grid range: docs/METHOD_NOTES.md.
"""
BUMP_DCLASS = {"D2v", "D2i", "DNb"}
SQRTS = 13600.0
FLOOR_DEFAULT = 150.0

# same-axis merges for budget counting (raw Dobs -> canonical observable)
# m(SH): X->S H with m_S free is the same two-scalar mass axis as X->HH; the extra
# m_S dimension is an event selection, counted at the selections level not here.
OBS_MERGE = {"m(HH) 4b": "m(HH)", "m(SH)": "m(HH)"}
def canon(obs): return OBS_MERGE.get(obs, obs)

# ---------------------------------------------------------------------------- lepton flavour
# LEAF CHANNELS. A generic-lepton axis is not one search: ee and mumu are scanned by separate
# analyses with different triggers and -- decisively -- different mass resolutions (EM calorimeter
# ~1.5% and roughly mass-independent; muon sagitta resolution degrades from ~2% at 200 GeV to
# ~15% at 3 TeV). Merging them forces one compromise `r` onto both and counts one look where the
# program makes two, so flavour is a first-class axis here: FLOOR / RESOLUTION / SCAN and
# CANON_ORDER are defined on LEAF channels only.
#
# FLAV_SPLIT maps a flavour-INCLUSIVE label (used when the flavour of a decay is genuinely
# unspecified, e.g. Z' -> ll) onto the leaves it spans. This is the mirror of OBS_MERGE: canon()
# merges same-axis decay splits DOWNWARD, flav_channels() expands a flavour-blind label UPWARD
# onto every channel it populates. Always compose them as
#     flav_channels(canon(Dobs))
# An inclusive label marks all its leaves as populated (a Z' -> ll prediction does motivate both
# the ee and the mumu scan); it never creates a channel of its own, so nothing is double-counted.
FLAV_SPLIT = {
    "m(ll)":         ["m(ee)", "m(mumu)"],
    "m(ll) (Zd)":    ["m(ee) (Zd)", "m(mumu) (Zd)"],
    "m(ll') LFV":    ["m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
    "m(l+-l+-) SS":  ["m(ee) SS", "m(emu) SS", "m(mumu) SS"],
    "m(lj)":         ["m(ej)", "m(muj)", "m(tauj)"],
    "mT(lv)":        ["mT(ev)", "mT(muv)"],
    "m(lljj)":       ["m(eejj)", "m(mumujj)"],
    "m(lgamma)":     ["m(egamma)", "m(mugamma)"],
    # a b-jet is a DIFFERENT OBJECT from a light jet (it is b-tagged, so a search selects on it
    # and scans its own spectrum), hence lepton+b is its own mass axis, not a b-tagged selection
    # of lepton+light-jet. Same rule that keeps m(bb) separate from m(jj).
    "m(lb)":         ["m(eb)", "m(mub)", "m(taub)"],
    "m(lZ)":         ["m(eZ)", "m(muZ)"],
    # NOT split: `multilepton` is an inclusive 3L/4L counting axis whose signal regions are
    # binned by flavour AND charge and then statistically combined into one spectrum -- the
    # flavour multiplicity there is an event selection (NSEL), not a separate mass axis.
}
def flav_channels(obs):
    """Leaf (flavour-resolved) channels an observable occupies. Identity for labels that are
    already leaves or are deliberately flavour-inclusive (multilepton)."""
    return FLAV_SPLIT.get(obs, [obs])

def is_inclusive(obs): return obs in FLAV_SPLIT

FLOOR = {
    "m(ee)":         (110, "Z-peak"),
    "m(mumu)":       (110, "Z-peak"),
    # sub-GeV dark photons are reachable as muon lepton-jets; an e+e- pair below ~1 GeV is a
    # single merged EM cluster, so the electron channel floor sits an order of magnitude higher.
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
    # ---- channels whose scan is published (or was published pre-13 TeV) ----
    "m(eejj)":       (300, "2e+2j kin (W_R)"),
    "m(mumujj)":     (300, "2mu+2j kin (W_R)"),
    "m(Vgamma)":     (200, "m_V + photon pT"),
    "m(egamma)":     (100, "electron+photon trig"),
    "m(mugamma)":    (100, "muon+photon trig"),
    "m(tW)":         (260, "m_t+m_W kin"),
    "m(Wb)":         (100, "m_W+b kin"),
    "m(Ht)":         (300, "m_h+m_t kin"),
}

def floor(obs, default=FLOOR_DEFAULT):
    """Analyzable lower scan edge [GeV] for a bump observable."""
    return float(FLOOR[obs][0]) if obs in FLOOR else float(default)

# Fractional mass resolution r = sigma_M / M per bump observable (rough central values; the one
# real physics input to the trials count). Leptons/photons sharp; jets/tau/MET broad. Default 0.06.
RESOLUTION = {
    # ---- lepton-flavour leaves. The e/mu difference is the whole reason these axes are split:
    # electrons are measured in the EM calorimeter (sigma/m ~1.5%, ~mass-independent), muons by
    # the sagitta of a curved track (sigma/m grows ~linearly with pT). The mumu values below are
    # the EFFECTIVE FLAT r that reproduces the trials integral of a rising resolution over that
    # channel's own scan window -- for m(mumu) over 150 GeV-8 TeV, r(M) ~ 2% at 200 GeV rising to
    # ~15% at 3 TeV integrates to n_s = 76, i.e. r_eff = 0.05. CAVEAT: a single flat r per muon
    # channel is a window-averaged stand-in; the framework assumes constant fractional resolution
    # (see n_s), so a mass-dependent r(M) would be the next refinement, not a free parameter here.
    "m(ee)": 0.015, "m(mumu)": 0.05,
    # at LOW mass the ordering REVERSES: a few-GeV muon pair is sharp (tracker, no sagitta
    # penalty at low pT) while the e+e- pair is a marginally-resolved cluster -- so here the MUON
    # channel is the better-resolved one, the opposite of the high-mass pattern above.
    "m(ee) (Zd)": 0.025, "m(mumu) (Zd)": 0.02,
    "m(emu) LFV": 0.03, "m(etau) LFV": 0.12, "m(mutau) LFV": 0.12,   # tau reco dominates
    "m(ee) SS": 0.015, "m(emu) SS": 0.025, "m(mumu) SS": 0.04,
    "m(ej)": 0.03, "m(muj)": 0.04, "m(tauj)": 0.10,
    # b-jets carry semileptonic neutrinos -> a b-tagged leg resolves slightly worse than a light
    # jet of the same energy, so each lepton+b axis sits one notch above its lepton+j twin
    "m(eb)": 0.04, "m(mub)": 0.05, "m(taub)": 0.12, "m(bj)": 0.06,
    # a leptonically-decaying Z is sharp; the trilepton mass inherits the third lepton's error
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
    # new channels
    "m(Vgamma)": 0.03,
    "m(tW)": 0.08, "m(Wb)": 0.08, "m(Ht)": 0.08,
}
RES_DEFAULT = 0.06

def res(obs): return RESOLUTION.get(obs, RES_DEFAULT)

# ---------------------------------------------------------------------------- scan windows
# Published-search windows (ATLAS Run-2/early-Run-3 papers; extended where the search family's
# published signal grid prepares a wider scan). One entry per CANONICAL observable (post-merge). The source
# string documents the anchor. These set n_s for the budget: n_s = (1/r) ln(hi/lo).
# An entry is (lo, hi, source) or ([(lo1,hi1),(lo2,hi2),...], source) when the published program
# scans DISJOINT segments of the axis (the gap between them is not scanned and not counted).
SCAN = {
    "m(gammagamma)": ([(60, 110), (150, 5000)],
                      "low-mass diphoton 66-110; spin-0/2 high-mass 150 GeV-5 TeV (incl. X->aa->4gamma)"),
    # Lepton-flavour leaves inherit their parent's PUBLISHED window: ATLAS scans the same mass
    # range in ee and mumu (the dilepton papers quote one range and combine at the limit stage),
    # so the split changes r and the channel count, NOT the window. The two exceptions below are
    # documented where they occur (Zd low edge; LFV tau channels).
    "m(ee)":         (150, 8000, "high-mass Drell-Yan dielectron scan (grid prepared to 8 TeV)"),
    "m(mumu)":       (150, 8000, "high-mass Drell-Yan dimuon scan (grid prepared to 8 TeV)"),
    # the electron channel cannot open below ~1 GeV (merged cluster), so its window starts there
    "m(ee) (Zd)":    (1.0, 400,  "e-LJ 1-10 + H->ZdZd->4e 1-60 + heavy ZdZd grid to 400"),
    "m(mumu) (Zd)":  (0.3, 400,  "prompt dark-photon muon LJ 0.3-10 + H->ZdZd->4mu 1-60 + grid to 400"),
    "m(emu) LFV":    (100, 8000, "Z'->emu LFV dilepton (grid 100 GeV-8 TeV)"),
    # the tau LFV channels are published only up to the ~3 TeV where the hadronic-tau limit dies
    "m(etau) LFV":   (100, 3000, "Z'->etau LFV (hadronic-tau channel, grid 100 GeV-3 TeV)"),
    "m(mutau) LFV":  (100, 3000, "Z'->mutau LFV (hadronic-tau channel, grid 100 GeV-3 TeV)"),
    "m(ee) SS":      (200, 1300, "H++ same-sign ee pair scan"),
    "m(emu) SS":     (200, 1300, "H++ same-sign e-mu pair scan"),
    "m(mumu) SS":    (200, 1300, "H++ same-sign mumu pair scan"),
    "multilepton":   (50, 10000, "type-III seesaw / VLL / prompt heavy-N mass hypotheses (grid to 10 TeV)"),
    "m(ej)":         (200, 5000, "leptoquark pair -> eq (2006.05872) + QBH jet-electron; grid to 5 TeV"),
    "m(muj)":        (200, 5000, "leptoquark pair -> muq (2006.05872) + QBH jet-muon; grid to 5 TeV"),
    "m(tauj)":       (200, 5000, "leptoquark pair -> tau q mass scan (LQ3; grid to 5 TeV)"),
    # lepton + b-tagged jet: RPV stop -> b l (2406.18367, 200 GeV-1.5 TeV) and the b/c-tagged
    # leg of the scalar-LQ pair search (2006.05872)
    "m(eb)":         (200, 2000, "RPV stop -> b e (2406.18367) + b/c-tagged LQ pair (2006.05872)"),
    "m(mub)":        (200, 2000, "RPV stop -> b mu (2406.18367) + b/c-tagged LQ pair (2006.05872)"),
    "m(taub)":       (200, 2000, "3rd-generation LQ -> tau b"),
    "m(bj)":         (200, 6000, "b-tagged + light jet: b* -> b g, b-tagged dijet resonance "
                                 "(two-body ML anomaly scan 2307.01612)"),
    # trilepton resonance: RPV chargino/neutralino pair, l + leptonic Z, 100-1100 GeV
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
    # ---- channels whose scan is published (or was published pre-13 TeV) ----
    "m(eejj)":       (400, 7000, "W_R -> eejj (Keung-Senjanovic), published scan 0.4-7 TeV"),
    "m(mumujj)":     (400, 7000, "W_R -> mumujj (Keung-Senjanovic), published scan 0.4-7 TeV"),
    "m(Vgamma)":     (220, 6000, "X -> Z gamma / W gamma (ll/qq + photon)"),
    "m(egamma)":     (200, 5000, "excited electron e* -> e gamma (7/8 TeV only, 1201.3293/1308.1364; no 13 TeV search)"),
    "m(mugamma)":    (200, 5000, "excited muon mu* -> mu gamma (7/8 TeV only, 1201.3293/1308.1364; no 13 TeV search)"),
    # m(Zl) is NOT merely the 8 TeV VLL scan (1506.01291) inside the multilepton axis:
    # arXiv:2011.10543 is a Run-2 13 TeV trilepton RESONANCE scan (RPV chargino/neutralino,
    # 100-1100 GeV) that scans e+Z and mu+Z as separate spectra -- the m(eZ)/m(muZ) axes above.
    # m(Zb) is still omitted: scanned once at 8 TeV (1409.5500), no Run-2+ successor, and
    # adding it would shift N by <1%.
    "m(tW)":         (500, 3000, "excited b* -> tW"),
    "m(Wb)":         (800, 3000, "single VLQ T/Y -> Wb"),
    "m(Ht)":         (1000, 3000,"single VLQ T -> Ht"),
}

def scan_segments(obs):
    """List of scanned (lo, hi) segments [GeV] for a CANONICAL observable."""
    v = SCAN[canon(obs)]
    if isinstance(v[0], (list, tuple)):
        return [(float(lo), float(hi)) for lo, hi in v[0]]
    return [(float(v[0]), float(v[1]))]

def scan(obs):
    """Full extent (lo of first segment, hi of last) of the published-search window."""
    seg = scan_segments(obs)
    return seg[0][0], seg[-1][1]

def scan_source(obs): return SCAN[canon(obs)][-1]

import math as _math
def n_s(m_lo, m_hi, r):
    """Effective independent looks across [m_lo, m_hi] at fractional resolution r
    (Gross-Vitells resolution elements: constant fractional resolution sigma_M = r*M)."""
    if m_hi <= m_lo: return 0.0
    return _math.log(m_hi / m_lo) / r

def ns_scan(obs):
    """Budget trials for one observable: published-search window(s) at its resolution."""
    return sum(n_s(lo, hi, res(obs)) for lo, hi in scan_segments(obs))

def ns_achievable(obs):
    """Reference: effective looks over the full accessible window [floor, sqrt(s)]."""
    return n_s(floor(obs), SQRTS, res(obs))

def z_local_for_global5(N):
    """Local significance needed for 5 sigma global: Z^2_glob = Z^2_loc - 2 ln N (one-sided
    Gaussian-tail / Bonferroni approximation across N quasi-independent looks)."""
    return _math.sqrt(25.0 + 2.0 * _math.log(N)) if N > 1 else 5.0

# Canonical row order shared by EVERY Dobs-keyed plot: descending achievable trial contribution
# (pure observable metadata -- it depends on nothing but this module). Import CANON_ORDER and
# iterate it so every observable-keyed plot reads top-to-bottom identically.
CANON_ORDER = sorted(FLOOR, key=lambda o: -ns_achievable(o))

# ---- import-time self-consistency: the flavour layer must be complete and non-overlapping.
# A leaf without a SCAN entry would silently vanish from the budget; an inclusive parent that
# kept its own SCAN entry would be counted alongside its own leaves (double counting).
for _parent, _leaves in FLAV_SPLIT.items():
    assert _parent not in SCAN, f"inclusive {_parent!r} must not keep a SCAN entry (double count)"
    assert _parent not in FLOOR, f"inclusive {_parent!r} must not keep a FLOOR entry"
    for _lf in _leaves:
        assert _lf in SCAN and _lf in FLOOR, f"leaf {_lf!r} missing a SCAN/FLOOR entry"
# FLOOR additionally carries the OBS_MERGE parents (m(HH) 4b, m(SH)): those are raw decay
# labels that plots draw a floor line for, and canon() merges them away before any SCAN lookup.
assert set(SCAN) == set(FLOOR) - set(OBS_MERGE), \
    f"SCAN/FLOOR key mismatch: {set(SCAN) ^ (set(FLOOR) - set(OBS_MERGE))}"
del _parent, _leaves, _lf

# Production-mode kinematic ceiling [GeV]. A pair-produced resonance needs two heavy objects in
# the final state -> 2M < sqrt(s) -> M < sqrt(s)/2. Single s-channel / VBF / EFT reach ~sqrt(s)
# (parton-luminosity limited, no hard cap). Associated costs the recoil mass. NOTE P4 (cascade)
# assumes a *BSM* parent that is itself produced near sqrt(s)/2; a cascade from an SM parent is
# capped by that parent instead (e.g. Zd from h(125) < 62.5 GeV, H+ from top < m_t). Every
# observable here also carries a non-cascade mode, so the envelope max() below is unaffected,
# but keep this in mind if a pure-SM-parent observable is added.
PCAP = {"P1": SQRTS, "P2": SQRTS/2, "P3": SQRTS-400, "P4": SQRTS/2, "P5": SQRTS, "P6": SQRTS}

def prod_ceiling(pmodes):
    """Best-case (least restrictive) kinematic mass ceiling given the production modes present.
    == sqrt(s) if any single/VBF/EFT mode exists; sqrt(s)/2 for pair/cascade-only observables."""
    return max((PCAP.get(p, SQRTS) for p in pmodes), default=SQRTS)

def fmt_mass(m):
    """Clean mass label using GeV below 1 TeV and TeV above -- never 'k GeV'."""
    m = float(m)
    if m < 1000:
        return (f"{m:.1f} GeV" if m < 10 else f"{m:.0f} GeV")
    t = m / 1000.0
    return (f"{t:.0f} TeV" if abs(t - round(t)) < 0.05 else f"{t:.1f} TeV")

def fmt_range(lo, hi):
    return f"{fmt_mass(lo)}–{fmt_mass(hi)}"
