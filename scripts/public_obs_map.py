#!/usr/bin/env python3
"""Public BSM model -> bump observable(s) mapping (which spectra does public theory
motivate). Factored out of search_budget.py so the coverage-matrix plot and the budget
share ONE mapping. Keyed on the same canonical observables as bump_observables.SCAN; models
whose only signature is non-bump (mono-X MET, dE/dx, displaced-only) are omitted."""

PUBLIC_OBS = {
    "Hidden Abelian Higgs (HAHM)":   ["m(ee) (Zd)", "m(mumu) (Zd)"],
    "Simplified DM (dijet mediator)":["m(jj)"],
    "2HDM (general/typeII/CPV)":     ["m(tautau)", "m(Vh)", "m(tt)", "m(cb) dijet", "m(tb)", "m(bb)"],
    # a->ll is Yukawa-ordered: the dimuon channel carries the search, a->ee is helicity-suppressed
    "N2HDM / 2HDM+S (h->aa)":        ["m(mumu) (Zd)", "m(bb)", "m(tautau)"],  # a->mumu/bb/tautau
    "Georgi-Machacek":               ["m(VV)", "m(ee) SS", "m(emu) SS", "m(mumu) SS"],
    "Singlet scalar / SM+Scalars":   ["m(HH)", "m(VV)", "m(Vgamma)"],
    "Higgs portal":                  ["m(HH)"],
    "TRSM (singlet/triplet)":        ["m(HH)", "m(VV)"],
    "HeavyHiggs THDM":               ["m(tt)", "m(VV)"],
    "Vector-like quark (VLQ)":       ["m(ttZ)/m(Zt)", "m(tb)", "m(Wb)", "m(Ht)"],
    "Vector-like lepton (VLL)":      ["multilepton"],
    "Leptophilic gauge boson":       ["m(ee)", "m(mumu)"],
    # the LQ pair search explicitly b/c-tags one leg (2006.05872), and a 3rd-generation LQ
    # decays to b+lepton -- a different object pairing from the light-jet leg
    "Scalar leptoquark (S1)":        ["m(ej)", "m(muj)", "m(tauj)",
                                      "m(eb)", "m(mub)", "m(taub)",
                                      "m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
    # U1 is the 3rd-generation-philic vector LQ: the tau+b leg is the one it really predicts
    "Vector leptoquark (U1)":        ["m(tauj)", "m(taub)", "m(tautau)",
                                      "m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
    "Leptoquark NLO (mix/nomix)":    ["m(ej)", "m(muj)", "m(tauj)", "m(eb)", "m(mub)", "m(taub)"],
    "Minimal Z' / U(1)":             ["m(ee)", "m(mumu)", "m(emu) LFV", "m(etau) LFV",
                                      "m(mutau) LFV", "m(jj)", "m(bb)", "m(tt)"],
    "W'":                            ["mT(ev)", "mT(muv)", "m(tb)", "mT(taunu)", "m(jj)"],
    "Vector triplet (HVT)":          ["m(VV)", "m(Vh)"],
    "Randall-Sundrum / Radion":      ["m(gammagamma)", "m(VV)", "m(HH)"],
    "KK graviton (Gstar)":           ["m(gammagamma)", "m(ee)", "m(mumu)", "m(VV)"],
    "KK gluon (RS)":                 ["m(tt)"],
    "Large ED / UED / HEIDI":        ["m(gammagamma)", "m(ee)", "m(mumu)", "m(multi)"],
    "Type-III seesaw":               ["multilepton"],
    "Type-II seesaw":                ["m(ee) SS", "m(emu) SS", "m(mumu) SS", "m(VV)"],
    "LRSM / Alt-LRSM":               ["m(ee)", "m(mumu)", "mT(ev)", "mT(muv)",
                                      "m(eejj)", "m(mumujj)"],
    "Heavy neutrino / HNL (prompt)": ["m(eejj)", "m(mumujj)", "multilepton"],
    "ALP":                           ["m(gammagamma)", "m(Vgamma)"],
    "Technicolor / TC2":             ["m(VV)", "m(ee)", "m(mumu)"],
    "Coloron / Axigluon":            ["m(jj)", "m(tt)", "m(tt)/m(jj)"],
    # b* -> b g is a b-tagged+light-jet resonance, distinct from the inclusive dijet axis
    "Excited quark (q*/b*)":         ["m(jgamma)", "m(jj)", "m(bj)", "m(tW)"],
    # QBH was missing from this map entirely despite being a catalogue category AND a benchmark
    # of the LFV dilepton search (2307.08567). Its jet+lepton grids are two-body resonances.
    "Quantum black hole":            ["m(multi)", "m(ej)", "m(muj)", "m(jgamma)",
                                      "m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
    # RPV chargino/neutralino pair -> l + leptonic Z: a genuine Run-2 trilepton RESONANCE scan
    # (2011.10543), not a counting experiment -- the reason m(eZ)/m(muZ) now exist as axes
    "RPV electroweakino (trilepton)":["m(eZ)", "m(muZ)"],
    "Excited lepton (l*)":           ["m(egamma)", "m(mugamma)"],
    "Excited boson (W*/Z*)":         ["m(VV)", "m(Vgamma)"],
    "Sgluon":                        ["m(tt)", "m(3j)"],
    "Top-philic":                    ["m(tt)"],
    "SILH / Little Higgs":           ["m(VV)", "m(ee)", "m(mumu)"],
    "MSSM/NMSSM/RPV SUSY":           ["m(tautau)", "m(bb)", "m(3j)"],
    # the RPV stop also has a LEPTONIC decay t~ -> b l (2406.18367), a lepton+b-jet resonance;
    # only its hadronic mode was mapped before
    "Stop/scharm (RPV)":             ["m(tt)/m(jj)", "m(eb)", "m(mub)"],
    "Composite / NJL":               ["m(jj)", "m(VV)"],
    "Diquark / color-sextet":        ["m(jj)", "m(tt)"],
    "Vector-like confinement":       ["m(jj)", "m(gammagamma)"],
    "Toponium":                      ["m(tt)"],
}

# ---------------------------------------------------------------- event selections
# NSEL[obs] = (n_channels, "the real analysis channels that scan this mass axis
# separately"), anchored to PUBLISHED ATLAS search designs -- public information.
# A channel = a distinct event selection producing its own bump spectrum (not merely a
# fit category that is statistically combined).
# NB LEPTON FLAVOUR IS NO LONGER COUNTED HERE. It used to enter as a multiplier
# (m(ll) -> 2 for ee/mumu, m(lj) -> 3 for e/mu/tau, ...); it is now a first-class mass axis
# (bump_observables.FLAV_SPLIT), so each flavour leaf already contributes its own n_s with its
# own resolution. Re-adding a flavour factor below would double-count it. What remains here is
# only the WITHIN-CHANNEL selection multiplicity: b-tag bins, boost regimes, sub-decay modes.
NSEL = {
    # --- sharp leptonic / photonic (flavour now lives in the axis, not in this factor) ---
    "m(ee)":         (1, "single dielectron scan (barrel/endcap are fit categories, combined)"),
    "m(mumu)":       (1, "single dimuon scan (combined charge/eta categories)"),
    "m(ee) (Zd)":    (1, "prompt e-LJ / 4e scan"),
    # the mixed 2e2mu pairing is folded into the two flavour leaves rather than counted as a
    # third selection -- keeping it here would re-introduce a flavour factor on top of the axis
    "m(mumu) (Zd)":  (1, "prompt muon LJ / 4mu scan"),
    "m(emu) LFV":    (1, "LFV e-mu"),
    "m(etau) LFV":   (1, "LFV e-tau (hadronic tau)"),
    "m(mutau) LFV":  (1, "LFV mu-tau (hadronic tau)"),
    "m(ee) SS":      (1, "H++ same-sign ee"),
    "m(emu) SS":     (1, "H++ same-sign e-mu"),
    "m(mumu) SS":    (1, "H++ same-sign mu-mu"),
    "multilepton":   (4, "3L / 4L signal regions by flavour & charge (axis deliberately unsplit)"),
    "m(gammagamma)": (2, "spin-0 (ggF) + spin-2/VBF (converted/unconverted categories)"),
    # LQ PAIR production gives 2 leptons + 2 jets with an AMBIGUOUS pairing, so the published
    # search histograms all four lepton-jet combinations (e0j0, e0j1, e1j0, e1j1; 2006.05872).
    # Each is its own scanned spectrum.
    "m(ej)":         (4, "LQ pair -> eq: all four lepton-jet pairings scanned (2006.05872)"),
    "m(muj)":        (4, "LQ pair -> muq: all four lepton-jet pairings scanned (2006.05872)"),
    "m(tauj)":       (2, "LQ->tau q: tau_had + tau_lep selections"),
    "m(eb)":         (2, "RPV stop -> b e + b-tagged LQ leg (both pairings)"),
    "m(mub)":        (2, "RPV stop -> b mu + b-tagged LQ leg (both pairings)"),
    "m(taub)":       (2, "LQ3 -> tau b: tau_had + tau_lep"),
    "m(bj)":         (2, "b* -> b g: b-tagged leading + subleading jet pairing"),
    "m(eZ)":         (1, "trilepton e+Z resonance"),
    "m(muZ)":        (1, "trilepton mu+Z resonance"),
    "mT(ev)":        (1, "W'->e nu"),
    "mT(muv)":       (1, "W'->mu nu"),
    "mT(taunu)":     (1, "W'->tau nu (single hadronic-tau channel)"),
    # --- hadronic / heavy-object: split by decay mode, b-tag, boost ---
    "m(jj)":         (2, "low-mass (TLA/ISR) + high-mass inclusive dijet"),
    "m(cb) dijet":   (2, "light H+->cb in tt: hadronic + semileptonic tag"),
    "m(3j)":         (1, "RPV/sgluon three-jet resonance (single SR)"),
    "m(bb)":         (3, "h->aa->4b low-mass + A/Y->bb resolved + boosted"),
    "m(jgamma)":     (1, "q*->q gamma (single photon+jet SR)"),
    "m(tt)":         (4, "ttbar res: l+jets resolved, l+jets boosted, all-had, dilepton"),
    "m(tt)/m(jj)":   (2, "coloron pair -> tt or jj decay legs"),
    "m(tb)":         (3, "W'/H+->tb: 0-lepton, 1-lepton x b-tag"),
    "m(VV)":         (9, "WW/WZ/ZZ x qqqq / lvqq / llqq / vvqq / lvlv, ggF+VBF"),
    "m(Vh)":         (6, "Wh/Zh x 0/1/2-lepton x h->bb resolved/boosted"),
    "m(HH)":         (6, "bbbb resolved + boosted, bbtautau, bbgammagamma, bbVV, multilepton"),
    "m(ttZ)/m(Zt)":  (2, "VLQ single/pair T->tZ selections"),
    "m(multi)":      (1, "QBH / multijet (single high-multiplicity SR)"),
    "m(tautau)":     (3, "H/Z'->tautau: tau_e tau_h, tau_mu tau_h, tau_h tau_h"),
    # --- channels whose scan is published (or was published pre-13 TeV) ---
    "m(eejj)":       (1, "W_R->eejj"),
    "m(mumujj)":     (1, "W_R->mumujj"),
    "m(Vgamma)":     (2, "X->Zgamma: Z->ll and boosted Z->qq (+Wgamma)"),
    "m(egamma)":     (1, "excited electron e*"),
    "m(mugamma)":    (1, "excited muon mu*"),
    "m(tW)":         (2, "b*->tW: leptonic + hadronic top"),
    "m(Wb)":         (1, "single VLQ T/Y->Wb (1-lepton SR)"),
    "m(Ht)":         (1, "single VLQ T->Ht (h->bb tagged)"),
}
NSEL_DEFAULT = 2

def nsel(obs):
    return NSEL.get(obs, (NSEL_DEFAULT, "default"))[0]

# ---- import-time guard: both maps must be written in LEAF (flavour-resolved) channels. A
# flavour-inclusive label here would silently be treated as a channel of its own and counted
# alongside its own leaves.
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from bump_observables import SCAN as _SCAN, FLAV_SPLIT as _FS

_bad = {o for objs in PUBLIC_OBS.values() for o in objs if o in _FS} | {o for o in NSEL if o in _FS}
assert not _bad, f"flavour-inclusive label used where a leaf is required: {sorted(_bad)}"
_unknown = {o for objs in PUBLIC_OBS.values() for o in objs if o not in _SCAN} | \
           {o for o in NSEL if o not in _SCAN}
assert not _unknown, f"observable with no SCAN window: {sorted(_unknown)}"
