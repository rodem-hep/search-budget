#!/usr/bin/env python3

PUBLIC_OBS = {
    "Hidden Abelian Higgs (HAHM)":   ["m(ee) (Zd)", "m(mumu) (Zd)"],
    "Simplified DM (dijet mediator)":["m(jj)"],
    "2HDM (general/typeII/CPV)":     ["m(tautau)", "m(Vh)", "m(tt)", "m(cb) dijet", "m(tb)", "m(bb)"],
    "N2HDM / 2HDM+S (h->aa)":        ["m(mumu) (Zd)", "m(bb)", "m(tautau)"],
    "Georgi-Machacek":               ["m(VV)", "m(ee) SS", "m(emu) SS", "m(mumu) SS"],
    "Singlet scalar / SM+Scalars":   ["m(HH)", "m(VV)", "m(Vgamma)"],
    "Higgs portal":                  ["m(HH)"],
    "TRSM (singlet/triplet)":        ["m(HH)", "m(VV)"],
    "HeavyHiggs THDM":               ["m(tt)", "m(VV)"],
    "Vector-like quark (VLQ)":       ["m(ttZ)/m(Zt)", "m(tb)", "m(Wb)", "m(Ht)"],
    "Vector-like lepton (VLL)":      ["multilepton"],
    "Leptophilic gauge boson":       ["m(ee)", "m(mumu)"],
    "Scalar leptoquark (S1)":        ["m(ej)", "m(muj)", "m(tauj)",
                                      "m(eb)", "m(mub)", "m(taub)",
                                      "m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
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
    "Excited quark (q*/b*)":         ["m(jgamma)", "m(jj)", "m(bj)", "m(tW)"],
    "Quantum black hole":            ["m(multi)", "m(ej)", "m(muj)", "m(jgamma)",
                                      "m(emu) LFV", "m(etau) LFV", "m(mutau) LFV"],
    "RPV electroweakino (trilepton)":["m(eZ)", "m(muZ)"],
    "Excited lepton (l*)":           ["m(egamma)", "m(mugamma)"],
    "Excited boson (W*/Z*)":         ["m(VV)", "m(Vgamma)"],
    "Sgluon":                        ["m(tt)", "m(3j)"],
    "Top-philic":                    ["m(tt)"],
    "SILH / Little Higgs":           ["m(VV)", "m(ee)", "m(mumu)"],
    "MSSM/NMSSM/RPV SUSY":           ["m(tautau)", "m(bb)", "m(3j)"],
    "Stop/scharm (RPV)":             ["m(tt)/m(jj)", "m(eb)", "m(mub)"],
    "Composite / NJL":               ["m(jj)", "m(VV)"],
    "Diquark / color-sextet":        ["m(jj)", "m(tt)"],
    "Vector-like confinement":       ["m(jj)", "m(gammagamma)"],
    "Toponium":                      ["m(tt)"],
}

WIDTH = {
    "Hidden Abelian Higgs (HAHM)":   ("narrow",    "Gamma/M ~ eps^2, <<1e-3"),
    "N2HDM / 2HDM+S (h->aa)":        ("narrow",    "light a, <<1%"),
    "ALP":                           ("narrow",    "Gamma/M ~ (M/f)^2, <<1%"),
    "Singlet scalar / SM+Scalars":   ("narrow",    "mixing-suppressed, <~1%"),
    "Higgs portal":                  ("narrow",    "mixing-suppressed, <~1%"),
    "TRSM (singlet/triplet)":        ("narrow",    "mixing-suppressed, <~1%"),
    "Leptophilic gauge boson":       ("narrow",    "small g, <~1%"),
    "Georgi-Machacek":               ("narrow",    "H++ 1e-3 to a few %"),
    "Type-II seesaw":                ("narrow",    "H++ <~1%"),
    "Excited quark (q*/b*)":         ("narrow",    "2-4% at f=1, Lambda=M"),
    "Diquark / color-sextet":        ("narrow",    "1-5% for perturbative coupling"),
    "Sgluon":                        ("narrow",    "<~2%"),
    "Excited boson (W*/Z*)":         ("narrow",    "a few %"),
    "Vector-like confinement":       ("narrow",    "bound states, <~1%"),
    "Stop/scharm (RPV)":             ("narrow",    "lambda''-suppressed, <<1%"),
    "LRSM / Alt-LRSM":               ("narrow",    "W_R 2-3%"),
    "Heavy neutrino / HNL (prompt)": ("narrow",    "N itself <<1%"),
    "RPV electroweakino (trilepton)":("narrow",    "electroweak width, <<1%"),
    "Scalar leptoquark (S1)":        ("narrow",    "lambda^2/16pi ~ 2% at lambda=1"),
    "Leptoquark NLO (mix/nomix)":    ("narrow",    "lambda^2/16pi ~ 2% at lambda=1"),
    "Randall-Sundrum / Radion":      ("narrow",    "radion <~1%"),
    "W'":                            ("narrow",    "SSM ~3%, below the mT and m(tb) resolutions"),

    "Minimal Z' / U(1)":             ("benchmark", "SSM ~3% exceeds r=0.015 on m(ee); E6 0.5-1.2%"),
    "Vector triplet (HVT)":          ("benchmark", "model A ~2%, B ~5%; broad scenarios published"),
    "KK graviton (Gstar)":           ("benchmark", "0.014% at k/Mpl=0.01, ~6% for bulk RS"),
    "Simplified DM (dijet mediator)":("benchmark", "1% to >30% over the coupling grid"),
    "Vector-like quark (VLQ)":       ("benchmark", "pair narrow, single production 10-50%"),
    "Vector leptoquark (U1)":        ("benchmark", "~20% at the lambda~3 flavour-anomaly point"),
    "Technicolor / TC2":             ("benchmark", "technirho a few %, walking variants broad"),
    "Excited lepton (l*)":           ("benchmark", "a few % against r=0.015-0.03 on the lgamma axes"),
    "SILH / Little Higgs":           ("benchmark", "a few %, marginal on m(ee)"),
    "MSSM/NMSSM/RPV SUSY":           ("benchmark", "A/H->tautau 10-20% at tanbeta~50, still under r"),
    "2HDM (general/typeII/CPV)":     ("benchmark", "narrow on the scalar axes; see NONPEAK_ON"),
    "HeavyHiggs THDM":               ("benchmark", "narrow on m(VV); see NONPEAK_ON"),
    "Top-philic":                    ("benchmark", "see NONPEAK_ON"),

    "KK gluon (RS)":                 ("broad",     "15-30%, 2-4x the m(tt) resolution"),
    "Coloron / Axigluon":            ("broad",     "5-20% against r=0.05 on m(jj)"),
    "Composite / NJL":               ("broad",     "10-50%, strongly coupled"),

    "Quantum black hole":            ("nonpeak",   "threshold turn-on plus continuum, no peak"),
    "Large ED / UED / HEIDI":        ("nonpeak",   "non-resonant high-mass tail / degenerate tower"),
    "Type-III seesaw":               ("nonpeak",   "pair-produced triplets, a counting signature"),
    "Vector-like lepton (VLL)":      ("nonpeak",   "pair-produced, a counting signature"),
    "Toponium":                      ("nonpeak",   "threshold effect fixed at 2m_t, not scannable"),
}

NONPEAK_ON = {
    "Heavy neutrino / HNL (prompt)": {"multilepton"},
    "2HDM (general/typeII/CPV)":     {"m(tt)"},
    "HeavyHiggs THDM":               {"m(tt)"},
    "Top-philic":                    {"m(tt)"},
}

NONPEAK = {m for m, (c, _) in WIDTH.items() if c == "nonpeak"}

def peaks(model, obs):
    return model not in NONPEAK and obs not in NONPEAK_ON.get(model, ())

def nonpeak_only(obs, models):
    return bool(models) and not any(peaks(m, obs) for m in models)

NSEL = {
    "m(ee)":         (1, "single dielectron scan (barrel/endcap are fit categories, combined)"),
    "m(mumu)":       (1, "single dimuon scan (combined charge/eta categories)"),
    "m(ee) (Zd)":    (1, "prompt e-LJ / 4e scan"),
    "m(mumu) (Zd)":  (1, "prompt muon LJ / 4mu scan"),
    "m(emu) LFV":    (1, "LFV e-mu"),
    "m(etau) LFV":   (1, "LFV e-tau (hadronic tau)"),
    "m(mutau) LFV":  (1, "LFV mu-tau (hadronic tau)"),
    "m(ee) SS":      (1, "H++ same-sign ee"),
    "m(emu) SS":     (1, "H++ same-sign e-mu"),
    "m(mumu) SS":    (1, "H++ same-sign mu-mu"),
    "multilepton":   (4, "3L / 4L signal regions by flavour & charge (axis deliberately unsplit)"),
    "m(gammagamma)": (2, "spin-0 (ggF) + spin-2/VBF (converted/unconverted categories)"),
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

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from bump_observables import SCAN as _SCAN, FLAV_SPLIT as _FS

_bad = {o for objs in PUBLIC_OBS.values() for o in objs if o in _FS} | {o for o in NSEL if o in _FS}
assert not _bad, f"flavour-inclusive label used where a leaf is required: {sorted(_bad)}"
_unknown = {o for objs in PUBLIC_OBS.values() for o in objs if o not in _SCAN} | \
           {o for o in NSEL if o not in _SCAN}
assert not _unknown, f"observable with no SCAN window: {sorted(_unknown)}"

assert set(WIDTH) == set(PUBLIC_OBS), \
    f"WIDTH and PUBLIC_OBS disagree: {sorted(set(WIDTH) ^ set(PUBLIC_OBS))}"
assert set(NONPEAK_ON) <= set(PUBLIC_OBS), f"NONPEAK_ON has no model: {sorted(NONPEAK_ON)}"
_stray = {(m, o) for m, os_ in NONPEAK_ON.items() for o in os_ if o not in PUBLIC_OBS[m]}
assert not _stray, f"NONPEAK_ON names an axis the model does not populate: {sorted(_stray)}"
