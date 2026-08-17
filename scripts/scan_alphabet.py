#!/usr/bin/env python3
"""The object alphabet a fully combinatorial scan is built from, and the rules that bound it.

Definitions only, no outputs: scaled_scan.py runs the scan over this alphabet and
budget_uncertainty.py reprices it under every declared input variation, so the two share one
statement of what the scan is.

  ORDER_WIDE / NMAX_WIDE / KMAX_WIDE / NOBJ_WIDE   ten object types, no per-type ceiling, at most
                                                   four objects per category and per mass
  SIGMA_WIDE                                       fractional sigma per object, sigma = 2 r of its
                                                   own symmetric published channel (no new input)
  MASS_WIDE / LEPTON_WIDE                          the mass floor each object contributes, and the
                                                   charged leptons that split a category OS/SS
  WIDE_ARGS / BASE / SIGMA_BASE_DERIVED            the argument sets combinatorial_budget runs on
  LENSES / lens_views                              the four event-level lenses laid over an
                                                   unchanged mass axis, one at a time

Missing energy is a category split and a yield factor, never an ingredient of a mass, so the
transverse-mass axes of the model-driven budget are outside this scan's reach.

Reads  scripts/bump_observables.py (published resolutions), scripts/combinatorial_budget.py (the
       five-object rules), scripts/yield_model.py (the symmetric channel per object type).
Pure standard library.
"""
import os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import res
import combinatorial_budget as CB
import yield_model as YM

# Integrated luminosity relative to the Run-2 dataset the yield anchor is set on, ignoring the rise in
# high-mass cross sections from 13 to 13.6 TeV, which acts in the same direction. The last entry is the
# dataset every headline is priced on.
DATASETS = [("Run 2, 140 fb-1", 1.0), ("Run 2+3, ~400 fb-1", 3.0)]

# key -> (label, mass-floor contribution in GeV). The symmetric channel each resolution is read from
# and the yield factor each object carries are yield_model.SYM and yield_model.F.
WIDE = {
    "e": ("e",           None),
    "m": ("mu",          None),
    "T": ("tau_had",     None),
    "g": ("gamma",       None),
    "j": ("light jet",   None),
    "b": ("b-jet",       None),
    "t": ("boosted top", 173.0),
    "V": ("boosted W/Z",  91.2),
    "H": ("boosted H",   125.0),
    "Z": ("leptonic Z",   91.2),
    "X": ("MET",         None),
}
MET_KEY    = "X"                                 # a category split and a yield factor, never a mass
ORDER_WIDE = "emTgjbtVHZ"                        # MET is not an entry: it is the category's met flag
KMAX_WIDE  = 4                                   # objects per mass combination, MET excluded
NOBJ_WIDE  = 4                                   # objects per category, strictly
NMAX_WIDE  = {k: NOBJ_WIDE for k in ORDER_WIDE}  # no per-type ceiling
MASS_WIDE  = {k: m for k, (_, m) in WIDE.items() if m is not None}
LEPTON_WIDE = "emT"                              # hadronic taus are charged: they split OS/SS too

# sigma = 2 r of the symmetric channel. Missing energy needs none: it never enters a mass.
SIGMA_WIDE = {_k: (CB.SIGMA["Z"] if _ch is None else 2.0 * res(_ch))
              for _k, _ch in YM.SYM.items() if _k != MET_KEY}

# five object types with resolutions from the same inversion: separates the alphabet from the
# prescription, since the five-object scan declares its sigma per object instead of deriving it
SIGMA_BASE_DERIVED = {k: SIGMA_WIDE[k] for k in CB.ORDER}

# the two argument sets enumerate_scan() is run on
BASE = dict(order=CB.ORDER, nmax=CB.NMAX, sigma=CB.SIGMA, mass=CB.MASS, trig=CB.TRIG,
            lepton=CB.LEPTON, nobj=CB.NOBJ, kmax=CB.KMAX)
WIDE_ARGS = dict(order=ORDER_WIDE, nmax=NMAX_WIDE, sigma=SIGMA_WIDE, mass=MASS_WIDE,
                 lepton=LEPTON_WIDE, kmax=KMAX_WIDE, nobj=NOBJ_WIDE, trig="")

# ------------------------------------------------------------------ selection lenses
# An extra event-level requirement on an unchanged mass axis is another view of the same spectrum, and
# another look. Four of the handles a wide search would reach for are already inside this enumeration
# and must not be counted twice: high MET is the category's met split, high jet or lepton multiplicity
# is what the exclusive categories are, and b-tag and tau enrichment are the b and T types of the
# alphabet. The four below are orthogonal to the object content and to the mass axis. Conservative on
# every count: one lens at a time and never a product of two, the lens leaves the axis, its resolution
# and its window alone, the objects a lens needs count against the same four-object ceiling, and the
# view has to pass the statistics requirement at its own efficiency (yield_model.LENS_EFF), so a lens
# on a thinly populated spectrum is not available at all. k is the objects in the mass, n those in the
# category.
ISR_MAX = 200.0                       # an ISR-recoil view buys acceptance only at the low-mass end
LENSES = [
    ("ht",   "high HT or Meff",    "activity outside the mass",
     lambda k, n, lo: n > k, None),
    ("disp", "displaced activity", "any reconstructed mass",
     lambda k, n, lo: True, None),
    ("vbf",  "forward jet pair",   "two free slots for the tag jets",
     lambda k, n, lo: n <= 2, None),
    ("isr",  "ISR jet",            "one free slot, and a low-mass end",
     lambda k, n, lo: n <= 3 and lo < ISR_MAX, ISR_MAX),
]


def lens_views(rows):
    """(row index, lens key, its index in LENSES, looks it costs, does it fit) per available view."""
    for i, r in enumerate(rows):
        for li, (lk, _label, _rule, ok, cap) in enumerate(LENSES, 1):
            if not ok(len(r.group), r.ncat, r.lo):
                continue
            hi = r.hi if cap is None else min(r.hi, cap)
            _hs, lns, _ev, fits = YM.gate(r.lo, hi, r.r, r.w * YM.LENS_EFF[lk])
            yield i, lk, li, lns * r.split, fits
