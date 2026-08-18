#!/usr/bin/env python3
import os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import res
import combinatorial_budget as CB
import yield_model as YM

DATASETS = [("Run 2, 140 fb-1", 1.0), ("Run 2+3, ~400 fb-1", 3.0)]

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
MET_KEY    = "X"
ORDER_WIDE = "emTgjbtVHZ"
KMAX_WIDE  = 4
NOBJ_WIDE  = 4
NMAX_WIDE  = {k: NOBJ_WIDE for k in ORDER_WIDE}
MASS_WIDE  = {k: m for k, (_, m) in WIDE.items() if m is not None}
LEPTON_WIDE = "emT"

SIGMA_WIDE = {_k: (CB.SIGMA["Z"] if _ch is None else 2.0 * res(_ch))
              for _k, _ch in YM.SYM.items() if _k != MET_KEY}

SIGMA_BASE_DERIVED = {k: SIGMA_WIDE[k] for k in CB.ORDER}

BASE = dict(order=CB.ORDER, nmax=CB.NMAX, sigma=CB.SIGMA, mass=CB.MASS, trig=CB.TRIG,
            lepton=CB.LEPTON, nobj=CB.NOBJ, kmax=CB.KMAX)
WIDE_ARGS = dict(order=ORDER_WIDE, nmax=NMAX_WIDE, sigma=SIGMA_WIDE, mass=MASS_WIDE,
                 lepton=LEPTON_WIDE, kmax=KMAX_WIDE, nobj=NOBJ_WIDE, trig="")

ISR_MAX = 200.0
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
    for i, r in enumerate(rows):
        for li, (lk, _label, _rule, ok, cap) in enumerate(LENSES, 1):
            if not ok(len(r.group), r.ncat, r.lo):
                continue
            hi = r.hi if cap is None else min(r.hi, cap)
            _hs, lns, _ev, fits = YM.gate(r.lo, hi, r.r, r.w * YM.LENS_EFF[lk])
            yield i, lk, li, lns * r.split, fits
