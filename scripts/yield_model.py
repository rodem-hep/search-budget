#!/usr/bin/env python3
"""Whether a hypothetical spectrum holds enough events to be fitted at all.

A trials count may only contain histograms a search could actually fit. The requirement, applied
wherever a spectrum is hypothetical rather than published:

    at least MIN_EVENTS events in the histogram, and at least MIN_BINS bins holding one event or
    more, with one bin = one resolution element, the unit the trials count is built from.

The yield model behind it is declared and order-of-magnitude. A histogram's content is its
background, and every background here is a steeply falling mass spectrum, so

    n(m) = N_REF * W * (m/M_REF)^(1-P) * (r/R_REF)          events in one resolution element,

anchored on the light-jet pair at M_REF for a Run-2 dataset, with W the product of a per-object factor
F relative to a light jet and r the fractional resolution (a sharper channel has narrower bins and so
fewer events in each). Because n(m) falls, the populated part of a window is a prefix of it: the
requirement acts by truncating each window at the mass where an element holds one event, and by
dropping the spectrum when fewer than MIN_BINS elements are left.

F is statistics, not signal cross section: a tagged hadronic object costs its mistag rate, and a lepton
or genuine missing energy costs the price of an electroweak process against QCD. The factors are set so
that the model reproduces the published symmetric channel of each object type to an order of magnitude
(run this file for that table). Two consequences of the single power law are worth stating: it
overestimates yields below a few hundred GeV, where real spectra turn over, which can only make the
MIN_EVENTS test easier to pass and so errs towards more looks; and it puts the one-event mass of the
dijet spectrum at 10 TeV against a published fit that stops near 8, in the same direction.

Published windows are never gated: a published search demonstrates its own feasibility. This module is
for the spectra nobody has scanned.

Reads  scripts/bump_observables.py (resolutions, for the calibration table only).
Imported by combinatorial_budget.py, scaled_scan.py, two_body_matrix.py. Pure standard library.
"""
import math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import n_s, res

MIN_EVENTS, MIN_BINS = 100.0, 25       # the requirement

N_REF = 1.0e6      # events in one resolution element of the light-jet pair spectrum at M_REF
M_REF = 1000.0     # GeV, the anchor mass
R_REF = 0.05       # the resolution that anchor is quoted at
P = 7.0            # dN/dm ~ m^-P for a hadronic background; band 6 to 8

# object -> yield factor relative to a light jet, and the symmetric channel it is calibrated on
F   = {"j": 1.0, "b": 0.1, "V": 0.02, "t": 0.02, "H": 0.01, "g": 4e-3, "T": 3.5e-3,
       "e": 3e-3, "m": 3e-3, "X": 3e-3, "Z": 1e-5}
SYM = {"j": "m(jj)", "b": "m(bb)", "V": "m(VV)", "t": "m(tt)", "H": "m(HH)",
       "g": "m(gammagamma)", "T": "m(tautau)", "e": "m(ee)", "m": "m(mumu)", "X": "mT(ev)",
       "Z": None}

# An event-level lens keeps the axis but costs statistics: the efficiency of the requirement.
LENS_EFF = {"ht": 0.1, "disp": 1e-3, "vbf": 0.02, "isr": 0.2}


def weight(objs, f=None):
    """Yield of a category relative to a light-jet pair: the product over its objects."""
    f = F if f is None else f
    return math.prod(f[t] for t in objs)


def per_element(m, r, w):
    """Events in the resolution element at mass m."""
    return N_REF * w * (m / M_REF) ** (1.0 - P) * (r / R_REF)


def one_event_mass(r, w):
    """Mass at which a resolution element holds a single event."""
    return M_REF * (N_REF * w * r / R_REF) ** (1.0 / (P - 1.0))


def events(lo, hi, w):
    """Events between lo and hi, summed over resolution elements (independent of r)."""
    if hi <= lo:
        return 0.0
    return (N_REF * w / R_REF) * ((lo / M_REF) ** (1.0 - P)
                                  - (hi / M_REF) ** (1.0 - P)) / (P - 1.0)


def gate(lo, hi, r, w):
    """(scannable high edge, elements over it, events in it, can it be fitted).

    The window is truncated where an element drops below one event; the verdict is the requirement.
    """
    hi_scan = min(hi, one_event_mass(r, w))
    n = n_s(lo, hi_scan, r)
    ev = events(lo, hi_scan, w)
    return hi_scan, n, ev, (n >= MIN_BINS and ev >= MIN_EVENTS)


if __name__ == "__main__":
    print(f"n(m) = {N_REF:.0e} * W * (m/{M_REF:.0f} GeV)^(1-{P:.0f}) * (r/{R_REF}) events per "
          f"resolution element")
    print(f"requirement: >= {MIN_BINS} elements holding >= 1 event, and >= {MIN_EVENTS:.0f} events\n")
    print("obj  factor  symmetric channel      r      events/element at 1 TeV   1-event mass")
    for k, ch in SYM.items():
        if ch is None:
            print(f"  {k}  {F[k]:6.1e}  {'(no symmetric channel)':22s}")
            continue
        r, w = res(ch), weight(k + k)
        print(f"  {k}  {F[k]:6.1e}  {ch:22s} {r:5.3f} {per_element(M_REF, r, w):17,.1f}"
              f"        {one_event_mass(r, w):8,.0f} GeV")
    print("\nlens efficiencies: " + ", ".join(f"{k} {v:g}" for k, v in LENS_EFF.items()))
