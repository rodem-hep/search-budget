import math

from .bump_observables import n_s, res

MIN_EVENTS, MIN_BINS = 100.0, 25

N_REF = 1.0e6
M_REF = 1000.0
R_REF = 0.05
P = 7.0

F   = {"j": 1.0, "b": 0.1, "V": 0.02, "t": 0.02, "H": 0.01, "g": 4e-3, "T": 3.5e-3,
       "e": 3e-3, "m": 3e-3, "X": 3e-3, "Z": 1e-5}
SYM = {"j": "m(jj)", "b": "m(bb)", "V": "m(VV)", "t": "m(tt)", "H": "m(HH)",
       "g": "m(gammagamma)", "T": "m(tautau)", "e": "m(ee)", "m": "m(mumu)", "X": "mT(ev)",
       "Z": None}

LENS_EFF = {"ht": 0.1, "disp": 1e-3, "vbf": 0.02, "isr": 0.2}


def weight(objs, f=None):
    f = F if f is None else f
    return math.prod(f[t] for t in objs)


def per_element(m, r, w):
    return N_REF * w * (m / M_REF) ** (1.0 - P) * (r / R_REF)


def one_event_mass(r, w):
    return M_REF * (N_REF * w * r / R_REF) ** (1.0 / (P - 1.0))


def events(lo, hi, w):
    if hi <= lo:
        return 0.0
    return (N_REF * w / R_REF) * ((lo / M_REF) ** (1.0 - P)
                                  - (hi / M_REF) ** (1.0 - P)) / (P - 1.0)


def gate(lo, hi, r, w):
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
