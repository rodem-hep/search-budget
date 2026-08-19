import collections
import itertools
import math

from .bump_observables import z_local_for_global5
from . import yield_model as YM

ORDER  = "emjbZ"
NMAX   = {"e": 2, "m": 2, "j": 3, "b": 3, "Z": 1}
NOBJ   = 4
KMIN, KMAX = 2, 4
TRIG   = "emZ"
LEPTON = "em"
MET_W  = "X"

SIGMA  = {"e": 0.04, "m": 0.04, "Z": 0.04, "j": 0.10, "b": 0.20}
WINDOW = {2: (100.0, 5000.0), 3: (150.0, 4000.0), 4: (200.0, 3000.0)}
MASS   = {"Z": 91.2}
LIGHT  = 40.0

Cat  = collections.namedtuple("Cat", "n met os_ss")
Row  = collections.namedtuple("Row", "cat group split r lo hi hi_scan n_s w ncat")
Scan = collections.namedtuple("Scan",
                              "rows N n_cat_mult n_cat n_hist n_thin by_size by_type looks")


def r_group(comp, sigma):
    s2 = [sigma[t] ** 2 for t in comp]
    return 0.5 * math.sqrt(sum(s2) / len(s2))


def window(comp, mass):
    lo, hi = WINDOW[len(comp)]
    return max(lo, sum(mass.get(t, LIGHT) for t in comp)), hi


def enumerate_categories(order, nmax, nobj=NOBJ, trig=TRIG, lepton=LEPTON):
    def vectors(i, left, acc):
        if i == len(order):
            yield tuple(acc)
            return
        for v in range(min(nmax[order[i]], left) + 1):
            acc.append(v)
            yield from vectors(i + 1, left - v, acc)
            acc.pop()

    cats = []
    for counts in vectors(0, nobj, []):
        n = dict(zip(order, counts))
        if sum(counts) < KMIN:
            continue
        if trig and not any(n[t] for t in trig):
            continue
        charge = 2 if sum(n[t] for t in lepton) == 2 else 1
        for met in (0, 1):
            cats.append(Cat(n, met, charge))
    return cats


def enumerate_scan(order=ORDER, nmax=NMAX, sigma=SIGMA, mass=MASS, kmax=KMAX, collect=True,
                   weight=None, **kw):
    weight = YM.F if weight is None else weight
    cats = enumerate_categories(order, nmax, **kw)
    rows, N_trials, n_hist, n_thin = [], 0.0, 0, 0
    by_size, by_type = collections.Counter(), collections.Counter()
    looks = collections.Counter()
    for c in cats:
        objs = [t for t in order for _ in range(c.n[t])]
        w = math.prod([weight[t] for t in objs] + ([weight[MET_W]] if c.met else []))
        for k in range(KMIN, min(kmax, len(objs)) + 1):
            for idx in itertools.combinations(range(len(objs)), k):
                group = tuple(objs[i] for i in idx)
                lo, hi = window(group, mass)
                r = r_group(group, sigma)
                hi_scan, ns, _ev, ok = YM.gate(lo, hi, r, w)
                ns = ns * c.os_ss if ok else 0.0
                key = "".join(sorted(group))
                N_trials += ns
                if ok:
                    n_hist += c.os_ss
                    by_size[len(group)] += c.os_ss
                    by_type[key] += c.os_ss
                    looks[key] += ns
                else:
                    n_thin += c.os_ss
                if collect:
                    rows.append(Row("".join(f"{t}{c.n[t]}" for t in order) + f"_{c.met}met",
                                    "".join(group), c.os_ss, r, lo, hi, hi_scan if ok else lo,
                                    ns, w, sum(c.n.values())))
    return Scan(rows, N_trials, len(cats), sum(c.os_ss for c in cats), n_hist, n_thin, by_size,
                by_type, looks)


def report(s):
    print(f"categories (exclusive, MET-split, >=2 objects)                       : {s.n_cat_mult:,}")
    print(f"  ... after splitting same-flavour dilepton cases into OS and SS     : {s.n_cat:,}")
    print(f"(category, mass-group) combinations                                  : {len(s.rows):,}")
    print(f"too thin to fit ({YM.MIN_BINS} elements of >=1 event, {YM.MIN_EVENTS:.0f} events)"
          f"{'':16s}: {s.n_thin:,} histograms")
    print(f"mass spectra (histograms, index-specific groups, k=2..4, OS/SS-split): {s.n_hist:,}")
    print(f"  by group size: " + ", ".join(f"k={k}: {v:,}" for k, v in sorted(s.by_size.items())))
    print(f"distinct object-type multisets (physics observables)                 : {len(s.by_type)}")
    print()
    print(f"N_trials = {s.N:,.0f}    ln N = {math.log(s.N):.2f}")
    print(f"Z_local for 5 sigma global = {z_local_for_global5(s.N):.2f}")
    print(f"  r band (x0.5..x2 -> N x2..x0.5): "
          f"{z_local_for_global5(s.N*0.5):.2f} - {z_local_for_global5(s.N*2):.2f}")
    print(f"  a local 5 sigma is worth Z_global = "
          f"{math.sqrt(max(25.0 - 2*math.log(s.N), 0)):.2f} sigma")
