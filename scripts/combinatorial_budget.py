#!/usr/bin/env python3
"""Search budget of a fully combinatorial scan: exclusive object-multiplicity categories x every
2-to-4-object invariant mass they contain.

The object budget below is the one this study fixes for a scan of that kind -- a single-lepton
trigger, at most four objects, MET as a category split rather than a mass ingredient. It is a
design choice, not a measurement: the point is to price a scan of this shape against the published
program, and the answer moves only through ln N. Rules documented in docs/METHOD_NOTES.md.

scaled_scan.py reuses enumerate_scan() for the wider object alphabet, so the two scans differ only
in their inputs.

Writes results/tables/combinatorial_budget.csv.
"""
import os, sys, math, itertools, collections, csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import n_s, z_local_for_global5

# ------------------------------------------------------------------ inputs
ORDER  = "emjbZ"                                    # object types, in category-label order
NMAX   = {"e": 2, "m": 2, "j": 3, "b": 3, "Z": 1}   # per-type multiplicity ceiling
NOBJ   = 4                                          # max objects per category (MET excluded)
KMIN, KMAX = 2, 4                                   # mass-group sizes
TRIG   = "emZ"                                      # one of these must be present (lepton trigger)
LEPTON = "em"                                       # charged leptons, for the OS/SS category split

SIGMA  = {"e": 0.04, "m": 0.04, "Z": 0.04, "j": 0.10, "b": 0.20}   # fractional pT resolution
# scan window per group size k: [lo, hi] GeV (single-lepton-triggered turn-on -> stats ceiling)
WINDOW = {2: (100.0, 5000.0), 3: (150.0, 4000.0), 4: (200.0, 3000.0)}
MASS   = {"Z": 91.2}       # object masses that set the low edge; everything else counts as LIGHT
LIGHT  = 40.0

Cat  = collections.namedtuple("Cat", "n met os_ss")
Scan = collections.namedtuple("Scan", "rows N n_cat_mult n_cat n_hist by_size by_type")


def r_group(comp, sigma):
    """Fractional mass resolution of a k-body group from its object composition."""
    s2 = [sigma[t] ** 2 for t in comp]
    return 0.5 * math.sqrt(sum(s2) / len(s2))


def window(comp, mass):
    lo, hi = WINDOW[len(comp)]
    return max(lo, sum(mass.get(t, LIGHT) for t in comp)), hi     # object-mass floor


def enumerate_categories(order, nmax, nobj=NOBJ, trig=TRIG, lepton=LEPTON):
    """Exclusive multiplicity categories, MET-split, with a trigger object required."""
    cats = []
    for counts in itertools.product(*(range(nmax[t] + 1) for t in order)):
        n = dict(zip(order, counts))
        if not KMIN <= sum(counts) <= nobj:
            continue
        if trig and not any(n[t] for t in trig):
            continue
        charge = 2 if sum(n[t] for t in lepton) == 2 else 1   # OS/SS split of the dilepton cases
        for met in (0, 1):
            cats.append(Cat(n, met, charge))
    return cats


def enumerate_scan(order=ORDER, nmax=NMAX, sigma=SIGMA, mass=MASS, **kw):
    """Every (category, mass group) of the scan, with the looks each group costs."""
    cats = enumerate_categories(order, nmax, **kw)
    rows, N_trials, n_hist = [], 0.0, 0
    by_size, by_type = collections.Counter(), collections.Counter()
    for c in cats:
        objs = [t for t in order for _ in range(c.n[t])]       # indexed objects of the category
        for k in range(KMIN, min(KMAX, len(objs)) + 1):
            for idx in itertools.combinations(range(len(objs)), k):
                comp = tuple(objs[i] for i in idx)
                lo, hi = window(comp, mass)
                r = r_group(comp, sigma)
                ns = n_s(lo, hi, r) * c.os_ss                  # OS and SS are two disjoint looks
                N_trials += ns
                n_hist += c.os_ss
                by_size[k] += c.os_ss
                by_type["".join(sorted(comp))] += c.os_ss
                rows.append(("".join(f"{t}{c.n[t]}" for t in order) + f"_{c.met}met",
                             "".join(comp), c.os_ss, r, lo, hi, ns))
    # Two category counts, and they differ by the OS/SS split of the same-flavour dilepton cases:
    # quote them together, since N_trials counts an OS and an SS look separately.
    return Scan(rows, N_trials, len(cats), sum(c.os_ss for c in cats), n_hist, by_size, by_type)


def report(s):
    print(f"categories (exclusive, MET-split, >=2 objects)                       : {s.n_cat_mult:,}")
    print(f"  ... after splitting same-flavour dilepton cases into OS and SS     : {s.n_cat:,}")
    print(f"(category, mass-group) combinations                                  : {len(s.rows):,}")
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


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "group", "charge_split", "r", "M_lo_GeV", "M_hi_GeV", "n_s"])
        for row in rows:
            w.writerow([row[0], row[1], row[2], f"{row[3]:.3f}", f"{row[4]:.0f}", f"{row[5]:.0f}",
                        f"{row[6]:.1f}"])


if __name__ == "__main__":
    s = enumerate_scan()
    report(s)
    print()
    print("top object-type multisets by #histograms:")
    for t, v in s.by_type.most_common(12):
        print(f"  {t:6s} {v:5d}")
    out = os.path.join(ROOT, "results", "tables", "combinatorial_budget.csv")
    write_csv(s.rows, out)
    print(f"\nwrote results/tables/combinatorial_budget.csv ({len(s.rows):,} group rows)")
