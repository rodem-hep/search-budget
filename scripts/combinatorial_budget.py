#!/usr/bin/env python3
"""Search budget of a fully combinatorial scan: exclusive object-multiplicity categories x every
2-to-4-object invariant mass they contain.

The object budget below is the one this study fixes for a scan of that kind -- a single-lepton
trigger, at most four objects, MET as a category split rather than a mass ingredient. It is a
design choice, not a measurement: the point is to price a scan of this shape against the published
program, and the answer moves only through ln N. Rules documented in docs/METHOD_NOTES.md.

Writes results/tables/combinatorial_budget.csv.
"""
import os, sys, math, itertools, collections, csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import n_s, z_local_for_global5

# ------------------------------------------------------------------ inputs
NMAX   = {"e": 2, "m": 2, "j": 3, "b": 3, "Z": 1}   # per-type multiplicity ceiling
NOBJ   = 4                                          # max objects per category (MET excluded)
KMIN, KMAX = 2, 4                                   # mass-group sizes

SIGMA  = {"e": 0.04, "m": 0.04, "Z": 0.04, "j": 0.10, "b": 0.20}   # fractional pT resolution
# scan window per group size k: [lo, hi] GeV (single-lepton-triggered turn-on -> stats ceiling)
WINDOW = {2: (100.0, 5000.0), 3: (150.0, 4000.0), 4: (200.0, 3000.0)}
MZ     = 91.2

def r_group(comp):
    """Fractional mass resolution of a k-body group from its object composition."""
    s2 = [SIGMA[t] ** 2 for t in comp]
    return 0.5 * math.sqrt(sum(s2) / len(s2))

def window(comp):
    lo, hi = WINDOW[len(comp)]
    lo = max(lo, MZ * comp.count("Z") + 40.0 * (len(comp) - comp.count("Z")))  # Z-mass floor
    return lo, hi

# ------------------------------------------------------------------ enumerate categories
Cat = collections.namedtuple("Cat", "n met os_ss")
cats = []
for ne in range(NMAX["e"] + 1):
    for nm in range(NMAX["m"] + 1):
        for nj in range(NMAX["j"] + 1):
            for nb in range(NMAX["b"] + 1):
                for nz in range(NMAX["Z"] + 1):
                    n = {"e": ne, "m": nm, "j": nj, "b": nb, "Z": nz}
                    tot = sum(n.values())
                    if tot > NOBJ or tot < KMIN:
                        continue
                    if ne + nm == 0 and nz == 0:      # single-lepton trigger
                        continue
                    charge = 2 if (ne + nm) == 2 else 1     # OS/SS split of the dilepton cases
                    for met in (0, 1):
                        cats.append(Cat(n, met, charge))

# ------------------------------------------------------------------ enumerate spectra
rows, N_trials = [], 0.0
n_hist = 0
by_size = collections.Counter()
by_type = collections.Counter()          # distinct object-type multisets (physics observables)
for c in cats:
    objs = [t for t in "emjbZ" for _ in range(c.n[t])]     # indexed objects of the category
    for k in range(KMIN, min(KMAX, len(objs)) + 1):
        for idx in itertools.combinations(range(len(objs)), k):
            comp = tuple(objs[i] for i in idx)
            lo, hi = window(comp)
            r = r_group(comp)
            ns = n_s(lo, hi, r) * c.os_ss                  # OS and SS are two disjoint looks
            N_trials += ns
            n_hist += c.os_ss
            by_size[k] += c.os_ss
            by_type["".join(sorted(comp))] += c.os_ss
            rows.append(("".join(f"{t}{c.n[t]}" for t in "emjbZ") + f"_{c.met}met",
                         "".join(comp), c.os_ss, r, lo, hi, ns))

# Two category counts, and they differ by the OS/SS split of the same-flavour dilepton cases:
# quote them together, since N_trials below counts an OS and an SS look separately.
n_cat_mult = len(cats)                                    # multiplicity x MET categories
n_cat = sum(c.os_ss for c in cats)                        # ... after the OS/SS split
print(f"categories (exclusive, MET-split, >=2 objects, >=1 lepton)           : {n_cat_mult:,}")
print(f"  ... after splitting same-flavour dilepton cases into OS and SS     : {n_cat:,}")
print(f"(category, mass-group) combinations                                  : {len(rows):,}")
print(f"mass spectra (histograms, index-specific groups, k=2..4, OS/SS-split): {n_hist:,}")
print(f"  by group size: " + ", ".join(f"k={k}: {v:,}" for k, v in sorted(by_size.items())))
print(f"distinct object-type multisets (physics observables)                 : {len(by_type)}")
print()
print(f"N_trials = {N_trials:,.0f}    ln N = {math.log(N_trials):.2f}")
print(f"Z_local for 5 sigma global = {z_local_for_global5(N_trials):.2f}")
print(f"  r band (x0.5..x2 -> N x2..x0.5): "
      f"{z_local_for_global5(N_trials*0.5):.2f} - {z_local_for_global5(N_trials*2):.2f}")
print(f"  a local 5 sigma is worth Z_global = "
      f"{math.sqrt(max(25.0 - 2*math.log(N_trials), 0)):.2f} sigma")
print()
print("top object-type multisets by #histograms:")
for t, v in by_type.most_common(12):
    print(f"  {t:6s} {v:5d}")

with open(os.path.join(ROOT, "results", "tables", "combinatorial_budget.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["category", "group", "charge_split", "r", "M_lo_GeV", "M_hi_GeV", "n_s"])
    for row in rows:
        w.writerow([row[0], row[1], row[2], f"{row[3]:.3f}", f"{row[4]:.0f}", f"{row[5]:.0f}",
                    f"{row[6]:.1f}"])
print(f"\nwrote results/tables/combinatorial_budget.csv ({len(rows):,} group rows)")
