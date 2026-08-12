#!/usr/bin/env python3
"""What a scaled-up scan costs: the combinatorial scan of combinatorial_budget.py over a wider
object alphabet.

The baseline scan reaches five object types (e, mu, light jet, b-jet, leptonic Z). A general search
built on today's reconstruction would also form masses from hadronic taus, photons and boosted
large-R W/Z, H and top candidates, which is ten types. Same rules as the baseline -- exclusive
multiplicity categories, at most four objects, every 2-to-4-object subset its own spectrum, MET a
category split rather than a mass ingredient -- so the alphabet is the only thing that changes.

Resolutions are not new inputs: each object's fractional sigma is read back out of the published
resolution of its own symmetric channel, sigma = 2 r, the same inversion two_body_matrix.py uses to
price its gaps. The leptonic Z has no symmetric published channel and keeps the baseline value.

Reads  scripts/bump_observables.py (published resolutions), scripts/combinatorial_budget.py (rules).
Writes results/tables/scaled_scan.csv          (one row per scan variant, the paper's table)
       results/tables/scaled_scan_groups.csv   (looks per object composition, headline variant)
Prints the report the Makefile captures as results/tables/scaled_scan.txt.
"""
import os, sys, math, csv, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import res, z_local_for_global5
import combinatorial_budget as CB

# ------------------------------------------------------------------ the wider alphabet
# key -> (label, symmetric channel the resolution comes from, mass floor contribution in GeV)
WIDE = {
    "e": ("e",              "m(ee)",             None),
    "m": ("mu",             "m(mumu)",           None),
    "T": ("tau_had",        "m(tautau)",         None),
    "g": ("gamma",          "m(gammagamma)",     None),
    "j": ("light jet",      "m(jj)",             None),
    "b": ("b-jet",          "m(bb)",             None),
    "t": ("boosted top",    "m(tt)",             173.0),
    "V": ("boosted W/Z",    "m(VV)",              91.2),
    "H": ("boosted H",      "m(HH)",             125.0),
    "Z": ("leptonic Z",     None,                 91.2),
}
ORDER_WIDE = "emTgjbtVHZ"
# The four bounds the mass combination, not the event: a category may hold more objects than any one
# mass is built from. No per-type ceiling either, unlike the five-object scan, whose ceilings describe
# one particular grid. A leptonic Z counts as one object here as it does there.
KMAX_WIDE  = 4                                   # objects per mass combination
NOBJ_WIDE  = (4, 6, 8)                           # objects per category, swept
NMAX_WIDE  = {k: max(NOBJ_WIDE) for k in ORDER_WIDE}
SIGMA_WIDE = {k: (CB.SIGMA["Z"] if ch is None else 2.0 * res(ch)) for k, (_, ch, _) in WIDE.items()}
MASS_WIDE  = {k: m for k, (_, _, m) in WIDE.items() if m is not None}
LEPTON_WIDE = "emT"                       # hadronic taus are charged: they split OS/SS too
TRIG_WIDE   = "emZ"                       # a hadronic tau is not a lepton trigger

# Same rules, five object types, resolutions from the same inversion: isolates the alphabet from
# the resolution prescription, since the baseline declares sigma per object instead of deriving it.
SIGMA_BASE_DERIVED = {k: SIGMA_WIDE[k] for k in CB.ORDER}

BASE = dict(order=CB.ORDER, nmax=CB.NMAX, sigma=CB.SIGMA, mass=CB.MASS, trig=CB.TRIG,
            lepton=CB.LEPTON, nobj=CB.NOBJ, kmax=CB.KMAX)
WIDE_ARGS = dict(order=ORDER_WIDE, nmax=NMAX_WIDE, sigma=SIGMA_WIDE, mass=MASS_WIDE,
                 lepton=LEPTON_WIDE, kmax=KMAX_WIDE)

VARIANTS = [
    ("five objects, declared resolutions", BASE),
    ("five objects, derived resolutions",  {**BASE, "sigma": SIGMA_BASE_DERIVED}),
    ("ten objects, lepton trigger",        {**WIDE_ARGS, "trig": TRIG_WIDE, "nobj": 4}),
]
# the sweep: the mass combination stays at four objects, the category ceiling rises
VARIANTS += [(f"ten objects, <={n} per category", {**WIDE_ARGS, "trig": "", "nobj": n,
                                                   "collect": n <= 6}) for n in NOBJ_WIDE]

print("object resolutions of the wider alphabet (sigma = 2 r of the symmetric channel)")
for k in ORDER_WIDE:
    label, ch, _ = WIDE[k]
    src = "baseline value (no symmetric published channel)" if ch is None else f"from {ch}"
    print(f"  {k}  {label:12s} sigma = {SIGMA_WIDE[k]:.3f}   {src}")
print()

# ------------------------------------------------------------------ price each variant
runs = []
for label, kw in VARIANTS:
    s = CB.enumerate_scan(**kw)
    runs.append((label, kw, s))
    print(f"=== {label}  ({len(kw['order'])} object types, K<={kw['kmax']} per mass, "
          f"<={kw['nobj']} per category, "
          f"{'lepton required' if kw.get('trig', CB.TRIG) else 'no trigger requirement'})")
    CB.report(s)
    print()

hdr = ("scan", "types", "K_max", "objects_per_category", "trigger", "categories", "spectra",
       "compositions", "N_trials", "Z_local", "Z_lo", "Z_hi")
row = lambda label, kw, s: [label, len(kw["order"]), kw["kmax"], kw["nobj"],
                            "lepton" if kw.get("trig", CB.TRIG) else "any", s.n_cat, s.n_hist,
                            len(s.by_type), f"{s.N:.0f}", f"{z_local_for_global5(s.N):.2f}",
                            f"{z_local_for_global5(s.N*0.5):.2f}",
                            f"{z_local_for_global5(s.N*2):.2f}"]

print("%-36s %5s %5s %5s %7s %10s %8s %6s %11s %7s" % hdr[:10])
for label, kw, s in runs:
    r = row(label, kw, s)
    print("%-36s %5s %5s %5s %7s %10s %8s %6s %11s %7s" % tuple(r[:10]))

base, wide4 = runs[0][2].N, runs[3][2].N
print(f"\nat the same ceiling, trigger and resolution prescription, the alphabet multiplies the "
      f"trials by {runs[2][2].N / runs[1][2].N:.1f}")
print(f"dropping the lepton requirement at that ceiling: "
      f"{z_local_for_global5(wide4) - z_local_for_global5(runs[2][2].N):+.2f} sigma")
print(f"raising the category ceiling from 4 to {NOBJ_WIDE[-1]} multiplies them by a further "
      f"{runs[-1][2].N / wide4:.1f}, i.e. "
      f"{z_local_for_global5(runs[-1][2].N) - z_local_for_global5(wide4):+.2f} sigma")
print(f"widest variant against the five-object scan: {runs[-1][2].N / base:.1f} times the trials, "
      f"{z_local_for_global5(runs[-1][2].N) - z_local_for_global5(base):+.2f} sigma on the bar")

with open(os.path.join(ROOT, "results", "tables", "scaled_scan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(hdr)
    for label, kw, s in runs:
        w.writerow(row(label, kw, s))

# ------------------------------------------------------------------ where the headline cost sits
# the four-objects-per-category variant, which is the one comparable with the five-object grid
head = runs[3][2]
looks = head.looks
with open(os.path.join(ROOT, "results", "tables", "scaled_scan_groups.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["composition", "K", "spectra", "looks"])
    for comp in sorted(looks, key=lambda c: (len(c), -looks[c])):
        w.writerow([comp, len(comp), head.by_type[comp], f"{looks[comp]:.0f}"])

byK = collections.Counter()
for comp, v in looks.items():
    byK[len(comp)] += v
NEW = set("TgtVH")
new_looks = sum(v for comp, v in looks.items() if NEW & set(comp))
new_spec = sum(head.by_type[c] for c in looks if NEW & set(c))
print(f"\nspectra reaching a new object type: {new_spec} of {head.n_hist} "
      f"({100*new_spec/head.n_hist:.0f} %), carrying {100*new_looks/head.N:.0f} % of N")
print("\nheadline scan, by group size:")
for k in sorted(byK):
    print(f"  K={k}: {head.by_size[k]:6d} spectra, {byK[k]:9.0f} looks "
          f"({100*byK[k]/head.N:4.1f} % of N)")
print("\ncostliest compositions of the headline scan:")
for comp in sorted(looks, key=lambda c: -looks[c])[:12]:
    print(f"  {comp:6s} {head.by_type[comp]:5d} spectra {looks[comp]:9.0f} looks")

print(f"\nwrote results/tables/scaled_scan.csv ({len(runs)} variants) and "
      f"scaled_scan_groups.csv ({len(looks)} compositions)", file=sys.stderr)
