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
NMAX_WIDE  = {"e": 2, "m": 2, "T": 2, "g": 2, "j": 3, "b": 3, "t": 2, "V": 2, "H": 2, "Z": 1}
SIGMA_WIDE = {k: (CB.SIGMA["Z"] if ch is None else 2.0 * res(ch)) for k, (_, ch, _) in WIDE.items()}
MASS_WIDE  = {k: m for k, (_, _, m) in WIDE.items() if m is not None}
LEPTON_WIDE = "emT"                       # hadronic taus are charged: they split OS/SS too
TRIG_WIDE   = "emZ"                       # a hadronic tau is not a lepton trigger

# Same rules, five object types, resolutions from the same inversion: isolates the alphabet from
# the resolution prescription, since the baseline declares sigma per object instead of deriving it.
SIGMA_BASE_DERIVED = {k: SIGMA_WIDE[k] for k in CB.ORDER}

VARIANTS = [
    ("baseline, declared resolutions", CB.ORDER, CB.NMAX, CB.SIGMA, CB.MASS, CB.TRIG, CB.LEPTON),
    ("baseline, derived resolutions",  CB.ORDER, CB.NMAX, SIGMA_BASE_DERIVED, CB.MASS,
                                       CB.TRIG, CB.LEPTON),
    ("wider alphabet, lepton trigger", ORDER_WIDE, NMAX_WIDE, SIGMA_WIDE, MASS_WIDE,
                                       TRIG_WIDE, LEPTON_WIDE),
    ("wider alphabet, any trigger",    ORDER_WIDE, NMAX_WIDE, SIGMA_WIDE, MASS_WIDE,
                                       "", LEPTON_WIDE),
]

print("object resolutions of the wider alphabet (sigma = 2 r of the symmetric channel)")
for k in ORDER_WIDE:
    label, ch, _ = WIDE[k]
    src = "baseline value (no symmetric published channel)" if ch is None else f"from {ch}"
    print(f"  {k}  {label:12s} sigma = {SIGMA_WIDE[k]:.3f}   {src}")
print()

# ------------------------------------------------------------------ price each variant
runs = []
for label, order, nmax, sigma, mass, trig, lepton in VARIANTS:
    s = CB.enumerate_scan(order=order, nmax=nmax, sigma=sigma, mass=mass, trig=trig, lepton=lepton)
    runs.append((label, len(order), trig, s))
    print(f"=== {label}  ({len(order)} object types, "
          f"{'lepton required' if trig else 'no trigger requirement'})")
    CB.report(s)
    print()

hdr = ("scan", "objects", "trigger", "categories", "spectra", "compositions", "N_trials",
       "Z_local", "Z_lo", "Z_hi")
print("%-32s %7s %8s %11s %8s %13s %10s %8s" % hdr[:8])
for label, nobj, trig, s in runs:
    print("%-32s %7d %8s %11d %8d %13d %10.4g %8.2f"
          % (label, nobj, "lepton" if trig else "any", s.n_cat, s.n_hist, len(s.by_type),
             s.N, z_local_for_global5(s.N)))

base = runs[0][3].N
print(f"\nthe alphabet alone multiplies the trials by {runs[2][3].N / runs[1][3].N:.1f}, "
      f"and the headline scan is {runs[3][3].N / base:.1f} times the baseline: "
      f"{z_local_for_global5(runs[3][3].N) - z_local_for_global5(base):+.2f} sigma on the bar")

with open(os.path.join(ROOT, "results", "tables", "scaled_scan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(hdr)
    for label, nobj, trig, s in runs:
        w.writerow([label, nobj, "lepton" if trig else "any", s.n_cat, s.n_hist, len(s.by_type),
                    f"{s.N:.0f}", f"{z_local_for_global5(s.N):.2f}",
                    f"{z_local_for_global5(s.N*0.5):.2f}", f"{z_local_for_global5(s.N*2):.2f}"])

# ------------------------------------------------------------------ where the headline cost sits
head = runs[-1][3]
looks = collections.defaultdict(float)
for r in head.rows:
    looks["".join(sorted(r[1]))] += r[6]
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
