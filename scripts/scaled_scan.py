#!/usr/bin/env python3
"""What a scaled-up scan costs, and what fits inside a fixed trials budget.

The five-object scan of combinatorial_budget.py reaches e, mu, light jets, b-jets and a leptonic Z. A
scan built on today's reconstruction would also form masses from hadronic taus, photons and the
boosted large-R candidates for W/Z, H and top. Rules:

  * ten object types, no per-type ceiling, and STRICTLY at most four objects per category;
  * missing energy splits every category and is never an ingredient of a mass, so the transverse-mass
    axes of the model-driven budget are outside this scan's reach;
  * every 2-to-4-object subset of a category is its own spectrum;
  * any trigger: a category needs no lepton;
  * same-flavour dilepton categories still split OS/SS.

That scan is far larger than anyone would run, so the second half of this script imposes a trials
budget of TRIALS_BUDGET and asks which spectra survive it. Priority is, in order:

  1. spectra whose object composition is one of the model-motivated axes of bump_observables (the 46
     of the model-driven budget), highest expected rate first;
  2. everything else, highest expected rate first.

Expected rate is a declared order-of-magnitude weight per object type, multiplied over a category's
content. Only the ordering of those weights matters for the ranking, not their values.

The last part adds selection lenses: an extra event-level requirement on an unchanged mass axis, which
is one more view of the same spectrum and one more look. Half the handles a wide search would use are
already in the enumeration (high MET, object multiplicity, b-tag and tau enrichment) and would be
double counted; the four that are not are priced one at a time, never in combination.

The statistics requirement every histogram has to pass is anchored on a Run-2 dataset, so the report
closes by rescaling that anchor to Run 2 plus Run 3.

Resolutions are not new inputs: each object's fractional sigma is read back out of the published
resolution of its own symmetric channel, sigma = 2 r, as two_body_matrix.py does. Missing energy has
no symmetric channel and is read out of mT(e,v) instead.

Reads  scripts/bump_observables.py (published resolutions), scripts/combinatorial_budget.py (rules).
Writes results/tables/scaled_scan.csv    (one row per scan variant)
       results/tables/priority_scan.csv  (per composition: tier, spectra and looks, kept or dropped)
       results/tables/lens_scan.csv      (per lens: what it adds, and what the budget keeps of it)
Prints the report the Makefile captures as results/tables/scaled_scan.txt.
"""
import os, sys, math, csv, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import res, n_s, scan_segments, z_local_for_global5
import combinatorial_budget as CB
import yield_model as YM

TRIALS_BUDGET = 5.0e5

# Integrated luminosity relative to the Run-2 dataset the yield anchor is set on, ignoring the rise in
# high-mass cross sections from 13 to 13.6 TeV, which acts in the same direction.
DATASETS = [("Run 2, 140 fb-1", 1.0), ("Run 2+3, ~400 fb-1", 3.0)]

# ------------------------------------------------------------------ the wider alphabet
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

# ------------------------------------------------------------------ the model-motivated axes
# Every observable of the model-driven budget, as the object composition(s) a scan would build it
# from. A tuple where the axis spans several compositions, None where this alphabet cannot form it.
MOTIVATED = {
    "m(gammagamma)": ("gg",),      "m(egamma)": ("eg",),     "m(mugamma)": ("gm",),
    "m(jgamma)": ("gj",),          "m(Vgamma)": ("Vg",),
    "m(ee)": ("ee",),              "m(ee) SS": ("ee",),      "m(ee) (Zd)": ("ee",),
    "m(mumu)": ("mm",),            "m(mumu) SS": ("mm",),    "m(mumu) (Zd)": ("mm",),
    "m(emu) LFV": ("em",),         "m(emu) SS": ("em",),
    "m(tautau)": ("TT",),          "m(etau) LFV": ("eT",),   "m(mutau) LFV": ("Tm",),
    "m(ej)": ("ej",),              "m(muj)": ("jm",),        "m(tauj)": ("Tj",),
    "m(eb)": ("be",),              "m(mub)": ("bm",),        "m(taub)": ("Tb",),
    "m(eZ)": ("eZ",),              "m(muZ)": ("mZ",),
    "m(jj)": ("jj",),              "m(bj)": ("bj",),         "m(bb)": ("bb",),
    "m(cb) dijet": ("bj",),        "m(3j)": ("jjj",),
    "m(eejj)": ("eejj",),          "m(mumujj)": ("jjmm",),
    "m(VV)": ("VV",),              "m(Vh)": ("HV",),         "m(HH)": ("HH",),
    "m(tt)": ("tt",),              "m(tt)/m(jj)": ("tt", "jj"),
    "m(tb)": ("bt",),              "m(tW)": ("Vt",),         "m(Wb)": ("Vb",),
    "m(Ht)": ("Ht",),              "m(ttZ)/m(Zt)": ("Vt", "Zt"),
    "multilepton": ("eee", "eem", "emm", "mmm"),
    "m(multi)": None,              # a many-object mass with no fixed composition
    "mT(ev)": None, "mT(muv)": None, "mT(taunu)": None,   # no mass carries missing energy
}
canon = lambda c: "".join(sorted(c))
MOTIVATED_COMPS = {canon(c) for cs in MOTIVATED.values() if cs for c in cs}

print("object alphabet: resolution and yield factor from the published symmetric channel")
for k in ORDER_WIDE:
    ch = YM.SYM[k]
    print(f"  {k}  {WIDE[k][0]:12s} sigma = {SIGMA_WIDE[k]:.3f}  yield {YM.F[k]:8.1e}   "
          f"({ch if ch else 'baseline value, no symmetric channel'})")
print(f"  {MET_KEY}  {WIDE[MET_KEY][0]:12s} {'':14s} yield {YM.F[MET_KEY]:8.1e}   "
      f"(a category split, never a mass)")
print(f"\nmodel-motivated compositions: {len(MOTIVATED_COMPS)} from "
      f"{sum(1 for v in MOTIVATED.values() if v)} of the {len(MOTIVATED)} axes "
      f"(no composition for {', '.join(k for k, v in MOTIVATED.items() if v is None)})")

# The size of one exemption: published windows are never gated, because a published search
# demonstrates its own feasibility. This is what the requirement would do to them if they were.
pub_axes = {a: cs for a, cs in MOTIVATED.items() if cs}
pub_N, pub_N_gated, pub_ok = 0.0, 0.0, 0
for _a, _cs in pub_axes.items():
    _r, _w = res(_a), max(YM.weight(c) for c in _cs)
    _full = sum(n_s(lo, hi, _r) for lo, hi in scan_segments(_a))
    _gated = 0.0
    for lo, hi in scan_segments(_a):
        _hs, _ns, _ev, _fits = YM.gate(lo, hi, _r, _w)
        _gated += _ns if _fits else 0.0
    pub_N += _full
    pub_N_gated += _gated
    pub_ok += _gated > 0
print(f"the requirement is never applied to a published window. On the {len(pub_axes)} published axes "
      f"this alphabet can form it would leave {pub_ok} of them and N = {pub_N_gated:,.0f} of "
      f"{pub_N:,.0f}, so exempting them is the conservative choice")
print()

# ------------------------------------------------------------------ price the variants
BASE = dict(order=CB.ORDER, nmax=CB.NMAX, sigma=CB.SIGMA, mass=CB.MASS, trig=CB.TRIG,
            lepton=CB.LEPTON, nobj=CB.NOBJ, kmax=CB.KMAX)
WIDE_ARGS = dict(order=ORDER_WIDE, nmax=NMAX_WIDE, sigma=SIGMA_WIDE, mass=MASS_WIDE,
                 lepton=LEPTON_WIDE, kmax=KMAX_WIDE, nobj=NOBJ_WIDE, trig="")

VARIANTS = [
    ("five objects, lepton trigger", BASE),
    ("five objects, derived sigma",  {**BASE, "sigma": SIGMA_BASE_DERIVED}),
    ("ten objects, any trigger",     WIDE_ARGS),
]

runs = []
for label, kw in VARIANTS:
    s = CB.enumerate_scan(**kw)
    runs.append((label, kw, s))
    print(f"=== {label}  ({len(kw['order'])} object types, K<={kw['kmax']} objects per mass, "
          f"<={kw['nobj']} per category, "
          f"{'lepton required' if kw.get('trig', CB.TRIG) else 'any trigger'})")
    CB.report(s)
    print()

full = runs[-1][2]

# ------------------------------------------------------------------ fit it into the budget
# Three priority tiers, and inside each one the highest expected rate first:
#   0  every model-motivated axis once, in the best-populated category it appears in, so that no
#      motivated axis can be lost to the budget;
#   1  the same axes in their remaining categories;
#   2  everything else.
# The budget then takes the priority-ordered prefix that fits. It stops at the first spectrum that
# does not, rather than topping up with whatever cheap spectrum happens to fit in the remainder.
# Rows that cannot be fitted score no looks and are not spectra: they never enter the priority order.
# An OS/SS-split row is two histograms and its looks already count both, so the two bases must not be
# mixed: everything below counts histograms.
spectra = [r for r in full.rows if r.n_s > 0]
best = {}
for i, r in enumerate(spectra):
    key = canon(r.group)
    if key in MOTIVATED_COMPS and (key not in best or r.w > spectra[best[key]].w):
        best[key] = i
once = set(best.values())
ranked = sorted(((0 if i in once else 1 if canon(r.group) in MOTIVATED_COMPS else 2,
                  -r.w, r.n_s, canon(r.group), r.cat, r.split) for i, r in enumerate(spectra)),
                key=lambda x: (x[0], x[1], x[2]))

kept_looks, kept_n = collections.Counter(), collections.Counter()
drop_looks, drop_n = collections.Counter(), collections.Counter()
tier_N = collections.Counter()
tier_n = collections.Counter()
N_sel, n_sel, cut_rate, stopped = 0.0, 0, None, False
kept_cats, tier0_cats, tier_axes = set(), set(), collections.Counter()
for tier, negrate, ns, key, cat, nh in ranked:
    tier_N[tier] += ns
    tier_n[tier] += nh
    tier_axes[tier] += 1
    if tier == 0:
        tier0_cats.add(cat)
    if not stopped and N_sel + ns <= TRIALS_BUDGET:
        N_sel += ns
        n_sel += nh
        kept_looks[key] += ns
        kept_n[key] += nh
        kept_cats.add(cat)
        cut_rate = -negrate
    else:
        stopped = True
        drop_looks[key] += ns
        drop_n[key] += nh
tierA_N = tier_N[0] + tier_N[1]
tierA_n = tier_n[0] + tier_n[1]

print(f"=== full ten-object scan")
print(f"spectra {full.n_hist:,}   N = {full.N:,.0f}   Z_local = {z_local_for_global5(full.N):.2f}")
print(f"tier 0, every motivated axis once : {tier_n[0]:6,d} spectra, N = {tier_N[0]:10,.0f} "
      f"({100*tier_N[0]/full.N:4.1f} % of the scan), Z_local = {z_local_for_global5(tier_N[0]):.2f}")
print(f"tier 1, those axes elsewhere      : {tier_n[1]:6,d} spectra, N = {tier_N[1]:10,.0f} "
      f"({100*tier_N[1]/full.N:4.1f} %)")
print(f"tier 2, the rest                  : {tier_n[2]:6,d} spectra, N = {tier_N[2]:10,.0f} "
      f"({100*tier_N[2]/full.N:4.1f} %)")
print(f"tiers 0+1 together: N = {tierA_N:,.0f}, which is "
      f"{tierA_N/TRIALS_BUDGET:.2f} times the budget on its own")
NEW_TYPES = "TgtVH"                                  # what the wider alphabet adds to the five
_nn = sum(v for c, v in full.by_type.items() if any(k in c for k in NEW_TYPES))
_nl = sum(v for c, v in full.looks.items() if any(k in c for k in NEW_TYPES))
print(f"spectra reaching a new object type: {_nn:,} of {full.n_hist:,} "
      f"({100*_nn/full.n_hist:.0f} %), carrying {100*_nl/full.N:.0f} % of N; "
      f"two-body groups are {100*full.by_size[2]/full.n_hist:.0f} % of the spectra")
_lost = sorted(MOTIVATED_COMPS - set(full.by_type))
print(f"motivated compositions with no fittable histogram in any category: "
      f"{', '.join(_lost) if _lost else 'none'}")
print()
print(f"=== priority prefix that fits N <= {TRIALS_BUDGET:,.0f}")
print(f"selected {n_sel:,} of {full.n_hist:,} spectra ({100*n_sel/full.n_hist:.1f} %) over "
      f"{len(kept_n)} of {len(full.by_type)} compositions and {len(kept_cats):,} of "
      f"{full.n_cat:,} categories")
print(f"tier 0 is {tier_axes[0]} axes in {len(tier0_cats)} categories, {tier_n[0]} histograms")
if stopped:
    print(f"of the {tier_n[1]:,} tier-1 spectra, {n_sel - tier_n[0]:,} fit "
          f"({100*(n_sel - tier_n[0])/tier_n[1]:.0f} %)")
    print(f"the cut lands at a category yield of {cut_rate:.1e}: thinner categories are dropped")
else:
    print(f"nothing is cut: every fittable spectrum fits inside the budget, so the statistics "
          f"requirement binds first and the priority order never has to be applied")
print(f"N = {N_sel:,.0f} ({100*N_sel/full.N:.1f} % of the full scan), "
      f"Z_local = {z_local_for_global5(N_sel):.2f} "
      f"(band {z_local_for_global5(N_sel*0.5):.2f}-{z_local_for_global5(N_sel*2):.2f})")
print(f"a local 5 sigma is then worth Z_global = "
      f"{math.sqrt(max(25.0 - 2*math.log(N_sel), 0)):.2f} sigma")
print()

# How hard the requirement bites depends on the yield anchor, which is the least certain input: a
# factorised per-object model prices every object at the full cost of its own production, and real
# objects arrive in pairs from one boson, so it under-counts high-multiplicity categories. Two orders
# of magnitude either way is the honest band.
print("yield anchor scaled (the model's own uncertainty): spectra, N, Z_local")
_n_ref = YM.N_REF
for _s in (0.01, 1.0, 100.0):
    YM.N_REF = _n_ref * _s
    _r = CB.enumerate_scan(**WIDE_ARGS, collect=False)
    print(f"  x{_s:<6g} {_r.n_hist:7,d} of {_r.n_hist + _r.n_thin:7,d} histograms, "
          f"N = {_r.N:11,.0f}, Z_local = {z_local_for_global5(_r.N):.2f}")
YM.N_REF = _n_ref
print()

# which object types survive, and through which compositions
print("object type    spectra in scan  selected   kept looks   what survives")
for k in ORDER_WIDE:
    tot = sum(v for c, v in full.by_type.items() if k in c)
    sel = sum(v for c, v in kept_n.items() if k in c)
    lk = sum(v for c, v in kept_looks.items() if k in c)
    comps = sorted((c for c in kept_n if k in c), key=lambda c: -kept_looks[c])
    extra = [c for c in comps if c not in MOTIVATED_COMPS]
    if not comps:
        what = "DROPPED ENTIRELY"
    else:
        what = f"{len(comps)} compositions: " + " ".join(comps[:6]) + \
               (" ..." if len(comps) > 6 else "") + \
               (f"  (+{len(extra)} unmotivated)" if extra else "  (all motivated)")
    print(f"  {k} {WIDE[k][0]:12s} {tot:9,d} {sel:9,d} {lk:12,.0f}   {what}")
print()

print("compositions dropped entirely, costliest first:")
gone = [(c, drop_looks[c], drop_n[c]) for c in drop_looks if c not in kept_n]
for c, lk, n in sorted(gone, key=lambda x: -x[1])[:15]:
    print(f"  {c:6s} {n:6,d} spectra, {lk:10,.0f} looks")
print(f"  ... {len(gone)} of {len(full.by_type)} compositions dropped entirely")
print()
print("costliest compositions that are kept:")
for c in sorted(kept_looks, key=lambda c: -kept_looks[c])[:15]:
    tag = "motivated" if c in MOTIVATED_COMPS else "by rate"
    print(f"  {c:6s} {kept_n[c]:6,d} of {full.by_type[c]:6,d} spectra, "
          f"{kept_looks[c]:9,.0f} looks   ({tag})")
print()

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


items, tier_of = [], {}
lens_n, lens_N = collections.Counter(), collections.Counter()
lens_thin = collections.Counter()
for i, r in enumerate(spectra):
    key = canon(r.group)
    tier = 0 if i in once else 1 if key in MOTIVATED_COMPS else 2
    tier_of[i] = (tier, key)
    items.append((tier, -r.w, 0, r.n_s, key, r.cat, r.split, ""))
for i, lk, li, lns, fits in lens_views(spectra):
    r, (tier, key) = spectra[i], tier_of[i]
    if not fits:
        lens_thin[lk] += r.split
        continue
    items.append((max(tier, 1), -r.w, li, lns, key, r.cat, r.split, lk))
    lens_n[lk] += r.split
    lens_N[lk] += lns
items.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
lensed_n, lensed_N = full.n_hist + sum(lens_n.values()), full.N + sum(lens_N.values())

L_N, L_views, L_spec, L_stop = 0.0, 0, 0, False
L_comps, L_cats = set(), set()
L_lens_n, L_lens_N = collections.Counter(), collections.Counter()
for tier, negrate, li, ns, key, cat, nh, lk in items:
    if L_stop or L_N + ns > TRIALS_BUDGET:
        L_stop = True
        continue
    L_N += ns
    L_views += nh
    L_comps.add(key)
    L_cats.add(cat)
    if li == 0:
        L_spec += nh
    else:
        L_lens_n[lk] += nh
        L_lens_N[lk] += ns

print("selection lenses on the same axes, one at a time:")
for key, label, rule, _ok, cap in LENSES:
    print(f"  {label:20s} {rule:34s} eff {YM.LENS_EFF[key]:5.3f} {lens_n[key]:7,d} views "
          f"{lens_N[key]:9,.0f} looks, {lens_thin[key]:6,d} too thin"
          f"{'' if cap is None else f' (window capped at {cap:.0f} GeV)'}")
print(f"  {'all four':20s} {'':34s} {'':9s} {sum(lens_n.values()):7,d} views "
      f"{sum(lens_N.values()):9,.0f} looks, {sum(lens_thin.values()):6,d} too thin")
print(f"the ten-object scan with lenses: {lensed_n:,} histograms "
      f"({lensed_n/full.n_hist:.1f} per spectrum), N = {lensed_N:,.0f}, "
      f"Z_local = {z_local_for_global5(lensed_N):.2f}")
print()
print(f"=== priority prefix with lenses that fits N <= {TRIALS_BUDGET:,.0f}")
print(f"selected {L_spec:,} of {full.n_hist:,} spectra ({100*L_spec/full.n_hist:.1f} %) plus "
      f"{L_views - L_spec:,} lens views of them: {L_views:,} histograms "
      f"({100*L_views/lensed_n:.1f} % of the lensed scan) over {len(L_comps)} compositions and "
      f"{len(L_cats):,} of {full.n_cat:,} categories")
print(f"N = {L_N:,.0f}, Z_local = {z_local_for_global5(L_N):.2f}")
for key, label, _rule, _ok, _cap in LENSES:
    frac = 100.0 * L_lens_n[key] / lens_n[key] if lens_n[key] else 0.0
    print(f"  {label:20s} {L_lens_n[key]:7,d} of {lens_n[key]:7,d} views kept "
          f"({frac:4.1f} %), {L_lens_N[key]:9,.0f} looks")
if L_spec < n_sel:
    print(f"the same budget with no lens at all reaches {n_sel:,} spectra, so the lenses are paid for "
          f"in coverage: {n_sel - L_spec:,} fewer axes-in-categories for {L_views - L_spec:,} lens "
          f"views of the ones that remain")
else:
    print(f"the lenses cost nothing in coverage: axes and lens views together still fit the budget, "
          f"and {sum(lens_thin.values()):,} further views are ruled out by statistics alone")
print()

# The yield anchor is a Run-2 dataset, so a larger one enters exactly as the band above does. It buys
# little reach per spectrum, the one-event mass going as luminosity^(1/(P-1)).
print(f"dataset scaled (the anchor is {DATASETS[0][0]}): spectra, N, Z_local, mass reach")
for _label, _s in DATASETS:
    YM.N_REF = _n_ref * _s
    _r = CB.enumerate_scan(**WIDE_ARGS)
    _sp = [x for x in _r.rows if x.n_s > 0]
    _lv = [(_sp[i].split, lns) for i, _lk, _li, lns, fits in lens_views(_sp) if fits]
    print(f"  {_label:20s} x{_s:<4g} {_r.n_hist:6,d} of {_r.n_hist + _r.n_thin:6,d} histograms, "
          f"N = {_r.N:9,.0f}, Z_local = {z_local_for_global5(_r.N):.2f}, "
          f"x{_s ** (1.0 / (YM.P - 1.0)):.2f} in mass; with lenses "
          f"{_r.n_hist + sum(n for n, _ in _lv):6,d} histograms, "
          f"N = {_r.N + sum(x for _, x in _lv):9,.0f}, "
          f"Z_local = {z_local_for_global5(_r.N + sum(x for _, x in _lv)):.2f}")
YM.N_REF = _n_ref
print()

# ------------------------------------------------------------------ tables
hdr = ("scan", "types", "K_max", "objects_per_category", "trigger", "categories",
       "spectra", "compositions", "N_trials", "Z_local", "Z_lo", "Z_hi")
row = lambda label, kw, cats, nsp, ncomp, N: [
    label, len(kw["order"]), kw["kmax"], kw["nobj"],
    "lepton" if kw.get("trig", CB.TRIG) else "any", cats, nsp, ncomp, f"{N:.0f}",
    f"{z_local_for_global5(N):.2f}", f"{z_local_for_global5(N*0.5):.2f}",
    f"{z_local_for_global5(N*2):.2f}"]

with open(os.path.join(ROOT, "results", "tables", "scaled_scan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(hdr)
    for label, kw, s in runs:
        w.writerow(row(label, kw, s.n_cat, s.n_hist, len(s.by_type), s.N))
    w.writerow(row("... model-motivated axes once each", VARIANTS[-1][1], len(tier0_cats),
                   tier_n[0], len(MOTIVATED_COMPS & set(full.by_type)), tier_N[0]))
    w.writerow(row(f"... prioritised to N <= {TRIALS_BUDGET:.0e}", VARIANTS[-1][1], len(kept_cats),
                   n_sel, len(kept_n), N_sel))
    w.writerow(row("ten objects, with selection lenses", VARIANTS[-1][1], full.n_cat, lensed_n,
                   len(full.by_type), lensed_N))
    w.writerow(row(f"... prioritised to N <= {TRIALS_BUDGET:.0e}, lenses included",
                   VARIANTS[-1][1], len(L_cats), L_views, len(L_comps), L_N))

with open(os.path.join(ROOT, "results", "tables", "priority_scan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["composition", "K", "tier", "spectra_total", "spectra_kept", "looks_total",
                "looks_kept"])
    for c in sorted(full.by_type, key=lambda c: (len(c), -full.looks[c])):
        w.writerow([c, len(c), "motivated" if c in MOTIVATED_COMPS else "rate",
                    full.by_type[c], kept_n.get(c, 0), f"{full.looks[c]:.0f}",
                    f"{kept_looks.get(c, 0.0):.0f}"])

with open(os.path.join(ROOT, "results", "tables", "lens_scan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["lens", "requirement", "window_cap_GeV", "views_added", "looks_added", "views_kept",
                "looks_kept"])
    for key, label, rule, _ok, cap in LENSES:
        w.writerow([label, rule, "" if cap is None else f"{cap:.0f}", lens_n[key],
                    f"{lens_N[key]:.0f}", L_lens_n[key], f"{L_lens_N[key]:.0f}"])

print(f"\nwrote results/tables/scaled_scan.csv ({len(runs) + 4} rows), priority_scan.csv "
      f"({len(full.by_type)} compositions) and lens_scan.csv ({len(LENSES)} lenses)",
      file=sys.stderr)
