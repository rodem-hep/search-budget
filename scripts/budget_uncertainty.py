#!/usr/bin/env python3
"""How much the trials count and the discovery bar are uncertain by.

Every number in this repository rests on a declared input: a fractional resolution, a scan window, a
yield anchor, a fittability threshold, or the convention that turns a resolution element into an
independent look. This script varies each of them over the range it is honestly known to and prices
the result on both counted bases -- the 46 published-window spectra of the model space and the
fully combinatorial scan -- so that no headline is quoted as a bare number.

Three things it does that a single band cannot:

  * it VARIES THE INPUT AND RECOMPUTES rather than scaling N. Scaling a resolution changes which
    hypothetical histograms can be fitted as well as how many looks each carries, and the two
    effects partly cancel;
  * it separates the CORRELATED from the uncorrelated part. Almost every input shifts both bases the
    same way, so the difference between two rows of the headline table is far better known than
    either row. The last column of the table below is that difference;
  * it prices WHAT COUNTS AS ONE LOOK, which is the one input with no measurement behind it at all.
    Counting elements of width sigma_M and calling each an independent look is a convention. Adjacent
    elements are correlated, since a resonance spans more than one, which argues for fewer looks;
    the up-crossing form of the Gross-Vitells estimate argues for more. For a unit-variance Gaussian
    process in x = ln M with correlation length r, Rice's formula gives an expected number of
    up-crossings of level Z of (1/2pi)(1/r) ln(M_hi/M_lo) exp(-Z^2/2), which is the element count N
    times Z/sqrt(2pi) ~ 2.6 at the Z of interest, once compared against the Gaussian tail the
    element count is multiplied by here. The band taken below, 0.5 to Z/sqrt(2pi), brackets both
    readings, and this term is the largest single uncertainty on the bar.

Conventions are NOT uncertainties and are kept out of the total: which granularity a spectrum is
counted at, which dataset the yields are priced on, and the design choices that fix the shape of the
hypothetical scan (ten object types, at most four of them, one lens at a time). They are reported
separately, as alternatives.

Reads  scripts/bump_observables.py, scripts/public_obs_map.py (the model space and its windows),
       scripts/scan_alphabet.py, scripts/combinatorial_budget.py, scripts/yield_model.py (the scan),
       results/tables/two_body_matrix.csv (the price of the two-body pairs nobody has scanned).
Writes results/tables/budget_uncertainty.csv  (one row per basis x source x direction)
       results/overviews/BUDGET_UNCERTAINTY.md
       results/tex/uncertainty_table.tex      (the appendix table, a bare tabular)
Pure standard library.
"""
import collections, contextlib, csv, math, os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import (canon, res, n_s, scan_segments, CANON_ORDER,
                             z_local_for_global5 as z5)
from public_obs_map import PUBLIC_OBS, nsel, nonpeak_only
import combinatorial_budget as CB
import scan_alphabet as SA
import yield_model as YM

def _p(*a): return os.path.join(ROOT, *a)

# ---------------------------------------------------------------- the ranges, and where each comes from
R_FACTOR     = 2.0            # r is known to a factor of a few; two is the band every headline carries
R_SCATTER    = 2.0            # the same factor, as an INDEPENDENT per-channel 1 sigma
WINDOW_F     = 1.4            # published-window edges: how far the curated family edge could move
GEN_WINDOW_F = 1.25           # the generic window a hypothetical spectrum is given
ANCHOR_F     = 100.0          # the yield anchor of a factorised per-object rate model
P_BAND       = (6.0, 8.0)     # the falling-background exponent
EV_BAND      = (30.0, 300.0)  # events required in a fittable histogram (nominal 100)
BIN_BAND     = (15, 50)       # resolution elements required of one event or more (nominal 25)
LENS_EFF_F   = 10.0           # the declared order-of-magnitude lens efficiencies
LOOK_LO      = 0.5            # a resonance spans about two elements: only every second is a new look
MC_DRAWS     = 20000
SEED         = 20260817

# The muon axes replace a resolution that rises with p_T by a window-averaged constant, which is the
# largest modelling approximation in the resolution set. The rising form is anchored on the two values
# bump_observables quotes for the sagitta-limited di-muon measurement, and the upper anchor is the one
# that is uncertain.
MU_AXES    = ("m(mumu)", "m(mumu) SS")
MU_ANCHORS = ((200.0, 0.02), (3000.0, 0.15))
MU_HI_BAND = (0.10, 0.20)

# The dataset every headline is priced on: the last entry of scan_alphabet.DATASETS.
LUMI_LABEL, LUMI = SA.DATASETS[-1]

# ================================================================ basis 1: the model space
pub_models = collections.defaultdict(set)
for _model, _obss in PUBLIC_OBS.items():
    for _o in _obss:
        pub_models[canon(_o)].add(_model)
AXES = [o for o in CANON_ORDER if o == canon(o) and o in pub_models]
AXES += [o for o in sorted(pub_models) if o not in AXES]


def N_model(rfun=res, rscale=1.0, lo_f=1.0, hi_f=1.0, axes=None, extra=0.0):
    """Trials over the published windows: 1/r ln(hi/lo) per scanned segment."""
    tot = extra
    for o in (AXES if axes is None else axes):
        r = rfun(o) * rscale
        for lo, hi in scan_segments(o):
            tot += n_s(lo * lo_f, hi * hi_f, r)
    return tot


def ns_rising(lo, hi, anchors, steps=4000):
    """Looks over [lo, hi] for a resolution rising linearly in M between two anchors."""
    (m0, r0), (m1, r1) = anchors
    slope = (r1 - r0) / (m1 - m0)
    x0, x1 = math.log(lo), math.log(hi)
    h = (x1 - x0) / steps
    tot = 0.0
    for i in range(steps + 1):
        r = max(r0 + slope * (math.exp(x0 + i * h) - m0), 1e-4)
        tot += (0.5 if i in (0, steps) else 1.0) * h / r
    return tot


def N_model_rising(r_hi):
    """The muon axes counted with r(M) rising to r_hi at the upper anchor, the rest unchanged."""
    anchors = (MU_ANCHORS[0], (MU_ANCHORS[1][0], r_hi))
    tot = N_model(axes=[o for o in AXES if o not in MU_AXES])
    for o in MU_AXES:
        tot += sum(ns_rising(lo, hi, anchors) for lo, hi in scan_segments(o))
    return tot


def r_effective(obs, r_hi):
    """The flat r that reproduces the rising-resolution trials integral over this axis's window."""
    anchors = (MU_ANCHORS[0], (MU_ANCHORS[1][0], r_hi))
    seg = scan_segments(obs)
    ns = sum(ns_rising(lo, hi, anchors) for lo, hi in seg)
    return math.log(seg[-1][1] / seg[0][0]) / ns


def scattered(seed=SEED, draws=MC_DRAWS, factor=R_SCATTER):
    """Percentiles of N when every channel's r is drawn independently, a factor `factor` per sigma."""
    rng = random.Random(seed)
    sig = math.log(factor)
    seg = {o: [(lo, hi) for lo, hi in scan_segments(o)] for o in AXES}
    out = []
    for _ in range(draws):
        tot = 0.0
        for o in AXES:
            r = res(o) * math.exp(rng.gauss(0.0, sig))
            tot += sum(n_s(lo, hi, r) for lo, hi in seg[o])
        out.append(tot)
    out.sort()
    return [out[int(q * (draws - 1))] for q in (0.16, 0.50, 0.84)]


def overlap_looks():
    """Looks the dark-photon axes double-count against the high-mass dilepton axes.

    The (Zd) windows run to 400 GeV and the Drell-Yan windows start at 150, so 150-400 GeV is
    scanned twice. Removing it from the finer-resolution axis is the de-duplicated count.
    """
    tot = 0.0
    for zd, dy in (("m(ee) (Zd)", "m(ee)"), ("m(mumu) (Zd)", "m(mumu)")):
        lo = max(scan_segments(zd)[0][0], scan_segments(dy)[0][0])
        hi = min(scan_segments(zd)[-1][1], scan_segments(dy)[-1][1])
        tot += n_s(lo, hi, res(zd))
    return tot


def two_body_gap():
    """Looks the object pairs with no published axis would cost, from two_body_matrix.py."""
    with open(_p("results", "tables", "two_body_matrix.csv")) as f:
        rows = [r for r in csv.DictReader(f) if r["status"].startswith("gap")]
    return len(rows), sum(float(r["ns"]) for r in rows)


N_MODEL = N_model()
N_SEL = sum(nsel(o) * sum(n_s(lo, hi, res(o)) for lo, hi in scan_segments(o)) for o in AXES)

# ================================================================ basis 2: the combinatorial scan
@contextlib.contextmanager
def yields(**over):
    """yield_model globals (and combinatorial_budget's generic windows) held at other values."""
    keep = {k: getattr(YM, k) for k in over if k != "WINDOW"}
    win = CB.WINDOW
    try:
        for k, v in over.items():
            if k != "WINDOW":
                setattr(YM, k, v)
        if "WINDOW" in over:
            CB.WINDOW = over["WINDOW"]
        yield
    finally:
        for k, v in keep.items():
            setattr(YM, k, v)
        CB.WINDOW = win


def scan(sigma_scale=1.0, r_group=None, window_f=None, lumi=LUMI, lens_f=1.0, **over):
    """(spectra, N, lensed histograms, lensed N) of the ten-object scan under one variation."""
    args = dict(SA.WIDE_ARGS)
    if sigma_scale != 1.0:
        args["sigma"] = {k: v * sigma_scale for k, v in SA.SIGMA_WIDE.items()}
    if window_f is not None:
        lo_f, hi_f = window_f
        over["WINDOW"] = {k: (lo * lo_f, hi * hi_f) for k, (lo, hi) in CB.WINDOW.items()}
    over["N_REF"] = over.get("N_REF", YM.N_REF) * lumi
    if lens_f != 1.0:
        over["LENS_EFF"] = {k: v * lens_f for k, v in YM.LENS_EFF.items()}
    keep_r = CB.r_group
    try:
        if r_group is not None:
            CB.r_group = r_group
        with yields(**over):
            s = CB.enumerate_scan(**args)
            fit = [r for r in s.rows if r.n_s > 0]
            n_l, N_l = s.n_hist, s.N
            for i, _lk, _li, lns, ok in SA.lens_views(fit):
                if ok:
                    n_l += fit[i].split
                    N_l += lns
    finally:
        CB.r_group = keep_r
    return s.n_hist, s.N, n_l, N_l


def r_worst(comp, sigma):
    """Model-space convention: a multi-object mass takes its worst leg."""
    return 0.5 * max(sigma[t] for t in comp)


def r_quadrature(comp, sigma):
    """Textbook propagation for a two-body mass, a factor sqrt(k) above the calibrated prescription."""
    return 0.5 * math.sqrt(sum(sigma[t] ** 2 for t in comp))


N_SPEC, N_SCAN, N_LSPEC, N_LSCAN = scan()

# ================================================================ what one look is
def z_exact(N, zglob=5.0):
    """Local Z solving N (1-Phi(Z)) = 1-Phi(zglob) exactly, instead of the closed-form LEE relation."""
    target = 0.5 * math.erfc(zglob / math.sqrt(2.0)) / N
    lo, hi = 1.0, 25.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 0.5 * math.erfc(mid / math.sqrt(2.0)) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def c_upcrossing(N):
    """Rice/Gross-Vitells up-crossings against element counting: the factor Z/sqrt(2 pi)."""
    return z5(N) / math.sqrt(2.0 * math.pi)

# ================================================================ the uncertainty budget
# One entry per source. `dirs` holds the variations of that source, aligned across the two bases so
# that the difference between them is taken direction by direction: the correlated part cancels there.
# `ref` overrides the nominal N a source's shifts are measured against, which the per-channel
# resolution scatter needs: its band belongs around its own median, not around the nominal count.
# `tex` is the same range written for the appendix table, since the plain form carries characters
# LaTeX would choke on.
Src = collections.namedtuple("Src", "label detail dirs kind ref tex", defaults=(None, None))
# N_scan is a (scan, lensed scan) pair, or None where the source does not touch the scan: the lens
# layer is a third basis and every source has to move it too, since a lens view inherits its parent's
# window and resolution and faces the same statistics requirement.
Dir = collections.namedtuple("Dir", "label N_model N_scan")


def sc(**kw):
    """(N, lensed N) of the ten-object scan under one variation."""
    _n, N, _nl, N_l = scan(**kw)
    return N, N_l

_mu_lo, _mu_hi = (N_model_rising(r) for r in MU_HI_BAND)
_p16, _p50, _p84 = scattered()
_ngap, _gap = two_body_gap()
_nonpeak = [o for o in AXES if nonpeak_only(o, pub_models[o])]

SOURCES = [
    Src("mass resolution, scale", f"every r x{R_FACTOR:g} either way",
        [Dir(f"r x{1/R_FACTOR:g}", N_model(rscale=1 / R_FACTOR), sc(sigma_scale=1 / R_FACTOR)),
         Dir(f"r x{R_FACTOR:g}", N_model(rscale=R_FACTOR), sc(sigma_scale=R_FACTOR))],
        "systematic", None, rf"every $r \times {R_FACTOR:g}$ either way"),
    # An alternative correlation assumption for the SAME input as the row above, so it is reported
    # against its own median and never added: the two would double-count the resolution.
    Src("mass resolution, per channel",
        f"each r independently, x{R_SCATTER:g} per sigma (16-84%)",
        [Dir("16%", _p16, None), Dir("84%", _p84, None)], "alternative", _p50,
        rf"each $r$ independently, $\times {R_SCATTER:g}$ per $\sigma$ (16--84\%)"),
    Src("mass resolution, shape",
        f"muon axes with r(M) rising to {MU_HI_BAND[0]:.2f}-{MU_HI_BAND[1]:.2f} at "
        f"{MU_ANCHORS[1][0]:.0f} GeV",
        [Dir("r(3 TeV) = 0.20", _mu_hi, None), Dir("r(3 TeV) = 0.10", _mu_lo, None)],
        "systematic", None,
        rf"muon axes, $r(M)$ rising to {MU_HI_BAND[0]:.2f}--{MU_HI_BAND[1]:.2f} at "
        rf"{MU_ANCHORS[1][0]/1000:g}\,TeV"),
    Src("mass resolution, prescription", "worst leg / quadrature sum instead of the calibrated mean",
        [Dir("worst leg", None, sc(r_group=r_worst)),
         Dir("quadrature", None, sc(r_group=r_quadrature))], "systematic", None,
        "worst leg or quadrature sum, not the calibrated mean"),
    Src("scan windows", f"every edge x{WINDOW_F:g} either way "
                        f"(published), x{GEN_WINDOW_F:g} (generic)",
        [Dir("narrower", N_model(lo_f=WINDOW_F, hi_f=1 / WINDOW_F),
             sc(window_f=(GEN_WINDOW_F, 1 / GEN_WINDOW_F))),
         Dir("wider", N_model(lo_f=1 / WINDOW_F, hi_f=WINDOW_F),
             sc(window_f=(1 / GEN_WINDOW_F, GEN_WINDOW_F)))], "systematic", None,
        rf"every edge $\times {WINDOW_F:g}$ (published), $\times {GEN_WINDOW_F:g}$ (generic)"),
    Src("yield anchor", f"N_ref x{1/ANCHOR_F:g} to x{ANCHOR_F:g}",
        [Dir(f"x{1/ANCHOR_F:g}", None, sc(N_REF=YM.N_REF / ANCHOR_F)),
         Dir(f"x{ANCHOR_F:g}", None, sc(N_REF=YM.N_REF * ANCHOR_F))], "systematic", None,
        rf"$N_{{\mathrm{{ref}}}} \times 10^{{\pm{math.log10(ANCHOR_F):.0f}}}$"),
    Src("background slope", f"P = {P_BAND[0]:g} to {P_BAND[1]:g}",
        [Dir(f"P = {P_BAND[0]:g}", None, sc(P=P_BAND[0])),
         Dir(f"P = {P_BAND[1]:g}", None, sc(P=P_BAND[1]))], "systematic", None,
        rf"$P = {P_BAND[0]:g}$ to ${P_BAND[1]:g}$"),
    Src("fittability requirement",
        f"{EV_BAND[0]:.0f}-{EV_BAND[1]:.0f} events, {BIN_BAND[0]}-{BIN_BAND[1]} elements",
        [Dir("loose", None, sc(MIN_EVENTS=EV_BAND[0], MIN_BINS=BIN_BAND[0])),
         Dir("tight", None, sc(MIN_EVENTS=EV_BAND[1], MIN_BINS=BIN_BAND[1]))],
        "systematic", None,
        rf"{EV_BAND[0]:.0f}--{EV_BAND[1]:.0f} events, {BIN_BAND[0]}--{BIN_BAND[1]} elements"),
    Src("the axis set", f"non-peaking axes and the dilepton overlap dropped; the {_ngap} unscanned "
                        f"two-body pairs added",
        [Dir("dropped", N_model(axes=[o for o in AXES if o not in _nonpeak],
                                extra=-overlap_looks()), None),
         Dir("added", N_model(extra=_gap), None)], "systematic", None,
        rf"non-peaking axes out, the {_ngap} unscanned pairs in"),
    Src("what counts as one look", f"N x{LOOK_LO:g} to x Z/sqrt(2 pi)",
        [Dir(f"x{LOOK_LO:g}", N_MODEL * LOOK_LO, (N_SCAN * LOOK_LO, N_LSCAN * LOOK_LO)),
         Dir("up-crossings", N_MODEL * c_upcrossing(N_MODEL),
             (N_SCAN * c_upcrossing(N_SCAN), N_LSCAN * c_upcrossing(N_LSCAN)))],
        "systematic", None,
        rf"$N \times {LOOK_LO:g}$ to $N Z/\sqrt{{2\pi}}$"),
    Src("the closed-form LEE relation", "exact Gaussian-tail solution instead",
        [Dir("exact", None, None)], "systematic", None,
        "exact Gaussian-tail solution"),
]

# conventions, reported but never added to the band
ALTERNATIVES = [
    ("counting granularity", "published event selections instead of inclusive spectra",
     N_SEL, None),
    ("the dataset", "Run 2 alone instead of Run 2 + Run 3", None, scan(lumi=1.0)[1]),
    ("the lens layer", "the four lenses added to the scan", None, N_LSCAN),
    ("lens efficiency", f"every declared efficiency x{1/LENS_EFF_F:g} to x{LENS_EFF_F:g}", None,
     (scan(lens_f=1 / LENS_EFF_F)[3], scan(lens_f=LENS_EFF_F)[3])),
    ("the costliest single axis", "the largest contributor removed",
     N_MODEL - max(sum(n_s(lo, hi, res(o)) for lo, hi in scan_segments(o)) for o in AXES), None),
]

# ---------------------------------------------------------------- assemble
def dz(N_nom, N_var):
    return None if N_var is None else z5(N_var) - z5(N_nom)


Row = collections.namedtuple("Row", "src kind detail direction N_model N_scan N_lens "
                                    "dz_model dz_scan dz_lens dz_gap")
rows = []
for s in SOURCES:
    for d in s.dirs:
        n_sc, n_ln = d.N_scan if d.N_scan is not None else (None, None)
        if d.label == "exact":                     # the one direction that is not a change of input
            dm, ds = z_exact(N_MODEL) - z5(N_MODEL), z_exact(N_SCAN) - z5(N_SCAN)
            dl = z_exact(N_LSCAN) - z5(N_LSCAN)
        else:
            dm, ds = dz(s.ref or N_MODEL, d.N_model), dz(N_SCAN, n_sc)
            dl = dz(N_LSCAN, n_ln)
        gap = (ds - dm) if (dm is not None and ds is not None) else (ds if dm is None else -dm)
        rows.append(Row(s.label, s.kind, s.detail, d.label, d.N_model, n_sc, n_ln,
                        dm, ds, dl, gap))


def envelope(vals):
    """(down, up) envelope of a set of shifts, each side including zero."""
    v = [x for x in vals if x is not None]
    return (min([0.0] + v), max([0.0] + v))


def quad(pairs):
    lo = math.sqrt(sum(p[0] ** 2 for p in pairs))
    hi = math.sqrt(sum(p[1] ** 2 for p in pairs))
    return -lo, hi


per_src = {}
for s in SOURCES:
    rs = [r for r in rows if r.src == s.label]
    per_src[s.label] = (envelope([r.dz_model for r in rs]),
                        envelope([r.dz_scan for r in rs]),
                        envelope([r.dz_gap for r in rs]),
                        envelope([r.dz_lens for r in rs]))

ADDED = [s for s in SOURCES if s.kind == "systematic"]
TOT_MODEL = quad([per_src[s.label][0] for s in ADDED])
TOT_SCAN = quad([per_src[s.label][1] for s in ADDED])
TOT_GAP = quad([per_src[s.label][2] for s in ADDED])
TOT_LENS = quad([per_src[s.label][3] for s in ADDED])

# ---------------------------------------------------------------- console
def pm(pair, w=6):
    return f"{pair[1]:+.2f}/{pair[0]:+.2f}".rjust(w + 6)


print(f"basis 1, the model space : {len(AXES)} spectra, N = {N_MODEL:,.0f}, "
      f"Z_local = {z5(N_MODEL):.2f}")
print(f"basis 2, the {LUMI_LABEL} scan: {N_SPEC:,} spectra, N = {N_SCAN:,.0f}, "
      f"Z_local = {z5(N_SCAN):.2f}   (with lenses {N_LSPEC:,} histograms, "
      f"N = {N_LSCAN:,.0f}, Z_local = {z5(N_LSCAN):.2f})")
print(f"per-channel resolution scatter: N 16/50/84% = {_p16:,.0f} / {_p50:,.0f} / {_p84:,.0f}")
print(f"up-crossing factor Z/sqrt(2pi): {c_upcrossing(N_MODEL):.2f} (model space), "
      f"{c_upcrossing(N_SCAN):.2f} (scan)")
print(f"exact tail solution: {z_exact(N_MODEL):.2f} against the closed form {z5(N_MODEL):.2f} "
      f"(model space), {z_exact(N_SCAN):.2f} against {z5(N_SCAN):.2f} (scan)\n")

print("the flat muon r against the rising form it stands in for:")
for o in MU_AXES:
    print(f"  {o:12s} declared {res(o):.3f}   r_eff = "
          + ", ".join(f"{r_effective(o, rh):.4f} at r(3 TeV) = {rh:.2f}"
                      for rh in (MU_HI_BAND[0], MU_ANCHORS[1][1], MU_HI_BAND[1])))
print()
print(f"{'source':32s} {'model space':>13} {'scan':>13} {'+ lenses':>13} "
      f"{'difference':>13}   varied over")
for s in SOURCES:
    m, c, g, l = per_src[s.label]
    mark = " " if s.kind == "systematic" else "*"
    print(f"{s.label:31s}{mark} {pm(m):>13} {pm(c):>13} {pm(l):>13} {pm(g):>13}   {s.detail}")
print(f"{'TOTAL (quadrature)':32s} {pm(TOT_MODEL):>13} {pm(TOT_SCAN):>13} {pm(TOT_LENS):>13} "
      f"{pm(TOT_GAP):>13}")
print("* an alternative reading of an input already counted above, never added to the total")
print(f"\nZ_local = {z5(N_MODEL):.2f} {TOT_MODEL[1]:+.2f}/{TOT_MODEL[0]:+.2f} (model space), "
      f"{z5(N_SCAN):.2f} {TOT_SCAN[1]:+.2f}/{TOT_SCAN[0]:+.2f} (scan), "
      f"{z5(N_LSCAN):.2f} {TOT_LENS[1]:+.2f}/{TOT_LENS[0]:+.2f} (with lenses), "
      f"gap {z5(N_SCAN)-z5(N_MODEL):.2f} {TOT_GAP[1]:+.2f}/{TOT_GAP[0]:+.2f}")
print("\nconventions, not uncertainties (never added to the band):")
for label, detail, nm, nc in ALTERNATIVES:
    def _fmt(v):
        if v is None:
            return "-"
        if isinstance(v, tuple):
            return " to ".join(f"{z5(x):.2f}" for x in v)
        return f"{z5(v):.2f}"
    print(f"  {label:26s} {_fmt(nm):>14} {_fmt(nc):>14}   {detail}")

# ---------------------------------------------------------------- CSV
with open(_p("results", "tables", "budget_uncertainty.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["source", "varied_over", "direction", "kind", "N_model_space",
                "N_combinatorial_scan", "N_scan_with_lenses", "dZ_model_space",
                "dZ_combinatorial_scan", "dZ_scan_with_lenses", "dZ_difference"])
    w.writerow(["nominal", f"published windows / {LUMI_LABEL}", "nominal", "nominal",
                f"{N_MODEL:.0f}", f"{N_SCAN:.0f}", f"{N_LSCAN:.0f}", f"{z5(N_MODEL):.3f}",
                f"{z5(N_SCAN):.3f}", f"{z5(N_LSCAN):.3f}", ""])
    for r in rows:
        w.writerow([r.src, r.detail, r.direction, r.kind,
                    "" if r.N_model is None else f"{r.N_model:.0f}",
                    "" if r.N_scan is None else f"{r.N_scan:.0f}",
                    "" if r.N_lens is None else f"{r.N_lens:.0f}",
                    "" if r.dz_model is None else f"{r.dz_model:+.3f}",
                    "" if r.dz_scan is None else f"{r.dz_scan:+.3f}",
                    "" if r.dz_lens is None else f"{r.dz_lens:+.3f}",
                    "" if r.dz_gap is None else f"{r.dz_gap:+.3f}"])
    for s in SOURCES:
        m, c, g, l = per_src[s.label]
        w.writerow([s.label, s.detail, "envelope", s.kind, "", "", "",
                    f"{m[1]:+.3f}/{m[0]:+.3f}", f"{c[1]:+.3f}/{c[0]:+.3f}",
                    f"{l[1]:+.3f}/{l[0]:+.3f}", f"{g[1]:+.3f}/{g[0]:+.3f}"])
    w.writerow(["total", "quadrature over the systematic sources above", "total", "systematic",
                "", "", "",
                f"{TOT_MODEL[1]:+.3f}/{TOT_MODEL[0]:+.3f}",
                f"{TOT_SCAN[1]:+.3f}/{TOT_SCAN[0]:+.3f}",
                f"{TOT_LENS[1]:+.3f}/{TOT_LENS[0]:+.3f}",
                f"{TOT_GAP[1]:+.3f}/{TOT_GAP[0]:+.3f}"])
    for label, detail, nm, nc in ALTERNATIVES:
        vals = []
        for v in (nm, nc):
            if v is None:
                vals.append("")
            elif isinstance(v, tuple):
                vals.append(" to ".join(f"{z5(x):.2f}" for x in v))
            else:
                vals.append(f"{z5(v):.2f}")
        w.writerow([label, detail, "alternative", "convention", "", "", "", vals[0], vals[1],
                    "", ""])

# ---------------------------------------------------------------- LaTeX
def tex_pm(pair):
    if pair == (0.0, 0.0):
        return "---"
    if max(abs(pair[0]), abs(pair[1])) < 0.005:
        return "$<0.01$"
    if abs(pair[0]) < 0.005:                       # one-sided: a single signed number reads better
        return f"${pair[1]:+.2f}$"
    if abs(pair[1]) < 0.005:
        return f"${pair[0]:+.2f}$"
    return f"$^{{{pair[1]:+.2f}}}_{{{pair[0]:+.2f}}}$"


with open(_p("results", "tex", "uncertainty_table.tex"), "w") as f:
    f.write("% Generated by scripts/budget_uncertainty.py. Do not edit: regenerate instead.\n")
    # a wrapping second column: the ranges are too long for an l-column at this text width
    f.write("\\begin{tabular}{@{}l>{\\raggedright\\arraybackslash}p{0.30\\textwidth}ccc@{}}\n"
            "\\toprule\n")
    f.write("source & varied over & 46 spectra & scan & difference \\\\\n\\midrule\n")
    for s in SOURCES:
        m, c, g, _l = per_src[s.label]
        lbl = s.label + ("" if s.kind == "systematic" else "$^{\\dagger}$")
        f.write(f"{lbl} & {s.tex or s.detail} & {tex_pm(m)} & {tex_pm(c)} & {tex_pm(g)} \\\\\n")
    f.write("\\midrule\n")
    f.write(f"total & in quadrature & {tex_pm(TOT_MODEL)} & {tex_pm(TOT_SCAN)} & "
            f"{tex_pm(TOT_GAP)} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------- report
def md_row(s):
    m, c, g, _l = per_src[s.label]
    lbl = s.label + ("" if s.kind == "systematic" else " (&dagger;)")
    return (f"| {lbl} | {s.detail} | {m[1]:+.2f}/{m[0]:+.2f} | {c[1]:+.2f}/{c[0]:+.2f} | "
            f"{g[1]:+.2f}/{g[0]:+.2f} |")


alt_lines = []
for label, detail, nm, nc in ALTERNATIVES:
    def _z(v):
        if v is None:
            return "-"
        if isinstance(v, tuple):
            return " to ".join(f"{z5(x):.2f}" for x in v)
        return f"{z5(v):.2f}"
    alt_lines.append(f"| {label} | {detail} | {_z(nm)} | {_z(nc)} |")

md = f"""# The uncertainty budget

What the two counted bases are uncertain by, source by source. Written by
`scripts/budget_uncertainty.py`, which varies each declared input over the range below and
**recomputes**: a resolution change moves both the looks a spectrum carries and, for a hypothetical
one, whether it can be fitted at all.

| basis | spectra | N_trials | Z_local(5s global) | band |
|---|--:|--:|--:|---|
| the public model space, published windows | {len(AXES)} | {N_MODEL:,.0f} | **{z5(N_MODEL):.2f}** | \
{TOT_MODEL[1]:+.2f}/{TOT_MODEL[0]:+.2f} |
| the combinatorial scan, {LUMI_LABEL} | {N_SPEC:,} | {N_SCAN:,.0f} | **{z5(N_SCAN):.2f}** | \
{TOT_SCAN[1]:+.2f}/{TOT_SCAN[0]:+.2f} |
| ... with the four selection lenses | {N_LSPEC:,} | {N_LSCAN:,.0f} | **{z5(N_LSCAN):.2f}** | \
{TOT_LENS[1]:+.2f}/{TOT_LENS[0]:+.2f} |

## Source by source

Each entry is the shift in `Z_local` [sigma]. The last column is the shift in the **difference**
between the model space and the scan, taken direction by direction, which is what the correlated part
of every source cancels in. The lens layer is not a fourth column because it moves with the scan it
sits on: its total band is {TOT_LENS[1]:+.2f}/{TOT_LENS[0]:+.2f} against the scan's
{TOT_SCAN[1]:+.2f}/{TOT_SCAN[0]:+.2f}, per source in `budget_uncertainty.csv`.

| source | varied over | model space | scan | difference |
|---|---|--:|--:|--:|
{chr(10).join(md_row(s) for s in SOURCES)}
| **total** | in quadrature | **{TOT_MODEL[1]:+.2f}/{TOT_MODEL[0]:+.2f}** | \
**{TOT_SCAN[1]:+.2f}/{TOT_SCAN[0]:+.2f}** | **{TOT_GAP[1]:+.2f}/{TOT_GAP[0]:+.2f}** |

(&dagger;) an alternative reading of an input already counted in the row above it, so it is quoted
against its own median and never added: the two would double-count the resolution.

So the headline numbers are

```
model space, published windows : Z_local = {z5(N_MODEL):.2f} {TOT_MODEL[1]:+.2f}/{TOT_MODEL[0]:+.2f}
combinatorial scan, {LUMI_LABEL}  : Z_local = {z5(N_SCAN):.2f} {TOT_SCAN[1]:+.2f}/{TOT_SCAN[0]:+.2f}
... with the four lenses       : Z_local = {z5(N_LSCAN):.2f} {TOT_LENS[1]:+.2f}/{TOT_LENS[0]:+.2f}
the gap between the first two  : {z5(N_SCAN)-z5(N_MODEL):.2f} {TOT_GAP[1]:+.2f}/{TOT_GAP[0]:+.2f}
```

**What the two bases share cancels in the difference.** The resolution scale and the look convention,
the two largest terms, are worth {per_src['mass resolution, scale'][2][1]:+.2f}/\
{per_src['mass resolution, scale'][2][0]:+.2f} and \
{per_src['what counts as one look'][2][1]:+.2f}/{per_src['what counts as one look'][2][0]:+.2f} on the
gap against up to {max(abs(x) for x in per_src['mass resolution, scale'][1]):.2f} and \
{max(abs(x) for x in per_src['what counts as one look'][1]):.2f} on the bars. What is left on the gap
is the yield model of the hypothetical scan, which the model space never uses: the difference is
therefore better determined than the scan's bar and rests on entirely different inputs from either.
Quote differences rather than bars wherever the argument allows it.

## What each range is

- **mass resolution, scale.** `r` is propagated from ATLAS object performance rather than quoted per
  bump hunt, so it is known to a factor of a few; x{R_FACTOR:g} either way is the band every headline
  carries. On the scan the same scaling also changes which histograms can be fitted, which is why the
  band there is wider and asymmetric: a coarser resolution costs looks twice, once per spectrum and
  again through the 25-element requirement that a coarse window now fails, while a finer one buys
  elements but empties each of them.
- **mass resolution, per channel.** The scale above is fully correlated, which is the conservative
  reading. Drawing each channel's `r` independently at the same factor per sigma leaves
  N = {_p16:,.0f} to {_p84:,.0f} (16-84%, median {_p50:,.0f}), a band of
  {per_src['mass resolution, per channel'][0][1]:+.2f}/\
{per_src['mass resolution, per channel'][0][0]:+.2f} sigma about that median rather than
  {per_src['mass resolution, scale'][0][1]:+.2f}/\
{per_src['mass resolution, scale'][0][0]:+.2f}: the errors average down over 46 spectra, so the
  correlated factor two is the pessimistic end of the range and not a 1 sigma. The median sits above
  the nominal count because 1/r is convex, which is a property of the log-normal draw rather than a
  statement about the resolutions.
- **mass resolution, shape.** The muon axes replace a resolution rising with p_T by a
  window-averaged constant, the one place the constant-`r` framework is a real approximation.
  Counting them with `r(M)` rising linearly from {MU_ANCHORS[0][1]:.2f} at
  {MU_ANCHORS[0][0]:.0f} GeV to {MU_HI_BAND[0]:.2f}-{MU_HI_BAND[1]:.2f} at
  {MU_ANCHORS[1][0]:.0f} GeV brackets it, and the approximation turns out to be a good one: over the
  `m(mumu)` window the rising form is worth an effective flat
  `r_eff` = {r_effective('m(mumu)', MU_ANCHORS[1][1]):.3f} at the central anchor
  ({r_effective('m(mumu)', MU_HI_BAND[1]):.3f} to {r_effective('m(mumu)', MU_HI_BAND[0]):.3f} across
  the band) against the {res('m(mumu)'):.3f} declared, so the whole question is worth less than
  0.01 sigma on the bar.
- **mass resolution, prescription.** The scan takes `r = 1/2 sqrt(mean sigma^2)` over the object
  group, calibrated on the published channels. The alternatives are the worst leg, which is the
  convention of the model-space budget, and textbook quadrature propagation.
- **scan windows.** Published edges are read off papers; what is uncertain is how far the curated
  search family extends, taken as x{WINDOW_F:g} on each edge. A hypothetical spectrum has a
  declared generic window instead, varied by x{GEN_WINDOW_F:g}.
- **yield anchor**, **background slope**, **fittability requirement.** The three declared inputs of
  the statistics requirement. All three act on the scan only, since published windows are exempt.
- **the axis set.** Down: the axes motivated only by non-peaking models
  ({', '.join(f'`{o}`' for o in _nonpeak)}) dropped, and the {overlap_looks():.0f} looks the
  dark-photon axes double-count against the high-mass dilepton axes removed. Up: the {_ngap} object
  pairs with no published axis added at {_gap:,.0f} looks.
- **what counts as one look.** The convention with no measurement behind it. A resonance spans more
  than one element, which argues for fewer independent looks (x{LOOK_LO:g}); the up-crossing form of
  the Gross-Vitells estimate argues for more, by Z/sqrt(2 pi) = {c_upcrossing(N_MODEL):.2f} at this
  Z, since Rice's formula counts (1/2 pi)(1/r) ln(M_hi/M_lo) exp(-Z^2/2) up-crossings against the
  Gaussian tail exp(-Z^2/2)/(Z sqrt(2 pi)) that element counting multiplies. This is the largest
  single term on the model space and, on the scan, second only to the resolution scale; being common
  to every basis, it cancels in every difference.
- **the closed-form LEE relation.** `Z_local = sqrt(25 + 2 ln N)` is the asymptotic solution of
  `N p_local = p(5 sigma)`; solving it exactly gives {z_exact(N_MODEL):.2f} rather than
  {z5(N_MODEL):.2f}, so the closed form is marginally strict.

## Conventions, not uncertainties

These are choices, reported as alternatives and never added to the band.

| choice | alternative | model space | scan |
|---|---|--:|--:|
{chr(10).join(alt_lines)}

The shape of the hypothetical scan is a design choice in the same sense: ten object types, at most
four per mass, one selection lens at a time. It fixes what is being priced rather than carrying an
error, and `results/tables/scaled_scan.csv` prices the variants side by side.

Source: `scripts/budget_uncertainty.py` -> `results/tables/budget_uncertainty.csv`,
`results/tex/uncertainty_table.tex`. The bases themselves: `results/overviews/SEARCH_BUDGET.md`,
`results/tables/scaled_scan.txt`.
"""
open(_p("results", "overviews", "BUDGET_UNCERTAINTY.md"), "w").write(md)
print("\nwrote results/tables/budget_uncertainty.csv, results/tex/uncertainty_table.tex, "
      "results/overviews/BUDGET_UNCERTAINTY.md", file=sys.stderr)
