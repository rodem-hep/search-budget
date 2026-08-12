#!/usr/bin/env python3
"""Does an imperfect estimator change WHICH stage-1 selection rule wins?

Re-runs the argmax / fixed-threshold / Benjamini-Hochberg comparison of
results/overviews/MAX_OF_GAUSSIANS.md Part IV (scripts/bh_fdr_ab.py) with the estimator-defect
model of scripts/ab_split_outliers.py bolted on.

The conjecture worth testing: argmax and BH are RANK-based, so an outlier anywhere in the scan
perturbs the decision about every other bin; the fixed threshold is ABSOLUTE, so an outlier at
1 TeV says nothing about the bin at 3 TeV. If that matters, contamination should hurt the two
rank rules more than the threshold and widen the gap the original study already found.

Same toy as Part IV: n bins, one may carry a signal of strength mu, stage 2 confirms a selected
bin if its INDEPENDENT second significance exceeds t_B = 3. Defects enter the reported z, and the
two classes differ only in whether the defect repeats in the confirmation measurement:
  GLITCH  redrawn per run  -> selected on an artefact, dies in stage 2
  BIAS    same pull twice  -> selected on an artefact, CONFIRMS
p-values are formed with the NOMINAL null, which is what an analyst without knowledge of the
contamination would do -- that is precisely why BH's guarantee fails.

Public inputs only (no repo data: this is a pure statistics toy). Writes
results/plots/max_of_gaussians/bh_outliers.png and prints the tables quoted in section 8.5 of the
note. Needs numpy/scipy. Run it on a batch node, not a shared login node.
"""
import os, sys, math
import numpy as np
from scipy.stats import norm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
OUT = os.path.join(ROOT, "results", "plots", "max_of_gaussians")
os.makedirs(OUT, exist_ok=True)

sf, cdf, pdf = norm.sf, norm.cdf, norm.pdf
trapz = getattr(np, "trapezoid", None) or np.trapz

n, tB = 30_000, 3.0
MUS = [3, 4, 5, 6]
M, T, CHUNK = 1000, 20_000, 5_000
QGRID = np.geomspace(1e-9, 0.5, 90)
TGRID = np.linspace(2.6, 14.0, 160)
BETA = 2.0
ZT = np.linspace(-6.0, 45.0, 30001)          # shared tabulation grid
rng = np.random.default_rng(20260810)

_GX, _GW = np.polynomial.legendre.leggauss(200)
def halfnormal(beta, zmax=10.0):
    hi = zmax * beta
    b = 0.5 * hi * (_GX + 1.0)
    w = 0.5 * hi * _GW * math.sqrt(2.0 / math.pi) / beta * np.exp(-0.5 * (b / beta) ** 2)
    return b, w / w.sum()


class Defect:
    """eps: rate per bin. kind='glitch' (redrawn per run) or 'bias' (same pull in both runs).

    Everything the Monte Carlo needs is tabulated on ZT once and interpolated afterwards: the
    quadrature is over a 200-node defect grid, which must not be evaluated on the 6e6-element
    arrays the MC produces."""
    def __init__(self, kind, eps, beta=BETA):
        self.kind, self.eps, self.beta = kind, eps, beta
        self.d, self.w = halfnormal(beta)
        self.Qt = self._Q(ZT)
        self.ft = self._f(ZT)
        self.pct = self._pconf(ZT)
        lo = -np.log(np.clip(self.Qt, 1e-320, 1.0))
        ok = np.concatenate([[True], np.diff(lo) > 0])           # strictly increasing branch
        self._ilo, self._iz = lo[ok], ZT[ok]

    def _Q(self, z):
        base = (1 - self.eps) * sf(z)
        if self.eps == 0: return base
        return base + self.eps * (self.w @ sf(z[None, :] - self.d[:, None]))

    def _f(self, z):
        base = (1 - self.eps) * pdf(z)
        if self.eps == 0: return base
        return base + self.eps * (self.w @ pdf(z[None, :] - self.d[:, None]))

    def _pconf(self, z):
        """P(second measurement > t_B | reported z) for a BACKGROUND bin.

        GLITCH: the artefact is redrawn, so the second look knows nothing -> sf(t_B).
        BIAS:   the same pull is there again, so a bin selected on a large artefact confirms.
                Posterior-average of sf(t_B - b) over the defect state the reported z implies."""
        if self.kind == "glitch" or self.eps == 0:
            return np.full(z.shape, sf(tB))
        num = (1 - self.eps) * pdf(z) * sf(tB)
        num = num + self.eps * (self.w @ (pdf(z[None, :] - self.d[:, None])
                                          * sf(tB - self.d)[:, None]))
        return num / self._f(z)

    def z_of_u(self, u):                  # inverse mixture tail, for the order statistics
        return np.interp(-np.log(np.clip(u, 1e-300, 1.0)), self._ilo, self._iz)

    def pconf(self, z):
        return np.interp(z, ZT, self.pct)


def smallest_uniforms(size, ntot, m):
    g = np.cumsum(rng.exponential(size=(size, m)), axis=1)
    return g / (g[:, -1] + rng.gamma(ntot + 1 - m, size=size))[:, None]


def bh_K(p_sorted):
    """K(q) for every q in QGRID. R_k = n p_(k)/k; the suffix minimum is non-decreasing in k,
    so K(q) is a searchsorted rather than a scan over all M ranks for each q."""
    k = np.arange(1, p_sorted.shape[1] + 1)
    c = np.minimum.accumulate((n * p_sorted / k)[:, ::-1], axis=1)[:, ::-1]
    return np.stack([np.searchsorted(row, QGRID, side="right") for row in c])


# ---------------------------------------------------------------- the three rules
def analytic(dfc):
    """Fixed threshold (all t) and argmax: 1-D quadrature, no MC needed."""
    fz, pc = dfc.ft, dfc.pct
    Fz = np.clip(1.0 - dfc.Qt, 0.0, 1.0)
    x_thr = np.array([trapz((fz * pc)[ZT >= t], ZT[ZT >= t]) * n for t in TGRID])
    y_thr = {mu: sf(TGRID - mu) * sf(tB - mu) for mu in MUS}
    lead = np.exp((n - 1) * np.log(np.clip(Fz, 1e-300, 1.0)))     # every other bin below z
    x_arg = float(trapz(n * fz * lead * pc, ZT))
    y_arg = {mu: float(trapz(pdf(ZT - mu) * lead, ZT)) * sf(tB - mu) for mu in MUS}
    sel_arg = {mu: float(trapz(pdf(ZT - mu) * lead, ZT)) for mu in MUS}
    return x_thr, y_thr, x_arg, y_arg, sel_arg


def bh(dfc):
    """MC over the M smallest reported p-values. E[false confirmations] and power per q."""
    xq, ek = np.zeros(len(QGRID)), np.zeros(len(QGRID))
    for done in range(0, T, CHUNK):
        s = min(CHUNK, T - done)
        z = -np.sort(-dfc.z_of_u(smallest_uniforms(s, n, M)), axis=1)   # descending in z
        K = bh_K(sf(z))
        cum = np.cumsum(dfc.pconf(z), axis=1)
        idx = np.clip(K - 1, 0, M - 1)
        got = np.take_along_axis(cum, idx, axis=1)
        xq += np.where(K > 0, got, 0.0).sum(axis=0)
        ek += K.sum(axis=0)
    xq /= T; ek /= T

    pw = {}
    for mu in MUS:
        psel = np.zeros(len(QGRID))
        for done in range(0, T, CHUNK):
            s = min(CHUNK, T - done)
            bz = dfc.z_of_u(smallest_uniforms(s, n - 1, M))
            zs = rng.normal(mu, 1.0, size=s)[:, None]
            rank = (bz > zs).sum(axis=1) + 1
            z = -np.sort(-np.hstack([bz, zs]), axis=1)[:, :M]
            psel += (rank[:, None] <= bh_K(sf(z))).sum(axis=0)
        pw[mu] = psel / T * sf(tB - mu)
    return xq, pw, ek


def at_budget(x, y, xt):
    """Power of a swept rule at false-confirmation budget xt. NaN if the rule cannot be tuned
    to that budget at all -- which is itself a result: a contaminated BH keeps rejecting bins
    however small q is made, so it has a FLOOR on its false-alarm rate."""
    x = np.asarray(x, float)
    if not (x.min() <= xt <= x.max()): return float("nan")
    o = np.argsort(x)
    return float(np.interp(np.log(xt), np.log(np.clip(x[o], 1e-300, None)), np.asarray(y)[o]))


def fmt(v):
    return "    n/a" if not np.isfinite(v) else f"{v:6.1f}%"


# ---------------------------------------------------------------- run
print(f"n = {n:,} bins   stage-2 bar t_B = {tB:g}   {T:,} experiments per point", flush=True)
print("power = P(signal selected in stage 1 AND confirmed in stage 2), in %", flush=True)
print("all three rules compared at the SAME expected number of false confirmations "
      "(the argmax's own budget under that contamination)\n", flush=True)

XREF = 1.35e-3          # the perfect-estimator argmax budget = 3.0 sigma global

CASES = [("PERFECT estimator (reproduces the original study)", Defect("glitch", 0.0)),
         ("GLITCH  eps=1e-4   (incoherent: dies in stage 2)", Defect("glitch", 1e-4)),
         ("GLITCH  eps=1e-3", Defect("glitch", 1e-3)),
         ("GLITCH  eps=1e-2", Defect("glitch", 1e-2)),
         ("BIAS    eps=1e-4   (coherent: confirms in stage 2)", Defect("bias", 1e-4)),
         ("BIAS    eps=1e-3", Defect("bias", 1e-3)),
         ("BIAS    eps=1e-2", Defect("bias", 1e-2))]

for label, dfc in ([] if "--figonly" in sys.argv else CASES):
    x_thr, y_thr, x_arg, y_arg, sel_arg = analytic(dfc)
    x_bh, y_bh, ek = bh(dfc)
    print("=" * 100, flush=True)
    print(label)
    print(f"  argmax's own budget E[false conf.] = {x_arg:.3e}  (it has no knob; perfect "
          f"estimator: {XREF:.2e})")
    print(f"  BH floor: smallest reachable budget = {x_bh.min():.2e} at q = {QGRID[0]:.0e}"
          f"   threshold floor = {x_thr.min():.2e}")
    print(f"  all three at the COMMON budget {XREF:.2e}:")
    print(f"  {'mu':>4} {'argmax':>9} {'threshold':>11} {'BH':>9}   "
          f"{'thr - argmax':>13} {'thr - BH':>9}")
    for mu in MUS:
        pa = y_arg[mu] * 100 if abs(math.log(x_arg / XREF)) < 0.05 else float("nan")
        pt = at_budget(x_thr, y_thr[mu], XREF) * 100
        pb = at_budget(x_bh, y_bh[mu], XREF) * 100
        d1 = pt - pa if np.isfinite(pa) else float("nan")
        d2 = pt - pb if np.isfinite(pb) else float("nan")
        print(f"  {mu:>3}s {fmt(pa):>9} {fmt(pt):>11} {fmt(pb):>9}   "
              f"{'    n/a' if not np.isfinite(d1) else f'{d1:+12.1f}':>13} "
              f"{'    n/a' if not np.isfinite(d2) else f'{d2:+8.1f}':>9}", flush=True)
    qa, ka = at_budget(x_bh, QGRID, XREF), at_budget(x_bh, ek, XREF)
    if np.isfinite(qa):
        chk = ka * sf(tB) / XREF
        print(f"  BH needs nominal q = {qa:.2e} to buy that budget; "
              f"E[bins handed to stage 2] = {ka:.3f}"
              + (f"   [closure E[K]*sf(t_B)/budget = {chk:.3f}, exact only for GLITCH]"
                 if dfc.kind == "glitch" else ""), flush=True)
    else:
        print(f"  BH CANNOT reach that budget at any q: even q -> {QGRID[0]:.0e} leaves "
              f"E[K] = {ek[0]:.2f} rejections and {x_bh.min():.2e} false confirmations",
              flush=True)

# ================================================================ scan for the figure
# One signal strength only (mu = 5), so the epsilon scan stays affordable; the table above
# carries the mu dependence.
EPS = np.geomspace(1e-5, 3e-2, 10)
MU_S = 5
CACHE = os.path.join(ROOT, "results", "tables", "bh_outliers_scan.npz")
if os.path.exists(CACHE) and "--refit" not in sys.argv:
    _d = np.load(CACHE)
    EPS, pw_arg, pw_thr, pw_bh, q_need, x_arg_g, x_arg_b = (
        _d["eps"], _d["arg"], _d["thr"], _d["bh"], _d["q"], _d["xg"], _d["xb"])
    print(f"\nloaded scan from {os.path.basename(CACHE)} (--refit to redo)", flush=True)
else:
  pw_arg, pw_thr, pw_bh, q_need, x_arg_g, x_arg_b = ([] for _ in range(6))
  for e in EPS:
    dg, db = Defect("glitch", e), Defect("bias", e)
    xt, yt, xa, ya, _ = analytic(dg)
    xb, yb, _ = bh(dg)
    pw_arg.append(ya[MU_S] * 100 if abs(math.log(xa / XREF)) < 0.05 else np.nan)
    pw_thr.append(at_budget(xt, yt[MU_S], XREF) * 100)
    pw_bh.append(at_budget(xb, yb[MU_S], XREF) * 100)
    q_need.append(at_budget(xb, QGRID, XREF))
    x_arg_g.append(xa)
    x_arg_b.append(analytic(db)[2])
    print(f"  scan eps={e:.1e}: argmax {pw_arg[-1]:.1f}  thr {pw_thr[-1]:.1f}  "
          f"BH {pw_bh[-1]:.1f}   q_needed {q_need[-1]:.2e}", flush=True)
  np.savez(CACHE, eps=EPS, arg=pw_arg, thr=pw_thr, bh=pw_bh, q=q_need,
           xg=x_arg_g, xb=x_arg_b)

from plot_style import (SURF as surf, INK as ink, INK2 as ink2, GRID as grid,
                        C_ARG as c_arg, C_THR as c_thr, C_BKG as c_bh, style)

fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(6.2, 2.75), facecolor=surf,
                                 gridspec_kw=dict(wspace=0.36))
for ax in (a0, a1, a2):
    style(ax); ax.grid(color=grid, lw=0.4, ls=":", alpha=0.8); ax.set_xscale("log")

# (a) the three rules as the estimator degrades
# bh() draws from the module-level `rng`, so without this reseed the perfect-estimator reference
# would depend on how many draws the code above consumed -- which differs between the cached and
# the --refit path, and differs again under --figonly. Reseeding to the ORIGINAL seed makes this
# call reproduce the first CASES iteration exactly, so the dotted reference lines and the q*
# annotation below are the same MC draw the per-case table prints, not an independent one.
rng = np.random.default_rng(20260810)
p0 = analytic(Defect("glitch", 0.0)); b0 = bh(Defect("glitch", 0.0))
ref = dict(arg=p0[3][MU_S] * 100, thr=at_budget(p0[0], p0[1][MU_S], XREF) * 100,
           bh=at_budget(b0[0], b0[1][MU_S], XREF) * 100)
for y, col, lbl, dash in ((pw_thr, c_thr, "threshold", "-"),
                          (pw_bh, c_bh, "BH", (0, (5, 2))),
                          (pw_arg, c_arg, "argmax", (0, (5, 1.5, 1.2, 1.5)))):
    a0.plot(EPS, y, lw=1.3, color=col, ls=dash, label=lbl)
    a0.axhline(ref["thr" if col is c_thr else "bh" if col is c_bh else "arg"],
               color=col, lw=1.0, ls=(0, (2, 3)), alpha=0.7)
a0.fill_between(EPS, pw_arg, pw_thr, color=c_thr, alpha=0.10)
a0.set_ylim(0, 95)
a0.set_xlabel(r"defect rate $\epsilon$", color=ink2)
a0.set_ylabel("signal confirmed [%]", color=ink2)

# (b) the nominal FDR stops meaning anything
a1.loglog(EPS, q_need, lw=1.3, color=c_bh, label="BH")
a1.axhline(0.05, color=ink2, lw=1.2, ls=(0, (4, 3)), label="$q = 0.05$")
a1.axhline(ref_q := at_budget(b0[0], QGRID, XREF), color=ink2, lw=1.0, ls=(0, (2, 3)),
           label=f"perfect: {ref_q:.2f}")
a1.set_xlabel(r"defect rate $\epsilon$", color=ink2)
a1.set_ylabel("nominal $q$", color=ink2)

# (c) the argmax has no dial at all
a2.loglog(EPS, x_arg_b, lw=1.3, color=c_arg, label="BIAS")
a2.loglog(EPS, x_arg_g, lw=1.3, ls=(0, (5, 2)), color=c_arg, alpha=0.65,
          label="GLITCH")
a2.axhline(XREF, color=ink2, lw=1.2, ls=(0, (4, 3)),
           label="target")
a2.set_ylim(3e-4, 3.0)
a2.set_xlabel(r"defect rate $\epsilon$", color=ink2)
a2.set_ylabel("false-confirmation rate", color=ink2)

# One frameless legend per panel and nothing floating: the dotted/dashed reference lines
# are legend entries rather than annotations parked in whichever corner was free.
from matplotlib.lines import Line2D
a0.legend(handles=a0.get_legend_handles_labels()[0]
                  + [Line2D([], [], color=ink2, lw=1.0, ls=(0, (2, 3)),
                            label="perfect estimator")],
          frameon=False, fontsize=7.5, loc="lower left", labelcolor=ink2,
          handlelength=1.7, borderaxespad=0.5)
a1.legend(frameon=False, fontsize=7.5, loc="lower left", labelcolor=ink2,
          handlelength=1.7, borderaxespad=0.5)
a2.legend(frameon=False, fontsize=7.5, loc="upper left", labelcolor=ink2,
          handlelength=1.7, borderaxespad=0.5)

fig.savefig(os.path.join(OUT, "bh_outliers.png"), dpi=400, facecolor=surf, bbox_inches="tight")
print("\nwrote max_of_gaussians/bh_outliers.png")
