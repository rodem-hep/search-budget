#!/usr/bin/env python3
"""Benjamini-Hochberg FDR control as a stage-1 selection rule, alongside the argmax and the fixed
threshold (results/overviews/MAX_OF_GAUSSIANS.md, Part IV).

Scans the nominal FDR q by Monte Carlo, cached in results/tables/bh_fdr_mc.npz (--refit to redo).
Writes results/tables/bh_fdr_scan.csv and three figures into results/plots/max_of_gaussians/.
Setup and the order-statistics trick behind the sampling: docs/METHOD_NOTES.md.
"""
import os, sys, math, csv
import numpy as np
from scipy.stats import norm
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(ROOT, "results", "tables")
PLOTS = os.path.join(ROOT, "results", "plots", "max_of_gaussians")
os.makedirs(PLOTS, exist_ok=True)
sf, cdf = norm.sf, norm.cdf

# ---------------------------------------------------------------- setup (identical to Part III)
n, nb, tB = 30_000, 29_999, 3.0
tstar = norm.isf(1.0 / nb)                     # threshold rule: lambda(t*) = 1
x_arg = sf(tB)                                 # argmax fake budget = 1.35e-3  (3.0 sigma global)
MUS = [3, 4, 5, 6]
M = 2500                                       # rank truncation
T = 300_000                                    # experiments per configuration
# q -> 1 is degenerate (BH rejects every bin, since p_(n) <= 1 always), so stop the scan at 0.8
QGRID = np.unique(np.concatenate([np.geomspace(1e-3, 0.8, 60), [0.05, 0.1, 0.2, 0.5]]))

p_win = lambda mu: quad(lambda s: norm.pdf(s - mu) * cdf(s) ** nb, mu - 12, mu + 12, limit=400)[0]
y_arg = lambda mu: p_win(mu) * sf(tB - mu) + (1 - p_win(mu)) * sf(tB)      # argmax power
x_thr = lambda t: nb * sf(t) * sf(tB)                                      # threshold ROC
y_thr = lambda mu, t: cdf(mu - t) * cdf(mu - tB)

rng = np.random.default_rng(20260713)


def smallest_uniforms(size, ntot, m):
    """Exact joint law of the m smallest order statistics of ntot iid U(0,1)."""
    g = np.cumsum(rng.exponential(size=(size, m)), axis=1)
    tot = g[:, -1] + rng.gamma(ntot + 1 - m, size=size)
    return g / tot[:, None]


def bh_counts(p_sorted, qgrid):
    """K(q) for the BH step-up, vectorized over experiments and over q.

    R_k = n p_(k) / k is the q at which rank k alone would be rejectable; the suffix minimum
    C_k = min_{j>=k} R_j is non-decreasing in k, and K(q) = #{k : C_k <= q}."""
    k = np.arange(1, p_sorted.shape[1] + 1)
    r = n * p_sorted / k
    c = np.minimum.accumulate(r[:, ::-1], axis=1)[:, ::-1]          # suffix min -> sorted
    return np.stack([(c <= q).sum(axis=1) for q in qgrid], axis=1)  # (T, nq)


def run_null(T, chunk=8_000):
    """H0: all n bins background. Returns E[K], P(K>=1), P(>=1 false confirmation)."""
    ek = np.zeros(len(QGRID)); p1 = np.zeros(len(QGRID)); pf = np.zeros(len(QGRID)); kmax = 0
    for done in range(0, T, chunk):
        s = min(chunk, T - done)
        K = bh_counts(smallest_uniforms(s, n, M), QGRID)
        ek += K.sum(axis=0); p1 += (K >= 1).sum(axis=0)
        pf += (1.0 - (1.0 - sf(tB)) ** K).sum(axis=0)               # each candidate fakes z_B>3
        kmax = max(kmax, K.max())
    assert kmax < M, f"rank truncation hit (K_max={kmax} >= M={M})"
    return ek / T, p1 / T, pf / T


def run_signal(mu, T, chunk=8_000):
    """n-1 background + one signal bin at mu. Returns P(signal selected), E[K_bkg selected]."""
    psel = np.zeros(len(QGRID)); ebkg = np.zeros(len(QGRID))
    for done in range(0, T, chunk):
        s = min(chunk, T - done)
        b = smallest_uniforms(s, nb, M)                             # background order stats
        ps = sf(rng.normal(mu, 1.0, size=s))[:, None]               # the signal's p-value
        rank = (b < ps).sum(axis=1) + 1                             # its rank among all n
        p = np.sort(np.hstack([b, ps]), axis=1)[:, :M]              # combined, truncated
        K = bh_counts(p, QGRID)
        sel = rank[:, None] <= K                                    # signal survives stage 1
        psel += sel.sum(axis=0)
        ebkg += (K - sel).sum(axis=0)
    return psel / T, ebkg / T


# ---------------------------------------------------------------- run (cached: the MC is the slow part)
CACHE = os.path.join(TABLES, "bh_fdr_mc.npz")
if os.path.exists(CACHE) and "--refit" not in sys.argv:
    d = np.load(CACHE)
    QGRID, EK0, P1_0, PF0 = d["q"], d["ek"], d["p1"], d["pf"]
    SEL = {mu: d[f"sel{mu}"] for mu in MUS}
    BKG_ALT = {mu: d[f"bkg{mu}"] for mu in MUS}
    EKs, PFs = d["eks"], d["pfs"]
    SELs = {mu: float(d[f"sels{mu}"]) for mu in MUS}
    BKGs = {mu: float(d[f"bkgs{mu}"]) for mu in MUS}
    q_star = float(d["q_star"])
    print(f"n = {n:,}   stage-B bar z_B > {tB:g}   loaded MC from {os.path.basename(CACHE)}\n")
else:
    print(f"n = {n:,}   stage-B bar z_B > {tB:g}   {T:,} experiments per point\n")
    EK0, P1_0, PF0 = run_null(T)
    SEL, BKG_ALT = {}, {}
    for mu in MUS:
        SEL[mu], BKG_ALT[mu] = run_signal(mu, T)

    q_star = float(np.interp(1.0, EK0, QGRID))        # E[K|H0] = 1  <=> argmax's fake budget

    # BH at q*, evaluated on its own (a dedicated run at the single q)
    QG_SAVE = QGRID.copy()
    QGRID = np.array([q_star])
    EKs, _, PFs = run_null(T)
    SELs, BKGs = {}, {}
    for mu in MUS:
        a, b = run_signal(mu, T)
        SELs[mu], BKGs[mu] = float(a[0]), float(b[0])
    QGRID = QG_SAVE

    np.savez(CACHE, q=QGRID, ek=EK0, p1=P1_0, pf=PF0, eks=EKs, pfs=PFs, q_star=q_star,
             **{f"sel{mu}": SEL[mu] for mu in MUS}, **{f"bkg{mu}": BKG_ALT[mu] for mu in MUS},
             **{f"sels{mu}": SELs[mu] for mu in MUS}, **{f"bkgs{mu}": BKGs[mu] for mu in MUS})

x_bh = EK0 * sf(tB)                                   # average false confirmations

print(f"Daniels' theorem check (H0): P(BH makes >=1 rejection) should equal q")
for q in (0.05, 0.2, 0.5):
    i = int(np.argmin(abs(QGRID - q)))
    print(f"   q = {QGRID[i]:.3f}   P(>=1 rejection) = {P1_0[i]:.4f}")
print()
print("BH under H0 -- the nominal FDR q as a fake-bin budget")
print(f"{'q':>7}  {'E[K|H0]':>9}  {'P(>=1 rej)':>11}  {'E[false conf.]':>15}  {'global Z':>9}")
for q in (0.01, 0.05, 0.1, 0.2, 0.5, 0.8):
    i = int(np.argmin(abs(QGRID - q)))
    xx = EK0[i] * sf(tB)
    print(f"{QGRID[i]:7.3f}  {EK0[i]:9.3f}  {P1_0[i]:11.4f}  {xx:15.2e}  "
          f"{norm.isf(PF0[i]):9.2f}")
print(f"\nq* (E[K|H0] = 1, the argmax's fake budget)        = {q_star:.3f}")
print(f"   effective single-bin bar  Phi^-1(1 - q*/n)     = {norm.isf(q_star/n):.2f}"
      f"   (threshold rule: t* = {tstar:.2f})")
print(f"   E[false confirmations] at q*                   = {EKs[0]*sf(tB):.2e}"
      f"   (argmax: {x_arg:.2e})")
print(f"   P(>=1 false confirmation) at q*                = {PFs[0]:.3e}"
      f"   -> {norm.isf(PFs[0]):.2f} sigma global (argmax: {x_arg:.3e} -> 3.00)")
print()
print("Power at a matched 3.0 sigma global false-alarm rate")
print(f"{'mu':>4}  {'argmax':>8}  {'threshold t*':>13}  {'BH at q*':>9}  {'BH selects (A)':>15}")
for mu in MUS:
    print(f"{mu:>3}s  {y_arg(mu)*100:7.1f}%  {y_thr(mu, tstar)*100:12.1f}%  "
          f"{SELs[mu]*cdf(mu-tB)*100:8.1f}%  {SELs[mu]*100:7.1f}% sig / "
          f"{BKGs[mu]:.2f} bkg")

with open(os.path.join(TABLES, "bh_fdr_scan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["q", "EK_null", "P_reject_null", "E_false_conf", "P_false_conf"]
               + [f"P_sel_mu{mu}" for mu in MUS] + [f"power_mu{mu}" for mu in MUS])
    for i, q in enumerate(QGRID):
        w.writerow([f"{q:.4f}", f"{EK0[i]:.4f}", f"{P1_0[i]:.4f}", f"{x_bh[i]:.3e}",
                    f"{PF0[i]:.3e}"] + [f"{SEL[mu][i]:.4f}" for mu in MUS]
                   + [f"{SEL[mu][i]*cdf(mu-tB):.4f}" for mu in MUS])
print(f"\nwrote results/tables/bh_fdr_scan.csv")

# ---------------------------------------------------------------- plots (style of Part III)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from plot_style import (SURF as surf, INK as ink, INK2 as ink2, GRID as grid, BLUE,
                        C_BKG as c_bkg, C_ARG as c_fake, C_ARG as c_arg, C_THR as c_thr, style)

# ============ FIGURE 1: the scan (analogue of threshold_scan.png) ============
fig, (a0, a1) = plt.subplots(1, 2, figsize=(13.4, 5.1), facecolor=surf,
                             gridspec_kw=dict(wspace=0.22))
for ax in (a0, a1): style(ax)

a0.grid(color=grid, lw=0.8, alpha=0.6)
a0.loglog(QGRID, EK0, color=c_bkg, lw=2.2, label="E[K]  (bins BH hands to B)")
a0.loglog(QGRID, P1_0, color=ink2, lw=1.4, ls=(0, (4, 3)),
          label="P(any rejection) $= q$  (Daniels)")
a0.axhline(1, color=ink2, lw=1, ls=(0, (2, 3)))
a0.axvline(q_star, color=ink2, lw=1, ls=(0, (2, 3)))
a0.plot([q_star], [1], "o", ms=9, color=c_bkg, mec=surf, mew=2, zorder=6)
a0.plot([0.05], [np.interp(0.05, QGRID, EK0)], "o", ms=7, color=c_bkg, mec=surf, mew=1.8, zorder=6)
a0.text(0.043, np.interp(0.05, QGRID, EK0) * 1.9, "$q = 0.05$", color=ink2, fontsize=10,
        ha="right")
a0.text(0.75, 1.35, "one expected fake", color=ink2, fontsize=10.5, ha="right")
a0.text(q_star * 1.05, 3e-3, f"$q^\\star = {q_star:.2f}$", color=ink, fontsize=12,
        fontweight="medium")
a0.set_xlim(9e-4, 1.05); a0.set_ylim(8e-4, 6)
a0.set_xlabel("nominal false discovery rate  $q$", color=ink2, fontsize=10.5)
a0.set_ylabel("expected background bins selected in A", color=ink2, fontsize=10.5)
a0.legend(frameon=False, fontsize=10, loc="lower right", labelcolor=ink2)
a0.set_title("Background bins BH selects",
             color=ink, fontsize=13, pad=12, loc="left")

a1.grid(color=grid, lw=0.8, alpha=0.6)
for m in MUS:
    a1.semilogx(QGRID, SEL[m] * 100, color=BLUE[m], lw=2.2)
    a1.semilogx(QGRID, cdf(m - norm.isf(QGRID / n)) * 100, color=BLUE[m], lw=1.0,
                ls=(0, (2, 2.5)), alpha=0.9)
    a1.plot([q_star], [SELs[m] * 100], "o", ms=7, color=BLUE[m], mec=surf, mew=1.8, zorder=6)
    a1.text(0.92, SELs[m] * 100, f"$\\mu={m}\\sigma$", color=BLUE[m], fontsize=11.5,
            va="center", fontweight="medium")
a1.axvline(q_star, color=ink2, lw=1, ls=(0, (2, 3)))
a1.axhline(50, color=ink2, lw=1, ls=(0, (4, 3)), alpha=0.7)
a1.text(1.1e-3, 92, "solid: BH step-up (MC)\ndashed: the rank-1 bar alone",
        color=ink2, fontsize=10.5, linespacing=1.5)
a1.text(q_star * 1.06, 4, "$q^\\star$", color=ink, fontsize=12, fontweight="medium")
a1.set_xlim(9e-4, 1.35); a1.set_ylim(0, 103)
a1.set_xlabel("nominal false discovery rate  $q$", color=ink2, fontsize=10.5)
a1.set_ylabel("signal bin passes stage 1   [%]", color=ink2, fontsize=10.5)
a1.set_title("Chance the signal bin passes BH",
             color=ink, fontsize=13, pad=12, loc="left")
fig.savefig(os.path.join(PLOTS, "bh_scan.png"), dpi=170, facecolor=surf, bbox_inches="tight")

# ============ FIGURE 2: power & false alarm vs q (analogue of threshold_vs_argmax.png) ========
fig, (b0, b1) = plt.subplots(1, 2, figsize=(13.4, 5.1), facecolor=surf,
                             gridspec_kw=dict(wspace=0.24))
for ax in (b0, b1): style(ax)

b0.grid(color=grid, lw=0.8, alpha=0.6)
for m in MUS:
    b0.semilogx(QGRID, SEL[m] * cdf(m - tB) * 100, color=BLUE[m], lw=2.2)
    pa = y_arg(m) * 100
    b0.plot([9e-4, 0.9], [pa, pa], color=BLUE[m], lw=1.1, ls=(0, (2, 2.5)), alpha=0.8)
    b0.text(0.95, pa, f"${m}\\sigma$", color=BLUE[m], fontsize=11.5, va="center",
            fontweight="medium")
    b0.plot([q_star], [SELs[m] * cdf(m - tB) * 100], "o", ms=7, color=BLUE[m], mec=surf,
            mew=1.8, zorder=6)
b0.axvline(q_star, color=ink2, lw=1, ls=(0, (2, 3)))
b0.text(1.1e-3, 110, "solid: BH at nominal FDR $q$        dashed: argmax rule",
        color=ink2, fontsize=10.5, linespacing=1.5)
b0.text(q_star * 1.06, 2.5, "$q^\\star$", color=ink, fontsize=12, fontweight="medium")
b0.set_xlim(9e-4, 1.5); b0.set_ylim(0, 124)
b0.set_xlabel("nominal false discovery rate  $q$", color=ink2, fontsize=10.5)
b0.set_ylabel("signal confirmed in B   [%]", color=ink2, fontsize=10.5)
b0.set_title("Confirmation power",
             color=ink, fontsize=13, pad=12, loc="left")

b1.grid(color=grid, lw=0.8, alpha=0.6)
b1.loglog(QGRID, PF0, color=c_fake, lw=2.2)
b1.axhline(sf(tB), color=ink2, lw=1.2, ls=(0, (4, 3)))
b1.axvline(q_star, color=ink2, lw=1, ls=(0, (2, 3)))
b1.plot([q_star], [np.interp(q_star, QGRID, PF0)], "o", ms=9, color=c_fake, mec=surf, mew=2,
        zorder=6)
b1.fill_between(QGRID, PF0, 1e-9, where=(QGRID <= q_star), color=c_fake, alpha=0.06)
b1.text(1.2e-3, 2e-4, "$P_{\\rm fake} \\simeq 1-e^{-E[K|H_0]\\,[1-\\Phi(3)]}$",
        color=c_fake, fontsize=12.5)
b1.text(0.78, 2.3e-3, r"argmax rule:  $1.35\times10^{-3}$", color=ink2, fontsize=10.5,
        ha="right")
b1.text(q_star * 1.06, 1.6e-5, f"$q^\\star = {q_star:.2f}$", color=ink, fontsize=11)
b1.set_xlim(9e-4, 0.85); b1.set_ylim(1e-6, 3e-3)
b1.set_xlabel("nominal false discovery rate  $q$", color=ink2, fontsize=10.5)
b1.set_ylabel("P(at least one false confirmation)", color=ink2, fontsize=10.5)
b1.set_title("Global false-alarm rate it buys", color=ink, fontsize=13,
             pad=12, loc="left")
fig.savefig(os.path.join(PLOTS, "bh_vs_argmax.png"), dpi=170, facecolor=surf, bbox_inches="tight")

# ============ FIGURE 3: the ROC (analogue of roc_threshold_vs_argmax.png) ============
fig, (c0, c1) = plt.subplots(1, 2, figsize=(14.2, 5.4), facecolor=surf,
                             gridspec_kw=dict(width_ratios=[1.18, 1], wspace=0.24))
for ax in (c0, c1): style(ax)

t = np.linspace(2.55, 7.0, 800)
c0.grid(color=grid, lw=0.8, alpha=0.6)
c0.axvline(x_arg, color=c_arg, lw=1.1, ls=(0, (3, 3)), alpha=0.85, zorder=1)
for m in MUS:
    c0.semilogx(x_thr(t), y_thr(m, t) * 100, color=c_thr, lw=1.4, ls=(0, (5, 2)), alpha=0.9,
                zorder=2)
    c0.semilogx(x_bh, SEL[m] * cdf(m - tB) * 100, color=BLUE[m], lw=2.3, zorder=3)
    c0.plot([x_arg], [y_arg(m) * 100], marker="D", ms=9, color=c_arg, mec=surf, mew=1.8, zorder=6)
    c0.plot([EKs[0] * sf(tB)], [SELs[m] * cdf(m - tB) * 100], "o", ms=8, color=BLUE[m], mec=surf,
            mew=1.8, zorder=6)
for qq in (0.05, 0.2, 0.8):
    xq = np.interp(qq, QGRID, x_bh); yq = np.interp(qq, QGRID, SEL[5] * cdf(5 - tB)) * 100
    c0.plot([xq], [yq], marker="|", ms=10, color=BLUE[5], mew=2, zorder=5)
    c0.text(xq, yq - 7, f"$q={qq:g}$", color=BLUE[5], fontsize=9.5, ha="center")
for m, ylab in zip(MUS, [30, 72, 92, 101]):
    c0.text(0.30, ylab, f"$\\mu={m}\\sigma$", color=BLUE[m], fontsize=11.5, va="center",
            fontweight="medium")
c0.text(1.6e-6, 122, "◆  argmax rule", color=c_arg, fontsize=11, fontweight="medium")
c0.text(1.6e-6, 112, "—  BH at nominal FDR $q$", color=BLUE[5], fontsize=11, fontweight="medium")
c0.text(1.6e-6, 102, "- -  fixed threshold $t$ (Part III)", color=c_thr, fontsize=11, fontweight="medium")
c0.text(1.0e-3, 4, "argmax budget", color=c_arg, fontsize=9.5,
        ha="right", linespacing=1.4)
c0.set_xlim(1.2e-6, 8e-1); c0.set_ylim(0, 132)
c0.set_xlabel("average number of false positives   $E[K_{\\rm fake}]$", color=ink2, fontsize=10.5)
c0.set_ylabel("signal confirmed in B   [%]", color=ink2, fontsize=10.5)
c0.set_title("Power against false-alarm rate",
             color=ink, fontsize=12.8, pad=12, loc="left")

# right panel: the realized FDR of the surviving candidates, before and after stage B
c1.grid(color=grid, lw=0.8, alpha=0.6)
for m in MUS:
    fdr_A = BKG_ALT[m] / np.maximum(BKG_ALT[m] + SEL[m], 1e-12)
    fdr_B = (BKG_ALT[m] * sf(tB)) / np.maximum(BKG_ALT[m] * sf(tB) + SEL[m] * cdf(m - tB), 1e-12)
    c1.semilogx(QGRID, fdr_A * 100, color=BLUE[m], lw=2.2)
    c1.semilogx(QGRID, fdr_B * 100, color=BLUE[m], lw=1.1, ls=(0, (2, 2.5)), alpha=0.9)
c1.semilogx(QGRID, QGRID * 100, color=ink2, lw=1.2, ls=(0, (4, 3)))
c1.axvline(q_star, color=ink2, lw=1, ls=(0, (2, 3)))
c1.text(0.012, 62, "nominal $q$", color=ink2, fontsize=10.5, rotation=32)
c1.text(1.15e-3, 92, "solid: after stage 1\ndashed: after the stage-B confirmation",
        color=ink2, fontsize=10.5, linespacing=1.5)
for m, qlab in zip(MUS, [0.006, 0.03, 0.09, 0.16]):
    fdr_A = BKG_ALT[m] / np.maximum(BKG_ALT[m] + SEL[m], 1e-12)
    c1.text(qlab, np.interp(qlab, QGRID, fdr_A) * 100 + 3.5, f"${m}\\sigma$", color=BLUE[m],
            fontsize=11, va="bottom", ha="center", fontweight="medium")
c1.text(q_star * 1.06, 4, "$q^\\star$", color=ink, fontsize=12, fontweight="medium")
c1.set_xlim(9e-4, 1.35); c1.set_ylim(0, 103)
c1.set_xlabel("nominal false discovery rate  $q$", color=ink2, fontsize=10.5)
c1.set_ylabel("realized false discovery rate   [%]", color=ink2, fontsize=10.5)
c1.set_title("Realized against nominal FDR", color=ink, fontsize=12.8,
             pad=12, loc="left")
fig.savefig(os.path.join(PLOTS, "roc_bh_vs_threshold.png"), dpi=170, facecolor=surf,
            bbox_inches="tight")
print("wrote bh_scan.png, bh_vs_argmax.png, roc_bh_vs_threshold.png")
