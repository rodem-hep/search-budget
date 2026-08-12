#!/usr/bin/env python3
"""The z-cut BH actually applies, pseudo-experiment by pseudo-experiment: z_cut = Phi^-1(1-p_(K)).

Samples its distribution under H0 and with a signal, against the nominal rank-1 bar and the fixed
threshold matched on the same fake budget. Writes results/tables/bh_zcut_per_pe.csv and
results/plots/max_of_gaussians/bh_zcut.png. Heavy -- run it on a batch node (docs/REPRODUCE.md).
"""
import os, csv
import numpy as np
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(ROOT, "results", "tables")
PLOTS = os.path.join(ROOT, "results", "plots", "max_of_gaussians")
os.makedirs(PLOTS, exist_ok=True)
sf, cdf = norm.sf, norm.cdf

n, nb, tB = 30_000, 29_999, 3.0
tstar = norm.isf(1.0 / nb)
M, T = 2500, 200_000
QS = [0.001, 0.01, 0.05, 0.1, 0.2, 0.381, 0.5, 0.8]
MUS = [4, 5]
rng = np.random.default_rng(20260714)


def smallest_uniforms(size, ntot, m):
    g = np.cumsum(rng.exponential(size=(size, m)), axis=1)
    tot = g[:, -1] + rng.gamma(ntot + 1 - m, size=size)
    return g / tot[:, None]


def bh_K(p_sorted, q):
    k = np.arange(1, p_sorted.shape[1] + 1)
    c = np.minimum.accumulate((n * p_sorted / k)[:, ::-1], axis=1)[:, ::-1]
    return (c <= q).sum(axis=1)


def zcuts(mu=None, chunk=8_000):
    """Realized z_cut per PE (nan where BH rejects nothing), plus K, for every q."""
    Z = {q: [] for q in QS}
    K_ = {q: [] for q in QS}
    for done in range(0, T, chunk):
        s = min(chunk, T - done)
        if mu is None:
            p = smallest_uniforms(s, n, M)
        else:
            b = smallest_uniforms(s, nb, M)
            ps = sf(rng.normal(mu, 1.0, size=s))[:, None]
            p = np.sort(np.hstack([b, ps]), axis=1)[:, :M]
        for q in QS:
            K = bh_K(p, q)
            z = np.full(s, np.nan)
            hit = K >= 1
            z[hit] = norm.isf(p[hit, K[hit] - 1])       # p_(K): the weakest accepted bin
            Z[q].append(z); K_[q].append(K)
    return {q: np.concatenate(Z[q]) for q in QS}, {q: np.concatenate(K_[q]) for q in QS}


print(f"n = {n:,}   {T:,} pseudo-experiments per point\n")
Z0, K0 = zcuts(None)
ZS = {mu: zcuts(mu) for mu in MUS}

rows = []
print("Realized BH cut per pseudo-experiment.  z_1(q) = Phi^-1(1-q/n) is the nominal (rank-1) bar;")
print("t_match is the FIXED threshold with the same fake budget E[K|H0].\n")
hdr = (f"{'q':>6} {'z_1(q)':>7} {'t_match':>8} | {'P(K>=1)':>8} {'E[K|K>=1]':>10} "
       f"{'med z_cut':>10} {'[16,84]%':>14} {'min z_cut':>10}")
print(hdr); print("-" * len(hdr))
for q in QS:
    z = Z0[q]; K = K0[q]
    hit = np.isfinite(z)
    zc = z[hit]
    ek = K0[q][K0[q] >= 1].mean() if hit.any() else np.nan
    ekn = K0[q].mean()
    tm = norm.isf(max(ekn, 1e-12) / n)
    z1 = norm.isf(q / n)
    lo, med, hi = np.percentile(zc, [16, 50, 84]) if hit.sum() > 20 else (np.nan,) * 3
    print(f"{q:6.3f} {z1:7.2f} {tm:8.2f} | {hit.mean():8.4f} {ek:10.2f} "
          f"{med:10.2f} {'['+f'{lo:.2f}, {hi:.2f}'+']':>14} {zc.min():10.2f}")
    row = dict(q=q, z1=z1, t_match=tm, EK_null=ekn, P_hit=hit.mean(), EK_given_hit=ek,
               z_med=med, z_p16=lo, z_p84=hi, z_min=zc.min())
    for mu in MUS:
        zs = ZS[mu][0][q]; hs = np.isfinite(zs)
        row[f"z_med_mu{mu}"] = np.median(zs[hs])
        row[f"P_hit_mu{mu}"] = hs.mean()
    rows.append(row)

print("\nSame, but with a signal in the sample (the extra small p-value loosens the bar):")
hdr = f"{'q':>6} {'z_1(q)':>7} | " + " ".join(f"{'med z_cut mu='+str(m):>16}" for m in MUS)
print(hdr); print("-" * len(hdr))
for q in QS:
    print(f"{q:6.3f} {norm.isf(q/n):7.2f} | "
          + " ".join(f"{np.nanmedian(ZS[m][0][q]):16.2f}" for m in MUS))

print("\nHow often the realized cut is LOOSER than the nominal rank-1 bar z_1(q)  (i.e. K >= 2):")
for q in QS:
    K = K0[q]
    print(f"   q = {q:5.3f}   P(K>=2 | K>=1) = {(K >= 2).sum() / max((K >= 1).sum(), 1):.3f}"
          f"     K=1 / 2 / 3+ : {(K==1).mean():.4f} / {(K==2).mean():.4f} / {(K>=3).mean():.4f}")

out = os.path.join(TABLES, "bh_zcut_per_pe.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
print(f"\nwrote {os.path.relpath(out, ROOT)}")

# ---------------------------------------------------------------- figure
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(ROOT, "scripts"))
from plot_style import (SURF as surf, INK as ink, INK2 as ink2, C_BKG as c_bh, C_THR as c_thr,
                        C_ARG as c_bar, style)

fig, ax = plt.subplots(figsize=(5.9, 3.4), facecolor=surf)
style(ax)

qq = np.array(QS)
z1 = norm.isf(qq / n)
tm = np.array([r["t_match"] for r in rows])
med = np.array([r["z_med"] for r in rows])
p16 = np.array([r["z_p16"] for r in rows])
p84 = np.array([r["z_p84"] for r in rows])

ax.fill_between(qq, p16, p84, color=c_bh, alpha=0.16, lw=0,
                label="BH realized cut, central 68% of PEs")
ax.semilogx(qq, med, color=c_bh, lw=1.3, marker="o", ms=5, mec=surf, mew=1.4,
            label="BH realized cut, median PE")
ax.semilogx(qq, z1, color=c_bar, lw=1.5, ls=(0, (5, 2)),
            label="nominal rank-1 bar  $\\Phi^{-1}(1-q/n)$")
ax.semilogx(qq, tm, color=c_thr, lw=1.5, ls=(0, (2, 2.5)),
            label="fixed threshold at the same fake budget")
ax.axvline(0.381, color=ink2, lw=1, ls=(0, (2, 3)))
ax.text(0.381 * 1.06, 3.15, "$q^\\star$", color=ink, fontweight="medium")
ax.set_xlabel("nominal false discovery rate  $q$", color=ink2)
ax.set_ylabel("$z$ of the weakest bin BH accepts", color=ink2)
ax.legend(frameon=False, loc="upper right", labelcolor=ink2)
ax.set_ylim(2.95, 6.3)
fig.savefig(os.path.join(PLOTS, "bh_zcut.png"), dpi=400, facecolor=surf, bbox_inches="tight")
print("wrote bh_zcut.png")
