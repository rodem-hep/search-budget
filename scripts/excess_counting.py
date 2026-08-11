#!/usr/bin/env python3
"""Expected vs observed 3-sigma / 5-sigma excesses across the ATLAS search program.

Under background-only, N quasi-independent tests give N * p_1sided(Z) upward fluctuations above Z.
Reads results/tables/search_budget.csv; writes results/plots/excess_counting.png.
Anchors and their sources: results/overviews/EXCESS_COUNTING.md.
"""
import os, csv, math, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from plot_style import style, title, labels, MARK, GRID, INK2, C_ARG
def p1(Z): return 0.5 * math.erfc(Z / math.sqrt(2.0))   # one-sided Gaussian tail
P3, P5 = p1(3.0), p1(5.0)

# our bump-hunt program: the budget headline, read live so it never goes stale
try:
    with open(os.path.join(ROOT, "results", "tables", "search_budget.csv")) as f:
        N_res = sum(float(r["ns_scan"]) for r in csv.DictReader(f))
except Exception:
    N_res = 3.7e3           # fallback: last computed headline
N_atlas = 5.0e4            # full ATLAS BSM program, central

N = np.logspace(3, 6, 400)
exp3 = N * P3
exp5 = N * P5

fig, ax = plt.subplots(figsize=(9.5, 6.0))
style(ax, grid=False)
ax.grid(True, which="major", color=GRID, lw=0.8, alpha=0.6)
ax.plot(N, exp3, color=C_ARG, lw=2.2, label=r"expected $\geq3\sigma$  ($N\cdot p_{3\sigma}$)")
ax.plot(N, exp5, color=MARK, lw=2.2, ls="--", label=r"expected $\geq5\sigma$  ($N\cdot p_{5\sigma}$)")
ax.axhspan(10, 70, color=C_ARG, alpha=0.09, lw=0,
           label=r"observed $\geq3\sigma$: order tens, none confirmed")

ax.axvline(N_res, color=INK2, ls=":", lw=1.2)
ax.plot([N_atlas], [N_atlas*P3], "o", ms=8, color=C_ARG, mec="white", mew=1.6, zorder=5)
ax.plot([N_atlas], [N_atlas*P5], "o", ms=8, color=MARK, mec="white", mew=1.6, zorder=5)
ax.annotate("this program", (N_res, 1.5e-4), textcoords="offset points", xytext=(6, 0),
            color=INK2, fontsize=9.5, va="bottom")
ax.annotate("full ATLAS program", (N_atlas, N_atlas*P5), textcoords="offset points",
            xytext=(10, -4), color=INK2, fontsize=9.5, va="top")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1e3, 1.3e6); ax.set_ylim(1e-4, 3e2)
labels(ax, x="trials factor  N  (independent looks)",
           y="expected number of upward fluctuations")
title(ax, "Background-only excesses the search program should produce", size=12.5)
ax.legend(loc="upper left", bbox_to_anchor=(0.55, 0.72), fontsize=9.5, frameon=False,
          labelcolor=INK2)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "results", "plots", "excess_counting.png"), dpi=130)
print("wrote results/plots/excess_counting.png")
print(f"N=5e4 -> expected >=3s = {N_atlas*P3:.0f}, >=5s = {N_atlas*P5:.3f}")
print(f"N=1e4 -> {1e4*P3:.0f} / {1e4*P5:.4f} ;  N=1e5 -> {1e5*P3:.0f} / {1e5*P5:.3f}")
print(f"resonance N={N_res:.0f} -> {N_res*P3:.1f} / {N_res*P5:.5f}")
