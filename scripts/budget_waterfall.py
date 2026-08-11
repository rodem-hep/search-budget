#!/usr/bin/env python3
"""The whole search-budget story in one figure: Z_local required for a 5-sigma-GLOBAL discovery
as a function of the trials factor N, with the granularity levels of this study placed on the
curve (inclusive spectra -> published event selections -> the full ATLAS BSM program).

The N values are read live from the CSVs search_budget.py writes, so the plot never goes stale.

Reads results/tables/search_budget{,_selections}.csv; writes results/plots/budget_waterfall.png.
"""
import os, csv, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import z_local_for_global5 as z5
from plot_style import style, title, labels, MARK, GRID, INK2
def _tab(n): return os.path.join(ROOT, "results", "tables", n)

# ---- read the levels live from the tables
sb = list(csv.DictReader(open(_tab("search_budget.csv"))))
N_incl = sum(float(r["ns_scan"]) for r in sb)
n_incl = len(sb)

se = list(csv.DictReader(open(_tab("search_budget_selections.csv"))))
N_sel = sum(float(r["ns_with_selections"]) for r in se)
n_sel = sum(int(r["n_event_selections"]) for r in se)

N_full = 5e4          # full ATLAS BSM program; its 1e4-1e5 range is in EXCESS_COUNTING.md

# the kinematic envelope is a reference bound rather than a granularity level, and it lands within
# 3% of the selections level -- it stays in search_budget.csv and out of this figure.
LEVELS = [  # (N, label, annotation offset (pts), ha)
    (N_incl, f"{n_incl} inclusive spectra",         (-12,  22), "right"),
    (N_sel,  f"{n_sel} published event selections", (14,  -26), "left"),
    (N_full, "full ATLAS program",                  (-12,  20), "right"),
]

N = np.logspace(2.8, 5.6, 400)
fig, ax = plt.subplots(figsize=(9.5, 5.8))
style(ax, grid=False)
ax.grid(True, which="major", color=GRID, lw=0.8, alpha=0.6)
ax.plot(N, [z5(v) for v in N], color=MARK, lw=2.2)
ax.set_xscale("log")

for Nv, name, off, ha in LEVELS:
    zv = z5(Nv)
    ax.plot([Nv], [zv], "o", ms=9, color=MARK, mec="white", mew=2, zorder=5)
    ax.annotate(f"{name}\nN = {Nv:,.0f},  Z = {zv:.2f}",
                (Nv, zv), textcoords="offset points", xytext=off, ha=ha,
                va="bottom" if off[1] > 0 else "top", fontsize=9.5, color=INK2,
                arrowprops=dict(arrowstyle="-", color=GRID, lw=0.8))

labels(ax, x="trials factor  N  (independent looks)",
           y=r"$Z_{local}$ needed for a $5\sigma$ global discovery")
title(ax, r"The look-elsewhere price grows as $\sqrt{25 + 2\ln N}$", size=12.5)
fig.tight_layout()
out = os.path.join(ROOT, "results", "plots", "budget_waterfall.png")
fig.savefig(out, dpi=130)
print(f"wrote {out}")
for Nv, name, *_ in LEVELS:
    print(f"  {name:28s} N = {Nv:8,.0f}   Z_local = {z5(Nv):.2f}")
