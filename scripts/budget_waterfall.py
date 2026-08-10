#!/usr/bin/env python3
"""The whole search-budget story in one figure: Z_local required for a 5-sigma-GLOBAL discovery
as a function of the trials factor N, with the granularity levels of this study placed on the
curve (inclusive spectra -> published event selections -> the full ATLAS BSM program).

The N values are read live from the CSVs search_budget.py writes, so the plot never goes stale.

Reads results/tables/search_budget{,_selections}.csv; writes results/plots/budget_waterfall.png.
"""
import os, csv, math, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import z_local_for_global5 as z5
def _tab(n): return os.path.join(ROOT, "results", "tables", n)

# ---- read the levels live from the tables
sb = list(csv.DictReader(open(_tab("search_budget.csv"))))
N_incl = sum(float(r["ns_scan"]) for r in sb)
n_incl = len(sb)

se = list(csv.DictReader(open(_tab("search_budget_selections.csv"))))
N_sel = sum(float(r["ns_with_selections"]) for r in se)
n_sel = sum(int(r["n_event_selections"]) for r in se)

N_full, N_full_lo, N_full_hi = 5e4, 1e4, 1e5   # full ATLAS BSM program (EXCESS_COUNTING.md)

# the kinematic envelope is a reference bound rather than a granularity level, and it lands within
# 3% of the selections level -- it stays in search_budget.csv and out of this figure.
LEVELS = [  # (N, label, description, annotation offset (pts), ha)
    (N_incl, "inclusive spectra",      f"{n_incl} spectra (one per bump observable)",  (-10,  30), "right"),
    (N_sel,  "published event selections", f"{n_sel} channels (flavour/b-tag/boost)",  (14, -40), "left"),
    (N_full, "full ATLAS program",     "~300-500 searches x O(100) SRs",               (-10,  25), "right"),
]

N = np.logspace(2.8, 5.6, 400)
fig, ax = plt.subplots(figsize=(10.5, 6.5))
ax.plot(N, [z5(v) for v in N], color="#4c78a8", lw=2.2)
ax.set_xscale("log")
ax.grid(True, which="both", ls=":", alpha=0.3)

for Nv, name, desc, off, ha in LEVELS:
    zv = z5(Nv)
    ax.plot([Nv], [zv], "o", ms=9, color="#2f4b6e", zorder=5)
    ax.annotate(f"{name}\nN = {Nv:,.0f} -> Z = {zv:.2f}\n{desc}",
                (Nv, zv), textcoords="offset points", xytext=off, ha=ha,
                va="bottom" if off[1] > 0 else "top", fontsize=8,
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8))
ax.axvspan(N_full_lo, N_full_hi, color="#f0e6b8", alpha=0.45, zorder=0)
ax.text(math.sqrt(N_full_lo*N_full_hi), ax.get_ylim()[1] - 0.02, "full-program range",
        ha="center", va="top", fontsize=8, color="#8a6d00")

ax.set_xlabel("trials factor  N  (effective independent looks)")
ax.set_ylabel(r"$Z_{local}$ required for a $5\sigma$ GLOBAL discovery")
ax.set_title("The Look-Elsewhere price of the search program is logarithmic\n"
             r"$Z_{local} = \sqrt{25 + 2\ln N}$ -- a 20x larger program costs only ~+0.5$\sigma$")
txt = (f"from the {n_incl}-spectrum program to the full ATLAS program (~x{N_full/N_incl:.0f} in N):\n"
       f"Z_local rises {z5(N_incl):.2f} -> {z5(N_full):.2f}  (+{z5(N_full)-z5(N_incl):.2f}$\\sigma$)\n"
       f"the budget is robust: any counting convention lands within ~0.2$\\sigma$")
ax.text(0.98, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
        bbox=dict(boxstyle="round", fc="#eef3fb", ec="#4c78a8"))
fig.tight_layout()
out = os.path.join(ROOT, "results", "plots", "budget_waterfall.png")
fig.savefig(out, dpi=130)
print(f"wrote {out}")
for Nv, name, *_ in LEVELS:
    print(f"  {name:28s} N = {Nv:8,.0f}   Z_local = {z5(Nv):.2f}")
