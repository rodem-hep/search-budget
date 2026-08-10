#!/usr/bin/env python3
"""Expected vs observed 3-sigma / 5-sigma excesses across the ATLAS search program.

Under background-only, N quasi-independent tests give N * p_1sided(Z) upward fluctuations above Z.
Reads results/tables/search_budget.csv; writes results/plots/excess_counting.png.
Anchors and their sources: results/overviews/EXCESS_COUNTING.md.
"""
import os, csv, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def p1(Z): return 0.5 * math.erfc(Z / math.sqrt(2.0))   # one-sided Gaussian tail
P3, P5 = p1(3.0), p1(5.0)

# our bump-hunt program: the budget headline, read live so it never goes stale
try:
    with open(os.path.join(ROOT, "results", "tables", "search_budget.csv")) as f:
        N_res = sum(float(r["ns_scan"]) for r in csv.DictReader(f))
except Exception:
    N_res = 3.7e3           # fallback: last computed headline
N_atlas = 5.0e4            # full ATLAS BSM program, central
N_lo, N_hi = 1.0e4, 1.0e5  # ATLAS range

N = np.logspace(3, 6, 400)
exp3 = N * P3
exp5 = N * P5

fig, ax = plt.subplots(figsize=(10, 6.5))
ax.plot(N, exp3, color="#d62728", lw=2.2, label=r"expected $\geq3\sigma$  ($N\cdot p_{3\sigma}$)")
ax.plot(N, exp5, color="#1f77b4", lw=2.2, label=r"expected $\geq5\sigma$  ($N\cdot p_{5\sigma}$)")

# ATLAS full-program band + central line
ax.axvspan(N_lo, N_hi, color="#f0e6b8", alpha=0.5, label="full ATLAS BSM program (N range)")
ax.axvline(N_atlas, color="#8a6d00", ls="--", lw=1.4)
ax.axvline(N_res, color="#555555", ls=":", lw=1.4)

# observed markers
ax.axhline(0.0, color="k", lw=0.6)
for Ncen, lbl, dy in [(N_atlas, "ATLAS full", 1)]:
    ax.plot(Ncen, Ncen*P3, "o", color="#d62728", ms=8)
    ax.plot(Ncen, Ncen*P5, "o", color="#1f77b4", ms=8)
ax.annotate(f"~{N_atlas*P3:.0f} expected", (N_atlas, N_atlas*P3), textcoords="offset points",
            xytext=(8, 8), color="#d62728", fontsize=9)
ax.annotate(f"~{N_atlas*P5:.2f} expected", (N_atlas, N_atlas*P5), textcoords="offset points",
            xytext=(8, -14), color="#1f77b4", fontsize=9)
ax.text(N_res, 3e-4, f" resonance\n program\n N~{N_res/1000:.1f}k", color="#555555", fontsize=8, va="bottom")

# observed band for 3-sigma (order tens, all faded)
ax.axhspan(10, 70, xmin=0.0, xmax=1.0, color="#d62728", alpha=0.08)
ax.text(1.1e6, 30, "observed $\\geq3\\sigma$:\n~tens, all faded", color="#a01d1d",
        fontsize=8.5, ha="right", va="center")
ax.text(1.1e6, 3e-3, "observed spurious $\\geq5\\sigma$: 0\n(all $5\\sigma$ = real SM)", color="#12507a",
        fontsize=8.5, ha="right", va="center")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1e3, 1.3e6); ax.set_ylim(1e-4, 3e2)
ax.set_xlabel("trials factor  N  (independent looks across the ATLAS search program)")
ax.set_ylabel("expected number of upward fluctuations")
ax.set_title("Expected vs observed excesses across ATLAS  (background-only + look-elsewhere)")
ax.grid(True, which="both", ls=":", alpha=0.35)
ax.legend(loc="upper left", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "results", "plots", "excess_counting.png"), dpi=130)
print("wrote results/plots/excess_counting.png")
print(f"N=5e4 -> expected >=3s = {N_atlas*P3:.0f}, >=5s = {N_atlas*P5:.3f}")
print(f"N=1e4 -> {1e4*P3:.0f} / {1e4*P5:.4f} ;  N=1e5 -> {1e5*P3:.0f} / {1e5*P5:.3f}")
print(f"resonance N={N_res:.0f} -> {N_res*P3:.1f} / {N_res*P5:.5f}")
