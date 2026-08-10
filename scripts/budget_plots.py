#!/usr/bin/env python3
"""Figures of the search budget: trials per spectrum, the scan-window map, and which public model
motivates which spectrum.

Inputs are the two data-free modules (bump_observables.py, public_obs_map.py) plus literature
constants. Writes results/plots/{search_budget,scan_windows,model_observable_matrix}.png.
The Z_local-vs-N waterfall is budget_waterfall.py; the tables and the report search_budget.py.
"""
import os, math, sys, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import (canon, CANON_ORDER, floor, res, scan_segments, ns_scan,
                              z_local_for_global5 as z5, fmt_range, SQRTS)
from public_obs_map import PUBLIC_OBS, nsel

OUT = os.path.join(ROOT, "results", "plots")
os.makedirs(OUT, exist_ok=True)

# ---- the public channel set (canonical observables predicted by >=1 public model)
pub_models = collections.defaultdict(set)
for m, obss in PUBLIC_OBS.items():
    for o in obss:
        pub_models[canon(o)].add(m)
order = [o for o in CANON_ORDER if o == canon(o) and o in pub_models]
order += [o for o in sorted(pub_models) if o not in order]

N_incl = sum(ns_scan(o) for o in order)
N_sel  = sum(nsel(o) * ns_scan(o) for o in order)
n_sel  = sum(nsel(o) for o in order)

# ================================================================ 1. budget bars
pv = [ns_scan(o) for o in order]
nm = [len(pub_models[o]) for o in order]
yp = np.arange(len(order))
fig, ax = plt.subplots(figsize=(9.5, 8.5))
cmap = plt.cm.viridis
norm = matplotlib.colors.Normalize(vmin=1, vmax=max(nm))
ax.barh(yp, pv, color=[cmap(norm(v)) for v in nm], edgecolor="#333333", linewidth=0.5)
for i, o in enumerate(order):
    segs = " + ".join(fmt_range(lo, hi) for lo, hi in scan_segments(o))
    ax.text(pv[i] + max(pv) * 0.01, i, segs, va="center", fontsize=7)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.045)
cb.set_label("# public models feeding the channel", fontsize=9)
ax.set_yticks(yp); ax.set_yticklabels(order, fontsize=8.5); ax.invert_yaxis()
ax.set_xlim(0, max(pv) * 1.42)
ax.set_xlabel(r"effective independent looks  $n_s = (1/r)\,\ln(M_{hi}/M_{lo})$  (published windows)")
ax.set_title("Search budget — PUBLIC information only\n"
             "channels predicted by public BSM models, published-search scan windows")
txt = (f"{len(order)} bump spectra from public model DBs\n"
       f"N_trials = {N_incl:,.0f}   (band {N_incl*0.5:,.0f}-{N_incl*2:,.0f})\n"
       f"Z$_{{local}}$ for 5$\\sigma$ global = {z5(N_incl):.2f}"
       f"   (band {z5(N_incl*2):.2f}-{z5(N_incl*0.5):.2f})\n"
       f"a local 5$\\sigma$ bump is only ~{math.sqrt(max(0.1, 25 - 2*math.log(N_incl))):.1f}$\\sigma$ global")
ax.text(0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
        bbox=dict(boxstyle="round", fc="#eef7ea", ec="#3d7a34"))
fig.tight_layout()
fig.savefig(os.path.join(OUT, "search_budget.png"), dpi=130)
print("wrote search_budget.png")

# ================================================================ 2. scan-window map
fig2, ax2 = plt.subplots(figsize=(11, 0.38 * len(order) + 1.8))
for i, o in enumerate(order):
    for lo, hi in scan_segments(o):
        ax2.plot([lo, hi], [i, i], lw=7, color="#4c78a8", solid_capstyle="butt", zorder=3)
    fl = floor(o)
    ax2.plot([fl, fl], [i - 0.32, i + 0.32], color="green", ls="--", lw=1.4,
             label="analyzable floor (trig/kin)" if i == 0 else None)
    ax2.text(SQRTS * 1.35, i, f"n_s={ns_scan(o):.0f}  r={res(o):g}", va="center", fontsize=7.5)
ax2.axvline(SQRTS, color="black", ls=":", lw=1.2, label=r"$\sqrt{s}$ = 13.6 TeV")
ax2.set_yticks(range(len(order)))
ax2.set_yticklabels([f"{o}  [{len(pub_models[o])} mod]" for o in order], fontsize=8.5)
ax2.invert_yaxis(); ax2.set_xscale("log")
ax2.set_xlim(0.2, SQRTS * 12)
ax2.grid(axis="x", ls=":", alpha=0.4)
ax2.set_xlabel("resonance mass [GeV]")
ax2.set_title("Published-search scan windows per bump spectrum  (PUBLIC information only)\n"
              "blue = scanned window (disjoint segments where the program has a gap)")
ax2.legend(loc="lower left", fontsize=8.5)
fig2.tight_layout()
fig2.savefig(os.path.join(OUT, "scan_windows.png"), dpi=130)
print("wrote scan_windows.png")

# ================================================================ 3. model x observable matrix
models = sorted(PUBLIC_OBS, key=lambda m: (-len(PUBLIC_OBS[m]), m.lower()))
M = np.zeros((len(models), len(order)))
for i, m in enumerate(models):
    for o in PUBLIC_OBS[m]:
        M[i, order.index(canon(o))] = 1
cmap3 = matplotlib.colors.ListedColormap(["#f4f4f4", "#4c78a8"])
fig3, ax3 = plt.subplots(figsize=(14.5, 11.5))
ax3.pcolormesh(M, cmap=cmap3, vmin=0, vmax=1, edgecolors="white", linewidth=2)
ax3.set_xticks(np.arange(len(order)) + 0.5)
ax3.set_xticklabels([f"{o}\n[{len(pub_models[o])} mod | n_s={ns_scan(o):.0f}]" for o in order],
                    rotation=90, fontsize=8)
ax3.set_yticks(np.arange(len(models)) + 0.5)
ax3.set_yticklabels([f"{m}  ({len(PUBLIC_OBS[m])})" for m in models], fontsize=8.5)
ax3.invert_yaxis(); ax3.tick_params(length=0)
for s in ax3.spines.values(): s.set_visible(False)
ax3.set_title("Which public BSM model motivates which bump spectrum  (PUBLIC information only)\n"
              f"{len(models)} public model classes x {len(order)} spectra;  "
              f"total budget N = {N_incl:,.0f} looks", fontsize=11)
ax3.legend(handles=[Patch(fc="#4c78a8", label="model predicts a resonance in this spectrum"),
                    Patch(fc="#f4f4f4", ec="#cccccc", label="not predicted")],
           loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=9)
fig3.tight_layout(rect=[0, 0, 0.84, 1])
fig3.savefig(os.path.join(OUT, "model_observable_matrix.png"), dpi=130)
print("wrote model_observable_matrix.png")

print(f"\nbudget: {len(order)} spectra, N_incl = {N_incl:,.0f} (Z {z5(N_incl):.2f}); "
      f"{n_sel} event selections, N = {N_sel:,.0f} (Z {z5(N_sel):.2f})")
