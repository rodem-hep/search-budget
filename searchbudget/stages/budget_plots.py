import matplotlib
import numpy as np
from matplotlib.patches import Patch

from .. import io, paths
from ..core.bump_observables import (SQRTS, canon, floor, ns_scan, scan_segments,
                                     z_local_for_global5 as z5)
from ..core.catalogue import N_inclusive, N_selections, canonical_order, models_by_spectrum
from ..core.public_obs_map import PUBLIC_OBS, nsel
from ..registry import stage
from ..viz.labels import mathify, textsafe
from ..viz.style import GRID, INK, MARK, MARK_PALE, style

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BUDGET_BAR = "#007E64"


@stage(
    name="budget-plots",
    group="budget",
    summary="the budget bar chart, the scan windows and the model-observable matrix",
    outputs=["plots/search_budget.png", "plots/scan_windows.png",
             "plots/model_observable_matrix.png"],
)
def main(options=None):
    pub_models = models_by_spectrum()
    order = canonical_order(pub_models)

    N_incl = N_inclusive(order)
    N_sel = N_selections(order)
    n_sel = sum(nsel(o) for o in order)

    ranked = sorted(order, key=lambda o: -ns_scan(o))
    fig, ax = plt.subplots(figsize=(5.7, 7.0))
    style(ax, grid=False, minor="x")
    ax.tick_params(axis="y", right=False)
    ax.grid(axis="x", color=GRID, lw=0.5, ls=":", alpha=0.9)
    ax.barh(np.arange(len(ranked)), [ns_scan(o) for o in ranked], color=BUDGET_BAR,
            edgecolor="none", height=0.72)
    ax.set_yticks(np.arange(len(ranked)))
    ax.set_yticklabels([mathify(o) for o in ranked], fontsize=7.5)
    ax.set_ylim(len(ranked) - 0.5, -0.5)
    ax.set_xlabel(r"independent looks $n_s = (1/r)\,\ln(M_{\mathrm{hi}}/M_{\mathrm{lo}})$")
    fig.tight_layout()
    io.save(fig, paths.plot("search_budget.png"), dpi=400)
    print("wrote search_budget.png")

    fig2, ax2 = plt.subplots(figsize=(5.9, 0.155 * len(order) + 1.0))
    style(ax2, grid=False, minor="x")
    ax2.tick_params(axis="y", right=False)
    ax2.grid(axis="x", color=GRID, lw=0.5, ls=":", alpha=0.9)
    for i, o in enumerate(order):
        for lo, hi in scan_segments(o):
            ax2.plot([lo, hi], [i, i], lw=4.0, color=MARK, solid_capstyle="butt", zorder=3)
        fl = floor(o)
        ax2.plot([fl, fl], [i - 0.32, i + 0.32], color=INK, ls="-", lw=1.0, zorder=4,
                 label="analysable mass floor" if i == 0 else None)
    ax2.axvline(SQRTS, color=INK, ls=":", lw=1.0, label=r"$\sqrt{s} = 13.6$ TeV")
    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels([mathify(o) for o in order], fontsize=7.5)
    ax2.set_ylim(len(order) - 0.5, -0.5)
    ax2.set_xscale("log")
    ax2.set_xlim(0.2, SQRTS * 1.6)
    ax2.set_xlabel("resonance mass [GeV]")
    ax2.legend(loc="lower left")
    fig2.tight_layout()
    io.save(fig2, paths.plot("scan_windows.png"), dpi=400)
    print("wrote scan_windows.png")

    models = sorted(PUBLIC_OBS, key=lambda m: (-len(PUBLIC_OBS[m]), m.lower()))
    M = np.zeros((len(models), len(order)))
    for i, m in enumerate(models):
        for o in PUBLIC_OBS[m]:
            M[i, order.index(canon(o))] = 1
    cmap3 = matplotlib.colors.ListedColormap([MARK_PALE, MARK])
    fig3, ax3 = plt.subplots(figsize=(6.0, 6.4))
    ax3.pcolormesh(M, cmap=cmap3, vmin=0, vmax=1, edgecolors="white", linewidth=0.7)
    ax3.set_xticks(np.arange(len(order)) + 0.5)
    ax3.set_xticklabels([mathify(o) for o in order], rotation=90, fontsize=6.0)
    ax3.set_yticks(np.arange(len(models)) + 0.5)
    ax3.set_yticklabels([textsafe(m) for m in models], fontsize=6.0)
    ax3.invert_yaxis()
    ax3.tick_params(which="both", length=0, colors=INK)
    for s in ax3.spines.values():
        s.set_visible(True)
        s.set_color(INK)
        s.set_linewidth(0.7)
    ax3.legend(handles=[Patch(fc=MARK, ec="none", label="model predicts a resonance here")],
               loc="lower left", bbox_to_anchor=(0.0, 1.005), borderaxespad=0)
    fig3.tight_layout()
    io.save(fig3, paths.plot("model_observable_matrix.png"), dpi=400)
    print("wrote model_observable_matrix.png")

    print(f"\nbudget: {len(order)} spectra, N_incl = {N_incl:,.0f} (Z {z5(N_incl):.2f}); "
          f"{n_sel} event selections, N = {N_sel:,.0f} (Z {z5(N_sel):.2f})")
