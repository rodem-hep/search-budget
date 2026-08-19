import matplotlib
import numpy as np

from .. import io, paths
from ..core.lee import p1
from ..registry import stage
from ..viz.style import C_ARG, GRID, INK2, MARK, labels, style, title

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

N_ATLAS = 5.0e4


@stage(
    name="excess-counting",
    group="budget",
    summary="the background-only excesses each trials count predicts",
    outputs=["plots/excess_counting.png"],
    needs=["tables/search_budget.csv", "tables/census_budget.csv"],
)
def main(options=None):
    P3, P5 = p1(3.0), p1(5.0)
    try:
        N_res = sum(float(r["ns_scan"]) for r in io.read_rows(paths.table("search_budget.csv")))
    except Exception:
        N_res = 3.7e3
    try:
        N_cen = sum(float(r["n_s"]) for r in io.read_rows(paths.table("census_budget.csv")))
    except Exception:
        N_cen = 7.7e3

    N = np.logspace(3, 6, 400)
    exp3 = N * P3
    exp5 = N * P5

    fig, ax = plt.subplots(figsize=(5.3, 3.35))
    style(ax, grid=False)
    ax.grid(True, which="major", color=GRID, lw=0.8, alpha=0.6)
    ax.plot(N, exp3, color=C_ARG, lw=1.3, label=r"expected $\geq3\sigma$  ($N\cdot p_{3\sigma}$)")
    ax.plot(N, exp5, color=MARK, lw=1.3, ls="--",
            label=r"expected $\geq5\sigma$  ($N\cdot p_{5\sigma}$)")
    ax.axhspan(10, 70, color=C_ARG, alpha=0.09, lw=0,
               label=r"observed $\geq3\sigma$: order tens, none confirmed")

    ax.axvline(N_res, color=INK2, ls=":", lw=1.2)
    ax.axvline(N_cen, color=INK2, ls=":", lw=1.2)
    for Nv in (N_cen, N_ATLAS):
        ax.plot([Nv], [Nv*P3], "o", ms=4, color=C_ARG, mec="white", mew=1.6, zorder=5)
        ax.plot([Nv], [Nv*P5], "o", ms=4, color=MARK, mec="white", mew=1.6, zorder=5)
    ax.annotate("the model space", (N_res, 1.5e-4), textcoords="offset points", xytext=(-6, 0),
                color=INK2, fontsize=9.5, va="bottom", ha="right")
    ax.annotate("the published record", (N_cen, 1.5e-4), textcoords="offset points", xytext=(6, 0),
                color=INK2, fontsize=9.5, va="bottom")
    ax.annotate("full ATLAS program", (N_ATLAS, N_ATLAS*P5), textcoords="offset points",
                xytext=(10, -4), color=INK2, fontsize=9.5, va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e3, 1.3e6)
    ax.set_ylim(1e-4, 3e2)
    labels(ax, x="trials factor  N  (independent looks)",
           y="expected number of upward fluctuations")
    title(ax, "Background-only excesses the search program should produce", size=12.5)
    ax.legend(loc="upper left", bbox_to_anchor=(0.55, 0.72), fontsize=9.5, frameon=False,
              labelcolor=INK2)
    fig.tight_layout()
    io.save(fig, paths.plot("excess_counting.png"), dpi=400)
    print("wrote results/plots/excess_counting.png")
    print(f"N=5e4 -> expected >=3s = {N_ATLAS*P3:.0f}, >=5s = {N_ATLAS*P5:.3f}")
    print(f"N=1e4 -> {1e4*P3:.0f} / {1e4*P5:.4f} ;  N=1e5 -> {1e5*P3:.0f} / {1e5*P5:.3f}")
    print(f"model space N={N_res:.0f} -> {N_res*P3:.1f} / {N_res*P5:.5f}")
    print(f"published record N={N_cen:.0f} -> {N_cen*P3:.1f} / {N_cen*P5:.5f}")
