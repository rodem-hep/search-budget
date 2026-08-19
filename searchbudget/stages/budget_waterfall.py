import matplotlib
import numpy as np

from .. import io, paths
from ..core.bump_observables import z_local_for_global5 as z5
from ..registry import stage
from ..viz.style import GRID, INK2, MARK, labels, style, title

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

N_FULL = 5e4


@stage(
    name="budget-waterfall",
    group="budget",
    summary="Z_local against N across the granularity levels",
    outputs=["plots/budget_waterfall.png"],
    needs=["tables/search_budget.csv", "tables/search_budget_selections.csv"],
)
def main(options=None):
    sb = io.read_rows(paths.table("search_budget.csv"))
    N_incl = sum(float(r["ns_scan"]) for r in sb)
    n_incl = len(sb)

    se = io.read_rows(paths.table("search_budget_selections.csv"))
    N_sel = sum(float(r["ns_with_selections"]) for r in se)
    n_sel = sum(int(r["n_event_selections"]) for r in se)

    LEVELS = [
        (N_incl, f"{n_incl} inclusive spectra",         (-12,  22), "right"),
        (N_sel,  f"{n_sel} published event selections", (14,  -26), "left"),
        (N_FULL, "full ATLAS program",                  (-12,  20), "right"),
    ]

    N = np.logspace(2.8, 5.6, 400)
    fig, ax = plt.subplots(figsize=(5.5, 3.35))
    style(ax, grid=False)
    ax.grid(True, which="major", color=GRID, lw=0.8, alpha=0.6)
    ax.plot(N, [z5(v) for v in N], color=MARK, lw=1.3)
    ax.set_xscale("log")

    for Nv, name, off, ha in LEVELS:
        zv = z5(Nv)
        ax.plot([Nv], [zv], "o", ms=4.5, color=MARK, mec="white", mew=0.7, zorder=5)
        ax.annotate(f"{name}\nN = {Nv:,.0f},  Z = {zv:.2f}",
                    (Nv, zv), textcoords="offset points", xytext=off, ha=ha,
                    va="bottom" if off[1] > 0 else "top", fontsize=9.5, color=INK2,
                    arrowprops=dict(arrowstyle="-", color=GRID, lw=0.8))

    labels(ax, x="trials factor  N  (independent looks)",
           y=r"$Z_{\mathrm{local}}$ needed for a $5\sigma$ global discovery")
    title(ax, r"The look-elsewhere price grows as $\sqrt{25 + 2\ln N}$", size=12.5)
    fig.tight_layout()
    out = io.save(fig, paths.plot("budget_waterfall.png"), dpi=400)
    print(f"wrote {out}")
    for Nv, name, *_ in LEVELS:
        print(f"  {name:28s} N = {Nv:8,.0f}   Z_local = {z5(Nv):.2f}")
