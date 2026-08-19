import math

import matplotlib
import matplotlib.ticker as mticker
import numpy as np

from .. import io, paths
from ..registry import stage
from ..stats import spectra as toys
from ..viz.style import C_ALT, C_ARG, C_BKG, size, style

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@stage(
    name="ab-spurious-guard",
    group="ab",
    summary="the split as a spurious-signal guard, on two toy spectra",
    outputs=["plots/ab_spurious_guard.png", "tables/ab_guard_toys.csv"],
)
def main(options=None):
    OUT = paths.plot("ab_spurious_guard.png")

    SIGMA_REL, F_ILL, Z_CUT = 0.05, 0.25, 3.0
    edges, ctr, width = toys.grid(200, 4000, 220)
    bkg_full = toys.background(ctr, width)

    seed, zA, zB = toys.first_toy_above(Z_CUT, 0.2, ctr, bkg_full, F_ILL, SIGMA_REL)
    zA_sig, zB_sig = toys.split_toy(998, ctr, bkg_full, F_ILL, SIGMA_REL,
                                    sig_mass=1200.0, sig_zfull=7.0)

    fig, axs = plt.subplots(1, 2, figsize=size(0.92, 0.40), sharey=True,
                            gridspec_kw=dict(wspace=0.10))
    panels = [(axs[0], zA, zB, "background-only", "dies in B", C_BKG),
              (axs[1], zA_sig, zB_sig, "injected signal", "confirmed in B", C_ARG)]
    for ax, za, zb, tag, verdict, vcol in panels:
        ax.set_xscale("log")
        style(ax, minor="y")
        sel = za >= Z_CUT
        first = True
        for j in np.where(sel)[0]:
            ax.axvspan(ctr[j] * (1 - 2 * SIGMA_REL), ctr[j] * (1 + 2 * SIGMA_REL),
                       color=C_ALT, alpha=0.22, lw=0, zorder=0,
                       label="pre-registered window" if first else None)
            first = False
        ax.plot(ctr, za, lw=1.1, color=C_BKG, label="stage A scan (25% of the data)")
        inwin = np.zeros_like(sel)
        for j in np.where(sel)[0]:
            inwin |= np.abs(np.log(ctr / ctr[j])) < 2 * SIGMA_REL
        ax.plot(ctr, np.where(inwin, zb, np.nan), lw=1.6, color=C_ARG,
                label="stage B, opened only here (75%)")
        k = max(1, np.count_nonzero(np.diff(np.where(sel)[0]) > 1) + 1) if sel.any() else 1
        thr = math.sqrt(25.0 + 2.0 * math.log(k))
        ax.axhline(Z_CUT, color=C_BKG, ls="--", lw=0.8)
        ax.axhline(thr, color=C_ARG, ls="--", lw=0.8)
        ax.text(0.03, 0.95, tag, transform=ax.transAxes, ha="left", va="top", fontsize=9)
        ax.text(0.97, 0.95, verdict, transform=ax.transAxes, ha="right", va="top",
                fontsize=9.5, style="italic", color=vcol)
        ax.set_xlim(200, 4000)
        ax.set_xlabel("mass [GeV]")
        ax.xaxis.set_major_locator(mticker.FixedLocator([300, 1000, 3000]))
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    axs[1].text(212, Z_CUT + 0.18, r"$Z_{\mathrm{cut}} = 3$", color=C_BKG, fontsize=8)
    axs[1].text(212, thr + 0.18, f"claim bar = {thr:.1f}", color=C_ARG, fontsize=8)
    axs[0].set_ylabel(r"local significance $Z(m)$")
    axs[0].set_ylim(-4.2, 7.9)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.01), ncol=3,
               fontsize=8.5, handlelength=1.9, columnspacing=2.2)
    io.save(fig, OUT)
    print(f"wrote {OUT}")
    print(f"  bkg toy seed {seed}: max Z_A = {zA.max():.2f}, "
          f"Z_B in window = {np.nanmax(np.where(zA >= Z_CUT, zB, np.nan)):.2f}")
    print(f"  signal toy: Z_A = {zA_sig.max():.2f} -> Z_B = {zB_sig.max():.2f}")

    zb_bkg = float(np.nanmax(np.where(zA >= Z_CUT, zB, np.nan)))
    io.write_rows(
        paths.table("ab_guard_toys.csv"),
        ["toy", "Z_cut", "claim_bar", "Z_A", "Z_B", "outcome"],
        [["background only", f"{Z_CUT:g}", f"{thr:.2f}", f"{zA.max():.2f}", f"{zb_bkg:.2f}",
          "pre-registered, dies in stage B"],
         ["injected Z_full = 7", f"{Z_CUT:g}", f"{thr:.2f}", f"{zA_sig.max():.2f}",
          f"{zB_sig.max():.2f}", "confirms" if zB_sig.max() >= thr else "fails to confirm"]])
    print("wrote results/tables/ab_guard_toys.csv")
