#!/usr/bin/env python3
"""Two toy spectra showing the A/B split acting as a spurious-signal guard.

Same falling-spectrum toy as ab_split_toys.py (spectrum section), shown as a
two-panel figure: a background-only fluctuation is pre-registered in stage A and
dies in stage B; an injected Z_full = 7 signal is pre-registered and confirmed.

Public inputs only. Writes results/plots/ab_spurious_guard.png.
"""
import os, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.insert(0, os.path.join(ROOT, "scripts"))
from plot_style import style, size, C_BKG, C_ARG, C_ALT, INK

OUT = os.path.join(ROOT, "results", "plots", "ab_spurious_guard.png")

SIGMA_REL, F_ILL, Z_CUT = 0.05, 0.25, 3.0
edges = np.geomspace(200, 4000, 220)
ctr = np.sqrt(edges[:-1] * edges[1:])
bkg_full = 3e6 * np.exp(-ctr / 350.0) * np.diff(edges) / ctr


def scan(counts, expect):
    z = np.zeros_like(ctr)
    for j, m in enumerate(ctr):
        w = np.abs(ctr - m) < SIGMA_REL * m
        n, b = counts[w].sum(), expect[w].sum()
        z[j] = (n - b) / math.sqrt(b) if b > 0 else 0.0
    return z


def split_toy(seed, sig_mass=None, sig_zfull=0.0):
    r = np.random.default_rng(seed)
    mu = np.zeros_like(ctr)
    if sig_mass:
        w = np.abs(ctr - sig_mass) < SIGMA_REL * sig_mass
        gauss = np.exp(-0.5 * ((ctr - sig_mass) / (SIGMA_REL * sig_mass)) ** 2)
        mu = gauss / gauss[w].sum() * sig_zfull * math.sqrt(bkg_full[w].sum())
    nA = r.poisson(F_ILL * (bkg_full + mu))
    nB = r.poisson((1 - F_ILL) * (bkg_full + mu))
    return scan(nA, F_ILL * bkg_full), scan(nB, (1 - F_ILL) * bkg_full)


seed = 0
while True:
    zA, zB = split_toy(seed)
    if zA.max() >= Z_CUT + 0.2:
        break
    seed += 1
zA_sig, zB_sig = split_toy(998, sig_mass=1200.0, sig_zfull=7.0)

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
fig.savefig(OUT)
print(f"wrote {OUT}")
print(f"  bkg toy seed {seed}: max Z_A = {zA.max():.2f}, "
      f"Z_B in window = {np.nanmax(np.where(zA >= Z_CUT, zB, np.nan)):.2f}")
print(f"  signal toy: Z_A = {zA_sig.max():.2f} -> Z_B = {zB_sig.max():.2f}")
