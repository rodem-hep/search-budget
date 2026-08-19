import math

import matplotlib
import numpy as np
from scipy.integrate import quad
from scipy.stats import gumbel_r, kstest, norm, skew

from .. import io, paths
from ..registry import stage
from ..viz.style import BLUE, C_ALT, C_ARG, C_THR, INK, INK2, SURF, labels, panels, title

matplotlib.use("Agg")


@stage(
    name="max-of-gaussians",
    group="stats",
    summary="the maximum of many Gaussians: argmax against a fixed threshold",
    outputs=["plots/max_of_gaussians/max_of_gaussians_light.png",
             "plots/max_of_gaussians/signal_wins_the_max.png",
             "plots/max_of_gaussians/ab_confirmation.png",
             "plots/max_of_gaussians/threshold_scan.png",
             "plots/max_of_gaussians/threshold_vs_argmax.png",
             "plots/max_of_gaussians/roc_threshold_vs_argmax.png"],
)
def main(options=None):
    PLOTS = paths.GAUSSIANS

    n, nb, tB = 30_000, 29_999, 3.0
    NPE, CHUNK = 10_000, 500
    MUS = [3, 4, 5, 6]
    sf, cdf, pdf = norm.sf, norm.cdf, norm.pdf
    tstar = norm.isf(1.0 / nb)

    pdf_max = lambda x: n * pdf(x) * cdf(x) ** (n - 1)
    p_win = lambda mu: quad(lambda s: pdf(s - mu) * cdf(s) ** nb, mu - 12, mu + 12, limit=400)[0]
    y_arg = lambda mu: p_win(mu) * sf(tB - mu) + (1 - p_win(mu)) * sf(tB)
    lam = lambda t: nb * sf(t)
    p_fake = lambda t: 1 - np.exp(-lam(t) * sf(tB))
    y_thr = lambda mu, t: cdf(mu - t) * cdf(mu - tB)
    x_thr = lambda t: nb * sf(t) * sf(tB)
    x_arg = sf(tB)

    rng = np.random.default_rng(0)
    maxima = np.empty(NPE)
    one_pe = None
    for i in range(0, NPE, CHUNK):
        x = rng.standard_normal((CHUNK, n))
        if one_pe is None:
            one_pe = x[0].copy()
        maxima[i:i + CHUNK] = x.max(axis=1)

    ln_n = math.log(n)
    b_n = 1 / math.sqrt(2 * ln_n)
    a_n = math.sqrt(2 * ln_n) - (math.log(ln_n) + math.log(4 * math.pi)) / (2 * math.sqrt(2 * ln_n))
    mu_f, sd_f = maxima.mean(), maxima.std(ddof=0)

    print(f"n = {n:,}   pseudo-experiments = {NPE:,}   Gumbel a_n = {a_n:.3f}  b_n = {b_n:.3f}\n")
    print(f"{'':>10}  {'simulated':>10}  {'Gumbel':>10}")
    g = gumbel_r(loc=a_n, scale=b_n)
    print(f"{'mean':>10}  {mu_f:10.3f}  {g.mean():10.3f}")
    print(f"{'std':>10}  {sd_f:10.3f}  {g.std():10.3f}")
    print(f"{'median':>10}  {np.median(maxima):10.3f}  {g.median():10.3f}")
    print(f"{'skewness':>10}  {skew(maxima):10.2f}  {1.1395:10.2f}")
    print("\nKS tests of the simulated maxima")
    for name, cdf_fn in (("exact Phi^n", lambda x: cdf(x) ** n),
                         ("Gumbel(a_n,b_n)", lambda x: g.cdf(x)),
                         ("best-fit Gaussian", lambda x: cdf((x - mu_f) / sd_f))):
        ks = kstest(maxima, cdf_fn)
        print(f"   {name:>18}   D = {ks.statistic:.4f}   p = {ks.pvalue:.2g}")

    fig, (a0, a1) = panels(figsize=(5.9, 2.6), wspace=0.20)

    a0.hist(one_pe, bins=90, color="#c9c8c4", lw=0)
    a0.axvline(one_pe.max(), color=BLUE[5], lw=1.3)
    a0.annotate(f"maximum = {one_pe.max():.2f}",
                xy=(one_pe.max(), a0.get_ylim()[1] * 0.42), xytext=(-4.6, a0.get_ylim()[1] * 0.66),
                color=BLUE[5], ha="left", linespacing=1.45,
                arrowprops=dict(arrowstyle="->", color=BLUE[5], lw=1.2,
                                connectionstyle="arc3,rad=-0.32"))
    labels(a0, f"$x$", "draws per bin")
    title(a0, f"One pseudo-experiment: {n:,} Gaussian draws")

    xs = np.linspace(3.2, 5.6, 600)
    a1.hist(maxima, bins=70, density=True, color="#c9c8c4", lw=0,
            label=f"{NPE:,} simulated maxima")
    a1.plot(xs, pdf_max(xs), color=BLUE[5], lw=1.3, label="exact:  $n\\,\\phi(x)\\,\\Phi(x)^{n-1}$")
    a1.plot(xs, g.pdf(xs), color=C_ALT, lw=1.3, ls=(0, (5, 2.5)),
            label=f"Gumbel asymptotic ($a_n$={a_n:.2f}, $b_n$={b_n:.2f})")
    a1.text(0.60, 0.68, f"mean {mu_f:.2f}    std {sd_f:.2f}    skew {skew(maxima):.2f}",
            transform=a1.transAxes, color=INK2)
    a1.legend(frameon=False, labelcolor=INK2, loc="upper left")
    labels(a1, f"maximum of {n:,} draws", "probability density")
    title(a1, "Distribution of the maxima")
    io.save(fig, paths.gaussians("max_of_gaussians_light.png"), dpi=400, facecolor=SURF,
                bbox_inches="tight")

    fig, (b0, b1) = panels(figsize=(5.9, 2.6), wspace=0.20)

    xs = np.linspace(1.0, 8.5, 900)
    b0.fill_between(xs, nb * pdf(xs) * cdf(xs) ** (nb - 1), color="#c9c8c4", lw=0, alpha=0.85)
    b0.plot(xs, nb * pdf(xs) * cdf(xs) ** (nb - 1), color=INK2, lw=1.2)
    for m in MUS:
        b0.plot(xs, pdf(xs - m), color=BLUE[m], lw=1.3)
        b0.text(m, pdf(0) * 1.14, f"${m}\\sigma$", color=BLUE[m], ha="center",
                fontweight="medium")
    b0.text(4.45, 1.30, "background maximum\nof the other bins", color=INK2,
            linespacing=1.45)
    b0.set_ylim(0, 1.72)
    labels(b0, "$z$", "probability density")
    title(b0, "Signal bin against the background")

    mus = np.linspace(1.5, 8.0, 90)
    pw = np.array([p_win(m) for m in mus])
    b1.plot(mus, pw * 100, color=BLUE[5], lw=1.3)
    b1.axhline(50, color=INK2, lw=1, ls=(0, (4, 3)), alpha=0.8)
    for m in MUS:
        v = p_win(m) * 100
        b1.plot([m], [v], "o", ms=4, color=BLUE[m], mec=SURF, mew=0.7, zorder=5)
        b1.text(m + 0.12, v - 4.5, f"{v:.0f}%", color=BLUE[m], fontweight="medium")
    b1.set_ylim(0, 104)
    labels(b1, "signal strength  $\\mu$   [$\\sigma$]", "signal bin is the maximum   [%]")
    title(b1, "$p(\\mu) = P(\\,$signal wins the scan$\\,)$")
    io.save(fig, paths.gaussians("signal_wins_the_max.png"), dpi=400, facecolor=SURF,
                bbox_inches="tight")

    fig, (c0, c1) = panels(figsize=(5.9, 2.6), wspace=0.20)

    w, xpos = 0.36, np.arange(len(MUS))
    s1 = np.array([p_win(m) for m in MUS]) * 100
    s2 = np.array([y_arg(m) for m in MUS]) * 100
    c0.bar(xpos - w / 2, s1, w, color=BLUE[3], lw=0, label="stage 1: signal bin wins the scan in A")
    c0.bar(xpos + w / 2, s2, w, color=BLUE[5], lw=0,
           label="stage 2: and $z_B > 3$  $\\rightarrow$  CONFIRMED")
    for i, (u, v) in enumerate(zip(s1, s2)):
        c0.text(i - w / 2, u + 2.5, f"{u:.0f}%", ha="center", color=INK)
        c0.text(i + w / 2, v + 2.5, f"{v:.0f}%", ha="center", color=BLUE[5])
    c0.set_xticks(xpos)
    c0.set_xticklabels([f"${m}\\sigma$" for m in MUS])
    c0.set_ylim(0, 132)
    c0.legend(frameon=False, labelcolor=INK2, loc="upper left")
    labels(c0, "signal strength  $\\mu$", "probability  [%]")
    title(c0, "Two-stage A/B confirmation")

    mus = np.linspace(1.5, 8.5, 100)
    c1.plot(mus, [p_win(m) * 100 for m in mus], color=BLUE[3], lw=1.3)
    c1.plot(mus, [y_arg(m) * 100 for m in mus], color=BLUE[5], lw=1.3)
    for m in MUS:
        v = y_arg(m) * 100
        c1.plot([m], [v], "o", ms=4, color=BLUE[5], mec=SURF, mew=0.7, zorder=5)
        c1.text(m - 0.15, v + 4, f"{v:.0f}%", color=BLUE[5], ha="right")
    c1.axhline(sf(tB) * 100, color=C_ARG, lw=1.2)
    c1.text(8.4, 22, r"false-alarm rate 0.135%  =  $3.0\sigma$ global", color=C_ARG,
            ha="right")
    c1.text(6.15, 84, "wins stage 1", color=BLUE[3])
    c1.text(6.15, 66, "confirmed in both stages", color=BLUE[5])
    c1.set_ylim(-4, 112)
    labels(c1, "signal strength  $\\mu$   [$\\sigma$]", "probability  [%]")
    title(c1, "Confirmation power")
    io.save(fig, paths.gaussians("ab_confirmation.png"), dpi=400, facecolor=SURF,
                bbox_inches="tight")

    fig, (d0, d1) = panels(figsize=(5.9, 2.6), wspace=0.22)

    ts = np.linspace(3.0, 6.0, 600)
    d0.semilogy(ts, lam(ts), color=C_THR, lw=1.3)
    d0.axhline(1, color=INK2, lw=1, ls=(0, (4, 3)))
    d0.axvline(tstar, color=INK2, lw=1, ls=(0, (2, 3)))
    d0.plot([tstar], [1], "o", ms=4.5, color=C_THR, mec=SURF, mew=0.7, zorder=6)
    d0.text(tstar + 0.07, 3.0, f"$t^\\star = {tstar:.3f}$\n$\\lambda(t^\\star) = 1$", color=INK,
            linespacing=1.45)
    d0.text(3.05, 2.5e-2, "$\\lambda(t) = (n-1)\\,[1-\\Phi(t)]$", color=C_THR)
    labels(d0, "selection threshold  $t$", "expected background bins above $t$")
    title(d0, "Expected background bins above $t$")

    for m in MUS:
        d1.plot(ts, cdf(m - ts) * 100, color=BLUE[m], lw=1.3)
        d1.plot([m], [50], "o", ms=7, color=BLUE[m], mec=SURF, mew=0.7, zorder=6)
        d1.text(m - 0.07, 54, f"$\\mu={m}\\sigma$", color=BLUE[m], ha="right",
                fontweight="medium")
    d1.axhline(50, color=INK2, lw=1, ls=(0, (4, 3)), alpha=0.8)
    d1.axvline(tstar, color=INK2, lw=1, ls=(0, (2, 3)))
    d1.text(tstar + 0.06, 25, "$t^\\star$", color=INK, fontweight="medium")
    d1.set_ylim(0, 119)
    labels(d1, "selection threshold  $t$", "signal bin passes stage 1   [%]")
    title(d1, "Chance the signal bin passes")
    io.save(fig, paths.gaussians("threshold_scan.png"), dpi=400, facecolor=SURF,
                bbox_inches="tight")

    fig, (e0, e1) = panels(figsize=(5.9, 2.6), wspace=0.24)

    for m in MUS:
        e0.plot(ts, y_thr(m, ts) * 100, color=BLUE[m], lw=1.3)
        pa = y_arg(m) * 100
        e0.plot([3.0, 6.0], [pa, pa], color=BLUE[m], lw=1.2, ls=(0, (2, 2.5)), alpha=0.85)
        e0.text(6.05, pa, f"${m}\\sigma$", color=BLUE[m], va="center",
                fontweight="medium")
    e0.axvline(tstar, color=INK2, lw=1, ls=(0, (2, 3)))
    e0.text(tstar + 0.06, 14, "$t^\\star$", color=INK, fontweight="medium")
    e0.text(3.05, 112, "solid: threshold at $t$        dashed: argmax rule",
            color=INK2)
    e0.set_xlim(3.0, 6.35)
    e0.set_ylim(0, 124)
    labels(e0, "selection threshold  $t$", "signal confirmed in B   [%]")
    title(e0, "Confirmation power")

    e1.semilogy(ts, p_fake(ts), color=C_ARG, lw=1.3)
    e1.axhline(sf(tB), color=INK2, lw=1.2, ls=(0, (4, 3)))
    e1.axvline(tstar, color=INK2, lw=1, ls=(0, (2, 3)))
    e1.plot([tstar], [p_fake(tstar)], "o", ms=4.5, color=C_ARG, mec=SURF, mew=0.7, zorder=6)
    e1.text(3.05, 2.1e-3, r"argmax rule:  $1.35\times10^{-3}$", color=INK2)
    e1.text(4.32, 8e-3, f"$t^\\star = {tstar:.2f}$", color=INK)
    e1.text(3.05, 3e-6, "$P_{\\rm fake} = 1 - e^{-\\lambda(t)\\,[1-\\Phi(3)]}$", color=C_ARG)
    e1.set_ylim(1e-6, 8e-2)
    labels(e1, "selection threshold  $t$", "P(at least one false confirmation)")
    title(e1, "Global false-alarm rate it buys")
    io.save(fig, paths.gaussians("threshold_vs_argmax.png"), dpi=400, facecolor=SURF,
                bbox_inches="tight")

    fig, (f0, f1) = panels(figsize=(5.9, 2.7), wspace=0.24, width_ratios=[1.12, 1])

    tr = np.linspace(2.6, 6.6, 800)
    f0.axvline(x_arg, color=C_ARG, lw=1.3, ls=(0, (3, 3)), alpha=0.9)
    for m in MUS:
        f0.semilogx(x_thr(tr), y_thr(m, tr) * 100, color=BLUE[m], lw=1.3)
        f0.plot([x_arg], [y_arg(m) * 100], marker="D", ms=4.5, color=C_ARG, mec=SURF, mew=0.7, zorder=6)
        f0.plot([x_thr(tstar)], [y_thr(m, tstar) * 100], "o", ms=4, color=BLUE[m], mec=SURF, mew=0.7,
                zorder=6)
        f0.text(0.42, {3: 32, 4: 76, 5: 93, 6: 106}[m], f"$\\mu={m}\\sigma$", color=BLUE[m],
                va="center", fontweight="medium")
    for tt in (3, 5, 6):
        f0.plot([x_thr(tt)], [y_thr(5, tt) * 100], marker="|", ms=11, color=BLUE[5], mew=0.7, zorder=5)
        f0.text(x_thr(tt), y_thr(5, tt) * 100 - 8, f"$t={tt}$", color=BLUE[5],
                ha="center")
    f0.text(2.2e-8, 122, "◆  argmax rule", color=C_ARG, fontweight="medium")
    f0.text(2.2e-8, 112, "—  threshold at $t$    ●  $t^\\star$", color=BLUE[5],
            fontweight="medium")
    f0.text(x_arg * 0.75, 4, "argmax budget", color=C_ARG, fontsize=9.5, ha="right")
    f0.set_xlim(1.5e-8, 1.0)
    f0.set_ylim(0, 132)
    labels(f0, "average number of false positives   $E[K_{\\rm fake}] = \\lambda(t)\\,[1-\\Phi(3)]$",
           "signal confirmed in B   [%]")
    title(f0, "Power against false positives", 12.8)

    for m in MUS:
        f1.semilogy(ts, y_thr(m, ts) / x_thr(ts), color=BLUE[m], lw=1.3)
        f1.axhline(y_arg(m) / x_arg, color=BLUE[m], lw=1.1, ls=(0, (2, 2.5)), alpha=0.9)
        f1.plot([tstar], [y_thr(m, tstar) / x_thr(tstar)], "o", ms=7, color=BLUE[m], mec=SURF,
                mew=1.6, zorder=6)
        f1.text(6.05, y_thr(m, 6.0) / x_thr(6.0), f"${m}\\sigma$", color=BLUE[m],
                va="center", fontweight="medium")
    f1.axvline(tstar, color=INK2, lw=1, ls=(0, (2, 3)))
    f1.text(tstar + 0.06, 2.2, "$t^\\star$", color=INK, fontweight="medium")

    f1.text(6.4, 5.0, "dashed: argmax", color=INK2, ha="right")
    f1.set_xlim(3.0, 6.5)
    f1.set_ylim(3.5, 7e7)
    labels(f1, "selection threshold  $t$", "power per false positive   $y/x$")
    title(f1, "Power per false positive", 12.8)
    io.save(fig, paths.gaussians("roc_threshold_vs_argmax.png"), dpi=400, facecolor=SURF,
                bbox_inches="tight")

    print("\nwrote 6 figures to results/plots/max_of_gaussians/")
