import math

import matplotlib
import numpy as np
from scipy.special import ndtri

from .. import io, paths
from ..core.bump_observables import ns_scan, z_local_for_global5 as z5
from ..core.catalogue import sorted_spectra
from ..core.lee import Phi, p1, phi_inv as PhiInv
from ..core.public_obs_map import nsel
from ..registry import stage
from ..viz.style import BLUE, C_ARG, INK, style

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@stage(
    name="ab-split-budget",
    group="ab",
    summary="the two-stage design priced on each basis: reach, claim bar and break-even",
    outputs=["plots/ab_split_reach.png", "plots/ab_split_crossover.png",
             "tables/ab_split_scan.csv"],
    needs=["tables/scaled_scan.csv"],
)
def main(options=None):
    obs = sorted_spectra()
    N_incl = sum(ns_scan(o) for o in obs)
    N_sel  = sum(nsel(o) * ns_scan(o) for o in obs)

    def k_eff(N, zcut): return N * p1(zcut) + 1.0

    def zB_req(N, zcut, widen=3.0, zglob=5.0):
        return math.sqrt(zglob * zglob + 2.0 * math.log(widen * k_eff(N, zcut)))

    def reach_median(f, N, zcut, widen=3.0):
        return max(zcut / math.sqrt(f), zB_req(N, zcut, widen) / math.sqrt(1.0 - f))

    def power(mu, f, N, zcut, widen=3.0, zglob=5.0):
        return Phi(math.sqrt(f) * mu - zcut) * Phi(math.sqrt(1 - f) * mu - zB_req(N, zcut, widen, zglob))

    def reach(f, N, zcut, target=0.5, widen=3.0, zglob=5.0):
        lo, hi = 1.0, 20.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if power(mid, f, N, zcut, widen, zglob) < target: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)

    N = N_sel
    Z_single = z5(N)
    print(f"budget (public, event-selection level): N = {N:,.0f}   "
          f"single-stage 50%-power reach = {Z_single:.2f}\n")
    print(f"{'Z_cut':>5} {'k_bkg':>6} {'Z_B req':>8} {'best f':>7} {'reach(50% pow)':>14} "
          f"{'median-arith':>12}  vs single-stage")
    opt = None
    for zcut in (2.0, 2.5, 3.0, 3.5, 4.0):
        fs = np.linspace(0.02, 0.98, 481)
        rs = [reach(f, N, zcut) for f in fs]
        i = int(np.argmin(rs))
        print(f"{zcut:5.1f} {N*p1(zcut):6.1f} {zB_req(N, zcut):8.2f} {fs[i]:7.2f} {rs[i]:14.2f} "
              f"{reach_median(fs[i], N, zcut):12.2f}  ({rs[i]-Z_single:+.2f})")
        if opt is None or rs[i] < opt[2]: opt = (zcut, fs[i], rs[i])
    print(f"\n50/50 split, Z_cut=3: 50%-power reach = {reach(0.5, N, 3.0):.2f}  "
          f"({reach(0.5, N, 3.0)-Z_single:+.2f} vs single-stage {Z_single:.2f})")
    print(f"optimum: Z_cut={opt[0]:.1f}, f={opt[1]:.2f} -> reach {opt[2]:.2f} "
          f"({opt[2]-Z_single:+.2f})  [two-stage always costs; the price of countable trials]")

    def reach_topk(k, widen=3.0):
        zcut = PhiInv(1.0 - k / N)
        zB = math.sqrt(25.0 + 2.0 * math.log(widen * k))
        best = None
        for fi in range(2, 97):
            f = fi / 100.0
            lo, hi = 3.0, 25.0
            for _ in range(60):
                mu = 0.5 * (lo + hi)
                if Phi(math.sqrt(f) * mu - zcut) * Phi(math.sqrt(1 - f) * mu - zB) < 0.5: lo = mu
                else: hi = mu
            if best is None or mu < best[0]: best = (mu, f, zcut, zB)
        return best

    print(f"\nreach vs number of unblinded B regions (top-k formulation, f re-optimized):")
    print(f"{'k':>6} {'eff Z_cut':>9} {'Z_B bar':>8} {'best f':>7} {'reach':>6}  vs single")
    for k in (1, 3, 10, 30, 100, 300, 1000):
        mu, f, zc, zb = reach_topk(k)
        print(f"{k:6d} {zc:9.2f} {zb:8.2f} {f:7.2f} {mu:6.2f}  {mu-Z_single:+.2f}")

    Z_single3 = math.sqrt(9.0 + 2.0 * math.log(N))
    print(f"\nB criterion 'globally > 3 sigma' (evidence trigger; same A selection, w = 3):")
    print(f"single-stage 3sigma-global reference: Z_local = {Z_single3:.2f}")
    print(f"{'Z_cut':>5} {'k_bkg':>6} {'Z_B local bar':>13} {'best f':>7} {'reach(50%)':>10}  vs single-3sig")
    for zcut in (2.5, 3.0, 3.5):
        fs = np.linspace(0.02, 0.98, 481)
        rs = [reach(f, N, zcut, zglob=3.0) for f in fs]
        i = int(np.argmin(rs))
        print(f"{zcut:5.1f} {N*p1(zcut):6.1f} {zB_req(N, zcut, zglob=3.0):13.2f} {fs[i]:7.2f} "
              f"{rs[i]:10.2f}  ({rs[i]-Z_single3:+.2f})")
    print(f"ladder (Z_cut = 3, k_eff = {k_eff(N,3.0):.1f}): evidence Z_glob>=3 -> Z_B,local >= "
          f"{zB_req(N,3.0,zglob=3.0):.2f};  discovery Z_glob>=5 -> Z_B,local >= {zB_req(N,3.0):.2f}")
    print(f"false-evidence rate per B opening (background-only, by construction): "
          f"p1(3) = {p1(3.0):.2e}")

    print(f"\nBH-FDR flagging in A (single isolated signal; B evidence bar Z_glob >= 3, w = 3):")
    r_fix = min(reach(f, N, 3.0, zglob=3.0) for f in np.linspace(0.02, 0.98, 481))
    print(f"fixed-Z_cut=3 reference reach: {r_fix:.2f}")
    print(f"{'q':>6} {'A bar z(q/N)':>12} {'P(open B|bkg)':>13} {'Z_B bar (k~1)':>13} {'best f':>7} "
          f"{'reach':>6}  vs fixed cut")
    for q in (0.01, 0.05, 0.10, 0.25):
        zA = PhiInv(1.0 - q / N)
        zb = math.sqrt(9.0 + 2.0 * math.log(3.0 * (1.0 + q)))
        best = None
        for f in np.linspace(0.02, 0.98, 481):
            lo, hi = 1.0, 20.0
            for _ in range(50):
                mu = 0.5 * (lo + hi)
                if Phi(math.sqrt(f) * mu - zA) * Phi(math.sqrt(1 - f) * mu - zb) < 0.5: lo = mu
                else: hi = mu
            if best is None or mu < best[0]: best = (mu, f)
        print(f"{q:6.2f} {zA:12.2f} {q:13.2f} {zb:13.2f} {best[1]:7.2f} {best[0]:6.2f}  "
              f"({best[0]-r_fix:+.2f})")
    print("multi-channel adaptivity: with c channels above the bar the BH threshold is z(q*c/N):")
    for c in (1, 2, 3, 5):
        print(f"  c = {c}: z_A >= {PhiInv(1.0 - 0.10 * c / N):.2f}   (q = 0.10)")


    def reach2_opt(Nv, widen=3.0):
        best = (1e9, None, None)
        for zc in np.arange(1.0, 5.751, 0.25):
            for f in np.arange(0.04, 0.965, 0.04):
                r = reach(f, Nv, zc, widen=widen)
                if r < best[0]: best = (r, zc, f)
        _, zc0, f0 = best
        for zc in np.arange(max(0.5, zc0 - 0.3), zc0 + 0.31, 0.05):
            for f in np.arange(max(0.02, f0 - 0.06), min(0.98, f0 + 0.061), 0.01):
                r = reach(f, Nv, zc, widen=widen)
                if r < best[0]: best = (r, zc, f)
        return best

    def reach_median_opt(Nv, widen=1.0):
        best = 1e9
        for zc in np.arange(0.5, 6.01, 0.05):
            zb = zB_req(Nv, zc, widen)
            best = min(best, math.sqrt(zc * zc + zb * zb))
        return best

    print(f"\ncrossover scan: fully optimized two-stage reach vs single-stage, as a function of N")
    print(f"{'N':>12} {'single Z1':>9} {'med2(w=1)':>9} {'joint(w=1)':>10} {'joint(w=3)':>10} "
          f"{'cost(w=3)':>9} {'R*':>8}")
    cross = []
    for Nv in [10.0**e for e in range(1, 11)]:
        Z1 = z5(Nv)
        m1 = reach_median_opt(Nv, 1.0)
        r1 = reach2_opt(Nv, 1.0)[0]
        r3 = reach2_opt(Nv, 3.0)[0]
        Rstar = math.exp(0.5 * (r3 * r3 - Z1 * Z1))
        cross.append((Nv, Z1, r1, r3, Rstar))
        print(f"{Nv:12,.0f} {Z1:9.2f} {m1:9.2f} {r1:10.2f} {r3:10.2f} {r3-Z1:9.2f} {Rstar:8.0f}")
    print("median arithmetic (w=1) tracks the single stage at every N (trials factorization);")
    print("the joint-power cost never crosses zero -> the split never wins on raw sensitivity.")
    print(f"it wins iff the defendable single-stage trials count exceeds R* x N_true "
          f"(R* = {math.exp(0.5*(reach2_opt(N,3.0)[0]**2 - z5(N)**2)):.0f} at N = {N:,.0f}).")


    def k_sym(Nv, zcut): return 2.0 * Nv * p1(zcut) + 1.0

    def power_sym(mu, f, Nv, zcut, widen=3.0, zglob=5.0):
        zr = math.sqrt(zglob * zglob + 2.0 * math.log(widen * k_sym(Nv, zcut)))
        ac, bc = Phi(math.sqrt(f) * mu - zcut), Phi(math.sqrt(1 - f) * mu - zcut)
        ar, br = Phi(math.sqrt(f) * mu - zr),   Phi(math.sqrt(1 - f) * mu - zr)
        return ac * br + bc * ar - ar * br

    def reach_sym(f, Nv, zcut, target=0.5, widen=3.0, zglob=5.0):
        lo, hi = 1.0, 25.0
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            if power_sym(mid, f, Nv, zcut, widen, zglob) < target: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)

    def opt_box(fn, Nv, widen=3.0, zlo=2.0, zhi=4.5):
        best = (1e9, None, None)
        for zc in np.arange(zlo, zhi + 1e-9, 0.05):
            for f in np.arange(0.04, 0.965, 0.01):
                r = fn(f, Nv, zc, widen=widen)
                if r < best[0]: best = (r, zc, f)
        return best

    fgrid = np.linspace(0.02, 0.98, 481)
    print(f"\nsymmetrized swap: run A->B AND B->A, claim if EITHER confirms (union rule)")
    print(f"{'Z_cut':>5} | {'1-way f*':>8} {'reach':>6} {'cost':>6} | {'swap f*':>8} {'reach':>6} "
          f"{'cost':>6} | swap gain")
    for zcut in (2.0, 2.5, 3.0, 3.5, 4.0):
        r1 = [reach(f, N, zcut) for f in fgrid];     i1 = int(np.argmin(r1))
        r2 = [reach_sym(f, N, zcut) for f in fgrid]; i2 = int(np.argmin(r2))
        print(f"{zcut:5.1f} | {fgrid[i1]:8.2f} {r1[i1]:6.2f} {r1[i1]-Z_single:+6.2f} | "
              f"{fgrid[i2]:8.2f} {r2[i2]:6.2f} {r2[i2]-Z_single:+6.2f} | {r1[i1]-r2[i2]:+.2f}")
    print(f"the swap moves the optimum back to f = 1/2: the second direction only pays when both "
          f"directions have comparable power.")
    for f, zc in ((0.50, 3.0), (0.30, 3.0)):
        print(f"  f = {f:.2f}, Z_cut = {zc:g}:  one-way {reach(f, N, zc):.2f} "
              f"({reach(f, N, zc)-Z_single:+.2f})   swapped {reach_sym(f, N, zc):.2f} "
              f"({reach_sym(f, N, zc)-Z_single:+.2f})")
    b1, b2 = opt_box(reach, N), opt_box(reach_sym, N)
    print(f"best design in the practical box (Z_cut = 2-4.5): one-way {b1[0]:.2f} "
          f"({b1[0]-Z_single:+.2f}) at f={b1[2]:.2f}, Z_cut={b1[1]:.2f}   swapped {b2[0]:.2f} "
          f"({b2[0]-Z_single:+.2f}) at f={b2[2]:.2f}, Z_cut={b2[1]:.2f}")

    print(f"\nhow well the price of the split is known: best design in the box, over the N band and w")
    print(f"{'N':>12} {'w':>4} {'single':>7} {'split':>7} {'cost':>7}   {'Z_cut':>5} {'f':>5}")
    costs = []
    for Nv, tag in ((N * 0.5, "N x0.5"), (N, "nominal"), (N * 2.0, "N x2")):
        for w in (1.0, 3.0):
            r, zc, f = opt_box(reach, Nv, widen=w)
            costs.append(r - z5(Nv))
            print(f"{Nv:12,.0f} {w:4.0f} {z5(Nv):7.2f} {r:7.2f} {r - z5(Nv):+7.2f}   {zc:5.2f} "
                  f"{f:5.2f}   {tag}")
    print(f"the price of the split is {min(costs):+.2f} to {max(costs):+.2f} sigma over that box: half a "
          f"sigma is the round number, and it is known to {0.5*(max(costs)-min(costs)):.2f}")
    print(f"\nbest-design gap vs N (both schemes optimized in the same box; neither beats single-stage)")
    print(f"{'N':>12} {'single':>7} {'1-way':>7} {'swap':>7} {'gap':>7}")
    gaps = []
    for Nv in [10.0**e for e in range(1, 11)]:
        Z1, a, b = z5(Nv), opt_box(reach, Nv), opt_box(reach_sym, Nv)
        gaps.append(b[0] - a[0])
        print(f"{Nv:12,.0f} {Z1:7.2f} {a[0]:7.2f} {b[0]:7.2f} {b[0]-a[0]:+7.2f}")
    print(f"gap stays within {max(abs(g) for g in gaps):.2f} sigma over nine decades and changes sign "
          f"with the design box: at best-design level the two")
    print("schemes are a WASH -- there is no clean crossover in N, only a swing in which (f, Z_cut) "
          "the optimizer prefers.")
    print("Both stay 0.3-0.5 above the single stage: the staircase fits the Neyman-Pearson half-plane "
          "better than the corner, never equals it.")

    figc, (axr, axg) = plt.subplots(2, 1, figsize=(5.6, 4.7), sharex=True,
                                    gridspec_kw={"height_ratios": [2.0, 1.15], "hspace": 0.07})
    style(axr); style(axg)
    Nx = [c[0] for c in cross]
    axr.plot(Nx, [c[1] for c in cross], "--", color=C_ARG,
             label=r"single-stage exact correction, $\sqrt{25+2\ln N}$")
    axr.plot(Nx, [c[2] for c in cross], "-o", color=BLUE[4], ms=3,
             label=r"two-stage, optimised, exact windows ($w=1$)")
    axr.plot(Nx, [c[3] for c in cross], "-o", color=BLUE[6], ms=3,
             label=r"two-stage, optimised, $\pm2\sigma_M$ freedom ($w=3$)")
    axr.axvline(N, color=INK, lw=0.6, ls=":")
    axr.text(N * 1.3, 5.6, "this program", color=INK)
    axr.set_xscale("log")
    axr.set_ylabel(r"reach: $Z_{\mathrm{local}}^{\mathrm{full}}$ at 50% power")
    axr.legend(loc="upper left")
    axg.plot(Nx, [c[4] for c in cross], "-o", color="#009e73", ms=3)
    axg.axvline(N, color=INK, lw=0.6, ls=":")
    axg.set_xscale("log")
    axg.set_ylim(0, 22)
    axg.set_xlabel(r"number of independent looks $N$")
    axg.set_ylabel(r"$R^{*}$ for the" "\n" "split to win")
    figc.tight_layout()
    outc = io.save(figc, paths.plot("ab_split_crossover.png"), dpi=400)
    print(f"wrote {outc}")

    fs = np.linspace(0.05, 0.95, 400)
    fig, ax = plt.subplots(figsize=(5.6, 3.65))
    style(ax)
    ramp = ["#c6dbef", "#9ecae1", "#4292c6", "#2171b5", "#08306b"]
    for zcut, col in zip((2.0, 2.5, 3.0, 3.5, 4.0), ramp):
        ax.plot(fs, [reach(f, N, zcut) for f in fs], color=col,
                label=rf"$Z_{{\mathrm{{cut}}}} = {zcut:g}$")
    ax.axhline(Z_single, color=C_ARG, ls="--")
    ax.text(0.965, Z_single - 0.10, rf"single-stage scan, $Z = {Z_single:.2f}$", color=C_ARG,
            ha="right", va="top")
    ax.plot([opt[1]], [opt[2]], "o", ms=4.5, color="#08306b", mec="white", mew=0.8, zorder=5)
    ax.annotate(rf"best split: $f = {opt[1]:.2f}$, reach ${opt[2]:.2f}$",
                (opt[1], opt[2]), textcoords="offset points", xytext=(14, -22), ha="left",
                color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=0.5))
    ax.set_xlabel(r"exploration fraction $f$")
    ax.set_ylabel(r"reach: full-dataset $Z_{\mathrm{local}}$ at 50% power")
    ax.set_ylim(5.6, 10)
    ax.legend(loc="upper center", ncol=2, title=r"stage-A threshold")
    fig.tight_layout()
    out = io.save(fig, paths.plot("ab_split_reach.png"), dpi=400)
    print(f"\nwrote {out}")

    scan_N = {r["scan"]: int(r["N_trials"])
              for r in io.read_rows(paths.table("scaled_scan.csv"))}
    N_SCAN = {"unlensed scan": scan_N["ten objects, any trigger (Run 2+3, ~400 fb-1)"],
              "lensed scan": scan_N["ten objects, with selection lenses (Run 2+3, ~400 fb-1)"]}

    rng_s = np.random.default_rng(20260818)

    ZCUTS = (2.0, 2.5, 3.0, 3.5, 4.0)
    rows = []

    print("\n" + "=" * 78)
    print("priced on the combinatorial scan instead of the model space")
    for tag, Nv in list(N_SCAN.items()) + [("model space", N_sel)]:
        Z1 = z5(Nv)
        r_opt, zc_opt, f_opt = opt_box(reach, Nv)
        r_naive = reach(0.5, Nv, 3.0)
        kb = Nv * p1(zc_opt)

        print(f"\n{tag}: N = {Nv:,.0f}   single-stage exactly corrected = {Z1:.2f}")
        print(f"  optimised split: Z_cut = {zc_opt:.2f}, f = {f_opt:.2f} -> reach {r_opt:.2f} "
              f"({r_opt-Z1:+.2f})   pre-registers {kb:.0f} windows, claim bar "
              f"{zB_req(Nv, zc_opt):.2f}")
        for zc in ZCUTS:
            fs2 = np.linspace(0.02, 0.98, 481)
            rs2 = [reach(f, Nv, zc) for f in fs2]
            i2 = int(np.argmin(rs2))
            print(f"    Z_cut = {zc:.1f}: k_bkg = {Nv*p1(zc):8.1f}  claim bar = "
                  f"{zB_req(Nv, zc):.2f}  best f = {fs2[i2]:.2f}  reach = {rs2[i2]:.2f} "
                  f"({rs2[i2]-Z1:+.2f})")
        print(f"  naive 50/50 at Z_cut = 3: reach {r_naive:.2f} ({r_naive-Z1:+.2f})")
        print(f"  break-even R* = {math.exp(0.5*(r_opt*r_opt - Z1*Z1)):.0f}")

        wide = {}
        for w in (1.0, 3.0):
            rg = min((min(reach(f, Nv, zc, widen=w) for f in np.linspace(0.02, 0.98, 481)), zc)
                     for zc in ZCUTS)
            rf, zf, ff = reach2_opt(Nv, widen=w)
            wide[w] = (rg[0], rf, math.exp(0.5 * (rf * rf - Z1 * Z1)))
            print(f"  [w = {w:.0f}] grid optimum {rg[0]:.2f} ({rg[0]-Z1:+.2f}) at Z_cut = {rg[1]:.1f}; "
                  f"fine optimum {rf:.2f} ({rf-Z1:+.2f}) at Z_cut = {zf:.2f}, f = {ff:.2f}; "
                  f"R* = {math.exp(0.5*(rf*rf - Z1*Z1)):.0f} (grid {math.exp(0.5*(rg[0]**2 - Z1*Z1)):.0f})")

        toys = {}
        for zc in (3.0, zc_opt):
            k = rng_s.binomial(Nv, p1(zc), size=20000)
            u = rng_s.random(20000)
            best = np.where(k > 0, ndtri(np.clip(u, 1e-300, 1) ** (1.0 / np.maximum(k, 1))), np.nan)
            toys[zc] = (k.mean(), np.nanmax(best), int(np.nansum(best >= zB_req(Nv, zc))))
            print(f"  2e4 background-only toys at Z_cut = {zc:.2f}: k = {k.mean():.0f} +- "
                  f"{k.std():.0f}, best confirmation max = {np.nanmax(best):.2f}, "
                  f"99% = {np.nanpercentile(best, 99):.2f}, claim bar = {zB_req(Nv, zc):.2f}, "
                  f"false claims = {int(np.nansum(best >= zB_req(Nv, zc)))}")

        rows.append({"basis": tag, "N_trials": round(Nv), "Z_single_stage": round(Z1, 2),
                     "Z_cut_opt": round(zc_opt, 2), "f_opt": round(f_opt, 2),
                     "reach_opt": round(r_opt, 2), "cost_opt": round(r_opt - Z1, 2),
                     "reach_5050_zcut3": round(r_naive, 2), "cost_5050": round(r_naive - Z1, 2),
                     "k_zcut3": round(Nv * p1(3.0)), "claim_bar_zcut3": round(zB_req(Nv, 3.0), 2),
                     "toy_best_zcut3": round(toys[3.0][1], 2), "toy_false_claims_zcut3": toys[3.0][2],
                     "reach_w1": round(wide[1.0][0], 2), "reach_w3": round(wide[3.0][0], 2),
                     "R_star_w1": round(wide[1.0][2]), "R_star_w3": round(wide[3.0][2])})

    out_csv = io.write_dicts(paths.table("ab_split_scan.csv"), rows)
    print(f"\nwrote {out_csv}")
