import math

import matplotlib
import numpy as np
from scipy.special import erfc

from .. import io, paths
from ..core.bump_observables import ns_scan, z_local_for_global5 as z5
from ..core.catalogue import sorted_spectra
from ..core.public_obs_map import nsel
from ..registry import stage
from ..stats import spectra as toys
from ..stats.defects import halfnormal as _halfnormal

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

QUADRATURE = 600


def halfnormal(beta, zmax=10.0):
    return _halfnormal(beta, QUADRATURE, zmax)


@stage(
    name="ab-split-outliers",
    group="ab",
    summary="what the split buys against an imperfect estimator, defect class by class",
    outputs=["plots/ab_split_outliers.png", "plots/ab_outliers_mechanism.png",
             "plots/ab_outliers_spectrum.png"],
)
def main(options=None):

    obs = sorted_spectra()
    N = int(round(sum(nsel(o) * ns_scan(o) for o in obs)))
    SQ2 = math.sqrt(2.0)

    Z_CUT, F_OPT, WIDEN = 3.0, 0.30, 3.0
    Z_SINGLE = z5(N)

    def Q(z):
        return 0.5 * erfc(np.asarray(z, dtype=float) / SQ2)


    class Clean:
        name = "perfectly calibrated"
        def p_single(self, z, mu):      return Q(z - mu)
        def p_selA(self, zc, mu, f):    return Q(zc - math.sqrt(f) * mu)
        def p_joint(self, zc, zr, mu, f):
            return Q(zc - math.sqrt(f) * mu) * Q(zr - math.sqrt(1 - f) * mu)

    class Bias:
        def __init__(self, eps, beta=2.0):
            self.eps, self.beta = eps, beta
            self.b, self.w = halfnormal(beta)
            self.name = f"BIAS   eps={eps:.1e} beta={beta:g}"
        def p_single(self, z, mu):
            e = self.eps
            return (1 - e) * Q(z - mu) + e * float(self.w @ Q(z - mu - self.b))
        def p_selA(self, zc, mu, f):
            e, sf = self.eps, math.sqrt(f)
            return (1 - e) * Q(zc - sf * mu) + e * float(self.w @ Q(zc - sf * mu - sf * self.b))
        def p_joint(self, zc, zr, mu, f):
            e, sf, sb = self.eps, math.sqrt(f), math.sqrt(1 - f)
            cl = Q(zc - sf * mu) * Q(zr - sb * mu)
            df = float(self.w @ (Q(zc - sf * mu - sf * self.b) * Q(zr - sb * mu - sb * self.b)))
            return (1 - e) * cl + e * df

    class Var:
        def __init__(self, eps, s=2.0):
            self.eps, self.s = eps, s
            self.name = f"VAR    eps={eps:.1e} s={s:g}"
        def p_single(self, z, mu):
            e, s = self.eps, self.s
            return (1 - e) * Q(z - mu) + e * Q(z / s - mu)
        def p_selA(self, zc, mu, f):
            e, s, sf = self.eps, self.s, math.sqrt(f)
            return (1 - e) * Q(zc - sf * mu) + e * Q(zc / s - sf * mu)
        def p_joint(self, zc, zr, mu, f):
            e, s, sf, sb = self.eps, self.s, math.sqrt(f), math.sqrt(1 - f)
            return ((1 - e) * Q(zc - sf * mu) * Q(zr - sb * mu)
                    + e * Q(zc / s - sf * mu) * Q(zr / s - sb * mu))

    class Glitch:
        def __init__(self, eps, beta=2.0):
            self.eps, self.beta = eps, beta
            self.t, self.w = halfnormal(beta)
            self.name = f"GLITCH eps={eps:.1e} beta={beta:g}"
        def _tail(self, z, shift):
            e = self.eps
            return (1 - e) * Q(z - shift) + e * float(self.w @ Q(z - shift - self.t))
        def p_single(self, z, mu):      return self._tail(z, mu)
        def p_selA(self, zc, mu, f):    return self._tail(zc, math.sqrt(f) * mu)
        def p_joint(self, zc, zr, mu, f):
            return self._tail(zc, math.sqrt(f) * mu) * self._tail(zr, math.sqrt(1 - f) * mu)

    class Cocktail:
        def __init__(self, *parts):
            self.parts = parts
            self.name = "COCKTAIL " + " + ".join(p.name.split()[0] for p in parts)
        def _mix(self, fn, clean, *a):
            return clean + sum(fn(p, *a) - clean for p in self.parts)
        def p_single(self, z, mu):
            return self._mix(lambda p, *a: p.p_single(*a), Clean().p_single(z, mu), z, mu)
        def p_selA(self, zc, mu, f):
            return self._mix(lambda p, *a: p.p_selA(*a), Clean().p_selA(zc, mu, f), zc, mu, f)
        def p_joint(self, zc, zr, mu, f):
            return self._mix(lambda p, *a: p.p_joint(*a), Clean().p_joint(zc, zr, mu, f),
                             zc, zr, mu, f)

    def zreq(model, zc=Z_CUT, f=F_OPT, widen=WIDEN):
        k = N * model.p_selA(zc, 0.0, f) + 1.0
        return math.sqrt(25.0 + 2.0 * math.log(widen * k)), k

    def claim_single(model, z, mu=0.0):
        pn = model.p_single(z, 0.0)
        ps = model.p_single(z, mu) if mu > 0 else pn
        return 1.0 - (1.0 - ps) * (1.0 - pn) ** (N - 1)

    def claim_split(model, mu=0.0, zc=Z_CUT, f=F_OPT):
        zr, _ = zreq(model, zc, f)
        pn = model.p_joint(zc, zr, 0.0, f)
        ps = model.p_joint(zc, zr, mu, f) if mu > 0 else pn
        return 1.0 - (1.0 - ps) * (1.0 - pn) ** (N - 1)

    def bisect(fn, target, lo, hi, n=80):
        for _ in range(n):
            mid = 0.5 * (lo + hi)
            if fn(mid) < target: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)

    def reach_single(model, z=None):
        z = Z_SINGLE if z is None else z
        return bisect(lambda m: claim_single(model, z, m), 0.5, 1.0, 30.0)

    def reach_split(model, zc=Z_CUT, f=F_OPT):
        return bisect(lambda m: claim_split(model, m, zc, f), 0.5, 1.0, 30.0)

    def matched_threshold(model, rate):
        return bisect(lambda z: -claim_single(model, z, 0.0), -rate, 3.0, 30.0)

    def gain(model, split_model=None, zc=Z_CUT, f=F_OPT):
        sm = split_model or model
        rate = claim_split(sm, 0.0, zc, f)
        if rate > 0.2: return float("nan")
        zm = matched_threshold(model, rate)
        return reach_single(model, zm) - reach_split(sm, zc, f)

    clean = Clean()
    zr0, k0 = zreq(clean)
    print(f"budget N = {N:,} looks    single-stage bar Z = {Z_SINGLE:.2f}")
    print(f"working point: f = {F_OPT:g}, Z_cut = {Z_CUT:g}, widen = {WIDEN:g}"
          f"  ->  k_bkg = {k0 - 1:.1f}, stage-B bar Z_B = {zr0:.2f}")
    print(f"perfect estimator:  single-stage reach {reach_single(clean):.2f}   "
          f"split reach {reach_split(clean):.2f}   (cost {reach_split(clean) - Z_SINGLE:+.2f})")
    print(f"                    spurious-claim prob: single {claim_single(clean, Z_SINGLE):.2e}   "
          f"split {claim_split(clean):.2e}\n")

    print("=" * 108)
    print("1. how much spurious-claim protection does the split buy, per defect class?")
    print("   'artefacts>3' = expected # of looks per scan where the defect alone pushes the "
          "reported Z above 3")
    print("=" * 108)
    hdr = (f"{'defect':<26} {'artefacts>3':>11} {'P(spurious) 1-stage':>19} {'split':>10} "
           f"{'suppr.':>8} {'1-stage bar to match':>21}")
    for beta_or_s, ctor, tag in ((2.0, Bias, "BIAS"), (2.0, Var, "VAR"), (2.0, Glitch, "GLITCH")):
        print("\n" + hdr)
        for eps in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
            m = ctor(eps)
            r1 = claim_single(m, Z_SINGLE)
            r2 = claim_split(m)
            n3 = N * eps * (Q(3.0 / beta_or_s) if tag == "VAR" else float(m.w @ Q(3.0 - m.b))
                            if tag == "BIAS" else float(m.w @ Q(3.0 - m.t)))
            zm = matched_threshold(m, r2)
            print(f"{m.name:<26} {n3:>11.2f} {r1:>19.2e} {r2:>10.2e} "
                  f"{r1 / max(r2, 1e-300):>8.1f}x {zm:>21.2f}")

    print("\n" + "=" * 108)
    print("2. THE COMPARISON THAT MATTERS: at EQUAL spurious-claim probability, who reaches lower?")
    print("   the single stage can always buy robustness by raising its bar -- the question is the "
           "exchange rate")
    print("=" * 108)
    print(f"{'defect':<26} {'split reach':>12} {'matched 1-stage bar':>20} "
          f"{'= its reach':>12} {'split gain':>12}")
    for ctor in (Bias, Var, Glitch):
        print()
        for eps in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
            m = ctor(eps)
            zm = matched_threshold(m, claim_split(m))
            print(f"{m.name:<26} {reach_split(m):>12.2f} {zm:>20.2f} "
                  f"{reach_single(m, zm):>12.2f} {gain(m):>+12.2f}")

    class CrossHalf:
        def __init__(self, eps, beta=2.0):
            self.g = Glitch(eps, beta)
            self.name = f"XHALF  eps={eps:.1e} beta={beta:g}"
        def p_single(self, z, mu):      return self.g.p_single(z, mu)
        def p_selA(self, zc, mu, f):    return self.g.p_selA(zc, mu, f)
        def p_joint(self, zc, zr, mu, f):
            return self.g.p_selA(zc, mu, f) * Q(zr - math.sqrt(1 - f) * mu)

    print("\n  training arrangements, when the background model is FIT TO DATA and overfits at rate eps")
    print("  (a single pass has no held-out half, so it always carries the defect):")
    print(f"  {'overfit rate':<26} {'1-stage bar':>12} {'k in A':>8} {'split reach':>12} "
          f"{'split gain':>12}   arrangement")
    for eps in (1e-5, 1e-4, 1e-3, 1e-2):
        for mdl, tag in ((Glitch(eps), "v1 self-contained (train A/score A, train B/score B)"),
                         (CrossHalf(eps), "v3 A-trained throughout (train A, score A then B)"),
                         (clean, "v2 full cross (train B/score A, train A/score B) -- blindness leak")):
            g = gain(Glitch(eps), split_model=mdl)
            print(f"  {Glitch(eps).name:<26} {matched_threshold(Glitch(eps), claim_split(mdl)):>12.2f} "
                  f"{zreq(mdl)[1] - 1:>8.1f} {reach_split(mdl):>12.2f} {g:>+12.2f}   {tag}")
        print()

    print("\n" + "=" * 108)
    print("3. THE DECIDING VARIABLE: what fraction of the outliers is COHERENT between the halves?")
    print("   total defect rate eps, of which rho is bias-type (survives the split) and 1-rho is")
    print("   glitch-type (redrawn per run). The split pays for itself below the break-even rho.")
    print("=" * 108)

    def mix(eps, rho, beta=2.0):
        return Cocktail(Bias(rho * eps, beta), Glitch((1 - rho) * eps, beta))

    rhos = np.linspace(0.0, 1.0, 41)
    curves = {}
    print(f"{'rho (coherent frac)':>20} " + " ".join(f"{f'eps={e:.0e}':>12}" for e in (1e-4, 1e-3, 1e-2)))
    for e in (1e-4, 1e-3, 1e-2):
        curves[e] = [gain(mix(e, r)) for r in rhos]
    for i in range(0, 41, 4):
        print(f"{rhos[i]:>20.2f} " + " ".join(f"{curves[e][i]:>+12.2f}" for e in (1e-4, 1e-3, 1e-2)))
    for e in (1e-4, 1e-3, 1e-2):
        lo, hi = 0.0, 1.0
        if gain(mix(e, 1.0)) > 0:
            print(f"  eps = {e:.0e}: the split wins at EVERY coherence")
            continue
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if gain(mix(e, mid)) > 0: lo = mid
            else: hi = mid
        n3 = N * e * float(Glitch(1.0).w @ Q(3.0 - Glitch(1.0).t))
        print(f"  eps = {e:.0e} ({n3:.1f} artefacts above 3 sigma per scan): "
              f"break-even at rho = {0.5 * (lo + hi):.2f}")

    print("\n" + "=" * 108)
    print("4. a realistic cocktail: rare hard bias + occasional bad fit + rare glitch")
    print("=" * 108)
    for label, ck in (("optimistic", Cocktail(Bias(1e-4), Var(1e-4, 1.5), Glitch(1e-4))),
                      ("plausible",  Cocktail(Bias(1e-3), Var(1e-3, 1.5), Glitch(1e-3))),
                      ("pessimistic", Cocktail(Bias(1e-2), Var(1e-2, 2.0), Glitch(1e-2)))):
        zm = matched_threshold(ck, claim_split(ck))
        print(f"{label:<12} P(spurious 5s): 1-stage@{Z_SINGLE:.2f} = {claim_single(ck, Z_SINGLE):.2e}"
              f"   split = {claim_split(ck):.2e}")
        print(f"{'':<12} split reach {reach_split(ck):.2f}  vs matched 1-stage bar {zm:.2f} -> reach "
              f"{reach_single(ck, zm):.2f}   split gain {gain(ck):+.2f} sigma")
        tot = claim_split(ck) - claim_split(clean)
        for p in ck.parts:
            frac = (claim_split(p) - claim_split(clean)) / tot if tot > 0 else 0.0
            print(f"{'':<16}{p.name:<26} 1-stage {claim_single(p, Z_SINGLE):.2e} -> split "
                  f"{claim_split(p):.2e}   = {100 * frac:5.1f}% of what the split still lets through")

    print("\n" + "=" * 108)
    print("5. toy MC check of the mean-field stage-B bar under contamination (GLITCH eps=1e-2)")
    print("=" * 108)
    rng = np.random.default_rng(20260810)
    m = Glitch(1e-2)
    NTOY, CHUNK = 60000, 3000
    hits_split = hits_single = 0
    kk = []
    for s in range(0, NTOY, CHUNK):
        n = min(CHUNK, NTOY - s)
        def draw(shape):
            g = rng.standard_normal(shape)
            d = rng.random(shape) < m.eps
            return g + np.where(d, np.abs(rng.standard_normal(shape)) * m.beta, 0.0)
        zA, zB, zF = draw((n, N)), draw((n, N)), draw((n, N))
        sel = zA >= Z_CUT
        k = sel.sum(axis=1); kk.append(k)
        bar = np.sqrt(25.0 + 2.0 * np.log(WIDEN * np.maximum(k, 1)))
        best = np.where(sel, zB, -np.inf).max(axis=1)
        hits_split += int(((k > 0) & (best >= bar)).sum())
        hits_single += int((zF.max(axis=1) >= Z_SINGLE).sum())
    kk = np.concatenate(kk)
    zr_m, k_m = zreq(m)
    print(f"  mean k in stage A: toys {kk.mean():.1f}  vs mean-field {k_m - 1:.1f}   "
          f"(clean would be {k0 - 1:.1f})")
    print(f"  P(spurious claim), split : toys {hits_split / NTOY:.2e} "
          f"[+-{math.sqrt(max(hits_split,1))/NTOY:.1e}]   analytic {claim_split(m):.2e}")
    print(f"  P(spurious claim), single: toys {hits_single / NTOY:.2e} "
          f"[+-{math.sqrt(max(hits_single,1))/NTOY:.1e}]   analytic {claim_single(m, Z_SINGLE):.2e}")

    fig, ax = plt.subplots(figsize=(5.7, 3.6))
    for e, col in ((1e-4, "#0072b2"), (1e-3, "#e69f00"), (1e-2, "#d55e00")):
        n3 = N * e * float(Glitch(1.0).w @ Q(3.0 - Glitch(1.0).t))
        ax.plot(rhos, curves[e], lw=1.3, color=col,
                label=f"defect rate $\\epsilon$ = {e:.0e}  ({n3:.1f} artefacts $>3\\sigma$ per scan)")
    ax.axhline(0.0, color="#444444", lw=1.2)
    ax.axhline(-(reach_split(clean) - Z_SINGLE), color="#000000", ls=":", lw=1.4)
    ax.text(0.012, -(reach_split(clean) - Z_SINGLE) + 0.03, "perfect estimator",
            fontsize=9.5, color="#000000")
    ax.fill_between([0, 1], 0, 3, color="#0072b2", alpha=0.07)
    ax.text(0.36, 1.9, "split ahead", color="#08306b")
    ax.text(0.62, -0.45, "single stage ahead", color="#8a5a2b")
    ax.set_xlim(0, 1); ax.set_ylim(-0.75, 2.6)
    ax.set_xlabel(r"coherent fraction $\rho$: outliers that repeat in both halves")
    ax.set_ylabel(r"advantage of the split  [$\sigma$ of reach]")
    ax.grid(ls=":", alpha=0.35); ax.legend(loc="upper right", fontsize=9.5, frameon=False)
    fig.tight_layout()
    io.save(fig, paths.plot("ab_split_outliers.png"), dpi=400)
    print("\nwrote ab_split_outliers.png")

    SF, SB = math.sqrt(F_OPT), math.sqrt(1 - F_OPT)
    BETA = 2.0
    rng2 = np.random.default_rng(4711)

    def selected(kind, n_target, beta=BETA, zsig=7.5):
        out_a, out_b = [], []
        while len(out_a) < n_target:
            m = 200000
            if kind == "glitch":
                za = rng2.standard_normal(m) + np.abs(rng2.standard_normal(m)) * beta
                zb = rng2.standard_normal(m)
            elif kind == "bias":
                b = np.abs(rng2.standard_normal(m)) * beta
                za = SF * b + rng2.standard_normal(m)
                zb = SB * b + rng2.standard_normal(m)
            elif kind == "signal":
                za = SF * zsig + rng2.standard_normal(m)
                zb = SB * zsig + rng2.standard_normal(m)
            else:
                za, zb = rng2.standard_normal(m), rng2.standard_normal(m)
            k = za >= Z_CUT
            out_a.append(za[k]); out_b.append(zb[k])
            if sum(len(x) for x in out_a) >= n_target: break
        return np.concatenate(out_a)[:n_target], np.concatenate(out_b)[:n_target]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.2, 2.9))

    zr_c, _ = zreq(clean)
    axL.axvspan(Z_CUT, 12, color="#dfe8f0", alpha=0.45, zorder=0)
    axL.fill_between([Z_CUT, 12], zr_c, 12, color="#0072b2", alpha=0.20, zorder=1)
    for kind, col, lbl, ms in (("bkg", "#b0b0b0", "background", 7),
                               ("glitch", "#d55e00", "GLITCH (per-run artefact)", 13),
                               ("bias", "#e08a1e", "BIAS (mismodelling)", 13),
                               ("signal", "#2f6f4e", r"real signal, $Z_{\mathrm{full}}$ = 7.5", 13)):
        a, b = selected(kind, 900 if kind == "bkg" else 350)
        axL.scatter(a, b, s=ms, alpha=0.45 if kind == "bkg" else 0.55, color=col, lw=0, label=lbl)
    xs = np.linspace(0, 12, 10)
    axL.plot(xs, SB / SF * xs, ls="--", lw=1.2, color="#2f6f4e")
    axL.text(7.9, 9.2, "signal locus", fontsize=9.5, color="#2f6f4e")
    axL.axvline(Z_CUT, color="#333333", lw=1.2)
    axL.axhline(zr_c, color="#333333", lw=1.2)
    axL.text(3.15, -1.6, f"$Z_{{cut}}$ = {Z_CUT:g}", fontsize=9)
    axL.text(8.4, zr_c + 0.2, f"claim bar $Z_B$ = {zr_c:.2f}", fontsize=9)
    axL.text(4.2, 7.6, "claim region", color="#08306b")
    axL.set_xlim(2.4, 12); axL.set_ylim(-2.5, 12)
    axL.set_xlabel("$Z_A$  (exploration half, 30% of the data)")
    axL.set_ylabel("$Z_B$  (confirmation half, 70%)")
    axL.set_title("Confirmation against exploration significance", loc="left")
    axL.legend(loc="lower right", fontsize=9, frameon=False)
    axL.grid(ls=":", alpha=0.3)

    eps_grid = np.logspace(-6, -1, 40)
    axR.plot(eps_grid, [claim_single(Glitch(e), Z_SINGLE) for e in eps_grid], lw=1.3, color="#333333",
             label=f"single stage @ {Z_SINGLE:.2f}: either defect")
    axR.plot(eps_grid, [claim_split(Bias(e)) for e in eps_grid], lw=1.3, ls="--", color="#e08a1e",
             label="two-stage split, BIAS")
    axR.plot(eps_grid, [claim_split(Glitch(e)) for e in eps_grid], lw=1.3, ls="--", color="#d55e00",
             label="two-stage split, GLITCH")
    axR.plot(eps_grid, 6e-2 * (eps_grid / 1e-3), color="#999999", lw=1.0, ls=":")
    axR.plot(eps_grid, 3e-4 * (eps_grid / 1e-3) ** 2, color="#999999", lw=1.0, ls=":")
    axR.text(2.2e-6, 4e-4, r"$\propto\epsilon$", color="#777777", rotation=17)
    axR.text(1.1e-4, 2e-7, r"$\propto\epsilon^2$", color="#777777")
    axR.set_xscale("log"); axR.set_yscale("log")
    axR.set_ylim(1e-9, 30.0)
    axR.set_xlabel(r"defect rate $\epsilon$ per look")
    axR.set_ylabel("P(spurious $5\\sigma$ claim) per scan")
    axR.set_title("Spurious-claim probability", loc="left")
    axR.legend(loc="upper left", fontsize=9, frameon=False)
    axR.grid(ls=":", alpha=0.35, which="both")
    fig.tight_layout()
    io.save(fig, paths.plot("ab_outliers_mechanism.png"), dpi=400)
    print("wrote ab_outliers_mechanism.png")

    SIGMA_REL, M0 = 0.05, 1200.0
    edges, ctr, wid = toys.grid(200, 4000, 260)
    bkg = toys.background(ctr, wid)
    inwin = np.abs(ctr - M0) < SIGMA_REL * M0
    gauss = np.exp(-0.5 * ((ctr - M0) / (SIGMA_REL * M0)) ** 2)

    def scan(counts, expect):
        return toys.window_scan(ctr, counts, expect, SIGMA_REL)

    BW, SW = bkg[inwin].sum(), (bkg * gauss)[inwin].sum()
    def delta_for(z_target, w=1.0): return z_target * math.sqrt(BW) / (math.sqrt(w) * SW)

    dA = delta_for(4.0, F_OPT)
    dC = delta_for(7.5)
    true = bkg * (1 + dC * gauss)

    def toy(seed):
        r = np.random.default_rng(seed)
        zA_g = scan(r.poisson(F_OPT * bkg), F_OPT * bkg * (1 - dA * gauss))
        zB_g = scan(r.poisson((1 - F_OPT) * bkg), (1 - F_OPT) * bkg)
        zA_c = scan(r.poisson(F_OPT * true), F_OPT * bkg)
        zB_c = scan(r.poisson((1 - F_OPT) * true), (1 - F_OPT) * bkg)
        return zA_g, zB_g, zA_c, zB_c

    for seed in range(300):
        zA_g, zB_g, zA_c, zB_c = toy(seed)
        if 3.7 <= zA_g.max() <= 4.5 and 3.7 <= zA_c.max() <= 4.6: break
    print(f"  spectrum toys: seed {seed}")

    fig2, axs = plt.subplots(1, 2, figsize=(6.2, 2.8), sharey=True)
    for ax, (za, zb), title, verdict, vcol in [
            (axs[0], (zA_g, zB_g),
             "GLITCH: the A-half background fit undershoots at 1.2 TeV\n"
             "(the events are not there -- the model is wrong, in A only)",
             "DIES in B", "#d55e00"),
            (axs[1], (zA_c, zB_c),
             "BIAS: the background really is mismodelled by "
             f"{100 * dC:.1f}% at 1.2 TeV\n(both halves see the same excess)",
             "CONFIRMS -- indistinguishable from a signal", "#2f6f4e")]:
        ax.plot(ctr, za, lw=1.2, color="#0072b2", label="stage A  $Z_A(m)$  (30% of data)")
        sel = za >= Z_CUT
        first = True
        for j in np.where(sel)[0]:
            ax.axvspan(ctr[j] * (1 - 2 * SIGMA_REL), ctr[j] * (1 + 2 * SIGMA_REL),
                       color="#f0e6b8", alpha=0.55, zorder=0,
                       label="pre-registered window" if first else None)
            first = False
        mask = np.zeros_like(sel)
        for j in np.where(sel)[0]:
            mask |= np.abs(np.log(ctr / ctr[j])) < 2 * SIGMA_REL
        ax.plot(ctr, np.where(mask, zb, np.nan), lw=1.3, color="#d55e00",
                label="stage B  $Z_B(m)$ -- unblinded ONLY here (70%)")
        ax.axhline(Z_CUT, color="#0072b2", ls="--", lw=1.1)
        ax.axhline(zr_c, color="#d55e00", ls="--", lw=1.1)
        ax.text(215, Z_CUT + 0.15, f"$Z_{{cut}}$ = {Z_CUT:g}", color="#08306b", fontsize=8.5)
        ax.text(215, zr_c + 0.15, f"claim bar $Z_B$ = {zr_c:.2f}", color="#d55e00", fontsize=8.5)
        ax.text(0.97, 0.055, verdict, transform=ax.transAxes, ha="right",
                color=vcol, weight="bold")
        ax.set_xscale("log"); ax.set_xlabel("mass [GeV]"); ax.set_title(title)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.tick_params(axis="x", which="minor", labelsize=7.5)
        ax.grid(ls=":", alpha=0.3)
    axs[0].set_ylabel("local significance in a $\\pm\\sigma_M$ window")
    axs[0].legend(fontsize=8.5, loc="upper left")
    fig2.tight_layout()
    io.save(fig2, paths.plot("ab_outliers_spectrum.png"), dpi=400)
    print(f"wrote ab_outliers_spectrum.png  (glitch: Z_A = {zA_g.max():.2f} -> "
          f"Z_B = {np.nanmax(np.where(zA_g >= Z_CUT, zB_g, np.nan)):.2f};  "
          f"bias: Z_A = {zA_c.max():.2f} -> Z_B = "
          f"{np.nanmax(np.where(zA_c >= Z_CUT, zB_c, np.nan)):.2f})")
