#!/usr/bin/env python3
"""Does an IMPERFECT significance estimator change the A/B verdict?

The two-stage study (scripts/ab_split_{budget,toys}.py, note section 9) assumes a perfectly
calibrated estimator: every look returns Z ~ N(0,1) under background. Then splitting can only
lose, because the look-elsewhere effect is conserved. This script relaxes that assumption and
adds three physically distinct classes of estimator defect, each with its own rate, and asks
how much ROBUSTNESS the split buys in exchange for its 0.5 sigma of reach.

The three defect classes differ only in how they correlate between the two halves, which is the
only thing the split can see:

  BIAS    a mismodelled background / detector artefact in that window. The pull grows with the
          integrated luminosity exactly like a signal: Z = delta*B/sqrt(B) = delta*sqrt(wB).
          Both halves see the SAME delta.                        -> split sees a signal
  VAR     the look's uncertainty is underestimated by a factor s (bad fit, underestimated
          systematic). The miscalibration is a property of the look and is shared, but what it
          amplifies is the independent statistical fluctuation of each half.
                                                                 -> split suppresses partially
  GLITCH  a per-RUN failure: a fit that did not converge on this particular subsample, a network
          artefact, a corrupted input. Redrawn independently every time the estimator is run.
                                                                 -> split suppresses quadratically

Everything is exact quadrature over the defect magnitude (per-look independence makes the claim
probability 1-(1-p)^N); a toy MC at the end validates the mean-field k used for the stage-B bar.

Public inputs only (bump_observables + public_obs_map). Writes
results/plots/ab_split_outliers.png and prints the tables quoted in section 9.7 of the
note and in results/overviews/TWO_STAGE_UNBLINDING.md. Needs numpy/scipy.
"""
import os, math, sys
import numpy as np
from scipy.special import erfc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import canon, ns_scan, z_local_for_global5 as z5
from public_obs_map import PUBLIC_OBS, nsel

OUT = os.path.join(ROOT, "results", "plots")
os.makedirs(OUT, exist_ok=True)

obs = sorted({canon(o) for objs in PUBLIC_OBS.values() for o in objs})
N = int(round(sum(nsel(o) * ns_scan(o) for o in obs)))
SQ2 = math.sqrt(2.0)

Z_CUT, F_OPT, WIDEN = 3.0, 0.30, 3.0     # recommended working point (TWO_STAGE_UNBLINDING.md)
Z_SINGLE = z5(N)

def Q(z):
    return 0.5 * erfc(np.asarray(z, dtype=float) / SQ2)

_GX, _GW = np.polynomial.legendre.leggauss(600)

def halfnormal(beta, zmax=10.0):
    """Nodes/weights for |b|, b ~ N(0, beta^2): sum(w) = 1."""
    hi = zmax * beta
    b = 0.5 * hi * (_GX + 1.0)
    w = 0.5 * hi * _GW * math.sqrt(2.0 / math.pi) / beta * np.exp(-0.5 * (b / beta) ** 2)
    return b, w / w.sum()

# ---------------------------------------------------------------- defect models
# Each returns three per-look probabilities for a look whose TRUE full-dataset signal is mu:
#   p_single(z)          reported full-data Z >= z
#   p_selA(zc)           reported Z on the A-half >= zc
#   p_joint(zc, zr)      selected in A AND confirmed in B (>= zr)

class Clean:
    name = "perfectly calibrated"
    def p_single(self, z, mu):      return Q(z - mu)
    def p_selA(self, zc, mu, f):    return Q(zc - math.sqrt(f) * mu)
    def p_joint(self, zc, zr, mu, f):
        return Q(zc - math.sqrt(f) * mu) * Q(zr - math.sqrt(1 - f) * mu)

class Bias:
    """Coherent: a luminosity-scaling pull b, identical in both halves."""
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
    """Shared miscalibration by factor s; the amplified fluctuation is per-half independent."""
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
    """Incoherent: an additive artefact redrawn independently at every estimator run."""
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
    """Independent superposition (rates are small, so the mixtures simply add)."""
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

# ---------------------------------------------------------------- procedures
def zreq(model, zc=Z_CUT, f=F_OPT, widen=WIDEN):
    """Stage-B bar. Contamination floods the A-selection, which raises it by itself."""
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
    """Single-stage threshold whose spurious-claim rate equals `rate`."""
    return bisect(lambda z: -claim_single(model, z, 0.0), -rate, 3.0, 30.0)

def gain(model, split_model=None, zc=Z_CUT, f=F_OPT):
    """Matched-robustness advantage of the split, in sigma of full-dataset reach.

    Positive = the split reaches a weaker signal than a single-stage scan whose bar has been
    raised until it is EQUALLY unlikely to make a spurious claim. `split_model` differs from
    `model` only for the cross-half-training case, where the split's estimator is clean by
    construction while the single pass cannot be.  Returns NaN once the contamination is so
    heavy that neither procedure controls the false-claim rate at all (spurious prob > 20%),
    where the reach comparison is meaningless."""
    sm = split_model or model
    rate = claim_split(sm, 0.0, zc, f)
    if rate > 0.2: return float("nan")
    zm = matched_threshold(model, rate)
    return reach_single(model, zm) - reach_split(sm, zc, f)

# ================================================================ 0. setup
clean = Clean()
zr0, k0 = zreq(clean)
print(f"budget N = {N:,} looks    single-stage bar Z = {Z_SINGLE:.2f}")
print(f"working point: f = {F_OPT:g}, Z_cut = {Z_CUT:g}, widen = {WIDEN:g}"
      f"  ->  k_bkg = {k0 - 1:.1f}, stage-B bar Z_B = {zr0:.2f}")
print(f"perfect estimator:  single-stage reach {reach_single(clean):.2f}   "
      f"split reach {reach_split(clean):.2f}   (cost {reach_split(clean) - Z_SINGLE:+.2f})")
print(f"                    spurious-claim prob: single {claim_single(clean, Z_SINGLE):.2e}   "
      f"split {claim_split(clean):.2e}\n")

# ================================================================ 1. per-class robustness
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

# ================================================================ 2. matched-robustness reach
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

# ---------------------------------------------------------------- cross-half training
# The sharpest version of an incoherent defect: a background estimator TRAINED on the data
# overfits its own fluctuations. A single pass cannot avoid it -- there is no held-out half.
# The two-stage scheme trains on the complementary half at each stage, so the confirmation Z is
# clean BY CONSTRUCTION: not eps^2, but zero.
class CrossHalf:
    """Variant 3: the background model is trained on A only. Stage A scores itself, so it still
    overfits at rate eps and its candidate list is inflated; stage B is scored with that same
    A-trained model, which never saw B, so the CONFIRMATION Z is overfit-free by construction --
    not eps^2-suppressed, zero. Unlike full cross-training (train B / score A), nothing about B
    enters the pre-registration, so the k-coin bookkeeping stays exact."""
    def __init__(self, eps, beta=2.0):
        self.g = Glitch(eps, beta)
        self.name = f"XHALF  eps={eps:.1e} beta={beta:g}"
    def p_single(self, z, mu):      return self.g.p_single(z, mu)   # no held-out half in one pass
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

# ================================================================ 3. coherence is the knob
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

# ================================================================ 4. a realistic cocktail
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

# ================================================================ 5. toy validation
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

# ================================================================ 6. figure
fig, ax = plt.subplots(figsize=(9.2, 5.8))
for e, col in ((1e-4, "#4c78a8"), (1e-3, "#d9a441"), (1e-2, "#d1495b")):
    n3 = N * e * float(Glitch(1.0).w @ Q(3.0 - Glitch(1.0).t))
    ax.plot(rhos, curves[e], lw=2.2, color=col,
            label=f"defect rate $\\epsilon$ = {e:.0e}  ({n3:.1f} artefacts $>3\\sigma$ per scan)")
ax.axhline(0.0, color="#444444", lw=1.2)
ax.axhline(-(reach_split(clean) - Z_SINGLE), color="#888888", ls=":", lw=1.4)
ax.text(0.012, -(reach_split(clean) - Z_SINGLE) + 0.03,
        f"perfect estimator: split costs {reach_split(clean) - Z_SINGLE:.2f}$\\sigma$",
        fontsize=9, color="#555555")
ax.fill_between([0, 1], 0, 3, color="#4c78a8", alpha=0.07)
ax.text(0.42, 1.75, "SPLIT WINS\nreaches a weaker signal at equal\nspurious-claim probability",
        fontsize=10.5, color="#2f4b6e")
ax.text(0.60, -0.45, "single stage wins\n(just raise the bar)", fontsize=10.5, color="#8a5a2b")
ax.annotate("break-even $\\rho \\approx 0.80$, essentially independent of the rate",
            xy=(0.80, 0.0), xytext=(0.30, 0.52), fontsize=9.5, color="#444444",
            arrowprops=dict(arrowstyle="->", color="#444444", lw=1.1))
ax.set_xlim(0, 1); ax.set_ylim(-0.75, 2.6)
ax.set_xlabel(r"coherent fraction $\rho$: outliers that repeat in BOTH halves"
              "\n" r"($\rho=0$ per-run glitches   $\to$   $\rho=1$ mismodelled background)")
ax.set_ylabel(r"matched-robustness advantage of the split  [$\sigma$ of reach]")
ax.set_title("An imperfect significance estimator flips the A/B verdict --\n"
             "but only for outliers that do not repeat in the confirmation half",
             fontsize=12)
ax.grid(ls=":", alpha=0.35); ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "ab_split_outliers.png"), dpi=130)
print("\nwrote ab_split_outliers.png")

# ================================================================ 7. mechanism + scaling figure
# Why the split separates one kind of outlier and not the other, in one picture: in the
# (Z_A, Z_B) plane a COHERENT pull sits on the locus Z_B = sqrt((1-f)/f) Z_A -- exactly where a
# real signal sits -- while an incoherent one sits on the Z_B ~ 0 axis. The split is a cut that
# separates the axis from the locus, and by construction cannot separate the locus from itself.
SF, SB = math.sqrt(F_OPT), math.sqrt(1 - F_OPT)
BETA = 2.0
rng2 = np.random.default_rng(4711)

def selected(kind, n_target, beta=BETA, zsig=7.5):
    """Draw (Z_A, Z_B) for looks that PASS the stage-A cut, for each population."""
    out_a, out_b = [], []
    while len(out_a) < n_target:
        m = 200000
        if kind == "glitch":                      # artefact in the A run only; B redrawn clean
            za = rng2.standard_normal(m) + np.abs(rng2.standard_normal(m)) * beta
            zb = rng2.standard_normal(m)
        elif kind == "bias":                      # one pull b, seen by both halves as sqrt(w) b
            b = np.abs(rng2.standard_normal(m)) * beta
            za = SF * b + rng2.standard_normal(m)
            zb = SB * b + rng2.standard_normal(m)
        elif kind == "signal":
            za = SF * zsig + rng2.standard_normal(m)
            zb = SB * zsig + rng2.standard_normal(m)
        else:                                     # clean background
            za, zb = rng2.standard_normal(m), rng2.standard_normal(m)
        k = za >= Z_CUT
        out_a.append(za[k]); out_b.append(zb[k])
        if sum(len(x) for x in out_a) >= n_target: break
    return np.concatenate(out_a)[:n_target], np.concatenate(out_b)[:n_target]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.2, 5.4))

zr_c, _ = zreq(clean)
axL.axvspan(Z_CUT, 12, color="#dfe8f0", alpha=0.45, zorder=0)
axL.fill_between([Z_CUT, 12], zr_c, 12, color="#4c78a8", alpha=0.20, zorder=1)
for kind, col, lbl, ms in (("bkg", "#b0b0b0", "background", 7),
                           ("glitch", "#d1495b", "GLITCH (per-run artefact)", 13),
                           ("bias", "#e08a1e", "BIAS (mismodelling)", 13),
                           ("signal", "#2f6f4e", "real signal, $Z_{full}$ = 7.5", 13)):
    a, b = selected(kind, 900 if kind == "bkg" else 350)
    axL.scatter(a, b, s=ms, alpha=0.45 if kind == "bkg" else 0.55, color=col, lw=0, label=lbl)
xs = np.linspace(0, 12, 10)
axL.plot(xs, SB / SF * xs, ls="--", lw=1.6, color="#2f6f4e")
axL.text(7.6, 8.6, "coherent locus\n$Z_B=\\sqrt{(1-f)/f}\\,Z_A$\n(where a real signal lives)",
         fontsize=9, color="#2f6f4e")
axL.axvline(Z_CUT, color="#333333", lw=1.2)
axL.axhline(zr_c, color="#333333", lw=1.2)
axL.text(3.15, -1.6, f"$Z_{{cut}}$ = {Z_CUT:g}", fontsize=9)
axL.text(8.4, zr_c + 0.2, f"claim bar $Z_B$ = {zr_c:.2f}", fontsize=9)
axL.text(4.6, 2.0, "pre-registered\nbut DIES in B", fontsize=9.5, color="#a02735")
axL.text(4.2, 7.6, "CLAIM", fontsize=12, color="#2f4b6e", weight="bold")
axL.set_xlim(2.4, 12); axL.set_ylim(-2.5, 12)
axL.set_xlabel("$Z_A$  (exploration half, 30% of the data)")
axL.set_ylabel("$Z_B$  (confirmation half, 70%)")
axL.set_title("Only looks passing $Z_A \\geq Z_{cut}$ are shown.\n"
              "A coherent defect is on the signal locus; an incoherent one is not.", fontsize=10.5)
axL.legend(loc="lower right", fontsize=8.5, framealpha=0.9)
axL.grid(ls=":", alpha=0.3)

eps_grid = np.logspace(-6, -1, 40)
# the two defects have IDENTICAL single-stage tails by construction -- one curve, and that is
# precisely the point: a single pass cannot tell them apart, the split is what separates them.
axR.plot(eps_grid, [claim_single(Glitch(e), Z_SINGLE) for e in eps_grid], lw=2.4, color="#333333",
         label=f"single stage @ {Z_SINGLE:.2f}: BIAS and GLITCH\nare indistinguishable")
axR.plot(eps_grid, [claim_split(Bias(e)) for e in eps_grid], lw=2.4, ls="--", color="#e08a1e",
         label="two-stage split, BIAS")
axR.plot(eps_grid, [claim_split(Glitch(e)) for e in eps_grid], lw=2.4, ls="--", color="#d1495b",
         label="two-stage split, GLITCH")
axR.plot(eps_grid, 6e-2 * (eps_grid / 1e-3), color="#999999", lw=1.0, ls=":")
axR.plot(eps_grid, 3e-4 * (eps_grid / 1e-3) ** 2, color="#999999", lw=1.0, ls=":")
axR.text(2.2e-6, 4e-4, "slope 1 ($\\propto\\epsilon$)", fontsize=9, color="#777777", rotation=17)
axR.text(1.1e-4, 2e-7, "slope 2 ($\\propto\\epsilon^2$):\nthe artefact must RECUR",
         fontsize=9, color="#777777")
axR.set_xscale("log"); axR.set_yscale("log")
axR.set_ylim(1e-9, 30.0)
axR.set_xlabel(r"defect rate $\epsilon$ per look")
axR.set_ylabel("P(spurious $5\\sigma$ claim) per scan")
axR.set_title("The split changes the SCALING, not just the rate --\n"
              "but only for defects that do not repeat", fontsize=10.5)
axR.legend(loc="upper left", fontsize=8.5, framealpha=0.92)
axR.grid(ls=":", alpha=0.35, which="both")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "ab_outliers_mechanism.png"), dpi=130)
print("wrote ab_outliers_mechanism.png")

# ================================================================ 8. toy spectra
# The same two failure modes on an actual falling spectrum, scored the way an analysis would.
SIGMA_REL, M0 = 0.05, 1200.0
edges = np.geomspace(200, 4000, 260)
ctr = np.sqrt(edges[:-1] * edges[1:]); wid = np.diff(edges)
bkg = 3e6 * np.exp(-ctr / 350.0) * wid / ctr
inwin = np.abs(ctr - M0) < SIGMA_REL * M0
gauss = np.exp(-0.5 * ((ctr - M0) / (SIGMA_REL * M0)) ** 2)

def scan(counts, expect):
    z = np.zeros_like(ctr)
    for j, m in enumerate(ctr):
        w = np.abs(ctr - m) < SIGMA_REL * m
        n, b = counts[w].sum(), expect[w].sum()
        z[j] = (n - b) / math.sqrt(b) if b > 0 else 0.0
    return z

# a fractional perturbation delta*gauss gives Z = delta * sqrt(w) * S/sqrt(B) inside the scan
# window, with S the bkg-weighted Gaussian integral -- not delta*sqrt(wB), since the window
# covers only +-1 sigma_M of a falling spectrum.
BW, SW = bkg[inwin].sum(), (bkg * gauss)[inwin].sum()
def delta_for(z_target, w=1.0): return z_target * math.sqrt(BW) / (math.sqrt(w) * SW)

dA = delta_for(4.0, F_OPT)          # (a) the A-half background FIT undershoots locally
dC = delta_for(7.5)                 # (b) the background is genuinely mismodelled, both halves
true = bkg * (1 + dC * gauss)

def toy(seed):
    r = np.random.default_rng(seed)
    zA_g = scan(r.poisson(F_OPT * bkg), F_OPT * bkg * (1 - dA * gauss))   # wrong model in A only
    zB_g = scan(r.poisson((1 - F_OPT) * bkg), (1 - F_OPT) * bkg)          # correct model in B
    zA_c = scan(r.poisson(F_OPT * true), F_OPT * bkg)
    zB_c = scan(r.poisson((1 - F_OPT) * true), (1 - F_OPT) * bkg)
    return zA_g, zB_g, zA_c, zB_c

# Show a TYPICAL realisation of each, not an extreme one: both A-half peaks within ~0.5 of their
# design values (4.0 and sqrt(f)*7.5 = 4.1). A single toy peak fluctuates by ~1 either way.
for seed in range(300):
    zA_g, zB_g, zA_c, zB_c = toy(seed)
    if 3.7 <= zA_g.max() <= 4.5 and 3.7 <= zA_c.max() <= 4.6: break
print(f"  spectrum toys: seed {seed}")

fig2, axs = plt.subplots(1, 2, figsize=(13.6, 5.3), sharey=True)
for ax, (za, zb), title, verdict, vcol in [
        (axs[0], (zA_g, zB_g),
         "GLITCH: the A-half background fit undershoots at 1.2 TeV\n"
         "(the events are not there -- the model is wrong, in A only)",
         "DIES in B", "#a02735"),
        (axs[1], (zA_c, zB_c),
         "BIAS: the background really is mismodelled by "
         f"{100 * dC:.1f}% at 1.2 TeV\n(both halves see the same excess)",
         "CONFIRMS -- indistinguishable from a signal", "#2f6f4e")]:
    ax.plot(ctr, za, lw=1.6, color="#4c78a8", label="stage A  $Z_A(m)$  (30% of data)")
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
    ax.plot(ctr, np.where(mask, zb, np.nan), lw=2.4, color="#d1495b",
            label="stage B  $Z_B(m)$ -- unblinded ONLY here (70%)")
    ax.axhline(Z_CUT, color="#4c78a8", ls="--", lw=1.1)
    ax.axhline(zr_c, color="#d1495b", ls="--", lw=1.1)
    ax.text(215, Z_CUT + 0.15, f"$Z_{{cut}}$ = {Z_CUT:g}", color="#2f4b6e", fontsize=8.5)
    ax.text(215, zr_c + 0.15, f"claim bar $Z_B$ = {zr_c:.2f}", color="#a02735", fontsize=8.5)
    ax.text(0.97, 0.055, verdict, transform=ax.transAxes, ha="right", fontsize=11.5,
            color=vcol, weight="bold")
    ax.set_xscale("log"); ax.set_xlabel("mass [GeV]"); ax.set_title(title, fontsize=10)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.tick_params(axis="x", which="minor", labelsize=7.5)
    ax.grid(ls=":", alpha=0.3)
axs[0].set_ylabel("local significance in a $\\pm\\sigma_M$ window")
axs[0].legend(fontsize=8.5, loc="upper left")
fig2.suptitle("The split can only reject what does not repeat: a bad fit in A versus a genuinely "
              "wrong background", fontsize=12)
fig2.tight_layout(rect=[0, 0, 1, 0.93])
fig2.savefig(os.path.join(OUT, "ab_outliers_spectrum.png"), dpi=130)
print(f"wrote ab_outliers_spectrum.png  (glitch: Z_A = {zA_g.max():.2f} -> "
      f"Z_B = {np.nanmax(np.where(zA_g >= Z_CUT, zB_g, np.nan)):.2f};  "
      f"bias: Z_A = {zA_c.max():.2f} -> Z_B = "
      f"{np.nanmax(np.where(zA_c >= Z_CUT, zB_c, np.nan)):.2f})")
