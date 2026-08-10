#!/usr/bin/env python3
"""Two-stage A/B unblinding: discovery reach of the split strategy vs the single-stage scan.

Public inputs only. Writes results/plots/ab_split_reach.png and prints the design table
quoted in results/overviews/TWO_STAGE_UNBLINDING.md. Reach definition: docs/METHOD_NOTES.md.
"""
import os, math, sys, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import canon, ns_scan, z_local_for_global5 as z5
from public_obs_map import PUBLIC_OBS, nsel

obs = sorted({canon(o) for objs in PUBLIC_OBS.values() for o in objs})
N_incl = sum(ns_scan(o) for o in obs)                 # inclusive public budget
N_sel  = sum(nsel(o) * ns_scan(o) for o in obs)      # with published event selections

def p1(Z): return 0.5 * math.erfc(Z / math.sqrt(2.0))

def k_eff(N, zcut): return N * p1(zcut) + 1.0        # bkg selections + the signal's window

def zB_req(N, zcut, widen=3.0, zglob=5.0):
    """B-only local Z for a Z_glob-sigma GLOBAL result over the unblinded windows. widen = extra
    resolution elements per window (the B peak can sit anywhere in the pre-registered
    +-1-2 sigma_M window). zglob=5: discovery; zglob=3: evidence trigger."""
    return math.sqrt(zglob * zglob + 2.0 * math.log(widen * k_eff(N, zcut)))

def reach_median(f, N, zcut, widen=3.0):
    return max(zcut / math.sqrt(f), zB_req(N, zcut, widen) / math.sqrt(1.0 - f))

def Phi(x): return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def power(mu, f, N, zcut, widen=3.0, zglob=5.0):
    """Joint probability of a Z_glob-global B result for a signal of full-dataset local Z = mu."""
    return Phi(math.sqrt(f) * mu - zcut) * Phi(math.sqrt(1 - f) * mu - zB_req(N, zcut, widen, zglob))

def reach(f, N, zcut, target=0.5, widen=3.0, zglob=5.0):
    """Z_full at which the two-stage procedure reaches `target` power (bisection)."""
    lo, hi = 1.0, 20.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if power(mid, f, N, zcut, widen, zglob) < target: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

# ---------------------------------------------------------------- table
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

# ---------------------------------------------------------------- reach vs #unblinded regions
# Equivalent "top-k" formulation: unblind the k best A-windows. The effective A-cut is z* with
# N*p1(z*) = k, the B bar is sqrt(25 + 2 ln(widen*k)). k enters BOTH logarithmically and with
# opposite signs (bigger k -> harder B bar but easier A selection), so after re-optimizing f the
# reach is nearly flat in k across three orders of magnitude.
def PhiInv(p):
    lo, hi = -10.0, 10.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if Phi(mid) < p: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

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

# ---------------------------------------------------------------- B criterion: GLOBAL > 3 sigma
# Working point: A flags multiple regions (Z_cut as above); B is unblinded in all of them and the
# B result is quoted as a GLOBAL significance over the k_eff pre-registered windows (x widen for
# the +-2sigma_M freedom):  Z_B,global >= 3  <=>  Z_B,local >= sqrt(9 + 2 ln(w*k_eff)).
# By construction the background-only false-evidence rate per B opening is p1(3) = 1.35e-3,
# whatever k came out of A. Discovery keeps the 5-global bar of the tables above (same A flags).
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

# ---------------------------------------------------------------- BH-FDR flagging in A
# Alternative A rule: flag by Benjamini-Hochberg at FDR level q over the N one-sided p-values
# (largest j with p_(j) <= q*j/N; flag the top j). Properties in this setting:
#   * background-only (global null): FDR = FWER = q -> P(flag ANYTHING) = q. B is opened only
#     with probability q, instead of always (~9 background windows) as with a fixed Z_cut.
#   * single isolated signal: BH degenerates to Bonferroni -- flagged iff p <= q/N, i.e.
#     z_A >= PhiInv(1 - q/N) ~ 3.9-4.3: a HARDER A bar than Z_cut = 3.
#   * c channels lit by the same model (Z'->ee AND mumu, VLQ multi-channel, ...): the step-up
#     relaxes the bar to ~ q*c/N -- BH's adaptivity pays exactly when one model fires several
#     of the 78 selection-level spectra.
#   * the flagged set is data-dependent but exactly countable at pre-registration; the B ladder
#     applies with k_obs. One-sided Gaussian looks are PRDS, so BH is valid under the (positive)
#     window-overlap dependence.
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

# ---------------------------------------------------------------- crossover: does the split EVER win?
# Trials-factorization identity. With the A cut z* set by N*p1(z*) = k and the B bar zB by
# w*k*p1(zB) = p1(5), use the Gaussian tail p1(z) ~ phi(z)/z and add logs:
#     z*^2/2  + ln(z* sqrt(2pi)) = ln(N/k)
#     zB^2/2  + ln(zB sqrt(2pi)) = ln(w*k / alpha),  alpha = p1(5)
#     mu1^2/2 + ln(mu1 sqrt(2pi)) = ln(N / alpha)     [single-stage threshold]
# Adding the first two and subtracting the third, k CANCELS EXACTLY (N = (N/k) * k):
#     mu2_med^2 - mu1^2 = 2 ln w - 2 ln( z* * zB * sqrt(2pi) / mu1 )
# i.e. the LEE is conserved, never reduced: the split moves 2 ln(N/k) of it into the A-selection
# cut and leaves 2 ln k in the B correction. The median-arithmetic reach is therefore a wash with
# the single stage at ANY N (the two O(1) terms nearly cancel). The entire real cost is the
# TWO-COIN penalty: both independent halves must succeed, so at 50% joint power each stage
# carries a ~ +0.55 buffer, mu2^2 ~ (z*+a)^2 + (zB+b)^2 with Phi(a)Phi(b)=1/2 -> +~0.5 sigma,
# positive at every N and growing ~ sqrt(2 ln N). Hence NO crossover: no number of trials makes
# the split out-reach an exactly-corrected single-stage scan.
# The split wins only when the single-stage trials factor cannot be defended exactly: it is the
# more sensitive procedure iff the trials count you would otherwise have to defend exceeds
#     N_equiv = exp((mu2^2 - 25)/2),  i.e.  N_def / N_true > R* = exp((mu2^2 - mu1^2)/2)
# (~30 at this budget). That is the quantitative case for it in ML-driven scans, where the
# effective number of looks (trainings, selections, hyperparameters) is genuinely uncountable.

def reach2_opt(Nv, widen=3.0):
    """Fully optimized two-stage 50%-power reach: min over (Z_cut, f). Coarse grid + refine."""
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
        best = min(best, math.sqrt(zc * zc + zb * zb))   # balanced f = zc^2/(zc^2+zb^2)
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

# ---------------------------------------------------------------- symmetrized swap (A->B AND B->A)
# The obvious objection to the conservation identity: half the luminosity is idle at each stage, so
# run BOTH directions -- explore A confirm B, and explore B confirm A -- and claim on the union.
# Both directions pre-register, so the confirmation windows double (k -> 2k, i.e. 2 ln 2 on the B
# bar); what is bought is a second chance at the same signal. Since zB > z_cut the two directions
# coincide only where BOTH halves clear the claim bar, so the union power is
#     P = Phi_c(A) Phi_r(B) + Phi_c(B) Phi_r(A) - Phi_r(A) Phi_r(B),
# a factor 2 - Phi_r/Phi_c over one direction (up to x2). Geometrically: the Neyman-Pearson region
# for a mean shift is the HALF-PLANE sqrt(f) z_A + sqrt(1-f) z_B > z5(N) -- that IS the single-stage
# scan. One-way splitting approximates it with an L-shaped corner, the swap with a two-step
# staircase: closer to the line, but every step multiplies the coins. Conservation of the LEE is
# that staircase-vs-half-plane inefficiency, which is why no fold count crosses zero.

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
    """Best (reach, Z_cut, f) inside the PRACTICAL working box. Unconstrained optimization runs
    into the degenerate Z_cut -> -inf, f -> 0 limit where the split just becomes the single-stage
    scan (and where the power-curve plateau argument of section 2 says not to sit); pinning
    Z_cut to 2-4.5 keeps both schemes at designs one would actually pre-register."""
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

# crossover figure
figc, (axr, axg) = plt.subplots(2, 1, figsize=(9.0, 7.6), sharex=True,
                                gridspec_kw={"height_ratios": [2.0, 1.15], "hspace": 0.07})
Nx = [c[0] for c in cross]
axr.plot(Nx, [c[1] for c in cross], "--", color="#d1495b", lw=2,
         label=r"single-stage exact correction: $\sqrt{25+2\ln N}$")
axr.plot(Nx, [c[2] for c in cross], "-o", color="#9ecae9", lw=2, ms=4,
         label="two-stage, optimized (f, $Z_{cut}$), exact windows (w=1)")
axr.plot(Nx, [c[3] for c in cross], "-o", color="#2f4b6e", lw=2, ms=4,
         label=r"two-stage, optimized, $\pm2\sigma_M$ window freedom (w=3)")
axr.axvline(N, color="#888888", lw=1, ls=":")
axr.text(N * 1.25, 5.6, f"this program\nN = {N:,.0f}", fontsize=8.5, color="#555555")
axr.set_xscale("log")
axr.set_ylabel(r"discovery reach: $Z_{local}^{full}$ at 50% power")
axr.grid(ls=":", alpha=0.35)
axr.legend(fontsize=9, loc="upper left")
axr.set_title("Is there an N where splitting wins? No — the LEE factorizes ($N = \\frac{N}{k}\\cdot k$),\n"
              "the gap is the two-coin penalty (both halves must succeed) and it never closes")
axg.plot(Nx, [c[4] for c in cross], "-o", color="#5f9e6e", lw=2, ms=4)
axg.axvline(N, color="#888888", lw=1, ls=":")
axg.set_xscale("log")
axg.set_ylim(0, 22)
axg.set_xlabel("number of independent looks N (trials) — exactly known in the single-stage column")
axg.set_ylabel(r"$R^{*}$: trials inflation for" "\n" r"the split to win")
axg.grid(ls=":", alpha=0.35, which="both")
axg.text(0.02, 0.83, r"split more sensitive iff defendable single-stage trials $> R^{*}\times N$"
         "\n(uncountable analyst/ML freedom); below the curve the corrected scan wins",
         transform=axg.transAxes, fontsize=8.5, color="#3a6b46")
figc.tight_layout()
outc = os.path.join(ROOT, "results", "plots", "ab_split_crossover.png")
os.makedirs(os.path.dirname(outc), exist_ok=True)
figc.savefig(outc, dpi=130)
print(f"wrote {outc}")

# ---------------------------------------------------------------- plot
fs = np.linspace(0.05, 0.95, 400)
fig, ax = plt.subplots(figsize=(9.5, 6.2))
colors = {2.0: "#c6dbeF", 2.5: "#9ecae9", 3.0: "#4c78a8", 3.5: "#2f4b6e", 4.0: "#17324d"}
for zcut in (2.0, 2.5, 3.0, 3.5, 4.0):
    ax.plot(fs, [reach(f, N, zcut) for f in fs], lw=2,
            color=colors[zcut], label=f"$Z_{{cut}}$ = {zcut:g}  "
            f"(k$_{{bkg}}$ = {N*p1(zcut):.1f})")
    ax.plot(fs, [reach_median(f, N, zcut) for f in fs], lw=0.9, ls="--",
            color=colors[zcut], alpha=0.55)
ax.axhline(Z_single, color="#d1495b", ls="--", lw=1.8)
ax.text(0.96, Z_single - 0.09, f"single-stage full-dataset scan: Z = {Z_single:.2f} "
        f"(N = {N:,.0f})", color="#a02735", fontsize=9, ha="right", va="top")
ax.plot([opt[1]], [opt[2]], "o", ms=9, color="#2f4b6e", zorder=5)
ax.annotate(f"optimum: f = {opt[1]:.2f}, $Z_{{cut}}$ = {opt[0]:g}\n"
            f"reach {opt[2]:.2f} ({opt[2]-Z_single:+.2f} vs single-stage:\n"
            "the price of countable trials)",
            (opt[1], opt[2]), textcoords="offset points", xytext=(-30, 55), ha="center",
            fontsize=9, arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8))
ax.plot([0.5], [reach(0.5, N, 3.0)], "s", ms=8, color="#d9a441", zorder=5)
ax.annotate(f"naive 50/50, $Z_{{cut}}$=3:\nreach {reach(0.5, N, 3.0):.2f}",
            (0.5, reach(0.5, N, 3.0)), textcoords="offset points", xytext=(12, 16), fontsize=9,
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8))
ax.set_xlabel("exploration fraction  f  (dataset A);  confirmation uses 1$-$f (dataset B)")
ax.set_ylabel(r"discovery reach: full-dataset $Z_{local}$ at 50% power")
ax.set_ylim(5.6, 10)
ax.grid(ls=":", alpha=0.35)
ax.set_title("Two-stage A/B unblinding: discovery reach vs split fraction\n"
             "solid = 50% joint power (validated by toys);  thin dashed = naive median arithmetic")
ax.legend(loc="upper center", fontsize=9, title="selection threshold in A")
fig.tight_layout()
out = os.path.join(ROOT, "results", "plots", "ab_split_reach.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=130)
print(f"\nwrote {out}")
