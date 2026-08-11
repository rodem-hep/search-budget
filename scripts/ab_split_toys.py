#!/usr/bin/env python3
"""Toy Monte Carlo validating the two-stage A/B unblinding strategy.

Public inputs only. Writes results/plots/ab_toys_{background,power,spectrum}.png (~1 min).
Toy model and the three procedures compared: docs/METHOD_NOTES.md; report:
results/overviews/TWO_STAGE_UNBLINDING.md.
"""
import os, math, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import canon, ns_scan, z_local_for_global5 as z5
from public_obs_map import PUBLIC_OBS, nsel

OUT = os.path.join(ROOT, "results", "plots")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260703)

obs = sorted({canon(o) for objs in PUBLIC_OBS.values() for o in objs})
N = int(round(sum(nsel(o) * ns_scan(o) for o in obs)))   # event-selection-level budget
P5 = 0.5 * math.erfc(5.0 / math.sqrt(2.0))               # one-sided 5-sigma tail

def p1(z): return 0.5 * math.erfc(z / math.sqrt(2.0))
def zthr(k): return np.sqrt(25.0 + 2.0 * np.log(np.maximum(k, 1)))

Z_CUT, F_OPT = 3.0, 0.22
Z_SINGLE = z5(N)
K_BKG = N * p1(Z_CUT)

# ================================================================ 1. background-only toys
NTOY, CHUNK = 20000, 2000
k_obs, zB_sel_max, p_two, p_one = [], [], [], []
for start in range(0, NTOY, CHUNK):
    n = min(CHUNK, NTOY - start)
    gA = rng.standard_normal((n, N)); gB = rng.standard_normal((n, N))
    selm = gA >= Z_CUT
    k = selm.sum(axis=1); k_obs.append(k)
    zB = np.where(selm, gB, -np.inf).max(axis=1)         # best confirmation significance
    zB_sel_max.append(np.where(k > 0, zB, np.nan))
    # corrected global p of each procedure (validity check: must be ~Uniform(0,1))
    p_two.append(np.where(k > 0, np.minimum(1.0, k * 0.5 *
                 np.vectorize(math.erfc)(zB / math.sqrt(2.0))), 1.0))
    zF = (math.sqrt(F_OPT) * gA + math.sqrt(1 - F_OPT) * gB).max(axis=1)
    p_one.append(np.minimum(1.0, N * 0.5 * np.vectorize(math.erfc)(zF / math.sqrt(2.0))))
k_obs = np.concatenate(k_obs); zB_sel_max = np.concatenate(zB_sel_max)
p_two = np.concatenate(p_two); p_one = np.concatenate(p_one)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
# (a) selections in A
kmax = int(k_obs.max()) + 2
axs[0].hist(k_obs, bins=np.arange(-0.5, kmax + 0.5), color="#4c78a8", edgecolor="white",
            density=True, label="toys")
from math import lgamma
kk = np.arange(0, kmax)
pois = np.exp(kk * np.log(K_BKG) - K_BKG - np.array([lgamma(v + 1) for v in kk]))
axs[0].plot(kk, pois, "o-", ms=4, lw=1.2, color="#d1495b",
            label=f"Poisson($N\\,p_1(Z_{{cut}})$ = {K_BKG:.1f})")
axs[0].set_xlabel(f"# windows selected in A  ($Z_A \\geq$ {Z_CUT:g})")
axs[0].set_ylabel("fraction of toys")
axs[0].set_title(f"Stage A under background-only\nN = {N:,} looks -> "
                 f"mean {k_obs.mean():.1f} selections (predicted {K_BKG:.1f})", fontsize=10)
axs[0].legend(fontsize=8.5)
# (b) best z_B among selected windows
v = zB_sel_max[np.isfinite(zB_sel_max)]
axs[1].hist(v, bins=60, color="#4c78a8", edgecolor="white", density=True,
            label="best $Z_B$ of the unblinded windows")
med_thr = float(np.median(zthr(k_obs[k_obs > 0])))
axs[1].axvline(med_thr, color="#d1495b", ls="--", lw=1.8,
               label=f"claim threshold $\\sqrt{{25+2\\ln k}}$ (median {med_thr:.2f})")
axs[1].set_xlabel("best confirmation significance $Z_B$")
axs[1].set_title(f"Stage B kills the fluctuations\nmax over {NTOY:,} toys: "
                 f"$Z_B$ = {np.nanmax(v):.2f}  (0 false claims)", fontsize=10)
axs[1].legend(fontsize=8.5, loc="upper right")
# (c) uniformity of the corrected global p
for p, lbl, col in [(p_two, "two-stage: $k\\,p_1(Z_B^{best})$", "#4c78a8"),
                    (p_one, "single-stage: $N\\,p_1(Z^{max})$", "#d9a441")]:
    q = np.sort(p)
    axs[2].plot(np.linspace(0, 1, len(q)), q, lw=1.8, color=col, label=lbl)
axs[2].plot([0, 1], [0, 1], color="#888888", ls=":", lw=1.2, label="exact Uniform(0,1)")
axs[2].set_xlabel("expected quantile"); axs[2].set_ylabel("corrected global p-value")
axs[2].set_title("The countable trials factor is CORRECT\n(corrected p ~ uniform under bkg-only)",
                 fontsize=10)
axs[2].legend(fontsize=8.5)
fig.suptitle("Two-stage A/B unblinding, background-only toys "
             f"(f = {F_OPT:g}, $Z_{{cut}}$ = {Z_CUT:g})", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "ab_toys_background.png"), dpi=130)
print(f"wrote ab_toys_background.png   (mean k = {k_obs.mean():.2f}, "
      f"predicted {K_BKG:.2f}; max z_B = {np.nanmax(v):.2f} vs threshold ~{med_thr:.2f})")

# ================================================================ 2. power curves
def toy_power(mu, f, zcut, ntoy=1500):
    """Discovery probability: single-stage and two-stage on the same toys, signal at look 0."""
    gA = rng.standard_normal((ntoy, N)); gB = rng.standard_normal((ntoy, N))
    gA[:, 0] += math.sqrt(f) * mu; gB[:, 0] += math.sqrt(1 - f) * mu
    zF = math.sqrt(f) * gA + math.sqrt(1 - f) * gB
    disc1 = (zF.max(axis=1) >= Z_SINGLE)
    selm = gA >= zcut
    k = selm.sum(axis=1)
    zB = np.where(selm, gB, -np.inf).max(axis=1)
    disc2 = (k > 0) & (zB >= zthr(k))
    return disc1.mean(), disc2.mean()

mus = np.arange(5.0, 9.51, 0.25)
pw_single, pw_5050, pw_opt = [], [], []
for mu in mus:
    s1, t5050 = toy_power(mu, 0.5, Z_CUT)
    s1b, topt = toy_power(mu, F_OPT, Z_CUT)
    pw_single.append(0.5 * (s1 + s1b))       # single-stage is f-independent; average both runs
    pw_5050.append(t5050); pw_opt.append(topt)

# analytic 50%-JOINT-power reach (same formula as ab_split_budget.py; exact-window toy world,
# i.e. widen=1 -- the real-life note uses widen=3 for the +-2 sigma_M window freedom)
k_eff = K_BKG + 1
def Phi(x): return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
def reach50(f):
    lo, hi = 3.0, 20.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if Phi(math.sqrt(f) * mid - Z_CUT) * Phi(math.sqrt(1 - f) * mid - float(zthr(k_eff))) < 0.5:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
reach = {"single": Z_SINGLE, "50/50": reach50(0.5), "opt": reach50(F_OPT)}

fig2, ax2 = plt.subplots(figsize=(9.5, 6))
for y, lbl, col in [(pw_single, f"single-stage full scan (needs $Z\\geq${Z_SINGLE:.2f})", "#d9a441"),
                    (pw_opt,   f"two-stage f = {F_OPT:g}, $Z_{{cut}}$ = {Z_CUT:g}", "#4c78a8"),
                    (pw_5050,  f"two-stage 50/50, $Z_{{cut}}$ = {Z_CUT:g}", "#d1495b")]:
    ax2.plot(mus, y, "o-", ms=4, lw=2, color=col, label=lbl)
for key, col in [("single", "#d9a441"), ("opt", "#4c78a8"), ("50/50", "#d1495b")]:
    ax2.axvline(reach[key], color=col, ls=":", lw=1.3)
ax2.axhline(0.5, color="#888888", ls="--", lw=0.9)
ax2.set_xlabel(r"injected signal strength: full-dataset local significance  $Z_{full}$")
ax2.set_ylabel(r"discovery power  P(claim $5\sigma$ global)")
ax2.set_ylim(0, 1.02)
ax2.grid(ls=":", alpha=0.35)
ax2.set_title("Toy power of the three procedures  "
              f"(N = {N:,} looks, {1500:,} toys/point)\n"
              "dotted verticals = analytic 50%-joint-power reach (Phi-product formula)")
ax2.legend(loc="upper left", fontsize=9)
fig2.tight_layout()
fig2.savefig(os.path.join(OUT, "ab_toys_power.png"), dpi=130)
z50 = {lbl: float(np.interp(0.5, y, mus)) for lbl, y in
       [("single", pw_single), ("opt", pw_opt), ("50/50", pw_5050)]}
print("wrote ab_toys_power.png")
print(f"  50%-power points (toys): single {z50['single']:.2f}, f=0.22 {z50['opt']:.2f}, "
      f"50/50 {z50['50/50']:.2f}   (analytic reach {reach['single']:.2f} / "
      f"{reach['opt']:.2f} / {reach['50/50']:.2f})")

# ================================================================ 2b. symmetrized swap (A<->B)
# Run BOTH directions -- explore A confirm B, and explore B confirm A -- and claim on the union.
# Both directions pre-register, so the bar is set by the TOTAL unblinded windows k_A + k_B: the
# second chance is paid for with twice the coins. Validates the union formula of
# ab_split_budget.py, which is the number TWO_STAGE_UNBLINDING.md quotes.
def toy_power_sym(mu, f, zcut, ntoy=1500):
    gA = rng.standard_normal((ntoy, N)); gB = rng.standard_normal((ntoy, N))
    gA[:, 0] += math.sqrt(f) * mu; gB[:, 0] += math.sqrt(1 - f) * mu
    selA, selB = gA >= zcut, gB >= zcut
    ktot = selA.sum(axis=1) + selB.sum(axis=1)
    best = np.maximum(np.where(selA, gB, -np.inf).max(axis=1),
                      np.where(selB, gA, -np.inf).max(axis=1))
    return ((ktot > 0) & (best >= zthr(ktot))).mean()

def reach50_sym(f, zcut=Z_CUT):
    zr = float(zthr(2 * K_BKG + 1))
    lo, hi = 3.0, 20.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        ac, bc = Phi(math.sqrt(f) * m - zcut), Phi(math.sqrt(1 - f) * m - zcut)
        ar, br = Phi(math.sqrt(f) * m - zr),   Phi(math.sqrt(1 - f) * m - zr)
        if ac * br + bc * ar - ar * br < 0.5: lo = m
        else: hi = m
    return 0.5 * (lo + hi)

pw_sym = [toy_power_sym(mu, 0.5, Z_CUT) for mu in mus]
z50_sym = float(np.interp(0.5, pw_sym, mus))
print(f"  symmetrized swap (50/50, union rule): toys {z50_sym:.2f} vs analytic "
      f"{reach50_sym(0.5):.2f}   [one-way 50/50: toys {z50['50/50']:.2f}, "
      f"analytic {reach['50/50']:.2f}]")
print(f"  -> the swap rescues the naive 50/50 split by "
      f"{reach['50/50'] - reach50_sym(0.5):.2f} sigma, still "
      f"{reach50_sym(0.5) - Z_SINGLE:+.2f} vs the single stage")

# ================================================================ 3. spectrum illustration
SIGMA_REL, F_ILL = 0.05, 0.25
edges = np.geomspace(200, 4000, 220)
ctr = np.sqrt(edges[:-1] * edges[1:]); width = np.diff(edges)
bkg_full = 3e6 * np.exp(-ctr / 350.0) * width / ctr      # falling toy background, ~2e5 events

def scan(counts, expect):
    """Sliding +-1 sigma_M window counting significance z(m)."""
    z = np.zeros_like(ctr)
    for j, m in enumerate(ctr):
        w = np.abs(ctr - m) < SIGMA_REL * m
        n, b = counts[w].sum(), expect[w].sum()
        z[j] = (n - b) / math.sqrt(b) if b > 0 else 0.0
    return z

def split_toy(seed, sig_mass=None, sig_zfull=0.0):
    r = np.random.default_rng(seed)
    exp_full = bkg_full.copy()
    mu = np.zeros_like(ctr)
    if sig_mass:
        w = np.abs(ctr - sig_mass) < SIGMA_REL * sig_mass
        gauss = np.exp(-0.5 * ((ctr - sig_mass) / (SIGMA_REL * sig_mass)) ** 2)
        mu = gauss / gauss[w].sum() * sig_zfull * math.sqrt(bkg_full[w].sum())
    nA = r.poisson(F_ILL * (bkg_full + mu))
    nB = r.poisson((1 - F_ILL) * (bkg_full + mu))
    return scan(nA, F_ILL * exp_full), scan(nB, (1 - F_ILL) * exp_full)

# find a background-only toy whose A-half fluctuates above Z_cut somewhere
seed = 0
while True:
    zA, zB = split_toy(seed)
    if zA.max() >= Z_CUT + 0.2: break
    seed += 1
zA_sig, zB_sig = split_toy(998, sig_mass=1200.0, sig_zfull=7.0)

fig3, axs3 = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
for ax, (za, zb), title in [
        (axs3[0], (zA, zB), f"background-only toy (seed {seed}): the A-selection DIES in B"),
        (axs3[1], (zA_sig, zB_sig), "signal-injected toy ($Z_{full}$ = 7 at 1.2 TeV): CONFIRMS")]:
    ax.plot(ctr, za, lw=1.6, color="#4c78a8", label="stage A scan  $Z_A(m)$  (25% of data)")
    sel = za >= Z_CUT
    first = True
    for j in np.where(sel)[0]:
        lo, hi = ctr[j] * (1 - 2 * SIGMA_REL), ctr[j] * (1 + 2 * SIGMA_REL)
        ax.axvspan(lo, hi, color="#f0e6b8", alpha=0.6, zorder=0,
                   label="pre-registered window" if first else None)
        first = False
    inwin = np.zeros_like(sel)
    for j in np.where(sel)[0]:
        inwin |= np.abs(np.log(ctr / ctr[j])) < 2 * SIGMA_REL
    zb_m = np.where(inwin, zb, np.nan)
    ax.plot(ctr, zb_m, lw=2.4, color="#d1495b",
            label="stage B  $Z_B(m)$ -- unblinded ONLY here (75%)")
    ax.axhline(Z_CUT, color="#4c78a8", ls="--", lw=1.1)
    ax.text(210, Z_CUT + 0.12, f"$Z_{{cut}}$ = {Z_CUT:g}", color="#2f4b6e", fontsize=8.5)
    thr = float(zthr(max(int(sel.any()), np.count_nonzero(np.diff(np.where(sel)[0]) > 1) + 1
                         if sel.any() else 1)))
    ax.axhline(thr, color="#d1495b", ls="--", lw=1.1)
    ax.text(210, thr + 0.12, f"B claim threshold $\\sqrt{{25+2\\ln k}}$ = {thr:.1f}",
            color="#a02735", fontsize=8.5)
    ax.set_xscale("log"); ax.set_xlabel("mass [GeV]"); ax.set_title(title, fontsize=10.5)
    ax.grid(ls=":", alpha=0.3)
axs3[0].set_ylabel("local significance in a $\\pm\\sigma_M$ window")
axs3[0].legend(fontsize=8.5, loc="upper right")
fig3.suptitle("Two-stage unblinding on a toy spectrum: only A-selected windows are opened in B",
              fontsize=12)
fig3.tight_layout(rect=[0, 0, 1, 0.94])
fig3.savefig(os.path.join(OUT, "ab_toys_spectrum.png"), dpi=130)
print(f"wrote ab_toys_spectrum.png  (bkg toy: max Z_A = {zA.max():.2f} -> "
      f"Z_B in window = {np.nanmax(np.where(zA>=Z_CUT, zB, np.nan)):.2f};  "
      f"signal toy: Z_A = {zA_sig.max():.2f} -> Z_B = {zB_sig.max():.2f})")
