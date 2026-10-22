import math
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf
from scipy.stats import norm

from .. import io, paths
from ..core.lee import z_local
from ..registry import stage
from .inverse_lee import BASES

# The toy bump hunt of Appendix C: a background falling as m^-5, 5e7 events over 200-3000 GeV,
# binned at a quarter of the resolution sigma_M/M = 2 %, scanned by a Gaussian template whose
# strength is fitted by Poisson maximum likelihood with the background fixed at its true shape.
M_LO, M_HI, EVENTS, POWER = 200.0, 3000.0, 5e7, 5
R = 0.02
BIN = R / 4
HALF = 16                       # template half-width in bins: 4 sigma, kept wholly inside the range
NTOY, N_GV, C0 = 2_000_000, 500, 1.0
CHUNK = 250
Z_GRID = np.round(np.arange(3.0, 6.01, 0.25), 2)
CACHE = paths.table("gv_toys.npz")


def background():
    edges = np.linspace(math.log(M_LO), math.log(M_HI), int(round(math.log(M_HI / M_LO) / BIN)) + 1)
    # dN/dln m ∝ m^(1-P): integrate per bin in ln m
    w = np.exp(-(POWER - 1) * edges)
    per_bin = w[:-1] - w[1:]
    return EVENTS * per_bin / per_bin.sum()


def template():
    k = np.arange(-HALF, HALF + 1)
    s = 0.5 * (erf((k + 0.5) * BIN / (R * math.sqrt(2))) - erf((k - 0.5) * BIN / (R * math.sqrt(2))))
    return s / s.sum()


def scan(n, b, s, iterations=8):
    """Z = sign(mu_hat) sqrt(q0) at every template position, for a (toys, bins) block of counts."""
    win = 2 * HALF + 1
    nw = np.lib.stride_tricks.sliding_window_view(n, win, axis=1)            # (toys, pos, win)
    bw = np.lib.stride_tricks.sliding_window_view(b, win)                    # (pos, win)
    sb = s / bw
    # Newton on the concave score from the Gaussian estimate; mu kept where every lambda > 0
    mu = ((nw - bw) * sb).sum(-1) / (s * sb).sum(-1)
    floor = -0.999 / sb.max(-1)
    for _ in range(iterations):
        ratio = nw / (bw + mu[..., None] * s)
        score = (ratio * s).sum(-1) - 1.0
        curv = (ratio * s * s / (bw + mu[..., None] * s)).sum(-1)
        mu = np.maximum(mu + score / np.maximum(curv, 1e-300), floor)
    t = mu[..., None] * sb
    q0 = 2.0 * (np.where(nw > 0, nw * np.log1p(t), 0.0).sum(-1) - mu)
    return np.sign(mu) * np.sqrt(np.maximum(q0, 0.0))


def _block(args):
    seed, size = args
    b, s = background(), template()
    z = scan(np.random.default_rng(seed).poisson(b, size=(size, b.size)), b, s)
    up = ((z[:, :-1] < math.sqrt(C0)) & (z[:, 1:] >= math.sqrt(C0))).sum(1)
    return z.max(1).astype(np.float32), up.astype(np.int16)


def run_toys():
    seeds = np.random.SeedSequence(20261022).spawn(NTOY // CHUNK)
    jobs = [(sq, CHUNK) for sq in seeds]
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    with ProcessPoolExecutor(workers) as pool:
        out = list(pool.map(_block, jobs, chunksize=4))
    zmax = np.concatenate([o[0] for o in out])
    ups = np.concatenate([o[1] for o in out])
    np.savez_compressed(io.ensure(CACHE), zmax=zmax, ups=ups)
    return zmax, ups


def z_gv(n):
    """Local Z for a 5 sigma global claim from Eq. (gv) with N(5) = n."""
    n0 = n / (5 * math.sqrt(2 * math.pi))
    return brentq(lambda z: norm.sf(z) + n0 * math.exp(-z * z / 2) - norm.sf(5), 3, 10)


def z_exact(n):
    return float(norm.isf(-math.expm1(math.log1p(-norm.sf(5)) / n)))


@stage(
    name="gross-vitells",
    group="stats",
    summary="toy bump hunt: the Gross-Vitells trials factor against the resolution count, and Eq. (1) against Gross-Vitells",
    outputs=["tables/gv_toys.csv", "tex/gv_table.tex"],
    caches=["tables/gv_toys.npz"],
)
def main(options=None):
    if os.path.exists(CACHE) and not getattr(options, "refit", False):
        d = np.load(CACHE)
        zmax, ups = d["zmax"], d["ups"]
        print(f"loaded {zmax.size:,} scans from {os.path.basename(CACHE)} (--refit to redo)")
    else:
        zmax, ups = run_toys()

    positions = background().size - 2 * HALF
    ns = (positions - 1) * BIN / R
    n0 = ups[:N_GV].mean() * math.exp(C0 / 2)
    n0_all = ups.mean() * math.exp(C0 / 2)
    n5 = 5 * math.sqrt(2 * math.pi) * n0
    rice = 5 / (2 * math.sqrt(math.pi))
    print(f"scanned range {M_LO * math.exp(HALF * BIN):.1f}-{M_HI * math.exp(-HALF * BIN):.1f} GeV, "
          f"{positions} template positions, resolution count n_s = {ns:.1f}")
    print(f"N0 = {n0:.2f} from the first {N_GV} scans ({n0_all:.2f} from all {ups.size:,}); "
          f"N(5) = {n5:.1f}; N(5)/n_s = {n5 / ns:.3f} against Rice 5/(2 sqrt pi) = {rice:.3f}")
    for label, n in (("model space", BASES[0][1]), ("combinatorial scan", BASES[2][1])):
        print(f"  shift of the {label} bar: {math.log(n5 / ns) / z_local(n):+.3f}")

    rows = []
    print(f"\n{'Z':>5}{'p toys':>12}{'+-':>10}{'p GV':>12}{'p N(5) looks':>14}")
    for z in Z_GRID:
        hits = int((zmax >= z).sum())
        p = hits / zmax.size
        gv = norm.sf(z) + n0 * math.exp(-z * z / 2)
        fixed = -math.expm1(n5 * math.log1p(-norm.sf(z)))
        rows.append([z, hits, f"{p:.4g}", f"{math.sqrt(p * (1 - p) / zmax.size):.2g}", f"{gv:.4g}", f"{fixed:.4g}"])
        print(f"{z:5.2f}{p:12.3e}{math.sqrt(p * (1 - p) / zmax.size):10.1e}{gv:12.3e}{fixed:14.3e}")
    io.write_rows(paths.table("gv_toys.csv"), ["z", "toys_above", "p_toys", "p_toys_err", "p_gv", "p_fixed_n5"], rows)

    fmt_n = lambda n: f"${n / 10 ** int(math.log10(n)):.1f}\\times10^{{{int(math.log10(n))}}}$" \
        if n >= 10 ** 5 else f"${n:,.0f}$".replace(",", "\\,")
    table = [(label, n) for label, n, _ in BASES[:3]] + \
        [("toy, resolution count", ns), ("toy, measured $N(5)$", n5)]
    body = "\n".join(
        f"{label:<24} & {fmt_n(n)} & {z_local(n):.3f} & {z_gv(n):.3f} & ${z_gv(n) - z_local(n):+.3f}$ & "
        f"{z_exact(n):.3f} \\\\" for label, n in table)
    io.write_text(paths.tex("gv_table.tex"),
                  "% Generated by searchbudget/stages/gross_vitells.py. Do not edit: regenerate instead.\n"
                  "\\begin{tabular}{lrrrrr}\n\\toprule\n"
                  "trials input & $N$ & Eq.~\\eqref{eq:lee} & Gross--Vitells & difference & exact $1-(1-p)^N$ \\\\\n"
                  "\\midrule\n" + body + "\n\\bottomrule\n\\end{tabular}\n")
