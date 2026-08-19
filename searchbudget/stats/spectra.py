import math

import numpy as np

NORM, SCALE = 3e6, 350.0


def grid(lo, hi, nedges):
    edges = np.geomspace(lo, hi, nedges)
    ctr = np.sqrt(edges[:-1] * edges[1:])
    return edges, ctr, np.diff(edges)


def background(ctr, width, norm=NORM, scale=SCALE):
    return norm * np.exp(-ctr / scale) * width / ctr


def window_scan(ctr, counts, expect, sigma_rel):
    z = np.zeros_like(ctr)
    for j, m in enumerate(ctr):
        w = np.abs(ctr - m) < sigma_rel * m
        n, b = counts[w].sum(), expect[w].sum()
        z[j] = (n - b) / math.sqrt(b) if b > 0 else 0.0
    return z


def bump(ctr, bkg, mass, z_full, sigma_rel):
    w = np.abs(ctr - mass) < sigma_rel * mass
    gauss = np.exp(-0.5 * ((ctr - mass) / (sigma_rel * mass)) ** 2)
    return gauss / gauss[w].sum() * z_full * math.sqrt(bkg[w].sum())


def split_toy(seed, ctr, bkg, frac, sigma_rel, sig_mass=None, sig_zfull=0.0):
    rng = np.random.default_rng(seed)
    mu = bump(ctr, bkg, sig_mass, sig_zfull, sigma_rel) if sig_mass else np.zeros_like(ctr)
    nA = rng.poisson(frac * (bkg + mu))
    nB = rng.poisson((1 - frac) * (bkg + mu))
    return (window_scan(ctr, nA, frac * bkg, sigma_rel),
            window_scan(ctr, nB, (1 - frac) * bkg, sigma_rel))


def first_toy_above(z_cut, margin, ctr, bkg, frac, sigma_rel, limit=10000):
    seed = 0
    while seed < limit:
        zA, zB = split_toy(seed, ctr, bkg, frac, sigma_rel)
        if zA.max() >= z_cut + margin:
            return seed, zA, zB
        seed += 1
    raise RuntimeError(f"no background toy reached Z_A >= {z_cut + margin}")


def in_window(ctr, selected, sigma_rel, halfwidth=2.0):
    mask = np.zeros_like(selected)
    for j in np.where(selected)[0]:
        mask |= np.abs(np.log(ctr / ctr[j])) < halfwidth * sigma_rel
    return mask


def island_count(selected):
    if not selected.any():
        return 1
    return max(1, np.count_nonzero(np.diff(np.where(selected)[0]) > 1) + 1)
