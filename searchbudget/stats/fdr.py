import numpy as np


def smallest_uniforms(rng, size, ntot, m):
    g = np.cumsum(rng.exponential(size=(size, m)), axis=1)
    tot = g[:, -1] + rng.gamma(ntot + 1 - m, size=size)
    return g / tot[:, None]


def cutoffs(p_sorted, n):
    k = np.arange(1, p_sorted.shape[1] + 1)
    return np.minimum.accumulate((n * p_sorted / k)[:, ::-1], axis=1)[:, ::-1]


def accepted(p_sorted, n, q):
    return (cutoffs(p_sorted, n) <= q).sum(axis=1)


def accepted_grid(p_sorted, n, qgrid):
    c = cutoffs(p_sorted, n)
    return np.stack([(c <= q).sum(axis=1) for q in qgrid], axis=1)


def accepted_searchsorted(p_sorted, n, qgrid):
    c = cutoffs(p_sorted, n)
    return np.stack([np.searchsorted(row, qgrid, side="right") for row in c])
