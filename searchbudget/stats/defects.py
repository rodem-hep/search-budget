import functools
import math

import numpy as np


@functools.lru_cache(maxsize=None)
def _quadrature(order):
    return np.polynomial.legendre.leggauss(order)


def halfnormal(beta, order, zmax=10.0):
    gx, gw = _quadrature(order)
    hi = zmax * beta
    b = 0.5 * hi * (gx + 1.0)
    w = 0.5 * hi * gw * math.sqrt(2.0 / math.pi) / beta * np.exp(-0.5 * (b / beta) ** 2)
    return b, w / w.sum()
