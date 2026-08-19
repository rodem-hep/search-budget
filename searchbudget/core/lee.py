import math

from .bump_observables import z_local_for_global5

SQRT2 = math.sqrt(2.0)
SQRT2PI = math.sqrt(2.0 * math.pi)
Z_GLOBAL = 5.0


def p1(z):
    return 0.5 * math.erfc(z / SQRT2)


def Phi(x):
    return 0.5 * (1.0 + math.erf(x / SQRT2))


def z_local(N):
    return z_local_for_global5(N)


def z_global(N, zglob=Z_GLOBAL):
    x = zglob * zglob - 2.0 * math.log(N)
    return math.sqrt(x) if x > 0 else None


def claim_bar(k, widen=1.0, zglob=Z_GLOBAL):
    return math.sqrt(zglob * zglob + 2.0 * math.log(widen * max(k, 1.0)))


def bisect(fn, target, lo, hi, steps=60):
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if fn(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def phi_inv(p, lo=-10.0, hi=10.0, steps=80):
    return bisect(Phi, p, lo, hi, steps)


def z_exact(N, zglob=Z_GLOBAL, steps=200):
    target = p1(zglob) / N
    return bisect(lambda z: -p1(z), -target, 1.0, 25.0, steps)


def upcrossing_factor(N):
    return z_local(N) / SQRT2PI


def band(N, factor=2.0, label="Z_local(5s global)"):
    lo, hi = N / factor, N * factor
    return (f"N = {N:,.0f}  (r x{1/factor:g}..x{factor:g} -> {lo:,.0f}-{hi:,.0f});  "
            f"{label} = {z_local(N):.2f}  ({z_local(lo):.2f}-{z_local(hi):.2f})")


def merge_segments(segments):
    out = []
    for lo, hi in sorted(segments):
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out
