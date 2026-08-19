import collections

from .bump_observables import CANON_ORDER, canon, ns_scan
from .public_obs_map import PUBLIC_OBS, nsel


def models_by_spectrum():
    out = collections.defaultdict(set)
    for model, spectra in PUBLIC_OBS.items():
        for obs in spectra:
            out[canon(obs)].add(model)
    return out


def canonical_order(models=None):
    models = models_by_spectrum() if models is None else models
    order = [o for o in CANON_ORDER if o == canon(o) and o in models]
    order += [o for o in sorted(models) if o not in order]
    return order


def sorted_spectra():
    return sorted({canon(o) for spectra in PUBLIC_OBS.values() for o in spectra})


def N_inclusive(order):
    return sum(ns_scan(o) for o in order)


def N_selections(order):
    return sum(nsel(o) * ns_scan(o) for o in order)


def n_channels(order):
    return sum(nsel(o) for o in order)
