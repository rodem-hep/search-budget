import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE = os.path.join(ROOT, "searchbudget")


def _env(name, default):
    return os.path.abspath(os.environ.get(name, default))


DATA = _env("SEARCH_BUDGET_DATA", os.path.join(ROOT, "data"))
RESULTS = _env("SEARCH_BUDGET_RESULTS", os.path.join(ROOT, "results"))
TABLES = os.path.join(RESULTS, "tables")
PLOTS = os.path.join(RESULTS, "plots")
GAUSSIANS = os.path.join(PLOTS, "max_of_gaussians")
TEX = os.path.join(RESULTS, "tex")
OVERVIEWS = os.path.join(RESULTS, "overviews")
README = os.path.join(ROOT, "README.md")

_TREES = {"tables": TABLES, "plots": PLOTS, "tex": TEX, "overviews": OVERVIEWS, "data": DATA}


def data(*a):
    return os.path.join(DATA, *a)


def table(*a):
    return os.path.join(TABLES, *a)


def plot(*a):
    return os.path.join(PLOTS, *a)


def gaussians(*a):
    return os.path.join(GAUSSIANS, *a)


def tex(*a):
    return os.path.join(TEX, *a)


def overview(*a):
    return os.path.join(OVERVIEWS, *a)


def resolve(spec):
    head, _, rest = spec.partition("/")
    if rest and head in _TREES:
        return os.path.join(_TREES[head], *rest.split("/"))
    return os.path.join(ROOT, *spec.split("/"))


def rel(path):
    return os.path.relpath(os.path.abspath(path), ROOT).replace(os.sep, "/")
