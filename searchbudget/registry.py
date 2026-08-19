import dataclasses
import importlib
import pkgutil

from . import paths

GROUPS = ("budget", "census", "scan", "ab", "stats", "paper", "fetch")


@dataclasses.dataclass(frozen=True)
class Stage:
    name: str
    group: str
    summary: str
    module: str
    function: object
    outputs: tuple
    inputs: tuple
    needs: tuple
    caches: tuple
    network: bool

    @property
    def output_paths(self):
        return tuple(paths.resolve(o) for o in self.outputs)

    @property
    def input_paths(self):
        return tuple(paths.resolve(i) for i in self.inputs)

    @property
    def cache_paths(self):
        return tuple(paths.resolve(c) for c in self.caches)

    def run(self, options):
        return self.function(options)


STAGES = {}


def stage(name, group, summary, outputs, inputs=(), needs=(), caches=(), network=False):
    if group not in GROUPS:
        raise ValueError(f"unknown group {group!r}; one of {GROUPS}")

    def register(fn):
        if name in STAGES:
            raise ValueError(f"stage {name!r} declared twice")
        STAGES[name] = Stage(name, group, summary, fn.__module__, fn, tuple(outputs),
                             tuple(inputs), tuple(needs), tuple(caches), network)
        return fn

    return register


def load():
    if STAGES:
        return STAGES
    from . import stages as package
    for info in sorted(pkgutil.iter_modules(package.__path__), key=lambda i: i.name):
        importlib.import_module(f"{package.__name__}.{info.name}")
    return STAGES


def get(name):
    load()
    if name not in STAGES:
        raise KeyError(name)
    return STAGES[name]


def all_stages(group=None, network=None):
    load()
    out = list(STAGES.values())
    if group is not None:
        out = [s for s in out if s.group == group]
    if network is not None:
        out = [s for s in out if s.network is network]
    return out


def producers():
    load()
    out = {}
    for s in STAGES.values():
        for spec in s.outputs + s.caches:
            if spec in out:
                raise ValueError(f"{spec!r} is written by both {out[spec]!r} and {s.name!r}")
            out[spec] = s.name
    return out


def check():
    load()
    made = producers()
    problems = []
    for s in STAGES.values():
        for spec in s.needs:
            if spec not in made:
                problems.append(f"{s.name}: needs {spec!r}, which no stage writes")
            elif made[spec] == s.name:
                problems.append(f"{s.name}: needs {spec!r}, its own output")
    try:
        order([s.name for s in STAGES.values()])
    except ValueError as exc:
        problems.append(str(exc))
    return problems


def dependencies(name):
    made = producers()
    return sorted({made[spec] for spec in get(name).needs if spec in made})


def order(names):
    load()
    made = producers()
    seen, out, path = set(), [], []

    def visit(n):
        if n in seen:
            return
        if n in path:
            raise ValueError("cycle through " + " -> ".join(path[path.index(n):] + [n]))
        path.append(n)
        for spec in STAGES[n].needs:
            if spec in made and made[spec] != n:
                visit(made[spec])
        path.pop()
        seen.add(n)
        out.append(n)

    for n in names:
        visit(n)
    return out


def closure(names):
    return order(list(names))
