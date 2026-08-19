import ast
import functools
import importlib.util
import os

PACKAGE = __name__.rsplit(".", 1)[0]


def _file(module):
    try:
        spec = importlib.util.find_spec(module)
    except (ImportError, ValueError):
        return None
    return spec.origin if spec and spec.origin and spec.origin.endswith(".py") else None


def _absolute(node, module):
    if isinstance(node, ast.Import):
        return [a.name for a in node.names]
    base = module.rsplit(".", node.level)[0] if node.level else ""
    root = f"{base}.{node.module}" if node.module else base
    return [root] + [f"{root}.{a.name}" for a in node.names]


@functools.lru_cache(maxsize=None)
def _imports(module):
    src = _file(module)
    if src is None:
        return ()
    with open(src) as fh:
        tree = ast.parse(fh.read(), src)
    found = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            found += [m for m in _absolute(node, module) if m.split(".")[0] == PACKAGE]
    return tuple(dict.fromkeys(found))


def sources(module):
    seen, out = set(), []
    stack = [module]
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.add(name)
        src = _file(name)
        if src is None:
            continue
        out.append(os.path.abspath(src))
        stack += [m for m in _imports(name) if m not in seen]
    return sorted(set(out))
