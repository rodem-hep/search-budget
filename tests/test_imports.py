import ast
import builtins
import os

import pytest

from searchbudget import paths

MODULES = sorted(os.path.join(root, f)
                 for root, dirs, files in os.walk(paths.PACKAGE)
                 if "__pycache__" not in root
                 for f in files if f.endswith(".py"))

MODULE_GLOBALS = {"__file__", "__name__", "__doc__", "__package__"}


def _bound(tree):
    names = set(dir(builtins)) | MODULE_GLOBALS
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.asname or a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            names.update(a.asname or a.name for a in node.names)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            names.add(node.id)
        elif isinstance(node, ast.arg):
            names.add(node.arg)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            names.add(node.name)
        elif isinstance(node, (ast.Global, ast.Nonlocal)):
            names.update(node.names)
    return names


@pytest.mark.parametrize("path", MODULES, ids=lambda p: os.path.relpath(p, paths.PACKAGE))
def test_no_name_is_used_without_being_bound(path):
    tree = ast.parse(open(path).read(), path)
    used = {n.id for n in ast.walk(tree)
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    assert sorted(used - _bound(tree)) == []


@pytest.mark.parametrize("path", MODULES, ids=lambda p: os.path.relpath(p, paths.PACKAGE))
def test_nothing_reaches_outside_the_package_for_a_path(path):
    source = open(path).read()
    if path.endswith(("paths.py", "io.py", "deps.py", "runner.py", "cli.py")):
        return
    assert "sys.path.insert" not in source
    assert '"results"' not in source and "'results'" not in source
