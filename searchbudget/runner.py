import concurrent.futures
import os
import subprocess
import sys
import time
import types

from . import deps, paths, registry

DEFAULT_OPTIONS = dict(refit=False, figonly=False)


def options(**kw):
    merged = dict(DEFAULT_OPTIONS)
    merged.update(kw)
    return types.SimpleNamespace(**merged)


def _newest(files):
    stamps = [os.path.getmtime(f) for f in files if os.path.exists(f)]
    return max(stamps) if stamps else 0.0


def prerequisites(stage):
    return sorted(set(deps.sources(stage.module))
                  | set(stage.input_paths)
                  | {paths.resolve(spec) for spec in stage.needs})


def is_stale(stage):
    produced = stage.output_paths
    if any(not os.path.exists(p) for p in produced):
        return True
    return min(os.path.getmtime(p) for p in produced) < _newest(prerequisites(stage))


def missing_inputs(stage):
    return [p for p in stage.input_paths if not os.path.exists(p)]


def plan(names, force=False):
    todo = []
    for name in registry.order(names):
        stage = registry.get(name)
        upstream = any(dep in todo for dep in registry.dependencies(name))
        if force or upstream or is_stale(stage):
            todo.append(name)
    return todo


def _argv(name, opts):
    argv = [sys.executable, "-m", "searchbudget", "run", name, "--only", "--force"]
    if getattr(opts, "refit", False):
        argv.append("--refit")
    if getattr(opts, "figonly", False):
        argv.append("--figonly")
    return argv


def execute(name, opts):
    stage = registry.get(name)
    lack = missing_inputs(stage)
    if lack:
        raise SystemExit(f"{name}: missing input {', '.join(paths.rel(p) for p in lack)}")
    for spec in stage.outputs + stage.caches:
        parent = os.path.dirname(paths.resolve(spec))
        if parent:
            os.makedirs(parent, exist_ok=True)
    return stage.run(opts)


def _spawn(name, opts, capture):
    started = time.time()
    if capture:
        proc = subprocess.run(_argv(name, opts), cwd=paths.ROOT, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return proc.returncode, proc.stdout, time.time() - started
    proc = subprocess.run(_argv(name, opts), cwd=paths.ROOT)
    return proc.returncode, "", time.time() - started


def _serial(todo, opts):
    for i, name in enumerate(todo, 1):
        head = f"[{i}/{len(todo)}] {name}"
        print(f"{head}\n{'-' * len(head)}", flush=True)
        code, _, took = _spawn(name, opts, capture=False)
        if code:
            raise SystemExit(f"{name}: failed with exit code {code}")
        print(f"... {name} done in {took:.1f}s\n", flush=True)


def blocking(todo):
    inside = set(todo)
    return {n: set(registry.dependencies(n)) & inside for n in todo}


def _parallel(todo, opts, jobs):
    waits_on = blocking(todo)
    done, started, running = set(), set(), {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        while len(done) < len(todo):
            for name in todo:
                if len(running) >= jobs:
                    break
                if name in started or not waits_on[name] <= done:
                    continue
                started.add(name)
                running[pool.submit(_spawn, name, opts, True)] = name
            if not running:
                raise SystemExit(f"deadlocked with {sorted(set(todo) - done)} left to run")
            batch = concurrent.futures.wait(list(running),
                                            return_when=concurrent.futures.FIRST_COMPLETED)
            for future in batch.done:
                name = running.pop(future)
                code, output, took = future.result()
                done.add(name)
                head = f"[{len(done)}/{len(todo)}] {name}  ({took:.1f}s)"
                print(f"{head}\n{'-' * len(head)}\n{output}".rstrip() + "\n", flush=True)
                if code:
                    raise SystemExit(f"{name}: failed with exit code {code}")


def run(names, force=False, dry=False, jobs=1, opts=None):
    opts = opts or options()
    todo = plan(names, force=force)
    if not todo:
        print("everything up to date")
        return []
    if dry:
        for name in todo:
            print(f"{name:24s} -> " + ", ".join(registry.get(name).outputs))
        return todo
    started = time.time()
    if jobs > 1 and len(todo) > 1:
        _parallel(todo, opts, jobs)
    else:
        _serial(todo, opts)
    print(f"{len(todo)} stage{'s' if len(todo) > 1 else ''} in {time.time() - started:.1f}s")
    return todo
