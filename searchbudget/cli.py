import argparse
import os
import shutil
import sys

from . import paths, registry, runner

DESCRIPTION = ("the search budget of the BSM resonance program: every table, figure, LaTeX "
               "fragment and generated report under results/")


def _stage_names(args):
    if getattr(args, "all", False) or not args.stages:
        return [s.name for s in registry.all_stages(network=False)]
    names = []
    for token in args.stages:
        if token in {s.group for s in registry.all_stages()}:
            names += [s.name for s in registry.all_stages(group=token)]
        else:
            try:
                registry.get(token)
            except KeyError:
                raise SystemExit(f"unknown stage or group {token!r}; try `search-budget list`")
            names.append(token)
    return list(dict.fromkeys(names))


def cmd_list(args):
    registry.load()
    width = max(len(s.name) for s in registry.all_stages())
    for group in registry.GROUPS:
        group_stages = sorted(registry.all_stages(group=group), key=lambda s: s.name)
        if not group_stages:
            continue
        print(f"\n{group}")
        for s in group_stages:
            mark = "*" if s.network else (" " if not runner.is_stale(s) else "~")
            print(f"  {mark} {s.name:{width}s}  {s.summary}")
    print("\n  ~ out of date    * needs the network, outside `run --all`")


def cmd_graph(args):
    registry.load()
    made = registry.producers()
    for name in registry.order([s.name for s in registry.all_stages()]):
        s = registry.get(name)
        upstream = sorted({made[spec] for spec in s.needs if spec in made} - {name})
        print(f"{name}")
        for spec in s.inputs:
            print(f"    <- {spec}")
        for dep in upstream:
            print(f"    <- {dep}")
        for spec in s.outputs:
            print(f"    -> {spec}")
        for spec in s.caches:
            print(f"    ~> {spec}   (cache: reused unless --refit)")


def cmd_check(args):
    problems = registry.check()
    for s in registry.all_stages():
        for spec in s.outputs + s.caches:
            if not os.path.exists(paths.resolve(spec)):
                problems.append(f"{s.name}: {spec} has not been built")
    for line in problems:
        print(line)
    print(f"{len(registry.STAGES)} stages, {len(registry.producers())} declared outputs, "
          f"{len(problems)} problem{'s' if len(problems) != 1 else ''}")
    return 1 if problems else 0


def cmd_run(args):
    names = _stage_names(args)
    if args.only:
        for name in names:
            runner.execute(name, runner.options(refit=args.refit, figonly=args.figonly))
        return 0
    runner.run(names, force=args.force, dry=args.dry_run, jobs=args.jobs,
               opts=runner.options(refit=args.refit, figonly=args.figonly))
    return 0


def cmd_clean(args):
    removed = 0
    for root, dirs, files in os.walk(paths.PACKAGE):
        for d in list(dirs):
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d))
                dirs.remove(d)
                removed += 1
    if args.outputs:
        for s in registry.all_stages():
            for spec in s.outputs + s.caches:
                target = paths.resolve(spec)
                if spec != "README.md" and os.path.exists(target):
                    os.remove(target)
                    removed += 1
    print(f"removed {removed} item{'s' if removed != 1 else ''}")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(prog="search-budget", description=DESCRIPTION)
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("list", help="every stage, its group and whether it is out of date")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("run", help="run stages and everything they depend on")
    p.add_argument("stages", nargs="*", help="stage or group names; none means all of them")
    p.add_argument("--all", action="store_true", help="every stage that needs no network")
    p.add_argument("-f", "--force", action="store_true", help="rebuild even when up to date")
    p.add_argument("-n", "--dry-run", action="store_true", help="print the plan and stop")
    p.add_argument("-j", "--jobs", type=int, default=1, help="stages to run concurrently")
    p.add_argument("--only", action="store_true",
                   help="run exactly these stages in this process, no dependencies")
    p.add_argument("--refit", action="store_true", help="redo the cached Monte Carlo")
    p.add_argument("--figonly", action="store_true", help="skip the tables, redraw the figures")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("graph", help="the dependency graph, stage by stage")
    p.set_defaults(func=cmd_graph)

    p = sub.add_parser("check", help="verify the registry and that every output exists")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("clean", help="drop __pycache__, and with --outputs the generated files")
    p.add_argument("--outputs", action="store_true", help="also delete everything under results/")
    p.set_defaults(func=cmd_clean)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    registry.load()
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
