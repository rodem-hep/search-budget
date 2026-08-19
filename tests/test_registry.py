import os

import pytest

from searchbudget import paths, registry, runner

STAGES = registry.load()


def test_every_group_has_stages():
    for group in registry.GROUPS:
        assert registry.all_stages(group=group), f"no stage in group {group!r}"


def test_registry_is_consistent():
    assert registry.check() == []


def test_every_stage_declares_an_output():
    for stage in STAGES.values():
        assert stage.outputs, f"{stage.name} declares no output"


def test_outputs_resolve_inside_the_repository():
    for stage in STAGES.values():
        for spec, path in zip(stage.outputs + stage.caches,
                              stage.output_paths + stage.cache_paths):
            assert os.path.abspath(path).startswith(paths.ROOT), f"{stage.name}: {spec}"


def test_declared_inputs_exist():
    for stage in STAGES.values():
        for spec, path in zip(stage.inputs, stage.input_paths):
            assert os.path.exists(path), f"{stage.name}: missing input {spec}"


def test_order_puts_producers_before_consumers():
    order = registry.order(sorted(STAGES))
    seen = set()
    for name in order:
        for dep in registry.dependencies(name):
            assert dep in seen, f"{name} runs before its dependency {dep}"
        seen.add(name)


def test_network_stages_are_outside_the_default_run():
    default = {s.name for s in registry.all_stages(network=False)}
    assert "fetch-census-meta" not in default
    assert "fetch-census-abstracts" not in default
    assert len(default) == len(STAGES) - 2


@pytest.mark.parametrize("name", sorted(STAGES))
def test_prerequisites_include_the_stage_source(name):
    stage = registry.get(name)
    sources = runner.prerequisites(stage)
    assert any(s.endswith(stage.module.rsplit(".", 1)[-1] + ".py") for s in sources)


def test_a_cache_never_makes_its_stage_look_stale():
    for stage in STAGES.values():
        for spec in stage.caches:
            assert spec not in stage.outputs, (
                f"{stage.name}: {spec} is a cache, so it must not gate freshness")


def test_a_cache_is_never_its_own_prerequisite():
    for stage in STAGES.values():
        pre = set(runner.prerequisites(stage))
        for path in stage.cache_paths + stage.output_paths:
            assert path not in pre, f"{stage.name} waits on a file it writes: {path}"


def test_the_parallel_schedule_never_starts_a_stage_before_its_dependencies():
    todo = registry.order(sorted(STAGES))
    waits_on = runner.blocking(todo)
    done = set()
    for name in todo:
        assert waits_on[name] <= done, f"{name} is scheduled before {waits_on[name] - done}"
        done.add(name)
    assert set(waits_on) == set(todo)
