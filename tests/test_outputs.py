import filecmp
import os
import shutil
import subprocess
import sys

import pytest

from searchbudget import io, paths, registry

STAGES = registry.load()

STDLIB_STAGES = ["search-budget", "combinatorial-budget", "composition-gap", "published-census",
                 "census-budget", "scaled-scan", "two-body-matrix", "model-map", "census-bib",
                 "reported-excesses"]


def _built(spec):
    return os.path.exists(paths.resolve(spec))


@pytest.mark.parametrize("name", sorted(STAGES))
def test_declared_outputs_are_committed(name):
    stage = registry.get(name)
    if stage.network:
        pytest.skip("network stage; its output is an input to everything else")
    for spec in stage.outputs + stage.caches:
        assert _built(spec), f"{name} has never written {spec}"


def test_the_headline_numbers_are_what_the_paper_quotes():
    rows = io.read_rows(paths.table("search_budget.csv"))
    assert len(rows) == 46
    N = sum(float(r["ns_scan"]) for r in rows)
    assert round(N) == 3685
    assert sum(int(r["n_event_selections"]) for r in rows) == 94

    scan = {r["dataset"]: r for r in io.read_rows(paths.table("scan_summary.csv"))}
    run23 = scan["Run 2+3, ~400 fb-1"]
    assert int(run23["fittable_spectra"]) == 4438
    assert round(float(run23["N_trials"])) == 201136
    assert round(float(run23["lensed_N_trials"])) == 362815

    census = io.read_rows(paths.table("census_budget.csv"))
    priced = [r for r in census if float(r["n_s"]) > 0]
    assert len(priced) == 100
    assert round(sum(float(r["n_s"]) for r in priced)) == 7710


def test_the_uncertainty_total_is_quoted_on_every_basis():
    rows = {(r["source"], r["direction"]): r
            for r in io.read_rows(paths.table("budget_uncertainty.csv"))}
    total = rows[("total", "total")]
    assert total["dZ_model_space"] == "+0.185/-0.163"
    assert total["dZ_combinatorial_scan"] == "+0.231/-0.308"
    assert total["dZ_difference"] == "+0.123/-0.232"


def test_the_readme_block_is_generated_not_typed():
    text = open(paths.README).read()
    start = text.index("<!-- paper-numbers:start -->")
    block = text[start:text.index("<!-- paper-numbers:end -->")]
    assert "N = 3,685" in block
    assert "searchbudget/stages/paper_numbers.py" in block


@pytest.mark.slow
def test_a_rebuild_of_the_standard_library_stages_is_byte_identical(tmp_path):
    env = dict(os.environ, SEARCH_BUDGET_RESULTS=str(tmp_path))
    proc = subprocess.run([sys.executable, "-m", "searchbudget", "run", *STDLIB_STAGES, "--force"],
                          cwd=paths.ROOT, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    assert proc.returncode == 0, proc.stdout

    rebuilt = []
    for root, _dirs, files in os.walk(tmp_path):
        for f in files:
            rebuilt.append(os.path.relpath(os.path.join(root, f), tmp_path))
    assert rebuilt, "the rebuild wrote nothing"

    differ = [f for f in sorted(rebuilt)
              if not filecmp.cmp(os.path.join(tmp_path, f),
                                 os.path.join(paths.RESULTS, f), shallow=False)]
    assert differ == [], f"a rebuild changed {differ}"


def test_every_output_is_named_in_the_output_index():
    index = open(os.path.join(paths.ROOT, "docs", "OUTPUTS.md")).read()
    for stage in registry.all_stages():
        for spec in stage.outputs + stage.caches:
            assert os.path.basename(spec) in index, f"{spec} is missing from docs/OUTPUTS.md"


def test_every_stage_is_named_in_the_output_index():
    index = open(os.path.join(paths.ROOT, "docs", "OUTPUTS.md")).read()
    for stage in registry.all_stages():
        assert stage.name in index, f"stage {stage.name} is missing from docs/OUTPUTS.md"


@pytest.mark.slow
def test_a_cached_monte_carlo_stage_settles_after_one_run(tmp_path):
    stage = registry.get("bh-fdr-ab")
    tables = tmp_path / "tables"
    tables.mkdir(parents=True)
    for spec in stage.caches:
        shutil.copy(paths.resolve(spec), tables / os.path.basename(spec))
    env = dict(os.environ, SEARCH_BUDGET_RESULTS=str(tmp_path))

    def run(*extra):
        return subprocess.run([sys.executable, "-m", "searchbudget", "run", stage.name, *extra],
                              cwd=paths.ROOT, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True)

    first = run("--force")
    assert first.returncode == 0, first.stdout
    second = run()
    assert second.returncode == 0, second.stdout
    assert "everything up to date" in second.stdout, second.stdout
