# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A data-analysis repo (no application, no test suite, no package): the **public half** of a
two-repository study of the BSM resonance search program. Three studies, all reproducible from a
clone with no external input of any kind:

1. **the search budget** — how many independent bump hunts the public model space implies
   (46 spectra, `N ≈ 3.7k`, `Z_local ≈ 6.44`), the combinatorial-scan variant (`N ≈ 143k`), and the
   never-scanned-composition audit (57 of 82);
2. **two-stage A/B unblinding** — does splitting the dataset make the trials factor cheaper, or only
   more defensible (it costs ~0.5σ and buys countability);
3. **selection rules** — argmax vs fixed threshold vs Benjamini-Hochberg on one ROC, with a perfect
   and an imperfect significance estimator.

The **ATLAS-internal half** is the companion `BSM_signals` repo: the verified MC signal catalogue and
which of these 46 spectra it populates. That repo vendors two modules from here (see below); this one
imports nothing from it and must never do so.

`README.md` documents the physics; `docs/METHOD_NOTES.md` carries every counting convention and the
reasoning behind it. Outputs in `results/` are **committed**, so a changed constant shows up as a
diff in the regenerated tables and figures — regenerate them in the same change.

## Environments

| task | interpreter |
|---|---|
| `search_budget.py`, `combinatorial_budget.py`, `composition_gap.py` | system `python3` (stdlib only) |
| everything else (numpy / scipy / matplotlib) | any interpreter with `requirements.txt` installed — pass it as `make all PY=...` |

All scripts anchor paths to `ROOT` derived from `__file__`, so the working directory does not matter;
they add `ROOT/"scripts"` to `sys.path` and import the two modules directly.

**Run the Monte Carlo on a compute node, never the login node** (`ab_split_toys.py`, `bh_zcut.py`,
`bh_fdr_ab.py --refit`, `bh_fdr_outliers.py`):
`srun --partition=shared-cpu --time=00:40:00 --mem=8G <py> scripts/bh_zcut.py`.

## Pipeline DAG

```
bump_observables.py  ─┬─► search_budget.py ──► tables/search_budget{,_selections}.csv
public_obs_map.py    ─┤                        overviews/SEARCH_BUDGET.md
                      │                              │
                      ├─► budget_plots.py            ├─► budget_waterfall.py ──► plots/budget_waterfall.png
                      │     plots/{search_budget,     └─► excess_counting.py ──► plots/excess_counting.png
                      │            scan_windows,
                      │            model_observable_matrix}.png
                      ├─► ab_split_budget.py · ab_split_toys.py · ab_split_outliers.py ──► plots/ab_*.png
                      └─► combinatorial_budget.py ──► tables/combinatorial_budget.csv
                                                            └─► composition_gap.py ──► composition_gap.txt

plot_style.py ──► max_of_gaussians_plots.py · bh_fdr_ab.py · bh_zcut.py · bh_fdr_outliers.py
                                                    └─► plots/max_of_gaussians/
```

## Invariants to preserve

- **`bump_observables.py` is the single source of truth** for bump observables, floors,
  resolutions, published scan windows, same-axis merges (`canon()`), the lepton-flavour split
  (`flav_channels()`) and the LEE math. Import it; never redefine a floor or a window in a
  consumer. Its only import-time side effect is a set of asserts that keep the flavour layer
  complete (every leaf has a `SCAN`+`FLOOR`, no inclusive parent keeps one).
- **`public_obs_map.py` owns the model → spectrum map and `NSEL`.** Lepton flavour is a mass axis,
  not a selection, so it must NOT also appear as a multiplier in `NSEL`. The invariant that catches
  a double count: the selections-level channel count stays at **94**.
- **These two modules are vendored by the internal repo.** Any edit here changes a number there, and
  its `make check-vendor` will fail until it re-syncs. That is intended — but keep the two files
  **free of internal information**, comments included: they must state published-search facts only,
  never which signal samples a collaboration has produced.
- **The budget window is the published-search `SCAN` window.** Never a range inferred from which
  signal samples exist: that understates sparsely sampled channels and makes the exercise circular.
- **Quote the combinatorial category and spectrum counts in matching units.** 150 multiplicity
  categories become 204 after the OS/SS split; 1094 (category, mass-group) rows become 1502
  histograms. `N = 143k` counts the split histograms, so pairing an unsplit category count with it
  is the easy mistake.
- **Monte Carlo is seeded** and the slow scans are cached in `results/tables/*.npz`, so the figures
  are bit-for-bit reproducible. A changed figure means a changed input, not a moved seed. Pass
  `--refit` to redo the cached scans.

## Commands

```bash
make all                 # budget + ab + stats  (make help for the stages)
make note                # docs/note/main.pdf   (needs pdflatex; figures from results/plots/)

python3 scripts/search_budget.py            # -> the tables + SEARCH_BUDGET.md (THE NUMBER)
python3 scripts/combinatorial_budget.py     # -> combinatorial_budget.csv
python3 scripts/composition_gap.py > results/tables/composition_gap.txt
$PY scripts/budget_plots.py                 # $PY = an interpreter with requirements.txt
```

`composition_gap.py` only prints and must be redirected to its `results/` path; every other script
writes its own outputs.

## Conventions

- Reports in `results/overviews/` are a mix of generated (`SEARCH_BUDGET.md`) and hand-written
  narrative (`MAX_OF_GAUSSIANS.md`, `TWO_STAGE_UNBLINDING.md`, `EXCESS_COUNTING.md`). Check the
  writing script before editing an overview by hand — see the `Writes:` line in each docstring.
- Script docstrings carry the authoritative `Reads:`/`Writes:` contract; keep them current when
  changing I/O.
- Output text is ASCII-ish by convention (`Z_local`, `5s global`, `sigma`) since tables are read in
  terminals; `README.md` and the LaTeX note use Unicode freely.
