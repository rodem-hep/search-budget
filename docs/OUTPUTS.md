# Output index

Every file under `results/`, with the script that writes it. "hand-written" means a report authored
against the generated tables, not a program output — edit it directly; nothing regenerates it.

## The note — `docs/note/`

`main.tex` plus `sections/*.tex` and `refs.bib`, built with `make note` (or `make` in that
directory). It pulls its figures straight from `results/plots/`, so rebuild those first if they are
stale. `main.pdf` is committed so the note can be read without a LaTeX installation.

## Reports — `results/overviews/`

| file | producer | content |
|---|---|---|
| `SEARCH_BUDGET.md` | `search_budget.py` | **the headline**: spectra, `N_trials`, `Z_local(5σ global)`, per-spectrum table |
| `EXCESS_COUNTING.md` | hand-written | expected vs observed 3σ/5σ excesses — the external check on `N` |
| `MAX_OF_GAUSSIANS.md` | hand-written | the statistics: argmax vs threshold vs Benjamini-Hochberg |
| `TWO_STAGE_UNBLINDING.md` | hand-written | the A/B split strategy, reach and caveats |

## Tables — `results/tables/`

| file | producer |
|---|---|
| `search_budget.csv` | `search_budget.py` — one row per spectrum: `r`, floor, window, `n_s`, envelope, models, selections |
| `search_budget_selections.csv` | `search_budget.py` — the event-selection multiplicity per spectrum |
| `combinatorial_budget.csv` | `combinatorial_budget.py` — 1094 (category, mass-group) rows |
| `composition_gap.txt` | `composition_gap.py` |
| `bh_fdr_scan.csv`, `bh_fdr_mc.npz` | `bh_fdr_ab.py` (the `.npz` is the Monte Carlo cache) |
| `bh_zcut_per_pe.csv` | `bh_zcut.py` |
| `bh_outliers_scan.npz` | `bh_fdr_outliers.py` (Monte Carlo cache) |

## Figures — `results/plots/`

| file | producer |
|---|---|
| `search_budget.png`, `scan_windows.png`, `model_observable_matrix.png` | `budget_plots.py` |
| `budget_waterfall.png` | `budget_waterfall.py` — `Z_local` vs `N` across granularity levels |
| `excess_counting.png` | `excess_counting.py` |
| `ab_split_reach.png`, `ab_split_crossover.png` | `ab_split_budget.py` |
| `ab_toys_background.png`, `ab_toys_power.png`, `ab_toys_spectrum.png` | `ab_split_toys.py` |
| `ab_split_outliers.png`, `ab_outliers_mechanism.png`, `ab_outliers_spectrum.png` | `ab_split_outliers.py` |

### `results/plots/max_of_gaussians/` — figures of `MAX_OF_GAUSSIANS.md`

| file | producer |
|---|---|
| `max_of_gaussians_light.png`, `signal_wins_the_max.png`, `ab_confirmation.png`, `threshold_scan.png`, `threshold_vs_argmax.png`, `roc_threshold_vs_argmax.png` | `max_of_gaussians_plots.py` (Parts I–III) |
| `bh_scan.png`, `bh_vs_argmax.png`, `roc_bh_vs_threshold.png` | `bh_fdr_ab.py` (Part IV) |
| `bh_zcut.png` | `bh_zcut.py` (Part IV) |
| `bh_outliers.png` | `bh_fdr_outliers.py` (Part IV, imperfect estimator; MC cached in `bh_outliers_scan.npz`) |

## Modules (no outputs of their own)

| file | what it owns |
|---|---|
| `bump_observables.py` | observables, resolutions, analyzable floors, published scan windows + sources, same-axis merges, the lepton-flavour split, the LEE math |
| `public_obs_map.py` | public model → spectrum map, published event-selection counts (`NSEL`) |
| `plot_style.py` | the shared palette and axis style of the statistics figures |
