# Output index

Every file under `results/`, with the script that writes it. "hand-written" means a report authored
against the generated tables, not a program output — edit it directly; nothing regenerates it.

Three files under `data/` are inputs rather than outputs, and are committed: `published_spectra.csv`
(the census, curated by hand), `census_papers.csv` (its bibliographic details) and
`census_abstracts.csv` (the paper abstracts). The latter two are refreshed from the arXiv API by
`fetch_census_meta.py` and `fetch_census_abstracts.py`, the only scripts that need the network and
the only ones that are not part of `make all`.

## Reports — `results/overviews/`

| file | producer | content |
|---|---|---|
| `SEARCH_BUDGET.md` | `search_budget.py` | **the headline**: spectra, `N_trials`, `Z_local(5σ global)`, per-spectrum table |
| `BUDGET_UNCERTAINTY.md` | `budget_uncertainty.py` | what every declared input is worth on the bar, on both counted bases and on the difference between them |
| `PUBLISHED_CENSUS.md` | `published_census.py` | the ATLAS resonance-search publication record: 86 spectra, 290 papers, recency and Run-3 coverage |
| `CENSUS_BUDGET.md` | `census_budget.py` | the same census priced in trials: `N` and `Z_local` from the publication record, on both the published-search and the axes-scanned basis |
| `CENSUS_REFERENCES.md` | `export_census_bib.py` | all 290 census papers written out in full, numbered, under the spectrum each is counted against |
| `EXCESS_COUNTING.md` | hand-written | expected vs observed 3σ/5σ excesses — the external check on `N` |
| `REPORTED_EXCESSES.md` | `reported_excesses.py` | the observed side made concrete: every excess ≥3σ reported in the 290 census abstracts |
| `ESTIMATOR_DEFECTS.md` | `estimator_defects.py` | the published defect rate of one network scanner, the trials inflation it implies, and what that does to each procedure |
| `MAX_OF_GAUSSIANS.md` | hand-written | the statistics: argmax vs threshold vs Benjamini-Hochberg |
| `TWO_STAGE_UNBLINDING.md` | hand-written | the A/B split strategy, reach and caveats, priced on both the model space and the combinatorial scan |

## Tables — `results/tables/`

| file | producer |
|---|---|
| `search_budget.csv` | `search_budget.py` — one row per spectrum: `r`, floor, window, `n_s`, envelope, models, selections |
| `budget_uncertainty.csv` | `budget_uncertainty.py` — one row per source × direction: the varied `N` on each basis and the shift in `Z_local`, plus the per-source envelopes and the quadrature total |
| `search_budget_selections.csv` | `search_budget.py` — the event-selection multiplicity per spectrum |
| `combinatorial_budget.csv` | `combinatorial_budget.py` — 1094 (category, mass-group) rows, with the window truncated at its one-event mass and `n_s = 0` where the histogram cannot be fitted |
| `published_census.csv` | `published_census.py` |
| `census_budget.csv` | `census_budget.py` — one row per (published search, axis): the range it scanned, where that range came from, `r` and `n_s` |
| `reported_excesses.csv` | `reported_excesses.py` — per paper: largest quoted local/global significance and the abstract sentence behind it |
| `composition_gap.txt` | `composition_gap.py` (also writes the LaTeX fragment below) |
| `scaled_scan.csv`, `priority_scan.csv`, `lens_scan.csv`, `scan_summary.csv`, `scaled_scan.txt` | `scaled_scan.py` — the ten-object scan with MET in the masses: one row per variant, per composition what the trials budget keeps or drops, per (dataset, selection lens) how many views it could reach, its efficiency, how many survive the statistics requirement and what they cost, one summary row per dataset (`scan_summary.csv`: combinations, categories, fittable spectra, tiers, motivated share, lens totals, the yield-anchor variations), and the report |
| `two_body_matrix.csv` | `two_body_matrix.py` — every pair of grid objects: scanned or not, and the looks it costs either way |
| `model_spectrum_map.csv` | `export_model_map_tex.py` — per spectrum: `n_s`, the model classes pointing at it, the event-selection multiplicity and what those selections are |
| `bh_fdr_scan.csv`, `bh_fdr_mc.npz` | `bh_fdr_ab.py` (the `.npz` is the Monte Carlo cache) |
| `bh_zcut_per_pe.csv` | `bh_zcut.py` |
| `bh_outliers_scan.npz` | `bh_fdr_outliers.py` (Monte Carlo cache) |
| `selection_rules.csv` | `bh_fdr_outliers.py` — per (estimator defect, signal strength): the confirmation probability of each rule at the common false-alarm budget, the nominal FDR that buys it, and how often the argmax *is* the signal |
| `ab_split_scan.csv` | `ab_split_budget.py` — the two-stage design priced on each basis: single-stage bar, optimised and 50/50 reach, the pre-registered list at `Z_cut = 3` with its claim bar and background-only toys, and the break-even trials inflation `R*` |
| `estimator_defects.csv` | `estimator_defects.py` — the published defect rates, the looks behind them and everything derived from them |
| `ab_split_toys.csv` | `ab_split_toys.py` — the 50%-power reach of each procedure, toys against the analytic formula |
| `ab_guard_toys.csv` | `ab_spurious_guard.py` — the two toy spectra of the guard figure: what stage A flagged, what stage B showed, and the claim bar |

## LaTeX fragments — `results/tex/`

Generated for a consuming document to `\input`, the same way `results/plots/` is generated for it to
`\includegraphics`.

| file | producer |
|---|---|
| `composition_appendix.tex` | `composition_gap.py` — all 82 reachable object compositions with the published search covering each, the row-by-row form of the coverage summary |
| `census_refs.bib` | `export_census_bib.py` — a BibTeX entry per census paper, keyed `census:<arxiv-id>` |
| `census_appendix.tex` | `export_census_bib.py` — the census spectrum by spectrum, citing every paper behind each |
| `two_body_matrix.tex` | `two_body_matrix.py` — the two-body object grid priced in trials (a costed version of arXiv:1907.06659 Table 14); a bare `tabular`, so the consuming document supplies the float and caption |
| `uncertainty_table.tex` | `budget_uncertainty.py` — the uncertainty budget as a bare `tabular` (needs `booktabs` and `array` in the consuming document) |
| `model_map_appendix.tex` | `export_model_map_tex.py` — one row per spectrum: the model classes behind it, and the event selections that make the inclusive spectra into 94 channels |

## Figures — `results/plots/`

| file | producer |
|---|---|
| `search_budget.png`, `scan_windows.png`, `model_observable_matrix.png` | `budget_plots.py` |
| `budget_waterfall.png` | `budget_waterfall.py` — `Z_local` vs `N` across granularity levels |
| `excess_counting.png` | `excess_counting.py` |
| `ab_split_reach.png`, `ab_split_crossover.png` | `ab_split_budget.py` |
| `ab_toys_background.png`, `ab_toys_power.png`, `ab_toys_spectrum.png` | `ab_split_toys.py` |
| `ab_split_outliers.png`, `ab_outliers_mechanism.png`, `ab_outliers_spectrum.png` | `ab_split_outliers.py` |
| `ab_spurious_guard.png` | `ab_spurious_guard.py` — the split acting as a spurious-signal guard, on two toy spectra |

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
| `scan_alphabet.py` | the ten object types of the hypothetical scan, their derived `sigma`, the datasets it is priced on, and the four selection lenses; shared by `scaled_scan.py` and `budget_uncertainty.py` so both price the same scan |
| `yield_model.py` | the fittability requirement (≥100 events, ≥25 elements of ≥1 event), the declared background-yield model behind it and the per-lens efficiencies. Run it for the calibration table |
| `plot_style.py` | the shared palette and axis style of the statistics figures |
| `obs_labels.py` | the observable keys rendered as physics, shared by the figures and the LaTeX tables (standard library only) |
