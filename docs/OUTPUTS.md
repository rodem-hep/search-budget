# Output index

Every file under `results/`, with the stage that writes it. One file outside `results/` is
generated too: the "Every number in the paper" section of `README.md`, written by
`paper-numbers` between the `<!-- paper-numbers -->` markers from the tables below. Edit the
tables' producers, never that block.

Every producer named here is a stage name: run one with `search-budget run <stage>`, see them all
with `search-budget list`, and see this table's machine-readable form, straight from the stage
declarations rather than from prose, with `search-budget graph`.

"hand-written" means a report authored against the generated tables, not a program output — edit it
directly; nothing regenerates it.

Four files under `data/` are inputs rather than outputs, and are committed: `published_spectra.csv`
(the census, curated by hand), `census_papers.csv` (its bibliographic details),
`census_abstracts.csv` (the paper abstracts) and `model_papers.csv` (the bibliographic details of
the literature-sourced model classes' theory papers, ids declared in `public_obs_map.LITERATURE`).
The last three are refreshed from the arXiv API by `fetch-census-meta`, `fetch-census-abstracts`
and `fetch-model-meta`, the only stages that need the network and the only ones outside `make all`.

## Reports — `results/overviews/`

| file | producer | content |
|---|---|---|
| `SEARCH_BUDGET.md` | `search-budget` | **the headline**: spectra, `N_trials`, `Z_local(5σ global)`, per-spectrum table |
| `BUDGET_UNCERTAINTY.md` | `budget-uncertainty` | what every declared input is worth on the bar, on both counted bases and on the difference between them |
| `PUBLISHED_CENSUS.md` | `published-census` | the ATLAS resonance-search publication record: 86 spectra, 290 papers, recency and Run-3 coverage |
| `CENSUS_BUDGET.md` | `census-budget` | the same census priced in trials: `N` and `Z_local` from the publication record, on both the published-search and the axes-scanned basis |
| `CENSUS_REFERENCES.md` | `census-bib` | all 290 census papers written out in full, numbered, under the spectrum each is counted against |
| `EXCESS_COUNTING.md` | hand-written | expected vs observed 3σ/5σ excesses — the external check on `N` |
| `REPORTED_EXCESSES.md` | `reported-excesses` | the observed side made concrete: every excess ≥3σ reported in the 290 census abstracts |
| `ESTIMATOR_DEFECTS.md` | `estimator-defects` | the published defect rate of one network scanner, the trials inflation it implies, and what that does to each procedure |
| `MAX_OF_GAUSSIANS.md` | hand-written | the statistics: argmax vs threshold vs Benjamini-Hochberg |
| `TWO_STAGE_UNBLINDING.md` | hand-written | the A/B split strategy, reach and caveats, priced on both the model space and the combinatorial scan |

## Tables — `results/tables/`

| file | producer |
|---|---|
| `search_budget.csv` | `search-budget` — one row per spectrum: `r`, floor, window, `n_s`, envelope, models, selections |
| `budget_uncertainty.csv` | `budget-uncertainty` — one row per source × direction: the varied `N` on each basis and the shift in `Z_local`, plus the per-source envelopes and the quadrature total |
| `search_budget_selections.csv` | `search-budget` — the event-selection multiplicity per spectrum |
| `combinatorial_budget.csv` | `combinatorial-budget` — 1094 (category, mass-group) rows, with the window truncated at its one-event mass and `n_s = 0` where the histogram cannot be fitted |
| `published_census.csv` | `published-census` |
| `census_budget.csv` | `census-budget` — one row per (published search, axis): the range it scanned, where that range came from, `r` and `n_s` |
| `model_independence.csv` | `model-independence` — one row per charged (search, axis) pair: whether any of the search's papers states a model-independent result in its abstract, with the arXiv ids and matched phrases as evidence |
| `reported_excesses.csv` | `reported-excesses` — per paper: largest quoted local/global significance and the abstract sentence behind it |
| `composition_gap.txt` | `composition-gap` (also writes the LaTeX fragment below) |
| `scaled_scan.csv`, `priority_scan.csv`, `lens_scan.csv`, `scan_summary.csv`, `scaled_scan.txt` | `scaled-scan` — the ten-object scan with MET in the masses: one row per variant, per composition what the trials budget keeps or drops, per (dataset, selection lens) how many views it could reach, its efficiency, how many survive the statistics requirement and what they cost, one summary row per dataset (`scan_summary.csv`: combinations, categories, fittable spectra, tiers, motivated share, lens totals, the yield-anchor variations), and the report |
| `two_body_matrix.csv` | `two-body-matrix` — every pair of grid objects: scanned by a published search, priced in the budget with no published search (`ns_unscanned`, `axes_unscanned`), or a gap and what closing it costs |
| `model_spectrum_map.csv` | `model-map` — per spectrum: `n_s`, the model classes pointing at it, the event-selection multiplicity and what those selections are, whether a published ATLAS search already scans the axis, and whether the combinatorial scan can form it |
| `model_classes.csv` | `model-map` — one row per model class: whether it comes from the FeynRules/UFO database or the literature sweep, its defining arXiv references, and the spectra it points at |
| `unscanned_spectra.csv` | `unscanned-spectra` — the model-motivated spectra no published ATLAS search scans, in the scan's (category, mass-group) units: the pair-produced axis resolved into its per-category legs, each flagged fittable or not |
| `unscanned_scan_units.csv` | `unscanned-spectra` — the same question in the combinatorial scan's own units: of its 4,438 fittable spectra (one axis in one category, Run 2+3), how many are a final state some public model predicts with the mass on its resonant sub-system (`scaled_scan.FINAL_STATES`), and how many of those sit on an axis no published ATLAS search scans |
| `new_spectra.csv` | `unscanned-spectra` — the headline list: every spectrum of the combinatorial scan that a public model predicts (final state and resonant sub-system) and no published ATLAS search scans, one row per (mass group, category) with its looks, the budget axes it sits on and the model classes behind it |
| `bh_fdr_scan.csv`, `bh_fdr_mc.npz` | `bh-fdr-ab` (the `.npz` is the Monte Carlo cache) |
| `bh_zcut_per_pe.csv` | `bh-zcut` |
| `bh_outliers_scan.npz` | `bh-fdr-outliers` (Monte Carlo cache) |
| `selection_rules.csv` | `bh-fdr-outliers` — per (estimator defect, signal strength): the confirmation probability of each rule at the common false-alarm budget, the nominal FDR that buys it, and how often the argmax *is* the signal |
| `ab_split_scan.csv` | `ab-split-budget` — the two-stage design priced on each basis: single-stage bar, optimised and 50/50 reach, the pre-registered list at `Z_cut = 3` with its claim bar and background-only toys, and the break-even trials inflation `R*` |
| `estimator_defects.csv` | `estimator-defects` — the published defect rates, the looks behind them and everything derived from them |
| `ab_split_toys.csv` | `ab-split-toys` — the 50%-power reach of each procedure, toys against the analytic formula |
| `ab_guard_toys.csv` | `ab-spurious-guard` — the two toy spectra of the guard figure: what stage A flagged, what stage B showed, and the claim bar |

## LaTeX fragments — `results/tex/`

Generated for a consuming document to `\input`, the same way `results/plots/` is generated for it to
`\includegraphics`.

| file | producer |
|---|---|
| `composition_appendix.tex` | `composition-gap` — all 82 reachable object compositions with the published search covering each, the row-by-row form of the coverage summary |
| `census_refs.bib` | `census-bib` — a BibTeX entry per census paper, keyed `census:<arxiv-id>` |
| `census_appendix.tex` | `census-bib` — the census spectrum by spectrum, citing every paper behind each |
| `two_body_matrix.tex` | `two-body-matrix` — the two-body object grid priced in trials (a costed version of arXiv:1907.06659 Table 14); a bare `tabular`, so the consuming document supplies the float and caption |
| `uncertainty_table.tex` | `budget-uncertainty` — the uncertainty budget as a bare `tabular` (needs `booktabs` and `array` in the consuming document) |
| `new_spectra.tex` | `unscanned-spectra` — the spectra of the combinatorial scan a public model predicts and no published ATLAS search scans, grouped by budget axis; a bare `tabular`, so the consuming document supplies the float and caption |
| `model_map_appendix.tex` | `model-map` — one row per spectrum: the model classes behind it (literature-sourced classes cited), and the event selections that make the inclusive spectra into 111 channels |
| `model_refs.bib` | `model-map` — a BibTeX entry per theory paper behind the literature-sourced model classes, keyed `lit:<arxiv-id>` |

## Figures — `results/plots/`

| file | producer |
|---|---|
| `search_budget.png`, `scan_windows.png`, `model_observable_matrix.png` | `budget-plots` |
| `budget_waterfall.png` | `budget-waterfall` — `Z_local` vs `N` across granularity levels |
| `excess_counting.png` | `excess-counting` |
| `ab_split_reach.png`, `ab_split_crossover.png` | `ab-split-budget` |
| `ab_toys_background.png`, `ab_toys_power.png`, `ab_toys_spectrum.png` | `ab-split-toys` |
| `ab_split_outliers.png`, `ab_outliers_mechanism.png`, `ab_outliers_spectrum.png` | `ab-split-outliers` |
| `ab_spurious_guard.png` | `ab-spurious-guard` — the split acting as a spurious-signal guard, on two toy spectra |

### `results/plots/max_of_gaussians/` — figures of `MAX_OF_GAUSSIANS.md`

| file | producer |
|---|---|
| `max_of_gaussians_light.png`, `signal_wins_the_max.png`, `ab_confirmation.png`, `threshold_scan.png`, `threshold_vs_argmax.png`, `roc_threshold_vs_argmax.png` | `max-of-gaussians` (Parts I–III) |
| `bh_scan.png`, `bh_vs_argmax.png`, `roc_bh_vs_threshold.png` | `bh-fdr-ab` (Part IV) |
| `bh_zcut.png` | `bh-zcut` (Part IV) |
| `bh_outliers.png` | `bh-fdr-outliers` (Part IV, imperfect estimator; MC cached in `bh_outliers_scan.npz`) |

## Modules (no outputs of their own)

| file | what it owns |
|---|---|
| `bump_observables.py` | observables, resolutions, analyzable floors, published scan windows + sources, same-axis merges, the lepton-flavour split, the LEE math |
| `public_obs_map.py` | public model → spectrum map, published event-selection counts (`NSEL`) |
| `scan_alphabet.py` | the ten object types of the hypothetical scan, their derived `sigma`, the datasets it is priced on, and the four selection lenses; shared by `scaled-scan` and `budget-uncertainty` so both price the same scan |
| `yield_model.py` | the fittability requirement (≥100 events, ≥25 elements of ≥1 event), the declared background-yield model behind it and the per-lens efficiencies. Run it for the calibration table |
| `plot_style.py` | the shared palette and axis style of the statistics figures |
| `obs_labels.py` | the observable keys rendered as physics, shared by the figures and the LaTeX tables (standard library only) |
| `paper-numbers` | no table of its own: it collects every number the paper quotes from the tables above and writes them into `README.md` (standard library only) |
