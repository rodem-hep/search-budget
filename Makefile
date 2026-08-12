# The search budget of the BSM resonance program -- reproducible pipeline.
#
#   make help       list the targets
#   make all        every table, report and figure in the repository
#
# PY must point at an interpreter with requirements.txt installed (numpy/scipy/matplotlib).
# The budget itself (tables + report) is pure standard library and runs under any python3.
#     make all PY=/path/to/venv/bin/python

PY ?= python3
S  := scripts
T  := results/tables
O  := results/overviews
P  := results/plots
G  := results/plots/max_of_gaussians
X  := results/tex

# the two data-free modules every budget number derives from
MOD := $(S)/bump_observables.py $(S)/public_obs_map.py

BUDGET := $(T)/search_budget.csv $(T)/search_budget_selections.csv $(O)/SEARCH_BUDGET.md \
          $(P)/search_budget.png $(P)/budget_waterfall.png $(P)/excess_counting.png \
          $(T)/combinatorial_budget.csv $(T)/composition_gap.txt \
          $(X)/composition_appendix.tex \
          $(T)/scaled_scan.csv $(T)/scaled_scan_groups.csv $(T)/scaled_scan.txt \
          $(T)/published_census.csv $(O)/PUBLISHED_CENSUS.md \
          $(X)/census_refs.bib $(X)/census_appendix.tex $(O)/CENSUS_REFERENCES.md \
          $(X)/model_map_appendix.tex $(T)/model_spectrum_map.csv \
          $(X)/two_body_matrix.tex $(T)/two_body_matrix.csv \
          $(T)/reported_excesses.csv $(O)/REPORTED_EXCESSES.md

AB := $(P)/ab_split_reach.png $(P)/ab_toys_power.png $(P)/ab_split_outliers.png

STATS := $(G)/max_of_gaussians_light.png $(G)/bh_scan.png $(G)/bh_zcut.png $(G)/bh_outliers.png

.PHONY: help all budget ab stats clean
.DELETE_ON_ERROR:

help:
	@echo "targets:"
	@echo "  all       budget + ab + stats  (everything; no external input of any kind)"
	@echo "  budget    trials factor N, Z_local, the combinatorial scan and the excess bookkeeping"
	@echo "  ab        two-stage A/B unblinding: reach, toy MC, imperfect-estimator robustness"
	@echo "  stats     selection rules: argmax vs threshold vs Benjamini-Hochberg (MC is cached)"
	@echo "  clean     drop __pycache__"
	@echo "PY = $(PY)"

all: budget ab stats

# ---------------------------------------------------------------- search budget
budget: $(BUDGET)

$(T)/search_budget.csv $(T)/search_budget_selections.csv $(O)/SEARCH_BUDGET.md &: $(S)/search_budget.py $(MOD)
	$(PY) $(S)/search_budget.py
$(P)/search_budget.png $(P)/scan_windows.png $(P)/model_observable_matrix.png &: $(S)/budget_plots.py $(MOD)
	$(PY) $(S)/budget_plots.py
$(P)/budget_waterfall.png: $(T)/search_budget.csv $(T)/search_budget_selections.csv $(S)/budget_waterfall.py
	$(PY) $(S)/budget_waterfall.py
$(P)/excess_counting.png: $(T)/search_budget.csv $(S)/excess_counting.py
	$(PY) $(S)/excess_counting.py
$(T)/published_census.csv $(O)/PUBLISHED_CENSUS.md &: $(S)/published_census.py data/published_spectra.csv
	$(PY) $(S)/published_census.py
# The bibliography of the census. fetch_census_meta.py refreshes data/census_papers.csv from the
# arXiv API and is deliberately NOT a prerequisite: its output is committed so this stays offline.
$(X)/census_refs.bib $(X)/census_appendix.tex $(O)/CENSUS_REFERENCES.md &: $(S)/export_census_bib.py \
                                                 data/published_spectra.csv data/census_papers.csv
	$(PY) $(S)/export_census_bib.py
# The model-to-spectrum map as an appendix table: the inverse of model_observable_matrix.png.
$(X)/model_map_appendix.tex $(T)/model_spectrum_map.csv &: $(S)/export_model_map_tex.py $(MOD) \
                                                 $(S)/obs_labels.py
	$(PY) $(S)/export_model_map_tex.py
# The two-body object grid, priced in trials: arXiv:1907.06659 Table 14 with prices for marks.
$(X)/two_body_matrix.tex $(T)/two_body_matrix.csv &: $(S)/two_body_matrix.py $(MOD)
	$(PY) $(S)/two_body_matrix.py
# The observed excesses, mined from the census abstracts. fetch_census_abstracts.py refreshes
# data/census_abstracts.csv from the arXiv API and, like the metadata fetch, is NOT a prerequisite.
$(T)/reported_excesses.csv $(O)/REPORTED_EXCESSES.md &: $(S)/reported_excesses.py $(MOD) \
                        data/census_abstracts.csv data/census_papers.csv \
                        data/published_spectra.csv $(T)/search_budget.csv
	$(PY) $(S)/reported_excesses.py

$(T)/combinatorial_budget.csv: $(S)/combinatorial_budget.py $(S)/bump_observables.py
	$(PY) $(S)/combinatorial_budget.py
$(T)/composition_gap.txt $(X)/composition_appendix.tex &: $(T)/combinatorial_budget.csv \
                                                 $(S)/composition_gap.py $(S)/obs_labels.py
	$(PY) $(S)/composition_gap.py > $(T)/composition_gap.txt
# The same scan over the wider object alphabet: hadronic taus, photons, boosted W/Z, H and top.
$(T)/scaled_scan.csv $(T)/scaled_scan_groups.csv $(T)/scaled_scan.txt &: $(S)/scaled_scan.py \
                                                 $(S)/combinatorial_budget.py $(MOD)
	$(PY) $(S)/scaled_scan.py > $(T)/scaled_scan.txt

# ---------------------------------------------------------------- two-stage A/B unblinding
ab: $(AB)

$(P)/ab_split_reach.png $(P)/ab_split_crossover.png &: $(S)/ab_split_budget.py $(MOD)
	$(PY) $(S)/ab_split_budget.py
$(P)/ab_toys_power.png: $(S)/ab_split_toys.py $(MOD)
	$(PY) $(S)/ab_split_toys.py
$(P)/ab_split_outliers.png: $(S)/ab_split_outliers.py $(MOD)
	$(PY) $(S)/ab_split_outliers.py

# ---------------------------------------------------------------- selection-rule statistics
# self-contained Monte Carlo; the slow scans are cached in results/tables/*.npz
stats: $(STATS)

$(G)/max_of_gaussians_light.png: $(S)/max_of_gaussians_plots.py $(S)/plot_style.py
	$(PY) $(S)/max_of_gaussians_plots.py
$(G)/bh_scan.png $(T)/bh_fdr_scan.csv &: $(S)/bh_fdr_ab.py $(S)/plot_style.py
	$(PY) $(S)/bh_fdr_ab.py
$(G)/bh_zcut.png $(T)/bh_zcut_per_pe.csv &: $(S)/bh_zcut.py $(S)/plot_style.py
	$(PY) $(S)/bh_zcut.py
$(G)/bh_outliers.png: $(S)/bh_fdr_outliers.py $(S)/plot_style.py
	$(PY) $(S)/bh_fdr_outliers.py

clean:
	rm -rf $(S)/__pycache__
