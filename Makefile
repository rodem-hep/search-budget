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

# the two data-free modules every budget number derives from
MOD := $(S)/bump_observables.py $(S)/public_obs_map.py

BUDGET := $(T)/search_budget.csv $(T)/search_budget_selections.csv $(O)/SEARCH_BUDGET.md \
          $(P)/search_budget.png $(P)/budget_waterfall.png $(P)/excess_counting.png \
          $(T)/combinatorial_budget.csv $(T)/composition_gap.txt

AB := $(P)/ab_split_reach.png $(P)/ab_toys_power.png $(P)/ab_split_outliers.png

STATS := $(G)/max_of_gaussians_light.png $(G)/bh_scan.png $(G)/bh_zcut.png $(G)/bh_outliers.png

.PHONY: help all budget ab stats note clean
.DELETE_ON_ERROR:

help:
	@echo "targets:"
	@echo "  all       budget + ab + stats  (everything; no external input of any kind)"
	@echo "  budget    trials factor N, Z_local, the combinatorial scan and the excess bookkeeping"
	@echo "  ab        two-stage A/B unblinding: reach, toy MC, imperfect-estimator robustness"
	@echo "  stats     selection rules: argmax vs threshold vs Benjamini-Hochberg (MC is cached)"
	@echo "  note      build docs/note/main.pdf (needs pdflatex; figures from results/plots/)"
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
$(T)/combinatorial_budget.csv: $(S)/combinatorial_budget.py $(S)/bump_observables.py
	$(PY) $(S)/combinatorial_budget.py
$(T)/composition_gap.txt: $(T)/combinatorial_budget.csv $(S)/composition_gap.py
	$(PY) $(S)/composition_gap.py > $@

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

note:
	$(MAKE) -C docs/note

clean:
	rm -rf $(S)/__pycache__
	$(MAKE) -C docs/note clean
