PY ?= python3
S  := scripts
T  := results/tables
O  := results/overviews
P  := results/plots
G  := results/plots/max_of_gaussians
X  := results/tex

MOD := $(S)/bump_observables.py $(S)/public_obs_map.py
YM  := $(S)/yield_model.py

BUDGET := $(T)/search_budget.csv $(T)/search_budget_selections.csv $(O)/SEARCH_BUDGET.md \
          $(P)/search_budget.png $(P)/scan_windows.png $(P)/model_observable_matrix.png \
          $(P)/budget_waterfall.png $(P)/excess_counting.png \
          $(T)/combinatorial_budget.csv $(T)/composition_gap.txt \
          $(X)/composition_appendix.tex \
          $(T)/scaled_scan.csv $(T)/priority_scan.csv $(T)/lens_scan.csv \
          $(T)/scan_summary.csv $(T)/scaled_scan.txt \
          $(T)/published_census.csv $(O)/PUBLISHED_CENSUS.md \
          $(T)/census_budget.csv $(O)/CENSUS_BUDGET.md \
          $(X)/census_refs.bib $(X)/census_appendix.tex $(O)/CENSUS_REFERENCES.md \
          $(X)/model_map_appendix.tex $(T)/model_spectrum_map.csv \
          $(X)/two_body_matrix.tex $(T)/two_body_matrix.csv \
          $(T)/reported_excesses.csv $(O)/REPORTED_EXCESSES.md \
          $(T)/budget_uncertainty.csv $(O)/BUDGET_UNCERTAINTY.md $(X)/uncertainty_table.tex

AB := $(P)/ab_split_reach.png $(P)/ab_toys_power.png $(P)/ab_split_outliers.png \
      $(P)/ab_spurious_guard.png $(T)/ab_split_scan.csv $(T)/ab_split_toys.csv \
      $(T)/ab_guard_toys.csv \
      $(T)/estimator_defects.csv $(O)/ESTIMATOR_DEFECTS.md

STATS := $(G)/max_of_gaussians_light.png $(G)/bh_scan.png $(G)/bh_zcut.png \
         $(G)/bh_outliers.png $(T)/selection_rules.csv

.PHONY: help all budget ab stats clean
.DELETE_ON_ERROR:

help:
	@echo "targets:"
	@echo "  all       budget + ab + stats + the README number block (everything)"
	@echo "  budget    trials factor N, Z_local, the combinatorial scan and the excess bookkeeping"
	@echo "  ab        two-stage A/B unblinding: reach, toy MC, the measured estimator defect rate"
	@echo "  stats     selection rules: argmax vs threshold vs Benjamini-Hochberg (MC is cached)"
	@echo "  clean     drop __pycache__"
	@echo "PY = $(PY)"

all: budget ab stats README.md

README.md: $(S)/paper_numbers.py $(YM) $(S)/bump_observables.py \
           $(T)/search_budget.csv $(T)/search_budget_selections.csv $(T)/published_census.csv \
           $(T)/census_budget.csv $(T)/scan_summary.csv $(T)/lens_scan.csv \
           $(T)/priority_scan.csv $(T)/budget_uncertainty.csv $(T)/two_body_matrix.csv \
           $(T)/ab_split_scan.csv $(T)/ab_split_toys.csv $(T)/ab_guard_toys.csv \
           $(T)/estimator_defects.csv $(T)/selection_rules.csv
	$(PY) $(S)/paper_numbers.py

budget: $(BUDGET)

$(T)/search_budget.csv $(T)/search_budget_selections.csv $(O)/SEARCH_BUDGET.md &: $(S)/search_budget.py $(MOD)
	$(PY) $(S)/search_budget.py
$(P)/search_budget.png $(P)/scan_windows.png $(P)/model_observable_matrix.png &: $(S)/budget_plots.py $(MOD)
	$(PY) $(S)/budget_plots.py
$(P)/budget_waterfall.png: $(T)/search_budget.csv $(T)/search_budget_selections.csv $(S)/budget_waterfall.py
	$(PY) $(S)/budget_waterfall.py
$(P)/excess_counting.png: $(T)/search_budget.csv $(T)/census_budget.csv $(S)/excess_counting.py
	$(PY) $(S)/excess_counting.py
$(T)/published_census.csv $(O)/PUBLISHED_CENSUS.md &: $(S)/published_census.py data/published_spectra.csv
	$(PY) $(S)/published_census.py
$(T)/census_budget.csv $(O)/CENSUS_BUDGET.md &: $(S)/census_budget.py data/published_spectra.csv \
                                                 $(S)/bump_observables.py $(T)/search_budget.csv
	$(PY) $(S)/census_budget.py
$(X)/census_refs.bib $(X)/census_appendix.tex $(O)/CENSUS_REFERENCES.md &: $(S)/export_census_bib.py \
                                                 data/published_spectra.csv data/census_papers.csv
	$(PY) $(S)/export_census_bib.py
$(X)/model_map_appendix.tex $(T)/model_spectrum_map.csv &: $(S)/export_model_map_tex.py $(MOD) \
                                                 $(S)/obs_labels.py
	$(PY) $(S)/export_model_map_tex.py
$(X)/two_body_matrix.tex $(T)/two_body_matrix.csv &: $(S)/two_body_matrix.py $(MOD) $(YM)
	$(PY) $(S)/two_body_matrix.py
$(T)/reported_excesses.csv $(O)/REPORTED_EXCESSES.md &: $(S)/reported_excesses.py $(MOD) \
                        data/census_abstracts.csv data/census_papers.csv \
                        data/published_spectra.csv $(T)/search_budget.csv
	$(PY) $(S)/reported_excesses.py

$(T)/budget_uncertainty.csv $(O)/BUDGET_UNCERTAINTY.md $(X)/uncertainty_table.tex &: \
                                                 $(S)/budget_uncertainty.py $(MOD) $(YM) \
                                                 $(S)/scan_alphabet.py $(S)/combinatorial_budget.py \
                                                 $(T)/two_body_matrix.csv
	$(PY) $(S)/budget_uncertainty.py

$(T)/combinatorial_budget.csv: $(S)/combinatorial_budget.py $(S)/bump_observables.py $(YM)
	$(PY) $(S)/combinatorial_budget.py
$(T)/composition_gap.txt $(X)/composition_appendix.tex &: $(T)/combinatorial_budget.csv \
                                                 $(S)/composition_gap.py $(S)/obs_labels.py
	$(PY) $(S)/composition_gap.py > $(T)/composition_gap.txt
$(T)/scaled_scan.csv $(T)/priority_scan.csv $(T)/lens_scan.csv $(T)/scan_summary.csv \
$(T)/scaled_scan.txt &: \
                                                 $(S)/scaled_scan.py $(S)/scan_alphabet.py \
                                                 $(S)/combinatorial_budget.py $(MOD) $(YM)
	$(PY) $(S)/scaled_scan.py > $(T)/scaled_scan.txt

ab: $(AB)

$(P)/ab_split_reach.png $(P)/ab_split_crossover.png $(T)/ab_split_scan.csv &: \
                                                 $(S)/ab_split_budget.py $(MOD) \
                                                 $(T)/scaled_scan.csv
	$(PY) $(S)/ab_split_budget.py
$(T)/estimator_defects.csv $(O)/ESTIMATOR_DEFECTS.md &: $(S)/estimator_defects.py \
                                                 $(T)/ab_split_scan.csv
	$(PY) $(S)/estimator_defects.py
$(P)/ab_toys_power.png $(P)/ab_toys_background.png $(P)/ab_toys_spectrum.png \
$(T)/ab_split_toys.csv &:                        $(S)/ab_split_toys.py $(MOD)
	$(PY) $(S)/ab_split_toys.py
$(P)/ab_split_outliers.png $(P)/ab_outliers_mechanism.png $(P)/ab_outliers_spectrum.png &: \
                                                 $(S)/ab_split_outliers.py $(MOD)
	$(PY) $(S)/ab_split_outliers.py
$(P)/ab_spurious_guard.png $(T)/ab_guard_toys.csv &: $(S)/ab_spurious_guard.py \
                                                 $(S)/plot_style.py
	$(PY) $(S)/ab_spurious_guard.py

stats: $(STATS)

$(G)/max_of_gaussians_light.png $(G)/signal_wins_the_max.png $(G)/ab_confirmation.png \
$(G)/threshold_scan.png $(G)/threshold_vs_argmax.png $(G)/roc_threshold_vs_argmax.png &: \
                                                 $(S)/max_of_gaussians_plots.py $(S)/plot_style.py
	$(PY) $(S)/max_of_gaussians_plots.py
$(G)/bh_scan.png $(G)/bh_vs_argmax.png $(G)/roc_bh_vs_threshold.png $(T)/bh_fdr_scan.csv &: \
                                                 $(S)/bh_fdr_ab.py $(S)/plot_style.py
	$(PY) $(S)/bh_fdr_ab.py
$(G)/bh_zcut.png $(T)/bh_zcut_per_pe.csv &: $(S)/bh_zcut.py $(S)/plot_style.py
	$(PY) $(S)/bh_zcut.py
$(G)/bh_outliers.png $(T)/selection_rules.csv &: $(S)/bh_fdr_outliers.py $(S)/plot_style.py
	$(PY) $(S)/bh_fdr_outliers.py

clean:
	rm -rf $(S)/__pycache__
