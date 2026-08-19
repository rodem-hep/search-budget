PY  ?= python3
CLI := $(PY) -m searchbudget
J   ?= 1

.PHONY: help all budget census scan ab stats paper list graph check test clean distclean

help:
	@echo "every target below wraps \`$(CLI)\`:"
	@echo "  all       every stage that needs no network, in dependency order"
	@echo "  budget    the model-space budget, its figures and its uncertainty band"
	@echo "  census    the ATLAS publication record: counted, priced and written out"
	@echo "  scan      the combinatorial scan, its tiers and its selection lenses"
	@echo "  ab        two-stage A/B unblinding, and the measured estimator defect rate"
	@echo "  stats     selection rules: argmax vs threshold vs Benjamini-Hochberg"
	@echo "  paper     the generated number block in README.md"
	@echo "  list      every stage, its group, and which are out of date"
	@echo "  graph     what each stage reads and writes"
	@echo "  check     the registry, and that every declared output exists"
	@echo "  test      the test suite (needs pytest)"
	@echo "  clean     drop __pycache__;  distclean also deletes results/"
	@echo ""
	@echo "  PY = $(PY)          run stages concurrently with:  make all J=4"
	@echo "  the CLI takes more than these:  $(CLI) --help"

all:
	$(CLI) run --all -j $(J)

budget census scan ab stats paper:
	$(CLI) run $@ -j $(J)

list graph check:
	$(CLI) $@

test:
	$(PY) -m pytest -q

clean:
	$(CLI) clean

distclean:
	$(CLI) clean --outputs
