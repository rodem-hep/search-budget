# Build the note.  Figures come from ../../results/plots/, so run `make all` in the repository
# root first if they are missing or stale.
#
#   make        -> main.pdf
#   make clean  -> drop the LaTeX intermediates

LATEX ?= pdflatex
TEX   := main.tex $(wildcard sections/*.tex) refs.bib

main.pdf: $(TEX)
	$(LATEX) -interaction=nonstopmode -halt-on-error main.tex
	-bibtex main
	$(LATEX) -interaction=nonstopmode -halt-on-error main.tex
	$(LATEX) -interaction=nonstopmode -halt-on-error main.tex

.PHONY: clean
clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.toc
