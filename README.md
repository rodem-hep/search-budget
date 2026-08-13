# The search budget of the BSM resonance program

A bump hunt scans many mass points, and every additional spectrum anyone scans is another place
where a fluctuation can turn up. This repository works out how large that penalty really is for the
BSM resonance program taken as a whole: how many independent bump hunts the published model space
implies, what local significance a discovery therefore needs before it deserves to be called a
global 5σ, and how a scan of that size is best run.

Every number comes from published searches and public model databases, so a clone reproduces all of
it: the tables, the reports and the figures.

| | | computed by | the list behind it |
|---|---|---|---|
| bump spectra in the public model space | **46** | [`search_budget.py`](scripts/search_budget.py) | [`search_budget.csv`](results/tables/search_budget.csv), one row per spectrum with its window and source |
| trials factor | `N ≈ 3.7k` inclusive, `6.6k` at published event-selection granularity | [`search_budget.py`](scripts/search_budget.py) | [`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md), [`search_budget_selections.csv`](results/tables/search_budget_selections.csv) |
| the discovery bar | `Z_local(5σ global) ≈ 6.44 → 6.53`, stable to ±0.1σ under every counting choice | [`search_budget.py`](scripts/search_budget.py) | [`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md) |
| fully combinatorial scan | 1502 mass spectra across 204 object categories, giving `N ≈ 143k` and `Z_local ≈ 6.98` | [`combinatorial_budget.py`](scripts/combinatorial_budget.py) | [`combinatorial_budget.csv`](results/tables/combinatorial_budget.csv) |
| scaled up to ten object types | adding hadronic taus, photons, boosted W/Z, H and top, and MET as a free fifth ingredient of the masses: 36906 spectra, `N ≈ 2.0M`, `Z_local ≈ 7.35`, past the point (`N ≈ 2.7e5`) where a local 5σ has any global significance left | [`scaled_scan.py`](scripts/scaled_scan.py) | [`scaled_scan.csv`](results/tables/scaled_scan.csv), [`scaled_scan.txt`](results/tables/scaled_scan.txt) |
| held to a trials budget | `N ≤ 5e5` forces a priority order (model-motivated axes first, then highest expected rate): 5842 of 36906 spectra survive, all of them on the 41 model-motivated compositions, `Z_local = 7.16` | [`scaled_scan.py`](scripts/scaled_scan.py) | [`priority_scan.csv`](results/tables/priority_scan.csv), per composition kept against dropped |
| selection lenses | HT/Meff, displaced, VBF and ISR views of the same axes, one lens at a time: 2.9 histograms per spectrum, `N ≈ 5.6M`. Under the same budget they cost coverage rather than significance: 2221 spectra with 4841 lens views instead of 5842 inclusive ones | [`scaled_scan.py`](scripts/scaled_scan.py) | [`lens_scan.csv`](results/tables/lens_scan.csv), per lens what it adds and what survives |
| never scanned | **57 of 82** reachable ≤4-object mass compositions have no published ATLAS bump hunt | [`composition_gap.py`](scripts/composition_gap.py) | [`composition_gap.txt`](results/tables/composition_gap.txt), which names all 82 |
| the published record | **86** catalogued spectra over **290** ATLAS papers; 19 not revisited since before 2019, of which 12 sit on an axis already counted; 6 with a Run-3 result | [`published_census.py`](scripts/published_census.py) | [`published_spectra.csv`](data/published_spectra.csv), one row per spectrum with its arXiv references, and [`census_papers.csv`](data/census_papers.csv) with the title, journal and DOI of every one |
| observed against expected excesses | the published 3σ and 5σ count is what `N` predicts | [`excess_counting.py`](scripts/excess_counting.py) | [`EXCESS_COUNTING.md`](results/overviews/EXCESS_COUNTING.md) |
| two-stage A/B unblinding | costs about 0.5σ in reach, buys a trials factor you can count exactly | [`ab_split_budget.py`](scripts/ab_split_budget.py) | [`TWO_STAGE_UNBLINDING.md`](results/overviews/TWO_STAGE_UNBLINDING.md) |
| selection rules under an imperfect estimator | a fixed threshold beats both argmax and Benjamini-Hochberg | [`bh_fdr_outliers.py`](scripts/bh_fdr_outliers.py) | [`MAX_OF_GAUSSIANS.md`](results/overviews/MAX_OF_GAUSSIANS.md), [`bh_fdr_scan.csv`](results/tables/bh_fdr_scan.csv) |

Every window carries the published search it was taken from, so any single row can be checked
against its source without rerunning anything.

## Quick start

```bash
pip install -r requirements.txt        # numpy / scipy / matplotlib
make all                               # every table, report and figure
make help                              # the individual stages
```

The budget itself (`search_budget.py`, `combinatorial_budget.py`, `composition_gap.py`) is pure
standard library and runs under any `python3`; only the figures and the Monte Carlo need the
dependencies. Targets are file-level, so a rerun redoes just what has gone stale.

One script is outside `make all` and needs the network: `fetch_census_meta.py` refreshes the
bibliographic details of the census papers from the arXiv API. Its output is committed, so
everything else, including the generated bibliography, builds offline. Run it only after adding an
arXiv ID to `published_spectra.csv`, and it will fetch just the new one.

## How the counting works

One bump hunt means one invariant-mass spectrum. Every model that peaks in the same spectrum is
tested by the same search, so the number of searches is the number of distinct mass axes rather than
the number of models.

Within one spectrum, a resonance is resolvable if it sits at least a width away from its neighbour.
Taking the fractional resolution `σ_M/M = r` as constant gives

```
n_s = (1/r) · ln(M_hi / M_lo)          independent looks in a window [M_lo, M_hi]
N   = Σ n_s                            summed over spectra
Z_local = sqrt(25 + 2 ln N)            the local bar for a 5σ global discovery
```

Two conventions matter more than anything else in the arithmetic.

The window is always the one a published search family actually scanned, recorded per spectrum with
a source note. Counting over the mass range covered by whatever signal samples happen to exist
instead would understate thinly sampled channels and make the whole exercise circular.

The resolution `r` is the one genuine physics input, and it is only known to a factor of a few. So
no result is quoted as a single number: each is given with the band that follows from scaling `r` by
two in either direction, which moves `N` by a factor of two and `Z_local` by about 0.1σ. That the
answer barely moves is the point. Because `N` enters through `ln N`, breadth is cheap, and the
budget is remarkably insensitive to how finely you choose to slice the program.

## What is here

```
├── Makefile                 one entry point for every step  (make help)
├── scripts/
│   ├── bump_observables.py  ** the observables: resolutions, floors, published scan windows **
│   ├── public_obs_map.py    ** public model to spectrum map, published event-selection counts **
│   └── ...                  one concern per file
├── data/
│   ├── published_spectra.csv  the ATLAS resonance-search record, one row per spectrum
│   └── census_papers.csv     title, journal and DOI of every paper it cites
├── results/
│   ├── tables/              machine-readable tables and the cached Monte Carlo (.npz)
│   ├── plots/               figures; plots/max_of_gaussians/ is the selection-rule study
│   ├── tex/                 the census bibliography and its per-spectrum appendix
│   └── overviews/           the written reports, incl. the full census reference list
└── docs/
    ├── METHOD_NOTES.md      the counting conventions and why each one is what it is
    └── OUTPUTS.md           every output file and the script that writes it
```

The reports in `results/overviews/` are the long form of the rows above:
[`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md) for the budget itself, with the full
per-spectrum table; [`EXCESS_COUNTING.md`](results/overviews/EXCESS_COUNTING.md) for the check
against the excesses ATLAS has actually published;
[`TWO_STAGE_UNBLINDING.md`](results/overviews/TWO_STAGE_UNBLINDING.md) for the A/B scheme and its
caveats; [`MAX_OF_GAUSSIANS.md`](results/overviews/MAX_OF_GAUSSIANS.md) for the selection rules; and
[`PUBLISHED_CENSUS.md`](results/overviews/PUBLISHED_CENSUS.md) for what the publication record itself
contains, including the spectra that have gone stale, with all 290 of its papers written out in full
in [`CENSUS_REFERENCES.md`](results/overviews/CENSUS_REFERENCES.md).
[`docs/OUTPUTS.md`](docs/OUTPUTS.md) maps every output file to the script that writes it, and
[`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md) gives the reasoning behind each counting convention.

## Where the definitions live

`scripts/bump_observables.py` and `scripts/public_obs_map.py` hold every definition the numbers rest
on: the fractional resolution and analysable mass floor of each spectrum, its published scan window
together with the source that window came from, the map from public model classes to the spectra they
populate, and how many event selections each published search family scans on a given axis. They
also settle the cases where two labels describe the same mass axis, so `m(HH) 4b` counts as one more
event selection of `m(HH)` rather than as a second spectrum.

Everything else imports them. If you need a window or a resolution, take it from there rather than
writing it down again locally, because a second copy of a number will eventually disagree with the
first. Both modules check themselves when imported, so an incomplete channel fails loudly instead of
quietly dropping out of a table.

## Scope

Invariant-mass and transverse-mass bump hunts only. Several dedicated programs are deliberately left
out rather than added in: displaced-vertex searches for heavy neutral leptons, dark photons, RPV
decays and axion-like particles, missing energy plus jet, `dE/dx` monopole searches, and off-shell
EFT interpretations. Each carries a look-elsewhere effect of its own, computed over a space that is
not a mass spectrum, and summing them into a single trials count would mix search strategies that
are not comparable.
