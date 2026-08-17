# The search budget of the BSM resonance program

A bump hunt scans many mass points, and every additional spectrum anyone scans is another place
where a fluctuation can turn up. This repository works out how large that penalty is for the BSM
resonance program taken as a whole: how many independent bump hunts the published model space
implies, what local significance a discovery therefore needs before it deserves to be called a
global 5σ, and how a scan of that size is best run.

Every input is public, a published search window or a public model database, so a clone reproduces
all of it. `make all` writes every table, figure, LaTeX fragment and generated report under
`results/`, and needs nothing else: no ATLAS data, no simulation, no network. These are the numbers
behind *The Search Budget of the BSM Resonance Program*; `results/plots/` and `results/tex/` are
generated for that write-up to include directly, so no figure or table in it is drawn by hand.

## The result

The trials count on every basis, narrowest first. Rows are counted differently and are never summed.

| what | the number | computed by | the list behind it |
|---|---|---|---|
| spectra the public model space populates | **46** | [`search_budget.py`](scripts/search_budget.py) | [`search_budget.csv`](results/tables/search_budget.csv), one row per spectrum with its window and source |
| the trials factor and the discovery bar | `N ≈ 3.7k` → `Z_local = 6.44`, or `6.6k` → `6.53` at published event-selection granularity | [`search_budget.py`](scripts/search_budget.py) | [`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md) |
| the published ATLAS program, all searches | `N ~ 5e4` → `Z_local = 6.83`, taken from the literature rather than from this enumeration | [`search_budget.py`](scripts/search_budget.py) | [`EXCESS_COUNTING.md`](results/overviews/EXCESS_COUNTING.md), where the anchor and its sources are set out |
| a fully combinatorial scan, Run 2+3 | every mass built from ≤4 of ten object types: **4438** fittable spectra of 21644, `N = 2.0e5` → `Z_local = 7.03` | [`scaled_scan.py`](scripts/scaled_scan.py) | [`scaled_scan.txt`](results/tables/scaled_scan.txt), the closing per-dataset block |
| how much of that scan theory motivates | **2411 of 4438** spectra (54%) sit on a model-motivated mass axis, over 37 of the 258 fittable object compositions | [`scaled_scan.py`](scripts/scaled_scan.py) | [`priority_scan.csv`](results/tables/priority_scan.csv), per composition |
| the published record | **86** catalogued spectra over **290** ATLAS papers; 19 not revisited since before 2019, of which 12 sit on an axis already counted | [`published_census.py`](scripts/published_census.py) | [`published_spectra.csv`](data/published_spectra.csv), one row per spectrum with its arXiv references |
| observed against expected excesses | the published 3σ and 5σ count is what `N` predicts, with no adjustment | [`excess_counting.py`](scripts/excess_counting.py) | [`REPORTED_EXCESSES.md`](results/overviews/REPORTED_EXCESSES.md) |
| two-stage A/B unblinding | costs about 0.5σ in reach, buys a trials factor you can count exactly | [`ab_split_budget.py`](scripts/ab_split_budget.py) | [`TWO_STAGE_UNBLINDING.md`](results/overviews/TWO_STAGE_UNBLINDING.md) |
| selection rules under an imperfect estimator | a fixed threshold beats both argmax and Benjamini-Hochberg, by more as the estimator degrades | [`bh_fdr_outliers.py`](scripts/bh_fdr_outliers.py) | [`MAX_OF_GAUSSIANS.md`](results/overviews/MAX_OF_GAUSSIANS.md) |

A factor 50 in `N` separates the known model space from the fully combinatorial scan, and moves the
bar by 0.6σ, because `N` enters only through its logarithm. That is the point of the whole exercise:
breadth is cheap, and the answer barely depends on how finely the program is sliced.

## Reproducing it

```bash
pip install -r requirements.txt        # numpy / scipy / matplotlib
make all                               # every table, report and figure
make help                              # the individual stages
```

A full rebuild is byte-identical, figures included, so `git status` stays clean afterwards and
anything that does differ is a real change. The Monte Carlo is seeded and cached in
`results/tables/*.npz` (pass `--refit` to redo it), and the figures reproduce exactly under the
library versions recorded in `requirements.txt`. Targets are file-level, so a rerun redoes only what
has gone stale. The three hand-written reports in `results/overviews/` are the only files under
`results/` that no script produces, and [`docs/OUTPUTS.md`](docs/OUTPUTS.md) marks them.

The budget itself, `search_budget.py`, `combinatorial_budget.py`, `scaled_scan.py` and
`composition_gap.py`, is pure standard library and runs under any `python3`; only the figures and
the Monte Carlo need the dependencies.

Two scripts are outside `make all` because they need the network: `fetch_census_meta.py` and
`fetch_census_abstracts.py` refresh the bibliographic details and abstracts of the census papers
from the arXiv API. Their output is committed, so everything else, including the generated
bibliography, builds offline. Run them only after adding an arXiv ID to `published_spectra.csv`,
and they will fetch just the new one.

## How the counting works

One bump hunt means one invariant-mass spectrum. Every model that peaks in the same spectrum is
tested by the same search, so the number of searches is the number of distinct mass axes rather than
the number of models, and reinterpretation costs nothing.

Within one spectrum, a resonance is resolvable if it sits at least a width away from its neighbour.
Taking the fractional resolution `σ_M/M = r` as constant gives

```
n_s = (1/r) · ln(M_hi / M_lo)          independent looks in a window [M_lo, M_hi]
N   = Σ n_s                            summed over spectra
Z_local = sqrt(25 + 2 ln N)            the local bar for a 5σ global discovery
```

Three conventions matter more than anything else in the arithmetic.

**The window is the published one.** Every spectrum is counted over the range a published search
family actually scanned, recorded per spectrum with a source note. Counting over the mass range
covered by whatever signal samples happen to exist instead would understate thinly sampled channels
and make the whole exercise circular.

**The resolution is the one physics input, and carries a band.** `r` is known only to a factor of a
few, so no result is quoted as a single number: each is given with the band that follows from
scaling `r` by two in either direction, which moves `N` by a factor of two and `Z_local` by about
0.1σ.

**A hypothetical spectrum has to hold enough events to fit.** Wherever a spectrum has not been
published, `yield_model.py` requires at least 100 events and at least 25 resolution elements holding
one event or more, on a declared power-law background calibrated per object type. It truncates every
window at its one-event mass and drops the spectrum if too little is left, which removes four of
every five combinations of the ten-object scan. Published windows are exempt, because a published
search demonstrates its own feasibility. Run the module for its calibration table.

## Where the definitions live

`scripts/bump_observables.py`, `scripts/public_obs_map.py` and `scripts/yield_model.py` hold every
definition the numbers rest on: the fractional resolution and analysable mass floor of each
spectrum, its published scan window together with the source that window came from, the map from
public model classes to the spectra they populate, how many event selections each published search
family scans on a given axis, and the yield model behind the statistics requirement. They also
settle the cases where two labels describe the same mass axis, so `m(HH) 4b` counts as one more
event selection of `m(HH)` rather than as a second spectrum.

Everything else imports them. If you need a window or a resolution, take it from there rather than
writing it down again locally, because a second copy of a number will eventually disagree with the
first. The modules check themselves when imported, so an incomplete channel fails loudly instead of
quietly dropping out of a table.

## What is here

```
├── Makefile                 one entry point for every step  (make help)
├── scripts/                 one concern per file; the three modules above own the definitions
├── data/                    the only hand-curated input: the ATLAS publication record, with the
│                            bibliographic details and abstracts fetched for it
├── results/                 written by `make all`
│   ├── tables/              machine-readable tables, and the cached Monte Carlo (.npz)
│   ├── plots/               figures; plots/max_of_gaussians/ is the selection-rule study
│   ├── tex/                 LaTeX fragments for a consuming document to \input
│   └── overviews/           the written reports: four generated, three authored by hand
└── docs/
    ├── METHOD_NOTES.md      the counting conventions and why each one is what it is
    └── OUTPUTS.md           every output file, and the script that writes it
```

The four generated reports are [`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md), the budget
with its full per-spectrum table; [`PUBLISHED_CENSUS.md`](results/overviews/PUBLISHED_CENSUS.md),
what the publication record contains, including the spectra that have gone stale;
[`REPORTED_EXCESSES.md`](results/overviews/REPORTED_EXCESSES.md), every excess ≥3σ reported in the
census abstracts; and [`CENSUS_REFERENCES.md`](results/overviews/CENSUS_REFERENCES.md), all 290
papers written out in full.

## Scope

Invariant-mass and transverse-mass bump hunts only. Several dedicated programs are deliberately left
out rather than added in: displaced-vertex searches for heavy neutral leptons, dark photons, RPV
decays and axion-like particles, missing energy plus jet, `dE/dx` monopole searches, and off-shell
EFT interpretations. Each carries a look-elsewhere effect of its own, computed over a space that is
not a mass spectrum, and summing them into a single trials count would mix search strategies that
are not comparable.
