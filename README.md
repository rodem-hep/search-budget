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
| spectra the public model space populates | **46** | [`search_budget.py`](searchbudget/stages/search_budget.py) | [`search_budget.csv`](results/tables/search_budget.csv), one row per spectrum with its window and source |
| the trials factor and the discovery bar | `N ≈ 3.7k` → `Z_local = 6.44 +0.18/−0.16`, or `6.6k` → `6.53` at published event-selection granularity | [`search_budget.py`](searchbudget/stages/search_budget.py) | [`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md) |
| what that bar is uncertain by | every declared input varied over its own range: **`+0.18/−0.16` σ** on the model space, `+0.23/−0.31` on the combinatorial scan, but only `+0.12/−0.23` on the difference between them | [`budget_uncertainty.py`](searchbudget/stages/budget_uncertainty.py) | [`BUDGET_UNCERTAINTY.md`](results/overviews/BUDGET_UNCERTAINTY.md), source by source |
| a fully combinatorial scan, Run 2+3 | every mass built from ≤4 of ten object types: **4438** fittable spectra of 21644, `N = 2.0e5` → `Z_local = 7.03 +0.23/−0.31` | [`scaled_scan.py`](searchbudget/stages/scaled_scan.py) | [`scaled_scan.txt`](results/tables/scaled_scan.txt), the closing per-dataset block |
| … with one event-level selection at a time | four selection lenses, never combined: **8211** histograms, `N = 3.6e5` → `Z_local = 7.11 +0.24/−0.32` | [`scaled_scan.py`](searchbudget/stages/scaled_scan.py) | [`lens_scan.csv`](results/tables/lens_scan.csv), per lens and per dataset |
| the same program, as the literature counts it | `N ~ 5e4` quasi-independent tests over all ATLAS searches, signal regions included — a broader object than these mass scans, quoted only as an order-of-magnitude check | — | [`EXCESS_COUNTING.md`](results/overviews/EXCESS_COUNTING.md), where the anchor and its sources are set out |
| how much of that scan theory motivates | **2411 of 4438** spectra (54%) sit on a model-motivated mass axis, over 37 of the 258 fittable object compositions | [`scaled_scan.py`](searchbudget/stages/scaled_scan.py) | [`priority_scan.csv`](results/tables/priority_scan.csv), per composition |
| the published record | **100** charged (search, axis) pairs, from **86** catalogued entries over **290** ATLAS papers on 62 distinct bump observables; 19 entries not revisited since before 2019, of which 12 sit on an axis already counted | [`published_census.py`](searchbudget/stages/published_census.py) | [`published_spectra.csv`](data/published_spectra.csv), one row per entry with its arXiv references |
| what the published record costs | the same census priced in trials: **`N = 7.7k` → `Z_local = 6.55`** over those 100 pairs, or `N = 3.7k` → `6.44` counting each axis once over the union of its published ranges — the model side, reached from the publication record instead | [`census_budget.py`](searchbudget/stages/census_budget.py) | [`CENSUS_BUDGET.md`](results/overviews/CENSUS_BUDGET.md), per published search and per axis |
| observed against expected excesses | the published 3σ and 5σ count is what `N` predicts, with no adjustment | [`excess_counting.py`](searchbudget/stages/excess_counting.py) | [`REPORTED_EXCESSES.md`](results/overviews/REPORTED_EXCESSES.md) |
| two-stage A/B unblinding | costs about 0.5σ in reach, buys a trials factor you can count exactly | [`ab_split_budget.py`](searchbudget/stages/ab_split_budget.py) | [`ab_split_scan.csv`](results/tables/ab_split_scan.csv), the design priced on each basis; [`TWO_STAGE_UNBLINDING.md`](results/overviews/TWO_STAGE_UNBLINDING.md) |
| what an imperfect estimator costs | a published scanner mis-estimates 1e-5 to 1e-4 of its looks, **78–310×** the Gaussian tail, against a break-even of **14×** for the split; **51 ± 7 %** of its spurious flags are coherent | [`estimator_defects.py`](searchbudget/stages/estimator_defects.py) | [`ESTIMATOR_DEFECTS.md`](results/overviews/ESTIMATOR_DEFECTS.md) |
| selection rules under an imperfect estimator | a fixed threshold beats both argmax and Benjamini-Hochberg, by more as the estimator degrades | [`bh_fdr_outliers.py`](searchbudget/stages/bh_fdr_outliers.py) | [`MAX_OF_GAUSSIANS.md`](results/overviews/MAX_OF_GAUSSIANS.md) |

A factor 47 in `N` separates the published record from a combinatorial scan with its selection
lenses, and moves the bar by 0.56σ, because `N` enters only through its logarithm. That is the point
of the whole exercise: breadth is cheap, and the answer barely depends on how finely the program is
sliced. The expensive part is scanning with an estimator whose mistakes cannot be counted, which is
what the last two rows are about.

## Reproducing it

```bash
pip install -r requirements.txt        # numpy / scipy / matplotlib, for the figures only
make all                               # every table, report, figure and the block below
make help                              # the stage groups, and what each one covers
```

Everything is one package, `searchbudget`, driven by one command. `make` is a thin wrapper over that
command and `pip install -e .` puts it on the path as `search-budget`, but neither is required:
`python -m searchbudget` works from a clone with nothing installed.

```bash
python -m searchbudget list               # every stage, its group, and which are out of date
python -m searchbudget graph              # what each stage reads and what it writes
python -m searchbudget run budget         # a whole group ...
python -m searchbudget run census-budget  # ... or one stage, with whatever it depends on
python -m searchbudget run --all -j 4     # everything, four stages at a time
python -m searchbudget check              # the registry, and that every declared output exists
```

Every stage declares the files it reads and the files it writes, so the dependency graph lives in
the code instead of being restated in a build file. `run` puts the stages you asked for in
dependency order, pulls in whatever they need, and skips any whose outputs are already newer than
every source feeding them — the declared inputs, the upstream tables, and the package modules that
stage imports. `--force` rebuilds regardless, `-n` prints the plan and stops, `-j` runs independent
stages concurrently, and each stage runs in its own process, so nothing one stage leaves behind can
reach the next.

A full rebuild is byte-identical, figures included, so `git status` stays clean afterwards and
anything that does differ is a real change; `make test` asserts exactly that for the stages that
need no dependencies, alongside the registry and the definitions (`pip install -e ".[dev]"` for
pytest). The Monte Carlo is seeded and cached in `results/tables/*.npz` (pass
`--refit` to redo it), and the figures reproduce exactly under the library versions recorded in
`requirements.txt`. The three hand-written reports in `results/overviews/` are the only files under
`results/` that no stage produces, and [`docs/OUTPUTS.md`](docs/OUTPUTS.md) marks them; the number
block below this section is generated into this README as well, so it cannot fall behind the tables.

The budget itself — the `budget`, `census` and `scan` groups — is pure standard library and runs
under any `python3`; only the figures and the Monte Carlo of the `ab` and `stats` groups need the
dependencies.

Two stages sit outside `run --all` because they need the network: `fetch-census-meta` and
`fetch-census-abstracts` refresh the bibliographic details and abstracts of the census papers from
the arXiv API. Their output is committed, so everything else, including the generated bibliography,
builds offline. Run them only after adding an arXiv ID to `published_spectra.csv`, and they will
fetch just the new one.

## Every number in the paper

Every figure quoted in *The Search Budget of the BSM Resonance Program*, with the file that
produces it. This section is generated by `make all`, so it cannot drift from the pipeline.

<!-- paper-numbers:start -->

### The ladder (Table 1)

| quantity | value | reproduced by |
|---|---|---|
| spectra the public model space populates | **56** | `search_budget.csv` |
| the model space, published windows | `N = 4,118` → `Z_local = 6.45` +0.18/-0.16 | `search_budget.csv`, `budget_uncertainty.csv` |
| the same at published event-selection granularity | 104 channels, `N = 7,030` → `6.54` | `search_budget_selections.csv` |
| the published ATLAS program | 104 charged (search, axis) pairs, `N = 7,875` → `6.55` | `census_budget.csv` |
| a fully combinatorial scan, Run 2+3 | 4,438 fittable spectra, `N = 201,136` → `7.03` +0.23/-0.31 | `scan_summary.csv` |
| ... with one event-level selection at a time | 8,211 histograms, `N = 362,815` → `7.11` +0.24/-0.32 | `scan_summary.csv` |
| published record → lensed scan | a factor 46 in `N`, worth +0.56σ | derived |
| the whole ladder | a factor 88 in `N`, worth +0.66σ | derived |
| model space → combinatorial scan | +0.58σ, +0.12/-0.23 as a difference | `budget_uncertainty.csv` |
| a second experiment scanning the same axes | doubles `N`, worth +0.11σ | derived |

### Counting the looks (Section 2)

| quantity | value | reproduced by |
|---|---|---|
| a histogram is fittable if it holds | ≥ 100 events and ≥ 25 elements of ≥ 1 event | `yield_model.py` |
| the declared background | `n(m) = 10^6·W·(m/1000 GeV)^(1−7)·(r/0.05)` per element | `yield_model.py` |
| the dataset | 140 fb⁻¹ at 13 TeV plus ~2× that at 13.6 TeV, i.e. ×3 the anchor and ×1.20 in mass reach | `scan_summary.csv` |
| what one look *is*, the leading uncertainty | `N` from 0.5·N to N·Z/√(2π): +0.14/-0.11 on the model space, +0.14/-0.10 on the scan, +0.01/+0.00 on the difference | `budget_uncertainty.csv` |
| the resolution scale | every `r` ×2 either way: +0.11/-0.11 / +0.13/-0.19 / +0.03/-0.08 | `budget_uncertainty.csv` |
| ... drawn per channel instead of in common | +0.02/-0.02 | `budget_uncertainty.csv` |
| Rice up-crossings against the element count | a factor `Z/√(2π)` = 2.6 at this bar | derived |

### The model-space budget (Section 3)

| quantity | value | reproduced by |
|---|---|---|
| cheapest and dearest spectrum | 14 looks (`m(Ht)`) to 402 (`m(gammagamma)`) | `search_budget.csv` |
| what the five sharpest axes carry | 36% of `N` | `search_budget.csv` |
| dropping the largest contributor | `Z_local` 6.45 → 6.44 | `search_budget.csv` |
| event selections instead of inclusive spectra | 104 channels, `6.54` | `search_budget_selections.csv` |
| the spectra the most model classes point at | `m(VV)` (13), `m(tt)` (9), `m(mumu)` (8) | `model_spectrum_map.csv` |
| where the model classes come from | 43 from the FeynRules/UFO database, 13 from the literature sweep, referenced class by class | `model_classes.csv` |
| model-motivated spectra with no ATLAS scan | **14** in the scan's (category, mass-group) units, the pair-produced axis resolved into its legs; 12 of them fittable | `unscanned_spectra.csv` |

### The published ATLAS program (Section 4)

| quantity | value | reproduced by |
|---|---|---|
| the record | **290** papers (2010–2026) in **86** catalogued searches over 64 bump observables | `published_census.csv`, `census_budget.csv` |
| (search, axis) pairs | 114, of which **104** carry a chargeable range, over 76 of the 86 searches | `census_budget.csv` |
| pairs charged a single look | 4 (a fixed mass, not a scan) | `census_budget.csv` |
| pairs with nothing chargeable | 10 | `census_budget.csv` |
| priced entry by entry | `N = 7,875` → `Z_local = 6.55` | `census_budget.csv` |
| priced once per axis, over the union of its ranges | `N = 3,837` → `6.44`, within 6.8% of the model space | `census_budget.csv` |
| recency | 19 of the 86 have no paper since 2019, 6 carry a published Run-3 result | `published_census.csv` |
| model independence | **29** of the 104 charged pairs carry a model-independent result, from 19 of the 76 searches | `model_independence.csv` |
| excesses a background-only sweep expects | 11 local ≥3σ and 0.0023 spurious ≥5σ over the census | `REPORTED_EXCESSES.md`, `EXCESS_COUNTING.md` |
| what a local 5σ is worth globally | 2.7σ over the published record, 0.8σ over the combinatorial scan, nothing at all over the lensed scan | derived |

### The combinatorial scan (Section 5)

| quantity | value | reproduced by |
|---|---|---|
| the alphabet | ten object types, masses of ≤4 objects: 21,644 combinations in 2,412 categories | `scan_summary.csv` |
| what statistics leaves | **4,438** fittable, `N = 201,136` → `7.03` | `scan_summary.csv` |
| two-body share of the survivors | 73% | `scan_summary.csv` |
| tier (i), every motivated composition once | 55 spectra, `N = 3,428` → `6.42` | `scan_summary.csv` |
| tier (ii), those compositions elsewhere | 2,702 spectra, `N = 135,875` | `scan_summary.csv` |
| tier (iii), nothing motivates it | 1,681 spectra, `N = 61,833` | `scan_summary.csv` |
| the share theory motivates | 2,757 of 4,438 spectra (62%), over 47 of the 258 fittable compositions | `scan_summary.csv`, per composition in `priority_scan.csv` |
| selection lenses | 9,407 conceivable views, 5,634 ruled out by statistics, 3,773 survive carrying 161,678 looks | `lens_scan.csv` |
| the lensed scan | 8,211 histograms (0.85 extra per spectrum), `N = 362,815` → `7.11`, so a lens costs +0.08σ | `scan_summary.csv` |
| where the relation loses meaning | `25 − 2 ln N < 0` at `N = 2.7·10^5`; the lensed scan is past it | derived |
| the yield anchor ×0.01 and ×100 | 1,729 to 8,923 spectra, `Z_local` 6.88 to 7.14 | `scan_summary.csv` |
| the same requirement applied to published axes | would keep 43 of 52 axes and 3,131 of 4,006 looks, so exempting them is conservative | `scan_summary.csv` |
| the two-body grid | 11 unscanned pairs cost 505 looks, from 30 (`tauH`) to 71 (`jV`); closing all of them takes `N` 4,118 → 4,623 and `Z_local` 6.45 → 6.47 | `two_body_matrix.csv` |
| the costliest compositions | `gg` (13,981), `eg` (13,838), `ee` (12,894), `ej` (12,018) | `scan_summary.csv`, per composition in `priority_scan.csv` |
| the most heavily scanned pair | `ee`, 630 looks over 3 axes | `two_body_matrix.csv` |

### The four selection lenses, Run 2+3 (Section 5)

One at a time, never a product of two (`lens_scan.csv`):

| lens | applies when | efficiency | spectra it could reach | views kept | looks |
|---|---|--:|--:|--:|--:|
| high HT or Meff | activity outside the mass | 0.1 | 4,046 | 2,615 | 114,012 |
| displaced activity | any reconstructed mass | 0.001 | 4,438 | 985 | 39,184 |
| forward jet pair | two free slots for the tag jets | 0.02 | 92 | 68 | 3,381 |
| ISR jet | one free slot, and a low-mass end, window capped at 200 GeV | 0.2 | 831 | 105 | 5,101 |

### Two-stage unblinding (Section 6)

| quantity | value | reproduced by |
|---|---|---|
| the basis it is priced on | the lensed scan, `N = 362,815`, single stage exactly corrected `7.11` | `ab_split_scan.csv` |
| the optimised split | reach `7.58` at `w = 3` (+0.47σ) and `7.43` with pinned windows (+0.32σ) | `ab_split_scan.csv` |
| the naive 50/50 split | reach `8.90` (+1.79σ) | `ab_split_scan.csv` |
| the pre-registered list at `Z_cut = 3` | 490 windows, claim bar `6.29` | `ab_split_scan.csv` |
| background-only toys | in 2·10⁴ toys the best confirmation reaches `5.30` against that bar: 0 false claims | `ab_split_scan.csv` |
| break-even trials inflation | `R* = 14` at `w = 3`, `5` with pinned windows | `ab_split_scan.csv` |
| the measured defect rate of a published scanner | 0.048% of analytical-function and 0.129% of simulated histograms flag above 5σ; 53 ± 7 spurious candidates on 39,746 histograms | `estimator_defects.csv` |
| ... per look, against the Gaussian tail | 2.2e-05 to 8.9e-05 against 2.87e-07, an inflation of **78 to 310×** | `estimator_defects.csv` |
| ... against the break-even | 6 to 22× past it, so the split is the more sensitive procedure | `estimator_defects.csv` |
| the coherent fraction the split cannot kill | 27 of 53 candidates, **51 ± 7 %** | `estimator_defects.csv` |
| toy validation of the reach | single stage on the full dataset 6.52 against 6.54 analytic; optimised split 7.14 against 7.10 analytic; naive 50/50 split 7.72 against 7.72 analytic (on the model space, `N = 7,030`) | `ab_split_toys.csv` |
| the procedure on one toy spectrum | background only: `Z_A = 3.27` → `Z_B = 0.62`; injected Z_full = 7: `Z_A = 3.03` → `Z_B = 6.35`, against a claim bar of 5.14 | `ab_guard_toys.csv` |
| the split's own exposure to that inflation | the list grows 490 → 49,000 and the bar 6.29 → 6.99 | `estimator_defects.csv` |

### Selecting candidates (Section 7)

Confirmation probability for a 5σ signal, every rule held to the same 1.35e-03 false-confirmation budget, in per cent (`selection_rules.csv`):

| estimator | argmax | threshold | BH | thr − argmax | nominal `q` | P(argmax *is* the signal) |
|---|--:|--:|--:|--:|--:|--:|
| PERFECT estimator | 78.5 | 82.5 | 80.9 | +4.0 | 3.82e-01 | 80.3 |
| GLITCH  eps=1e-4 | 73.2 | 81.2 | 80.1 | +8.0 | 3.22e-01 | 74.9 |
| GLITCH  eps=1e-3 | 43.8 | 57.3 | 57.2 | +13.6 | 1.39e-02 | 44.8 |
| GLITCH  eps=1e-2 | 6.2 | 5.8 | 6.0 | -0.4 | 4.24e-07 | 6.3 |
| BIAS    eps=1e-4 | *untunable* | 0.2 | *unreachable* | -- | -- | 74.9 |
| BIAS    eps=1e-3 | *untunable* | 0.0 | *unreachable* | -- | -- | 44.8 |
| BIAS    eps=1e-2 | *untunable* | 0.0 | *unreachable* | -- | -- | 6.3 |

### The uncertainty budget (Appendix)

Every declared input, moved over the range it is known to (`budget_uncertainty.csv`):

| source | varied over | model space | scan | difference |
|---|---|--:|--:|--:|
| mass resolution, scale | every r x2 either way | +0.107/-0.108 | +0.133/-0.191 | +0.026/-0.083 |
| mass resolution, per channel | each r independently, x2 per sigma (16-84%) | +0.023/-0.021 | +0.000/+0.000 | +0.021/-0.023 |
| mass resolution, shape | muon axes with r(M) rising to 0.10-0.20 at 3000 GeV | +0.001/-0.000 | +0.000/+0.000 | +0.000/-0.001 |
| mass resolution, prescription | worst leg / quadrature sum instead of the calibrated mean | +0.000/+0.000 | +0.000/-0.102 | +0.000/-0.102 |
| scan windows | every edge x1.4 either way (published), x1.25 (generic) | +0.031/-0.039 | +0.023/-0.024 | +0.015/-0.008 |
| yield anchor | N_ref x0.01 to x100 | +0.000/+0.000 | +0.109/-0.148 | +0.109/-0.148 |
| background slope | P = 6 to 8 | +0.000/+0.000 | +0.013/-0.012 | +0.013/-0.012 |
| fittability requirement | 30-300 events, 15-50 elements | +0.000/+0.000 | +0.044/-0.116 | +0.044/-0.116 |
| the axis set | non-peaking axes and the dilepton overlap dropped; the 11 unscanned two-body pairs added | +0.018/-0.004 | +0.000/+0.000 | +0.004/-0.018 |
| the definition of one look | N x0.5 to x Z/sqrt(2 pi) | +0.145/-0.108 | +0.145/-0.099 | +0.009/+0.000 |
| the closed-form LEE relation | exact Gaussian-tail solution instead | +0.000/-0.037 | +0.000/-0.045 | +0.000/-0.009 |
| total | quadrature over the systematic sources above | +0.183/-0.162 | +0.231/-0.308 | +0.123/-0.231 |

### The yield model (Appendix)

Calibrated object by object on the Run-2 anchor (`yield_model.py`, run it for this table):

| object | `F` | symmetric channel | `r` | events/element at 1 TeV | one-event mass |
|---|--:|---|--:|--:|--:|
| `j` | 1 | `m(jj)` | 0.050 | 1,000,000.0 | 10,000 GeV |
| `b` | 0.1 | `m(bb)` | 0.100 | 20,000.0 | 5,210 GeV |
| `V` | 0.02 | `m(VV)` | 0.060 | 480.0 | 2,798 GeV |
| `t` | 0.02 | `m(tt)` | 0.080 | 640.0 | 2,936 GeV |
| `H` | 0.01 | `m(HH)` | 0.080 | 160.0 | 2,330 GeV |
| `g` | 0.004 | `m(gammagamma)` | 0.010 | 3.2 | 1,214 GeV |
| `T` | 0.0035 | `m(tautau)` | 0.120 | 29.4 | 1,757 GeV |
| `e` | 0.003 | `m(ee)` | 0.015 | 2.7 | 1,180 GeV |
| `m` | 0.003 | `m(mumu)` | 0.050 | 9.0 | 1,442 GeV |
| `X` | 0.003 | `mT(ev)` | 0.100 | 18.0 | 1,619 GeV |
| `Z` | 1e-05 | -- | -- | -- | -- |

Every row above is written by [`paper_numbers.py`](searchbudget/stages/paper_numbers.py) from the committed tables, so `make all` keeps it and the paper in step; nothing here is typed by hand.

<!-- paper-numbers:end -->

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
0.1σ. That is one line of `budget_uncertainty.py`, which does the same for every other declared
input — the windows, the yield model, the fittability thresholds — and for the one input with no
measurement behind it, the convention that makes a resolution element an independent look. That last
one is the largest term, and because it moves every basis together it cancels in every difference
this study quotes.

**A hypothetical spectrum has to hold enough events to fit.** Wherever a spectrum has not been
published, `yield_model.py` requires at least 100 events and at least 25 resolution elements holding
one event or more, on a declared power-law background calibrated per object type. It truncates every
window at its one-event mass and drops the spectrum if too little is left, which removes four of
every five combinations of the ten-object scan. Published windows are exempt, because a published
search demonstrates its own feasibility. Run the module for its calibration table.

## Where the definitions live

`searchbudget/core/` holds every definition the numbers rest on: `bump_observables.py` the
fractional resolution, analysable mass floor and published scan window of each spectrum together
with the source that window came from; `public_obs_map.py` the map from public model classes to the
spectra they populate and how many event selections each published search family scans on a given
axis; `yield_model.py` the yield model behind the statistics requirement; `scan_alphabet.py` and
`combinatorial_budget.py` the ten object types of the hypothetical scan and the enumeration over
them. They also settle the cases where two labels describe the same mass axis, so `m(HH) 4b` counts
as one more event selection of `m(HH)` rather than as a second spectrum. `lee.py` and
`catalogue.py` hold the arithmetic and the canonical spectrum ordering every stage shares.

Every stage imports them. If you need a window or a resolution, take it from there rather than
writing it down again locally, because a second copy of a number will eventually disagree with the
first. The modules check themselves when imported, so an incomplete channel fails loudly instead of
quietly dropping out of a table.

## What is here

```
├── Makefile                 a wrapper over the CLI, for the usual entry points  (make help)
├── pyproject.toml           `pip install -e .` -> the `search-budget` command
├── searchbudget/            the package: one stage per concern, one place per definition
│   ├── cli.py               list / run / graph / check / clean
│   ├── registry.py          what a stage is, and what each one declares it reads and writes
│   ├── runner.py            the dependency order, what is stale, and running it
│   ├── deps.py              the package modules each stage imports, for staleness
│   ├── paths.py io.py       every path and every read and write in the repository
│   ├── arxiv.py             the arXiv API client the two fetch stages share
│   ├── core/                the definitions the numbers rest on, and the shared arithmetic
│   ├── stats/               machinery shared by the toy studies: FDR, defects, toy spectra
│   ├── viz/                 the shared palette, axis style and observable labels
│   └── stages/              one module per stage; nothing else writes to results/
├── tests/                   the registry, the definitions, and a byte-identical rebuild
├── data/                    the only hand-curated input: the ATLAS publication record, with the
│                            bibliographic details and abstracts fetched for it
├── results/                 written by `make all`
│   ├── tables/              machine-readable tables, and the cached Monte Carlo (.npz)
│   ├── plots/               figures; plots/max_of_gaussians/ is the selection-rule study
│   ├── tex/                 LaTeX fragments for a consuming document to \input
│   └── overviews/           the written reports: seven generated, three authored by hand
└── docs/
    ├── METHOD_NOTES.md      the counting conventions and why each one is what it is
    └── OUTPUTS.md           every output file, and the stage that writes it
```

Adding a step means adding one module under `searchbudget/stages/` with a `@stage(...)` declaring
its group, its inputs and its outputs; the CLI, the dependency order, the staleness check and the
tests pick it up with nothing else to edit.

The seven generated reports are [`SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md), the budget
with its full per-spectrum table; [`BUDGET_UNCERTAINTY.md`](results/overviews/BUDGET_UNCERTAINTY.md),
what each declared input is worth on the bar;
[`PUBLISHED_CENSUS.md`](results/overviews/PUBLISHED_CENSUS.md),
what the publication record contains, including the spectra that have gone stale;
[`CENSUS_BUDGET.md`](results/overviews/CENSUS_BUDGET.md), that record priced in trials;
[`REPORTED_EXCESSES.md`](results/overviews/REPORTED_EXCESSES.md), every excess ≥3σ reported in the
census abstracts; [`CENSUS_REFERENCES.md`](results/overviews/CENSUS_REFERENCES.md), all 290
papers written out in full; and
[`ESTIMATOR_DEFECTS.md`](results/overviews/ESTIMATOR_DEFECTS.md), the measured defect rate of a
published network scanner and what it does to both procedures.

## Scope

Invariant-mass and transverse-mass bump hunts only. Several dedicated programs are deliberately left
out rather than added in: displaced-vertex searches for heavy neutral leptons, dark photons, RPV
decays and axion-like particles, missing energy plus jet, `dE/dx` monopole searches, and off-shell
EFT interpretations. Each carries a look-elsewhere effect of its own, computed over a space that is
not a mass spectrum, and summing them into a single trials count would mix search strategies that
are not comparable.
