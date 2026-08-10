# The search budget of the BSM resonance program

How many **independent bump hunts** does the BSM resonance model space imply, what local
significance does a 5σ *global* discovery in that program therefore cost, and how should such a
scan be run? Everything here derives from **published searches and public model databases** — no
experiment-internal input of any kind, so a clone reproduces every number.

| | |
|---|---|
| bump spectra in the public model space | **46** (lepton flavour and b-jet content are part of the mass axis) |
| trials factor | `N ≈ 3.7k` inclusive, `6.6k` at published event-selection granularity |
| the discovery bar | `Z_local(5σ global) ≈ 6.44 → 6.53`, stable to ±0.1σ under every counting choice |
| fully combinatorial scan | 150 object categories, 1094 mass groups → `N ≈ 143k`, `Z_local ≈ 6.98` |
| never scanned | **57 of 82** reachable ≤4-object mass compositions have no published ATLAS bump hunt |
| two-stage A/B unblinding | costs ~0.5σ in reach, buys an exactly countable trials factor |

**Start here:** [`docs/note/main.pdf`](docs/note/main.pdf) — the note, which brings the studies
together for a reader who wants the physics rather than the code.

## Quick start

```bash
pip install -r requirements.txt        # numpy / scipy / matplotlib
make all                               # every table, report and figure
make help                              # the individual stages
```

The budget itself (`search_budget.py`, `combinatorial_budget.py`, `composition_gap.py`) is pure
standard library and runs under any `python3`; only the figures and the Monte Carlo need the
dependencies. Targets are file-level, so a rerun redoes only what is stale.

## Layout

```
├── Makefile                 one entry point for every step  (make help)
├── scripts/
│   ├── bump_observables.py  ** the observables: resolutions, floors, published scan windows **
│   ├── public_obs_map.py    ** public model → spectrum map, published event-selection counts **
│   └── ...                  one concern per file
├── results/
│   ├── tables/              machine-readable tables + the cached Monte Carlo (.npz)
│   ├── plots/               figures; plots/max_of_gaussians/ is the selection-rule study
│   └── overviews/           the written reports
└── docs/
    ├── METHOD_NOTES.md      the counting conventions and why each one is what it is
    ├── OUTPUTS.md           every output file → the script that writes it
    └── note/                the LaTeX note (main.pdf is committed)
```

## Key results

| file | what |
|---|---|
| [`results/overviews/SEARCH_BUDGET.md`](results/overviews/SEARCH_BUDGET.md) | **the number**: 46 spectra, `N ≈ 3.7k`, `Z_local ≈ 6.44` |
| [`results/plots/budget_waterfall.png`](results/plots/budget_waterfall.png) | one-figure summary: `Z_local` vs `N` across granularity levels |
| [`results/overviews/EXCESS_COUNTING.md`](results/overviews/EXCESS_COUNTING.md) | expected vs observed 3σ/5σ excesses across ATLAS — the external check on `N` |
| [`results/overviews/MAX_OF_GAUSSIANS.md`](results/overviews/MAX_OF_GAUSSIANS.md) | selection rules: argmax vs threshold vs Benjamini-Hochberg, with and without a perfect estimator |
| [`results/overviews/TWO_STAGE_UNBLINDING.md`](results/overviews/TWO_STAGE_UNBLINDING.md) | A/B split unblinding: logic, reach vs split fraction, caveats |
| [`results/tables/composition_gap.txt`](results/tables/composition_gap.txt) | which reachable mass compositions no published ATLAS search has ever scanned |

A complete index of outputs is in [`docs/OUTPUTS.md`](docs/OUTPUTS.md).

## The two modules everything rests on

`scripts/bump_observables.py` and `scripts/public_obs_map.py` are the **single source of truth**.
Import them; never redefine a floor, a window or a resolution in a consumer. Between them they own
every choice that decides what counts as one search:

- the published scan window per spectrum, each with a source note (the budget is counted over
  **published** windows, never over the range some set of signal samples happens to cover — that
  would understate sparsely sampled channels and make the exercise circular);
- the same-axis merges (`m(HH) 4b` is an event selection of `m(HH)`, not a second spectrum);
- the lepton-flavour split (`m(ll)` → `m(ee)`, `m(mumu)`), which is a mass axis and therefore
  deliberately **not** also a multiplier in `NSEL`. The invariant that catches a double count: the
  selections-level channel count stays at 94.

Both carry import-time assertions keeping the flavour layer complete: every leaf channel has a
window and a floor, and no flavour-inclusive parent keeps one.

## Scope

Invariant- and transverse-mass bump hunts only. The dedicated non-bump programs (displaced
HNL/Zd/RPV/ALP, MET+jet, dE/dx monopole, off-shell EFT) carry their own look-elsewhere effect and
are excluded rather than summed.

One question is deliberately out of scope: which of these 46 spectra a given collaboration has
actually produced signal samples for. That is a statement about a production program, it cannot be
made from public information, and it changes no number here — the windows are published-search
windows either way. It is treated in a companion ATLAS-internal note and repository, which vendors
the two modules above unchanged so the two cannot drift apart on what a spectrum is.
