# Expected vs observed excesses across the ATLAS search program

**Question:** from a trials-factor estimate of the whole ATLAS search program, how many 3-sigma and
5-sigma excesses *should* we have seen under background-only, and how many did we *actually* see?

**Answer:** ~tens of >=3-sigma local excesses (all faded) and essentially zero spurious >=5-sigma --
which is exactly what background-only + the look-elsewhere effect predicts. The observed rate is
consistent with pure statistics; nothing anomalous, and no fake discovery.

Script: `scripts/excess_counting.py` -> `results/plots/excess_counting.png`.

## 1. Trials factor N (independent looks)
The relevant number is the count of quasi-independent statistical tests (mass-resolution elements +
signal-region bins) across all ATLAS searches.

| level | independent looks N | basis |
|---|---:|---|
| our bump-hunt budget (public channels, published windows) | ~2.5e3 | this study, `SEARCH_BUDGET.md` |
| ATLAS two-body invariant-mass searches | ~2-3e4 | lit. quotes **>50,000** two-body-mass looks LHC-wide (ATLAS ~ half) |
| **full ATLAS BSM program (+SUSY/exotics SRs)** | **~5e4** | ~300-500 searches x ~O(100) effective looks each ("a typical search now has hundreds of signal regions") |

Central **N ~ 5e4**, defensible range **1e4 - 1e5**. N is the dominant uncertainty (~x10), but the
conclusion is robust across the whole span.

## 2. Expected under background-only
One-sided Gaussian tail `p(>=Z) = 0.5*erfc(Z/sqrt2)`: `p(3s)=1.35e-3`, `p(5s)=2.87e-7`.
Expected count = `N * p`:

| N_trials | expected >=3s | expected >=5s (spurious) |
|---:|---:|---:|
| 1e4 | 14 | 0.003 |
| **5e4 (central)** | **~67** | **~0.014** |
| 1e5 | 135 | 0.029 |
| resonance subset (2.5e3) | ~3.4 | 7e-4 |

**Prediction: tens (~15-135) of >=3-sigma local excesses, and ~0.01 (i.e. none) spurious >=5-sigma.**

## 3. Observed
**For the resonance subset the observed side is now measured rather than anecdotal:**
`REPORTED_EXCESSES.md` mines the abstracts of the 288 census papers and finds **6 reported local
excesses >= 3 sigma** (every quoted global <= 2.1s), against 5.0 expected per background-only
sweep of the 46-spectrum budget (8.9 at selection granularity).

**Program-wide, take the number from the curated catalogue** rather than from a list of remembered
excesses: the LHC BSM Working Group maintains the list of excesses currently open, grouped by
signature, admitting anything above **2.4 sigma local** whatever the search type:
<https://lhc-bsm-wg.docs.cern.ch/excesses/> (steering committee, updated ~6-monthly). Of order 100
entries, **none reaches 5 sigma** and **no quoted global exceeds ~2.8 sigma**.

Scope differences to respect before comparing it with our numbers:

- it spans **ATLAS + CMS + LHCb** and mixes bump hunts with shape fits, SR counting and
  unconventional signatures, so it is not the 46-axis basis;
- it lists **open** excesses only, while the census keeps every paper;
- it reads the **paper bodies**, so it carries ATLAS resonance excesses >=3 sigma that our abstract
  mining misses (trigger-level dijet at 650 GeV, `A -> ZH -> llbb`, charged Higgs -> cb at 130 GeV).
  Our 6 is therefore a floor.

**>=5-sigma: every one is a genuine Standard-Model process**, none BSM. Spurious BSM >=5-sigma:
exactly zero.

## 4. Verdict -- consistent with background
| | expected (N~5e4) | observed |
|---|---:|---|
| >=3s local | ~15-135 | ~tens, all faded  OK |
| >=5s spurious | ~0.01 | 0  OK |

The parade of 3-sigma "hints" that come and go is **not a mystery** -- it is the expected rate of
statistical noise once ~1e4-1e5 looks are multiplied by the 3-sigma tail. The complete absence of a
fake 5-sigma matches `N*p(5s) ~ 0.01` -- which is *why* 5 sigma is the discovery threshold. The same
logic (arXiv 2605.24441) notes that with >50,000 looks even a local 5 sigma is only ~1.5-2.5 sigma
global, and argues for raising the fully-agnostic-scan bar toward 7 sigma -- the direct extension of
the trials-factor argument in `SEARCH_BUDGET.md`.

## Caveats
- N spans a factor ~10; the *conclusion* survives the whole 1e4-1e5 range (>=3s stays "tens", >=5s
  stays "<<1").
- Counts mix ATLAS-only and ATLAS+CMS excesses; for ATLAS alone the numbers are smaller but same order.
- "Observed ~tens" is the count of *notable/highlighted* >=3s excesses; the raw number of individual
  >=3s bins is larger and unremarked (they are the expected noise), also consistent.
- One-sided (upward) tail used, appropriate for excess searches; deficits are a separate (also
  background-consistent) population.

## Sources
- On the Statistical Interpretation of Discoveries in LHC Data -- arXiv:2605.24441 (>50,000 two-body
  looks; 5s local -> ~1.5-2.5s global; recommend 7s)
- Digging Deeper for New Physics in the LHC Data -- arXiv:1707.05783 (hundreds of SRs/search; hidden ~3s)
- Gross & Vitells, Trial factors for the look-elsewhere effect -- arXiv:1005.1891
- 750 GeV diphoton excess -- en.wikipedia.org/wiki/750_GeV_diphoton_excess
- ATLAS Searches public results -- twiki.cern.ch/twiki/bin/view/AtlasPublic/SearchesPublicResults
