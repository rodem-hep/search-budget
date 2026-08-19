# What a published estimator gets wrong, and what that costs

Every trials count in this repository is enumerable because a published window and a resolution can
be looked up. A network-driven scan offers neither, and it adds a failure the Gaussian arithmetic
does not cover: the significance estimator is itself imperfect. For one such estimator the failure
rate has been measured and published, so the size of the problem is not a guess.

Source: BumpNet, JHEP 02 (2025) 122, arXiv:2501.05603. Written by `searchbudget/stages/estimator_defects.py`; the two-stage numbers it compares
against come from `ab_split_scan.csv` (`searchbudget/stages/ab_split_budget.py`).

## The measured rate

BumpNet predicts per-bin significances across families of mass histograms. On background-only
inputs it flags above 5 sigma in **0.048%** of histograms built from smooth
analytical functions and **0.129%** of those built from simulated Standard Model
samples, and applied to its **39,746** application histograms with no signal injected
it returns **53 +- 7** spurious 5 sigma candidates, which is
0.133% per histogram.

Those are not fluctuations. At 15 to 60 elements per histogram the
set holds 6.0e+05 to 2.4e+06 looks, for which the Gaussian expectation at 5 sigma is
0.17 to 0.68 -- below one. The measured rate is

| | |
|---|--:|
| defect rate per look | 2.2e-05 to 8.9e-05 |
| Gaussian tail at 5 sigma | 2.87e-07 |
| **inflation** | **78 to 310** |

one to two orders of magnitude above the statistical one. It is a **lower bound** on the fraction
`eps` of looks whose significance is simply wrong, since only mis-estimates that reach 5 sigma are
counted; at the 3 sigma level where candidates are selected, `eps` may well reach 1e-3. Those are
the two values the selection-rule study runs its GLITCH rows at (`searchbudget/stages/bh_fdr_outliers.py`,
`selection_rules.csv`): `eps = 1e-4` is the rate measured here at the 5 sigma flag level and
`eps = 1e-3` the value it plausibly reaches where candidates are actually selected, so the pair spans
the measured regime rather than being chosen for convenience.

## What it does to a single-stage scan

`Z_local = sqrt(25 + 2 ln N)` calibrates the tail of `N` Gaussian looks. If the measured tail is
78 to 310 times heavier, the effective trials factor is a property of the
network that nobody can enumerate afterwards, and the correction is not conservative but wrong by an
unknown factor.

## What it does to the two-stage split

The split converts that inflation into a countable list: whatever produced the entries, the
confirmation stage carries the trials factor of the frozen list. Two numbers price it, both on the
lensed combinatorial scan of `scaled_scan.txt` (N = 362,815).

| | nominal | if the inflation persists at Z_cut |
|---|--:|--:|
| pre-registered windows at `Z_cut = 3` | 490 | 49000 |
| claim bar `sqrt(25 + 2 ln(3k))` | 6.29 | 6.99 |

A hundredfold inflation of the list costs 0.69 sigma on the bar, because the
list enters logarithmically -- and even that assumes the length is modelled rather than read off the
frozen document.

Against that, the split is the more sensitive procedure whenever the trials factor a single-stage
analysis would honestly have to defend exceeds **14x** the counted one
(`ab_split_scan.csv`). The measured inflation clears that break-even by a factor
6 to 22, 7 at the round 1e2, so for an estimator of this
class the split is not a concession to auditability but the stronger procedure outright.

That comparison charges the single stage the full measured inflation, which is the pessimistic
reading. Suppressing defects upstream narrows the margin, and the veto below is exactly such a
measure -- but it has to demonstrate a factor **7** before the verdict changes
sign.

## The defect the split cannot kill

Splitting suppresses what fails to repeat. A coherent defect -- a mismodelled background, a detector
artefact, a network that sculpts the same shape in both halves -- passes confirmation as readily as
a real signal. The published application measures how large that class is: after the veto on
cross-histogram correlations, **27 of the 53** candidates survive, a
coherent fraction of **51% +- 7%**, clustered at the start of their histograms and
so shape-driven rather than statistical.

That fraction is the one number to measure before unblinding: evaluate the estimator on two
independent background-only samples under identical histogram definitions and count how often a flag
lands in the same window twice. It is also the BIAS row of the selection-rule study, which no
selection rule survives -- it has to be removed upstream.
