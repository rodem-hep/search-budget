# Reported excesses, mined from the census abstracts

The observed side of `EXCESS_COUNTING.md`, made concrete for the resonance subset: of the 290 census papers, **6 report a local excess >= 3 sigma in their abstract**, and 5 more quantify a largest deviation below 3 sigma; the rest report no significant excess or none at all. Extraction: `searchbudget/stages/reported_excesses.py`; the full rows, with the abstract sentence each number came from, are in `results/tables/reported_excesses.csv`.

| arXiv | year | spectrum | Z_local | Z_global |
|---|--:|---|--:|--:|
| 1606.03833 | 2016 | Diphoton (high mass) | 3.9 | 2.1 |
| 2404.12915 | 2024 | X -> S H (scalar + Higgs) | 3.5 | 2.0 |
| 1707.06958 | 2017 | VH (llbb, lvbb, vvbb) | 3.3 | 2.1 |
| 2110.00313 | 2021 | H -> aa (exotic Higgs decays) | 3.3 | 1.7 |
| 2209.10910 | 2022 | HH -> bb tautau | 3.1 | 2.0 |
| 1503.03290 | 2015 | SFOS dilepton + MET (edge) | 3.0 | - |

Every quoted global significance above is <= 2.1 sigma, and none has since grown into a discovery.

## Against the background-only expectation

A single background-only sweep of the 56-spectrum budget expects `N x p(>=3 sigma)` = 4118 x 1.35e-3 = **5.6** reports (**9.5** at selection granularity, N = 7030). The observed 6 sits between the two.

Two corrections pull in opposite directions and are left uncorrected:

- **Abstracts under-report.** A sub-3 sigma maximum is routinely left to the body
  (the 2 TeV diboson excess of 1506.00962 appears in its abstract only as a 2.5
  sigma *global*), so the count is a lower bound on what the papers contain.
- **The reports are not independent trials.** Successive papers on one axis re-scan
  partially shared data, and the record sweeps most axes more than once across
  Run 1 and Run 2, so the effective number of independent sweeps is above one.

Neither changes the conclusion: the reported population is what a counted trials
factor predicts under background-only, at the subset level and not just for the
program-wide anchor of `EXCESS_COUNTING.md`.
