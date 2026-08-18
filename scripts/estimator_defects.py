#!/usr/bin/env python3
import csv, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
OVR = os.path.join(ROOT, "results", "overviews")

SOURCE = "BumpNet, JHEP 02 (2025) 122, arXiv:2501.05603"

PUBLISHED = {
    "rate_analytic": (4.8e-4, "fraction of background-only histograms built from smooth "
                              "analytical functions flagged above 5 sigma"),
    "rate_simulated": (1.29e-3, "the same fraction for histograms built from simulated Standard "
                                "Model samples"),
    "histograms": (39746, "application histograms scanned with no signal injected"),
    "flagged": (53, "spurious 5 sigma candidates returned on them"),
    "flagged_err": (7, "quoted uncertainty on that count"),
    "coherent": (27, "candidates surviving the veto on cross-histogram correlations"),
    "elements_lo": (15, "resolution elements per histogram, low end"),
    "elements_hi": (60, "resolution elements per histogram, high end"),
}

WIDEN = 3.0
INFLATION_ROUND = 100.0


def p1(z):
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def claim_bar(k, widen=WIDEN):
    return math.sqrt(25.0 + 2.0 * math.log(widen * max(k, 1.0)))


v = {k: val for k, (val, _) in PUBLISHED.items()}
P5, P3 = p1(5.0), p1(3.0)

looks_lo = v["histograms"] * v["elements_lo"]
looks_hi = v["histograms"] * v["elements_hi"]
exp_lo, exp_hi = looks_lo * P5, looks_hi * P5
rate_lo, rate_hi = v["flagged"] / looks_hi, v["flagged"] / looks_lo
infl_lo, infl_hi = rate_lo / P5, rate_hi / P5
per_hist = v["flagged"] / v["histograms"]

coh = v["coherent"] / v["flagged"]
coh_err = math.sqrt(coh * (1.0 - coh) / v["flagged"])

scan = {r["basis"]: r for r in csv.DictReader(open(os.path.join(TAB, "ab_split_scan.csv")))}
lensed = scan["lensed scan"]
N_scan = float(lensed["N_trials"])
k_nom = float(lensed["k_zcut3"])
k_infl = k_nom * INFLATION_ROUND
bar_nom, bar_infl = claim_bar(k_nom), claim_bar(k_infl)
R_star = float(lensed["R_star_w3"])
margin_lo, margin_hi = infl_lo / R_star, infl_hi / R_star
margin_round = INFLATION_ROUND / R_star

rows = [
    ("published rate, analytical-function histograms", f"{v['rate_analytic']:.3%}", SOURCE),
    ("published rate, simulated Standard Model histograms", f"{v['rate_simulated']:.3%}", SOURCE),
    ("application histograms, no signal injected", f"{v['histograms']:,}", SOURCE),
    ("spurious 5 sigma candidates returned", f"{v['flagged']} +/- {v['flagged_err']}", SOURCE),
    ("... implied rate per histogram", f"{per_hist:.3%}", "derived"),
    ("resolution elements per histogram", f"{v['elements_lo']} to {v['elements_hi']}", SOURCE),
    ("looks the set holds", f"{looks_lo:.1e} to {looks_hi:.1e}", "derived"),
    ("Gaussian expectation at 5 sigma", f"{exp_lo:.2f} to {exp_hi:.2f}", "derived"),
    ("measured defect rate per look", f"{rate_lo:.1e} to {rate_hi:.1e}", "derived"),
    ("inflation over the Gaussian tail", f"{infl_lo:.0f} to {infl_hi:.0f}", "derived"),
    ("candidates surviving the correlation veto", f"{v['coherent']} of {v['flagged']}", SOURCE),
    ("coherent fraction", f"{100*coh:.0f} +/- {100*coh_err:.0f} %", "derived"),
    ("break-even inflation for the two-stage split", f"{R_star:.0f}", "ab_split_scan.csv"),
    ("margin of the measured inflation over break-even",
     f"{margin_lo:.0f} to {margin_hi:.0f}", "derived"),
    ("pre-registered list at Z_cut = 3, nominal", f"{k_nom:.0f}", "ab_split_scan.csv"),
    (f"... if the inflation persists at Z_cut", f"{k_infl:.0f}", "derived"),
    ("claim bar, nominal list", f"{bar_nom:.2f}", "derived"),
    ("claim bar, inflated list", f"{bar_infl:.2f}", "derived"),
]

os.makedirs(TAB, exist_ok=True)
with open(os.path.join(TAB, "estimator_defects.csv"), "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["quantity", "value", "source"])
    wr.writerows(rows)

md = f"""# What a published estimator gets wrong, and what that costs

Every trials count in this repository is enumerable because a published window and a resolution can
be looked up. A network-driven scan offers neither, and it adds a failure the Gaussian arithmetic
does not cover: the significance estimator is itself imperfect. For one such estimator the failure
rate has been measured and published, so the size of the problem is not a guess.

Source: {SOURCE}. Written by `scripts/estimator_defects.py`; the two-stage numbers it compares
against come from `ab_split_scan.csv` (`scripts/ab_split_budget.py`).

## The measured rate

BumpNet predicts per-bin significances across families of mass histograms. On background-only
inputs it flags above 5 sigma in **{v['rate_analytic']:.3%}** of histograms built from smooth
analytical functions and **{v['rate_simulated']:.3%}** of those built from simulated Standard Model
samples, and applied to its **{v['histograms']:,}** application histograms with no signal injected
it returns **{v['flagged']} +- {v['flagged_err']}** spurious 5 sigma candidates, which is
{per_hist:.3%} per histogram.

Those are not fluctuations. At {v['elements_lo']} to {v['elements_hi']} elements per histogram the
set holds {looks_lo:.1e} to {looks_hi:.1e} looks, for which the Gaussian expectation at 5 sigma is
{exp_lo:.2f} to {exp_hi:.2f} -- below one. The measured rate is

| | |
|---|--:|
| defect rate per look | {rate_lo:.1e} to {rate_hi:.1e} |
| Gaussian tail at 5 sigma | {P5:.2e} |
| **inflation** | **{infl_lo:.0f} to {infl_hi:.0f}** |

one to two orders of magnitude above the statistical one. It is a **lower bound** on the fraction
`eps` of looks whose significance is simply wrong, since only mis-estimates that reach 5 sigma are
counted; at the 3 sigma level where candidates are selected, `eps` may well reach 1e-3. Those are
the two values the selection-rule study runs its GLITCH rows at (`scripts/bh_fdr_outliers.py`,
`selection_rules.csv`): `eps = 1e-4` is the rate measured here at the 5 sigma flag level and
`eps = 1e-3` the value it plausibly reaches where candidates are actually selected, so the pair spans
the measured regime rather than being chosen for convenience.

## What it does to a single-stage scan

`Z_local = sqrt(25 + 2 ln N)` calibrates the tail of `N` Gaussian looks. If the measured tail is
{infl_lo:.0f} to {infl_hi:.0f} times heavier, the effective trials factor is a property of the
network that nobody can enumerate afterwards, and the correction is not conservative but wrong by an
unknown factor.

## What it does to the two-stage split

The split converts that inflation into a countable list: whatever produced the entries, the
confirmation stage carries the trials factor of the frozen list. Two numbers price it, both on the
lensed combinatorial scan of `scaled_scan.txt` (N = {N_scan:,.0f}).

| | nominal | if the inflation persists at Z_cut |
|---|--:|--:|
| pre-registered windows at `Z_cut = 3` | {k_nom:.0f} | {k_infl:.0f} |
| claim bar `sqrt(25 + 2 ln({WIDEN:.0f}k))` | {bar_nom:.2f} | {bar_infl:.2f} |

A hundredfold inflation of the list costs {bar_infl - bar_nom:.2f} sigma on the bar, because the
list enters logarithmically -- and even that assumes the length is modelled rather than read off the
frozen document.

Against that, the split is the more sensitive procedure whenever the trials factor a single-stage
analysis would honestly have to defend exceeds **{R_star:.0f}x** the counted one
(`ab_split_scan.csv`). The measured inflation clears that break-even by a factor
{margin_lo:.0f} to {margin_hi:.0f}, {margin_round:.0f} at the round 1e2, so for an estimator of this
class the split is not a concession to auditability but the stronger procedure outright.

That comparison charges the single stage the full measured inflation, which is the pessimistic
reading. Suppressing defects upstream narrows the margin, and the veto below is exactly such a
measure -- but it has to demonstrate a factor **{margin_round:.0f}** before the verdict changes
sign.

## The defect the split cannot kill

Splitting suppresses what fails to repeat. A coherent defect -- a mismodelled background, a detector
artefact, a network that sculpts the same shape in both halves -- passes confirmation as readily as
a real signal. The published application measures how large that class is: after the veto on
cross-histogram correlations, **{v['coherent']} of the {v['flagged']}** candidates survive, a
coherent fraction of **{coh:.0%} +- {coh_err:.0%}**, clustered at the start of their histograms and
so shape-driven rather than statistical.

That fraction is the one number to measure before unblinding: evaluate the estimator on two
independent background-only samples under identical histogram definitions and count how often a flag
lands in the same window twice. It is also the BIAS row of the selection-rule study, which no
selection rule survives -- it has to be removed upstream.
"""

os.makedirs(OVR, exist_ok=True)
open(os.path.join(OVR, "ESTIMATOR_DEFECTS.md"), "w").write(md)

print(f"looks {looks_lo:.2e}-{looks_hi:.2e}   Gaussian expectation {exp_lo:.2f}-{exp_hi:.2f}")
print(f"measured {rate_lo:.2e}-{rate_hi:.2e} per look   inflation {infl_lo:.0f}-{infl_hi:.0f}x")
print(f"coherent fraction {coh:.1%} +- {coh_err:.1%}")
print(f"break-even inflation {R_star:.0f}x (margin {margin_lo:.0f}-{margin_hi:.0f}x); "
      f"list {k_nom:.0f} -> {k_infl:.0f}, "
      f"bar {bar_nom:.2f} -> {bar_infl:.2f}")
print("wrote results/tables/estimator_defects.csv, results/overviews/ESTIMATOR_DEFECTS.md")
