# Method notes

Conventions and the reasoning behind them. Anything here that a script depends on is *only* here —
the code keeps its docstrings to a line or two.

## Search budget

**One bump hunt = one invariant/transverse-mass spectrum.** Every model peaking in the same
spectrum is tested by the same search, so the number of searches is the number of distinct bump
observables (after merging same-axis decay splits via `OBS_MERGE`: `m(HH) 4b` is an event selection
of the `m(HH)` axis, not a second spectrum).

Two rules decide what "distinct" means, and both are physical rather than bookkeeping:

* **A b-jet is a different object from a light jet.** A b-tagged leg is selected on and scanned as
  its own spectrum, so `m(eb)`/`m(mub)`/`m(taub)` are axes distinct from `m(ej)`/`m(muj)`/`m(tauj)`,
  and `m(bj)` (b\* → bg, b-tagged dijet) is distinct from `m(jj)` — the same rule that always kept
  `m(bb)` separate from `m(jj)`.
* **Lepton flavour is a mass axis, not an event selection.** `m(ee)` and `m(mumu)` are separate
  spectra with their own resolutions (`bump_observables.FLAV_SPLIT`): electrons are measured in the
  EM calorimeter (`r` ~ 0.015, ~mass-independent), muons by track sagitta (`r` rising to ~15% at
  3 TeV, entered as a window-averaged `r` ~ 0.05). At low mass the ordering reverses — a few-GeV
  `mumu` pair is sharp while `ee` merges into one cluster, which is why the `(Zd)` electron channel
  starts at 1 GeV and the muon channel at 0.3 GeV. Because flavour is an axis, it is **not** also a
  multiplier in `NSEL` (`public_obs_map.py`); adding it back there would double-count it.
  Cross-check: the selections-level channel count must stay at **94**.

With a constant fractional resolution `sigma_M = r·M`, the effective independent looks over a
scanned window are the Gross-Vitells resolution elements

```
n_s = (1/r) ln(M_hi / M_lo)      summed over disjoint scan segments
N   = sum over spectra of n_s
Z_local = sqrt(25 + 2 ln N)      for a 5 sigma GLOBAL discovery
```

`r` is the one real physics input, so every headline is quoted as a band: `r` ×2 (→ `N` ×0.5) and
×0.5 (→ `N` ×2). Because `N` enters through `ln N`, that factor-4 spread in `N` moves `Z_local` by
about ±0.1σ.

**Window convention (uniform, non-circular).** Every observable is counted over its
**published-search scan window** (`bump_observables.SCAN`, each with a per-channel source), extended
where that search family's published signal grid prepares a wider scan. The tempting alternative —
counting over the mass range for which signal samples happen to exist — is wrong twice over: it
understates sparsely sampled channels, since a search does not become cheaper because few samples
were made for it, and it ties a model-space quantity to one collaboration's production choices. The
model-independent envelope `[floor, sqrt(s)]` is kept as a reference upper bound, where the floor is
`max(lowest trigger, kinematic 2·m_daughter, Z-peak)`.

### The narrow-resonance assumption

`n_s` counts **detector** resolution elements. That is the right step size only while the natural
width of the signal stays below `r`; a resonance broader than the resolution correlates
neighbouring mass points, so fewer of them are independent looks. The direction matters more than
the size: counting resolution elements for a wide signal **over**-counts `N`, which makes `Z_local`
too strict rather than too loose. Nothing in the budget can be made *easier* by a width correction.

`public_obs_map.WIDTH` places all 43 public model classes in four bands, each with the Γ/M its
published benchmark carries:

| band | classes | what it means |
|---|--:|---|
| `narrow` | 22 | Γ/M below `r` on every axis the class populates (HAHM ~ε², LQ `λ²/16π` ≈ 2%, W_R 2–3%, q\* 2–4%) |
| `benchmark` | 13 | narrow at the point ATLAS publishes, broad elsewhere in its own parameter space (Z′_SSM 3% against `r`=0.015 on `m(ee)`; the DM-mediator coupling grid spans 1% to >30%; single-VLQ 10–50%; U1 ~20% at the flavour-anomaly point) |
| `broad` | 3 | already wider than `r` at the standard benchmark: KK gluon 15–30% against `r`=0.08 on `m(tt)`, coloron/axigluon, composite/NJL |
| `nonpeak` | 5 | no Breit-Wigner peak at all: QBH is a threshold turn-on, ADD/HEIDI a non-resonant tail, Type-III seesaw and VLL are pair-produced counting signatures, toponium is a threshold effect pinned at 2·m_t |

Peaking-ness belongs to the (model, axis) pair rather than to the model, so `NONPEAK_ON` records
the exceptions: prompt HNL is a genuine `m(lljj)` resonance but a counting signature in
`multilepton`, and 2HDM / heavy-Higgs / top-philic scalars interfere with the SM `tt̄` continuum,
giving a peak-dip lineshape rather than a bump on `m(tt)`.

**What it costs.** Only two axes are motivated *exclusively* by non-peaking models: `m(multi)`
(ADD/HEIDI + QBH) and `multilepton` (HNL + Type-III + VLL). Dropping both takes `N` from 3685 to
3620 and `Z_local` from 6.44 to 6.43. `search_budget.py` recomputes that sensitivity rather than
quoting it, so it follows the map. The combinatorial scan is unaffected by construction: the four
axes `scaled_scan.MOTIVATED` maps to `None` (`m(multi)` and the three transverse masses) are
exactly the non-peak observables, so `m(multi)` never entered the 3603. `multilepton` is the one
axis that counts as motivated in the scan while none of its models peaks there.

### Splitting the budget finer

**Event selections** (`search_budget.py`, the `NSEL` column): real searches slice a final state much
finer than the inclusive spectrum — b-tag categories, boosted vs resolved, sub-decay modes,
0/1/2-lepton bins, ambiguous pairings (a leptoquark pair search scans all four lepton-jet
pairings). `public_obs_map.NSEL[obs]` is a hand-curated, publication-anchored count of the spectra a
real search family scans, quoted as a band, and `N = Σ NSEL(obs) · n_s(obs)`.

Channels within one search family are partly correlated and are often statistically *combined* into
one limit, so treating each as an independent scan is a mild **over-count**: the true `N` sits
between the inclusive and the selections level. That is the direction one wants for a discovery
threshold.

A third level exists but needs a sample inventory to evaluate, so it is not computed here: splitting
by *production mechanism* on the same mass axis (resolved vs boosted, inclusive vs VBF-tagged,
prompt vs associated). Two effects compete there — multiplicity against each split scanning a
narrower window (pair/cascade production caps at √s/2, single/VBF ≈ √s) — and because `n_s` is
logarithmic, splitting adds trials **sub-linearly** in the number of modes.

### The combinatorial budget

`combinatorial_budget.py` counts exclusive object categories × all 2-to-4-object invariant masses.
The object budget is a **design choice of this study**, chosen as the shape a data-directed scan of
this kind naturally takes, and the point is to price it rather than to describe any particular
implementation:

* categories are **exclusive** multiplicity bins over `(n_e, n_mu, n_j, n_b, n_Z, MET)` — disjoint
  event sets, hence statistically independent looks;
* object budget `n_e ≤ 2, n_mu ≤ 2, n_j ≤ 3, n_b ≤ 3, n_Z ≤ 1`, with `2 ≤ Σn ≤ 4` (MET excluded);
* single-lepton trigger: a category is populated only if it contains a lepton — a loose e/μ, or a
  leptonic Z whose two OS-SF leptons fire the trigger and leave the loose collections;
* MET splits every category in two (0met/1met) but is never used in a mass and does not count
  against the four-object limit;
* a category with exactly two loose leptons (ee, eμ, μμ) splits into OS and SS;
* spectra ("groups"): every size-2…4 subset of the *indexed* objects is its own histogram —
  `j0j1`, `j0j2`, `j1j2` are three groups.

That gives **150** multiplicity categories (204 after the OS/SS split) and **1094**
(category, mass-group) combinations, 1502 histograms after the same split, of which **712** hold
enough events to be fitted (next section). `N` counts the OS/SS-split histograms, so quote the two
together — the pairing of an unsplit category count with a split `N` is the easy mistake here.

Per-group resolution comes from the object composition, `r = ½·sqrt(⟨sigma_i²⟩)` with fractional
p_T resolutions `sigma(e/mu/Z) = 0.04`, `sigma(j) = 0.10`, `sigma(b) = 0.20`; this reproduces the
per-channel values in `bump_observables.RESOLUTION` (jj → 0.05, ll → 0.02, bb → 0.10, lj → 0.04).
Note the ½ prefactor is an empirical calibration to those channels, not textbook propagation
(which for a 2-body mass gives `½·sqrt(Σ sigma_i²)`, a factor √2 larger), and it makes `r`
independent of group size — a 4-body mass is assumed as sharp as a 2-body one. Both choices only
move `ln N`.

### Only fittable histograms are counted

A trials count may contain only histograms a search could actually fit, so `yield_model.py` imposes
one requirement wherever a spectrum is hypothetical: **at least 100 events, and at least 25 bins
holding one event or more**, with one bin = one resolution element, the unit `n_s` is built from.

The yield model behind it is declared and order-of-magnitude. Every background here is a steeply
falling mass spectrum, so

```
n(m) = N_REF · W · (m/M_REF)^(1-P) · (r/R_REF)     events in one resolution element
```

with `N_REF = 1e6` events per element in the light-jet pair spectrum at `M_REF = 1 TeV` for a Run-2
dataset, `R_REF = 0.05` the resolution that anchor is quoted at, `P = 7`, and `W` the product of a
per-object factor `F` relative to a light jet. A sharper channel has narrower bins and therefore
fewer events in each, which is the `r/R_REF`.

`F` is **statistics, not signal cross section**: a histogram's content is its background, so a tagged
hadronic object costs its mistag rate and a lepton or genuine MET costs the price of an electroweak
process against QCD. The factors (`j` 1, `b` 0.1, `V`/`t` 0.02, `H` 0.01, `γ` 4e-3, `τ_h` 3.5e-3,
`e`/`μ`/MET 3e-3, leptonic `Z` 1e-5) are set so the model reproduces the published symmetric channel
of each type to an order of magnitude: run `python3 scripts/yield_model.py` for that table.

Because `n(m)` falls, the populated part of a window is a prefix of it, so the requirement acts by
**truncating each window at the one-event mass** and dropping the spectrum when fewer than 25 elements
survive. The 25-element test is what binds; at this slope a spectrum with 25 populated elements holds
far more than 100 events.

Three honest caveats, and all three err towards more looks rather than fewer:

* a single power law overestimates yields below a few hundred GeV, where real spectra turn over, which
  can only make the 100-event test easier;
* it puts the one-event mass of the dijet spectrum at 10 TeV against a published fit stopping near 8;
* the factorised form charges every object the full price of its own production, while real objects
  arrive in pairs from one boson, so it *under*-counts high-multiplicity leptonic categories. This is
  the one caveat in the other direction, and the yield-anchor band below is the answer to it.

**Published windows are never gated**: a published search demonstrates its own feasibility. The size
of that exemption is printed — applied to the 42 published axes this alphabet can form, the
requirement would leave 31 of them and `N = 2566` of 3573.

The requirement is a large cut on the enumeration and a small one on the bar. Scaling the yield anchor
by ×100 and ×0.01 moves the ten-object scan over `1284 … 7768` fittable spectra out of 21644 and
`Z_local` over `6.84 … 7.12`, i.e. four orders of magnitude in assumed statistics for ±0.14σ. Breadth
enters through `ln N`, and that is why the answer survives a model this crude.

**Which dataset the anchor describes: Run 2, 140 fb⁻¹.** That is the consistent choice, because every
scan window in the budget comes from a published Run-2 search family, so gate and windows describe the
same data. A larger dataset enters as a rescaling of `N_REF` (`scaled_scan.DATASETS`), and buys reach
slowly, the one-event mass going as `lumi^(1/(P-1))`, the sixth root. Run 2 plus Run 3, taken as ×3 and
ignoring the rise in high-mass cross sections from 13 to 13.6 TeV (which acts in the same direction),
extends each window by a fifth in mass and gives

| | fittable spectra | `N` | `Z_local` |
|---|---|---|---|
| ten objects, Run 2 | 3603 of 21644 | 1.6e5 | 7.00 |
| ten objects, Run 2+3 | 4438 of 21644 | 2.0e5 | 7.03 |
| … with the four lenses, Run 2 | 6611 histograms | 2.9e5 | 7.08 |
| … with the four lenses, Run 2+3 | 8211 histograms | 3.6e5 | 7.11 |

so the full dataset does not change the conclusion: the lensed scan still fits the `5e5` budget, and
statistics still binds before the trials factor does. On the five-object grid the same ×3 leaves 755 of
1502 histograms rather than 712, `N = 3.9e4`, `Z_local = 6.79`.

### The scaled-up scan and its trials budget

`scaled_scan.py` runs the same enumeration over the alphabet a general search would actually have,
and then asks what fits inside a fixed trials budget. The rules that differ from the five-object
scan, all of them design choices:

* **ten object types** — e, mu, hadronic tau, photon, light jet, b-jet, boosted top, boosted W/Z,
  boosted H, leptonic Z — with **no per-type ceiling**, since those describe one particular grid;
* **strictly at most four objects per category**, MET excluded, since MET splits every category but is
  **never an ingredient of a mass**: no transverse mass is formed, so the `mT` axes of the model-driven
  budget are outside this scan's reach. A category requiring MET carries its yield factor;
* **any trigger**: a category needs no lepton;
* resolutions are **derived, not declared**: `sigma = 2r` of each object's own symmetric published
  channel (`m(ee)` → 0.030, `m(tautau)` → 0.24, `m(gammagamma)` → 0.02, ...), which is the inversion
  `two_body_matrix.py` uses.

The enumeration is 2412 categories and 21644 possible histograms; **3603** of them can be fitted, at
`N = 1.6e5` and `Z_local = 7.00`. The priority order below is the policy for cutting a scan down to
`TRIALS_BUDGET = 5e5`: spectra are ordered and the scan is the longest prefix that fits, stopping at the
first spectrum that does not (rather than topping up with whatever cheap spectrum still fits the
remainder).

0. **every model-motivated axis once**, in the best-populated category it appears in, so no motivated
   axis can be lost to the budget. `MOTIVATED` maps each of the 46 observables of the model-driven
   budget onto the composition(s) a scan would build it from; 42 of them have one (`m(multi)` has no
   fixed composition, and the three `mT` axes need MET in the mass), and they collapse onto **38
   distinct compositions**, of which 36 can be fitted somewhere (`m(tautau)` and `m(Zt)` cannot);
1. **those same axes in their remaining categories**, highest yield first;
2. **everything else**, highest yield first.

Yield is `yield_model.F` multiplied over the category's content, i.e. the same model that decides
fittability, so "best-populated first" means most background events. Tier 0 costs 2.7e3 looks against
the 3.7e3 of the model-driven budget, which is the one place the two prescriptions can be checked
against each other.

The same split also measures how much of a combinatorial scan theory motivates at all, so the report
prints it per dataset: the tier boundaries move with fittability. On Run 2 the tiers hold 43 / 2023 /
1537 spectra, i.e. 2066 of 3603 (57 %) on a model-motivated axis over 36 of 217 fittable compositions;
on Run 2 + Run 3, 45 / 2366 / 2027, i.e. 2411 of 4438 (54 %) over 37 of 258. The motivated set is
bounded above by the 46 axes and the rest of the scan is not bounded at all, so the share falls as the
dataset grows.

Counting basis: a *spectrum* here is one axis in one category, i.e. one fitted histogram, and an
OS/SS-split category holds two of them. `charge_split` is therefore the weight in every spectrum count
below, exactly as in `N` — mixing the row count (17600 `(category, mass-group)` pairs) with the
histogram count is the easy mistake.

**The outcome is that the budget never binds.** The whole fittable scan is 1.6e5 looks, and with the
lenses of the next section 2.9e5, both inside 5e5, so no spectrum has to be dropped and the priority
order is never applied. What limits a scaled-up scan is statistics, not the trials factor: 18041 of the
21644 histograms cannot be fitted, while the 3603 that can cost a bar of `Z_local = 7.00`.

### Selection lenses

A lens is an extra event-level requirement laid over an unchanged mass axis: one more view of the same
spectrum, one more histogram, one more look. Of the eight handles a wide search would reach for, four
are **already inside the enumeration** and would be double counted: high MET is the category's met
split, high jet or lepton multiplicity is exactly what the exclusive categories are, and b-tag and tau
enrichment are the `b` and `T` types of the alphabet. The four that are orthogonal to both the object
content and the mass axis are priced:

| lens | applies when | efficiency | views | looks | ruled out by statistics |
|---|---|---|---|---|---|
| high HT or Meff | the category holds an object outside the mass | 0.1 | 2134 | 9.0e4 | 1145 |
| displaced activity | any reconstructed mass | 1e-3 | 743 | 2.9e4 | 2860 |
| forward jet pair (VBF) | two of the four slots are free for the tag jets | 0.02 | 53 | 2.7e3 | 29 |
| ISR jet | one slot free, and the window reaches below 200 GeV | 0.2 | 78 | 4.0e3 | 652 |

Every rule is deliberately conservative:

* **one lens at a time**, never a product of two. Pairs would square the count;
* a lens **leaves the axis, its resolution and its window alone**, so a lensed view costs the same
  looks as its inclusive parent. The one exception is the ISR lens, which buys acceptance only at the
  low-mass end and is therefore capped at 200 GeV, i.e. costs a fraction of a look;
* the objects a lens needs **count against the same four-object ceiling**, which is what makes the VBF
  lens rare: a forward tag pair fits only in a two-object category;
* an HT or Meff threshold on a mass with nothing else in the event is a cut on the resonance mass
  itself, not an independent look, so the lens requires activity outside the mass;
* a lens **costs the statistics of its own requirement** (`yield_model.LENS_EFF`) and the view then has
  to pass the same fittability test, which is what removes most of them: 4686 of the 7694 possible views
  are ruled out by statistics, the displaced lens worst of all at an efficiency of 1e-3.

What survives is **1.8 histograms per spectrum**, taking the ten-object scan to 6611 histograms and
`N = 2.9e5` (`Z_local = 7.08`), still inside the 5e5 budget. Lens views inherit their parent's tier and
yield and rank immediately behind it, so with nothing to cut the ordering never matters here either.

`composition_gap.py` then asks which of those flavour compositions any published ATLAS bump hunt
has ever scanned. The coverage mapping is deliberately **generous** — a composition counts as
covered if any published search scans a mass built from those object types, allowing a published
"V" to be a Z, a jj pair or an ll pair, a published "t" to contain b+jets, and so on — so the
uncovered list is a lower bound on the gap.

### The publication census is a different base

`published_census.py` counts from the other end: not the spectra models motivate, but the searches
ATLAS has actually published. **Its 86 entries and the budget's 46 spectra must never be summed**, and
a raw entry count is not a trials factor: the publication record keeps two analyses on a shared mass
axis as separate entries when they are separate papers, so the 86 entries carry only 62 distinct
observables, while the budget merges everything onto one axis and counts resolution elements along it.

Scope is matched to the budget's on purpose: bump hunts for new states, so hadron-spectroscopy
measurements are excluded even when they are literally bump hunts in a mass spectrum.

The crossing between the two bases is recorded in the data rather than asserted: every census row
carries a `budget_axis` naming the canonical axis (or axes) it scans, `-` when it falls on none of
the 46, and a `scan_GeV` transcribing the range that entry actually scanned (`fixed` for a
single-mass search, empty where the census does not record one). Two things read it.

**The stale list.** Revisiting a stale spectrum is free in trials **only** when its axis is already
in `N` (12 of the 17); the other 5 extend the axis count and are priced like any new spectrum.

**The census priced in trials** (`census_budget.py`). With an axis and a range per entry, the
publication record can be run through the same rule as the model space, and the comparison stops
being a category error: charging every published search for the range it scanned gives `N = 7,710`
(`Z_local = 6.55`) over 100 (search, axis) looks, and counting each axis once over the union of every
published range on it gives `N = 3,672` (`6.44`) — against the model side's 3,685 and 6.44, i.e. the
two enumerations land on the same bar from opposite directions. Where the census records no range,
the axis' own published window is the fallback, which is what separates the two bases; entries on no
axis are priced at `RES_DEFAULT`, and the 12 that carry neither an axis nor a range (the generic
multi-spectrum anomaly-detection scans and the displaced programs, whose trials belong in the
combinatorial count instead) stay unpriced, so `N` is a lower bound by that much.

## Uncertainties

`budget_uncertainty.py` is the one place that says how well any of this is known. It moves each
declared input over a stated range and **recomputes** the count rather than scaling `N`, which
matters for every hypothetical spectrum: a coarser `r` costs looks per spectrum *and* fails the
25-element test on windows that used to pass, so the two effects do not factorise. Its ranges live at
the top of that file, one constant per input, and nothing else in the repository declares an
uncertainty.

Three rules keep the accounting honest.

* **What counts as one look is priced, and it leads.** Counting elements of width `sigma_M` and
  calling each an independent test is a convention, not a measurement. Adjacent elements are
  correlated — a resonance lifts more than one — which argues for fewer looks; the up-crossing form of
  Gross-Vitells argues for more. For a smooth unit-variance Gaussian process in `x = ln M` with
  correlation length `r`, Rice's formula gives `<N_Z> = (1/2pi)(1/r) ln(M_hi/M_lo) exp(-Z^2/2)`
  up-crossings of level `Z`, which is the element count times `Z/sqrt(2pi) ~ 2.6` once compared
  against the Gaussian tail `p_local` that the Bonferroni step multiplies. The band taken is
  `0.5 N` to `N Z/sqrt(2pi)`, worth `+0.15/-0.11 sigma`: more than the factor-two resolution band.
* **Correlated sources cancel in differences.** Every input except the yield model enters both counted
  bases the same way, so the report carries a third column, the shift in the *difference* between
  them, computed variation by variation. The bars are known to `+0.18/-0.16` (model space) and
  `+0.23/-0.31` (scan); the 0.59σ gap between them to `+0.12/-0.23`, and that residual is the yield
  model, which the model space never uses. Quote differences, not bars, wherever the argument allows.
* **Conventions are not uncertainties.** Granularity (46 spectra vs 94 channels), the dataset the
  yields are priced on, the lens layer, and the shape of the hypothetical scan change *what* is
  counted. They are reported as alternatives and never added to the band. The published-program row is
  a literature count, so its band is the literature range `1e4` to `1e5`.

Two by-products worth keeping in mind. The per-channel resolution scatter (each `r` drawn
independently at a factor two per sigma) gives `±0.03σ`, so the correlated factor-two band is the
pessimistic end and not a 1σ; and the window-averaged constant that stands in for the muon axes'
rising `r(M)` reproduces the rising-resolution integral to 5% in effective `r`, so that
approximation — once described as the largest modelling approximation here — is worth under 0.01σ.

Monte Carlo results carry their own errors and are quoted with them: the BH power in
`bh_fdr_outliers.py` is a binomial fraction of `T = 2e4` toy experiments (`±0.3` percentage points at
`mu = 5`), while the threshold and argmax rules in the same table are quadrature integrals with no MC
error at all. The A/B reach is analytic, reproduced by toys to 0.03σ, and `ab_split_budget.py` prints
the price of the split over the `N` band and the window-widening factor `w`: `+0.32` to `+0.49σ`,
essentially flat in `N`.

## Two-stage A/B unblinding

Split the dataset into fractions `f` (A, exploration) and `1-f` (B, confirmation). Scan **all**
spectra in A with **no** global correction — that is selection, not inference; pre-register every
window with local `Z_A ≥ Z_cut`; unblind only those in B; claim on the B-only p-value corrected for
the `k` unblinded windows, whose trials factor is exactly countable.

Reach at 50% **joint** power (the honest figure of merit, validated by `ab_split_toys.py`): the two
stages fluctuate independently, so for a signal of full-dataset local significance `Z_full`

```
P(disc) = Phi(sqrt(f)·Z_full - Z_cut) · Phi(sqrt(1-f)·Z_full - Z_B_req)
Z_B_req = sqrt(25 + 2 ln(widen · k_eff)),   k_eff = N·p1(Z_cut) + 1
```

and the reach is the `Z_full` at which `P = 50%`. The median arithmetic
`Z_med = max(Z_cut/sqrt(f), Z_B_req/sqrt(1-f))` *underestimates* the reach wherever both
constraints bind (each stage at ~50% → joint ~25%); it is kept only as a dashed reference. The
single-stage 50%-power reach is `sqrt(25 + 2 ln N)`.

The toy model treats the budget's `N` effective looks as iid standard normals under
background-only, with `z_A = sqrt(f)·mu + g_A`, `z_B = sqrt(1-f)·mu + g_B`, recovering the
full-dataset scan exactly as `z_full = sqrt(f)·z_A + sqrt(1-f)·z_B ~ N(mu, 1)`.

## The statistics study

`results/overviews/MAX_OF_GAUSSIANS.md` derives, for `n = 30,000` bins, the exact law of the
background maximum and compares three stage-1 selection rules on one ROC — argmax, fixed threshold,
and Benjamini-Hochberg FDR. Two implementation notes:

* BH depends on the p-vector only through its smallest few order statistics, and the `m` smallest
  of `n` iid uniforms have the exact representation `U_(k) = Gamma_k / Gamma_(n+1)` with `Gamma_k`
  a sum of `k` iid Exp(1). One gamma draw plus `m` exponentials therefore gives the exact joint law
  without ever generating `n` values. The scripts assert at runtime that the rank truncation `m`
  was never reached.
* `K(q)` vectorises over experiments and over the `q`-grid by a suffix minimum: with
  `R_k = n·p_(k)/k` and `C_k = min_{j≥k} R_j` (non-decreasing in `k`), `K(q) = #{k : C_k ≤ q}`.

`bh_fdr_outliers.py` and `ab_split_outliers.py` repeat the comparison with a deliberately
*imperfect* significance estimator, under two defect classes: **glitch** (an incoherent artefact,
independent between dataset halves) and **bias** (a coherent one that reproduces in both). The
ordering of the three rules is set by whether a rule is rank-based or absolute, and it does not
depend on the defect magnitude; the size of the gaps does.

Monte Carlo in this repository is seeded explicitly, so the figures and the quoted validation
numbers are reproducible bit-for-bit. `scripts/plot_style.py` holds the shared palette; import it
rather than redefining colours per figure. `scripts/obs_labels.py` does the same for the observable
labels, so a figure axis and an appendix row spell a spectrum the same way; it is standard library
only, so the document generators can use it without matplotlib installed.
