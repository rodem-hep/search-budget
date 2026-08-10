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
(category, mass-group) combinations, 1502 histograms after the same split. `N` counts the
OS/SS-split histograms, so quote the two together — the pairing of an unsplit category count with a
split `N` is the easy mistake here.

Per-group resolution comes from the object composition, `r = ½·sqrt(⟨sigma_i²⟩)` with fractional
p_T resolutions `sigma(e/mu/Z) = 0.04`, `sigma(j) = 0.10`, `sigma(b) = 0.20`; this reproduces the
per-channel values in `bump_observables.RESOLUTION` (jj → 0.05, ll → 0.02, bb → 0.10, lj → 0.04).

`composition_gap.py` then asks which of those flavour compositions any published ATLAS bump hunt
has ever scanned. The coverage mapping is deliberately **generous** — a composition counts as
covered if any published search scans a mass built from those object types, allowing a published
"V" to be a Z, a jj pair or an ll pair, a published "t" to contain b+jets, and so on — so the
uncovered list is a lower bound on the gap.

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
rather than redefining colours per figure.
