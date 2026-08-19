# The uncertainty budget

What the two counted bases are uncertain by, source by source. Written by
`searchbudget/stages/budget_uncertainty.py`, which varies each declared input over the range
below and **recomputes**: a resolution change moves both the looks a spectrum carries and, for a
hypothetical one, whether it can be fitted at all.

| basis | spectra | N_trials | Z_local(5s global) | band |
|---|--:|--:|--:|---|
| the public model space, published windows | 46 | 3,685 | **6.44** | +0.18/-0.16 |
| the combinatorial scan, Run 2+3, ~400 fb-1 | 4,438 | 201,136 | **7.03** | +0.23/-0.31 |
| ... with the four selection lenses | 8,211 | 362,815 | **7.11** | +0.24/-0.32 |

## Source by source

Each entry is the shift in `Z_local` [sigma]. The last column is the shift in the **difference**
between the model space and the scan, taken direction by direction, which is what the correlated part
of every source cancels in. The lens layer is not a fourth column because it moves with the scan it
sits on: its total band is +0.24/-0.32 against the scan's
+0.23/-0.31, per source in `budget_uncertainty.csv`.

| source | varied over | model space | scan | difference |
|---|---|--:|--:|--:|
| mass resolution, scale | every r x2 either way | +0.11/-0.11 | +0.13/-0.19 | +0.03/-0.08 |
| mass resolution, per channel (&dagger;) | each r independently, x2 per sigma (16-84%) | +0.03/-0.02 | +0.00/+0.00 | +0.02/-0.03 |
| mass resolution, shape | muon axes with r(M) rising to 0.10-0.20 at 3000 GeV | +0.00/-0.00 | +0.00/+0.00 | +0.00/-0.00 |
| mass resolution, prescription | worst leg / quadrature sum instead of the calibrated mean | +0.00/+0.00 | +0.00/-0.10 | +0.00/-0.10 |
| scan windows | every edge x1.4 either way (published), x1.25 (generic) | +0.03/-0.04 | +0.02/-0.02 | +0.01/-0.01 |
| yield anchor | N_ref x0.01 to x100 | +0.00/+0.00 | +0.11/-0.15 | +0.11/-0.15 |
| background slope | P = 6 to 8 | +0.00/+0.00 | +0.01/-0.01 | +0.01/-0.01 |
| fittability requirement | 30-300 events, 15-50 elements | +0.00/+0.00 | +0.04/-0.12 | +0.04/-0.12 |
| the axis set | non-peaking axes and the dilepton overlap dropped; the 15 unscanned two-body pairs added | +0.03/-0.01 | +0.00/+0.00 | +0.01/-0.03 |
| the definition of one look | N x0.5 to x Z/sqrt(2 pi) | +0.14/-0.11 | +0.15/-0.10 | +0.01/+0.00 |
| the closed-form LEE relation | exact Gaussian-tail solution instead | +0.00/-0.04 | +0.00/-0.05 | +0.00/-0.01 |
| **total** | in quadrature | **+0.18/-0.16** | **+0.23/-0.31** | **+0.12/-0.23** |

(&dagger;) an alternative reading of an input already counted in the row above it, so it is quoted
against its own median and never added: the two would double-count the resolution.

So the headline numbers are

```
model space, published windows : Z_local = 6.44 +0.18/-0.16
combinatorial scan, Run 2+3, ~400 fb-1  : Z_local = 7.03 +0.23/-0.31
... with the four lenses       : Z_local = 7.11 +0.24/-0.32
the gap between the first two  : 0.59 +0.12/-0.23
```

**The inputs the two bases share cancel in the difference.** The resolution scale and the look convention,
the two largest terms, are worth +0.03/-0.08 and +0.01/+0.00 on the
gap against up to 0.19 and 0.15 on the bars. All the gap is left
carrying is the yield model of the hypothetical scan, which the model space never uses: the difference is
therefore better determined than the scan's bar and rests on entirely different inputs from either.
Quote differences rather than bars wherever the argument allows it.

## The declared ranges

- **mass resolution, scale.** `r` is propagated from ATLAS object performance rather than quoted per
  bump hunt, so it is known to a factor of a few; x2 either way is the band every headline
  carries. On the scan the same scaling also changes which histograms can be fitted, which is why the
  band there is wider and asymmetric: a coarser resolution costs looks twice, once per spectrum and
  again through the 25-element requirement that a coarse window now fails, while a finer one buys
  elements but empties each of them.
- **mass resolution, per channel.** The scale above is fully correlated, which is the conservative
  reading. Drawing each channel's `r` independently at the same factor per sigma leaves
  N = 3,940 to 5,428 (16-84%, median 4,584), a band of
  +0.03/-0.02 sigma about that median rather than
  +0.11/-0.11: the errors average down over 46 spectra, so the
  correlated factor two is the pessimistic end of the range and not a 1 sigma. The median sits above
  the nominal count because 1/r is convex, which is a property of the log-normal draw rather than a
  statement about the resolutions.
- **mass resolution, shape.** The muon axes replace a resolution rising with p_T by a
  window-averaged constant, the one place the constant-`r` framework is a real approximation.
  Counting them with `r(M)` rising linearly from 0.02 at
  200 GeV to 0.10-0.20 at
  3000 GeV brackets it, and the approximation turns out to be a good one: over the
  `m(mumu)` window the rising form is worth an effective flat
  `r_eff` = 0.047 at the central anchor
  (0.053 to 0.040 across
  the band) against the 0.050 declared, so the whole question is worth less than
  0.01 sigma on the bar.
- **mass resolution, prescription.** The scan takes `r = 1/2 sqrt(mean sigma^2)` over the object
  group, calibrated on the published channels. The alternatives are the worst leg, which is the
  convention of the model-space budget, and textbook quadrature propagation.
- **scan windows.** Published edges are read off papers; what is uncertain is how far the curated
  search family extends, taken as x1.4 on each edge. A hypothetical spectrum has a
  declared generic window instead, varied by x1.25.
- **yield anchor**, **background slope**, **fittability requirement.** The three declared inputs of
  the statistics requirement. All three act on the scan only, since published windows are exempt.
- **the axis set.** Down: the axes motivated only by non-peaking models
  (`m(multi)`, `multilepton`) dropped, and the 88 looks the
  dark-photon axes double-count against the high-mass dilepton axes removed. Up: the 15 object
  pairs with no published axis added at 698 looks.
- **the definition of one look.** The convention with no measurement behind it. A resonance spans more
  than one element, which argues for fewer independent looks (x0.5); the up-crossing form of
  the Gross-Vitells estimate argues for more, by Z/sqrt(2 pi) = 2.57 at this
  Z, since Rice's formula counts (1/2 pi)(1/r) ln(M_hi/M_lo) exp(-Z^2/2) up-crossings against the
  Gaussian tail exp(-Z^2/2)/(Z sqrt(2 pi)) that element counting multiplies. This is the largest
  single term on the model space and, on the scan, second only to the resolution scale; being common
  to every basis, it cancels in every difference.
- **the closed-form LEE relation.** `Z_local = sqrt(25 + 2 ln N)` is the asymptotic solution of
  `N p_local = p(5 sigma)`; solving it exactly gives 6.40 rather than
  6.44, so the closed form is marginally strict.

## Conventions, not uncertainties

These are choices, reported as alternatives and never added to the band.

| choice | alternative | model space | scan |
|---|---|--:|--:|
| counting granularity | published event selections instead of inclusive spectra | 6.53 | - |
| the dataset | Run 2 alone instead of Run 2 + Run 3 | - | 7.00 |
| the lens layer | the four lenses added to the scan | - | 7.11 |
| lens efficiency | every declared efficiency x0.1 to x10 | - | 7.08 to 7.15 |
| the costliest single axis | the largest contributor removed | 6.42 | - |

The shape of the hypothetical scan is a design choice in the same sense: ten object types, at most
four per mass, one selection lens at a time. It fixes what is being priced rather than carrying an
error, and `results/tables/scaled_scan.csv` prices the variants side by side.

Source: `searchbudget/stages/budget_uncertainty.py` ->
`results/tables/budget_uncertainty.csv`,
`results/tex/uncertainty_table.tex`. The bases themselves: `results/overviews/SEARCH_BUDGET.md`,
`results/tables/scaled_scan.txt`.
