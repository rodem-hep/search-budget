# Two-stage A/B unblinding: managing the LEE by data splitting — review & logic

**Proposal under review:** split the dataset in two — half **A for exploration** (scan everything),
half **B for confirmation** (unblind *only* the regions selected in A). The Look-Elsewhere penalty
is then paid only on the pre-registered B windows instead of on the whole program.

**Which program.** The design is priced here on two bases, both in `ab_split_scan.csv`
(`ab_split_budget.py`): the model space at event-selection granularity (**N = 7,030**), which is
where the derivations below are worked out, and the **lensed combinatorial scan** of
`scaled_scan.txt` (**N = 362,815**), which is the scan a network-driven hunt would actually run and
therefore the headline. Every window is allowed `w = 3` resolution elements of mass freedom
(≈ ±2σ_M) unless stated; `w = 1` pins it.

| basis | N | single stage | optimised split | naive 50/50, Z_cut = 3 |
|---|--:|--:|--:|--:|
| model space, event selections | 7,030 | 6.54 | 7.02 (+0.49) | 7.99 (+1.46) |
| **combinatorial scan + lenses** | **362,815** | **7.11** | **~7.6 (+0.47)** | **8.90 (+1.79)** |

**Verdict up front.** The idea is sound and its *main* value is real — the trials factor in B
becomes **exactly countable** (a pre-registered list, no Gross–Vitells modelling of
a machine-generated scan) and the procedure is immune to human/procedural biases. But two things
must be understood before adopting it:

1. **The split costs sensitivity, and 50/50 costs a lot.** The LEE you avoid is logarithmic
   (cheap) while the luminosity you give up costs √2 on Z (expensive). At 50% discovery power
   (toy-validated), on the combinatorial scan: a single stage exactly corrected needs a **7.11σ**
   signal; the optimized asymmetric split (**f ≈ 0.2–0.3** exploration) needs **~7.6σ**
   (+0.47σ at `w = 3`, +0.32σ with pinned windows); the naive **50/50** split needs **8.90σ**
   (+1.79σ). The ~0.5σ at the optimum is the honest, irreducible price of the countable trials
   factor, and it barely moves with N: a factor two in the trials count changes it by 0.01σ.
2. **The split only kills *statistical* look-elsewhere.** A spurious bump from correlated
   systematics — background mis-modelling, and specifically a network that sculpts a fake
   peak — reproduces in **both** halves and will happily "confirm". A/B splitting must be
   combined with decorrelated ML trainings (train-on-A→predict-B and vice versa) to also attack
   that failure mode. That class is not hypothetical and its size is measured:
   **51 ± 7 %** of one published scanner's spurious candidates are coherent
   (`ESTIMATOR_DEFECTS.md`).

**Why bother at all.** Because for a network-driven scan the single-stage correction is not merely
awkward but *wrong by an unmeasured factor*. The one estimator whose failure rate has been published
mis-estimates 10⁻⁵–10⁻⁴ of its looks, **78–310×** the Gaussian tail it is corrected against, while
the split's break-even sits at **14×** (`ESTIMATOR_DEFECTS.md`, section 2b here). At that measured
rate the split is not a concession to auditability — it is the more sensitive procedure outright.

Numbers: `searchbudget/stages/ab_split_budget.py` → `results/plots/ab_split_reach.png` and
`ab_split_crossover.png` (reach vs N: no crossover at any trials count — section 2b); toy Monte
Carlo validation: `searchbudget/stages/ab_split_toys.py` → `results/plots/ab_toys_*.png`
(all public-information-only; budget N from `SEARCH_BUDGET.md`).

## 1. The logic, step by step

Let the full dataset have luminosity L, split A : B = f : (1−f), statistically independent
(random event-level split). Significances scale as Z ∝ √L, so a signal worth `Z_full` on the
full dataset gives median `Z_A = √f·Z_full` and `Z_B = √(1−f)·Z_full`.

**Stage A — exploration (no LEE bookkeeping at all).**
Scan every spectrum, every selection, with the full machinery. Because A is used only to *select*
— never to claim — no global correction is needed here. This is the big conceptual
simplification: the entire N ≈ 6,400-look budget of `SEARCH_BUDGET.md` is spent in a half where
trials don't have to be corrected for, only *counted* to predict how many selections background
will produce. Every window with local `Z_A ≥ Z_cut` is **pre-registered**: exact mass window,
spectrum, event selection, test statistic — frozen in writing before B is touched.

Expected number of background selections (pure statistics, one-sided):

    k_bkg = N · p1(Z_cut),   p1(Z) = ½ erfc(Z/√2)

| Z_cut | k_bkg (N = 7,030) |
|---:|---:|
| 2.0 | 150.1 |
| 2.5 | 41.0 |
| 3.0 | **8.9** |
| 3.5 | 1.5 |
| 4.0 | 0.2 |

`Z_cut = 3` is the natural working point: ~9 windows to confirm — few enough to scrutinise
individually, many enough not to strangle the signal efficiency of stage A. (This is the same
arithmetic as `EXCESS_COUNTING.md`: tens of 3σ wiggles are *guaranteed* under background-only;
stage B is precisely the machine that disposes of them.)

**Stage B — confirmation with a countable trials factor.**
Unblind only the `k_eff ≈ k_bkg + 1` pre-registered windows (the +1 is the signal's own).
Because the peak position jitters between halves by ~σ_M, each window must be widened to
±(1–2)σ_M, i.e. ~3 resolution elements of freedom per window. The B-stage global correction is
then **exact and tiny**:

    Z_B required for 5σ global  =  √(25 + 2·ln(3·k_eff))
    Z_cut = 3  →  k_eff ≈ 9.5  →  Z_B ≥ 5.65     (vs 6.54 for the single-stage full scan)

No Gross–Vitells toys, no up-crossing counting, no argument about what "independent look" means
— the trials factor is the number of lines in a frozen document. That auditability is the
strongest argument for the scheme.

**Discovery claim = B-only p-value** (corrected by 2·ln(3·k_eff)). A+B combined is then used for
*characterisation* (mass, cross-section) but not for the discovery claim itself — re-using A in
the claim re-imports the selection bias (you picked A's upward fluctuations, so A's contribution
is biased high; removing that bias needs conditional/selective inference, doable but hard —
cleanest to just not do it).

## 2. The cost: reach at 50% power (why 50/50 is the wrong split — and why there's no free lunch)

The two stages fluctuate **independently**, so the discovery probability for a signal of
full-dataset local significance `Z_full` is a *product*:

    P(discovery) = Φ(√f·Z_full − Z_cut) · Φ(√(1−f)·Z_full − Z_B_req),
    Z_B_req = √(25 + 2·ln(3·k_eff))

and the reach is the `Z_full` where P = 50%. (A first draft of this note used the simpler
"median arithmetic" `max(Z_cut/√f, Z_B_req/√(1−f))`, which claimed the optimized split even
*beats* the single-stage scan by 0.1σ — **the toys exposed that as wrong**: near the optimum
*both* stages sit at ~50% pass probability, so the joint power there is only ~25%. The product
formula below reproduces the toys to 0.03σ.)

| strategy | Z_full at 50% power (analytic) | toys | vs single-stage |
|---|---:|---:|---:|
| single-stage full-dataset scan (N = 7,030) | **6.54** | 6.52 | — |
| A/B **50/50**, Z_cut = 3, B-only claim | **7.99** | 7.72* | **+1.46** |
| A/B optimized: **f ≈ 0.2–0.3**, Z_cut = 2–3 | **7.02–7.10** | 7.14* | **+0.5** |

*toy numbers use exact one-look windows (no ±2σ_M widening), worth ~0.1–0.3σ less than the
analytic column which includes the factor-3 window freedom.

**How much does the number of unblinded regions matter? Almost not at all.** The assumption
above is "unblind everything over Z_cut", so k is Poisson with mean `N·p1(Z_cut)` and the B bar
uses the actual k. In the equivalent "unblind the top-k A-windows" formulation, k enters twice,
both logarithmically and with opposite signs — a larger k *raises* the B bar
`√(25+2·ln(3k))` but *lowers* the effective A-selection cut `z*` (where `N·p1(z*) = k`) — and
after re-optimizing f the two nearly cancel:

| k unblinded | eff. Z_cut | Z_B bar | best f | reach | vs single-stage |
|--:|--:|--:|--:|--:|--:|
| 1 | 3.61 | 5.22 | 0.37 | 7.09 | +0.56 |
| 10 | 2.96 | 5.64 | 0.28 | 7.07 | +0.54 |
| 100 | 2.17 | 6.03 | 0.19 | 7.02 | +0.50 |
| 1000 | 1.03 | 6.40 | 0.10 | 6.94 | +0.41 |

Three orders of magnitude in k move the reach by ~0.15σ (and as k → N the procedure smoothly
degenerates into the single-stage scan: cut → −∞, f → 0, bar → √(25+2 ln 3N)). So k should be
chosen on **practical** grounds, not statistical ones: ~10 regions (Z_cut ≈ 3) are few enough to
scrutinise one by one, and a *very* small k (tight Z_cut) is actually riskier because of the
selection-efficiency plateau discussed below.

The intuition for the overall cost: the LEE penalty enters as `√(25+2 ln N)` — going from N ≈ 10
to N ≈ 6,400 costs less than 1σ — while halving luminosity costs a factor √2 ≈ 1.41 on Z, i.e.
~2.3σ at the 5σ working point. **You cannot buy back √2 of luminosity with a logarithm.** An asymmetric split
(A just large enough to flag candidates, B keeping the luminosity for the expensive
confirmation) reduces the damage from +1.4σ to ~+0.5σ; it cannot eliminate it, because the
selection stage always burns some luminosity and some efficiency. The optimum is broad in both
f and Z_cut (see `ab_split_reach.png`; the misleading median arithmetic is the third column of
the table `ab_split_budget.py` prints). One more subtlety visible in the toy power curves: with small f the *plateau* of
the power curve sits below 1 (even a huge signal fails the A-selection with probability
1−Φ(√f·Z_full−Z_cut)) — another reason to keep Z_cut moderate (2.5–3) rather than tight.

## 2b. How many trials until the split actually *wins*? — the trials-factorization identity

Natural question: the LEE penalty grows with N while the B trials factor stays ~ln k, so surely
at *some* number of looks the split becomes a net sensitivity gain? **No — at no finite N**, and
the reason is an exact identity: **the LEE is conserved under splitting; it factorizes, it does
not shrink.**

**Derivation.** Write all three thresholds with the same Gaussian tail `p1(z) ≈ φ(z)/z` and take
logs (α = p1(5) = 2.87·10⁻⁷; w ≈ 3 is the per-window ±2σ_M freedom in B; k = number of
unblinded windows; z* = effective A-cut, defined by N·p1(z*) = k):

    A-selection:      z*²/2  + ln(z*·√2π)  = ln(N/k)
    B-confirmation:   z_B²/2 + ln(z_B·√2π) = ln(w·k/α),   z_B ≈ √(25 + 2·ln(w·k))
    single-stage:     μ₁²/2  + ln(μ₁·√2π)  = ln(N/α),     μ₁  ≈ √(25 + 2·ln N)

Add the first two, subtract the third — **k cancels exactly** because `N = (N/k) · k`:

    μ₂,med² − μ₁²  =  2·ln w  −  2·ln( z*·z_B·√2π / μ₁ )

where `μ₂,med² = z*² + z_B²` is the median-arithmetic two-stage reach (each stage exactly at its
bar, split balanced at f\* = z\*²/(z\*²+z_B²)). Interpretation: of the total trials cost `2·ln N`,
the split pays `2·ln(N/k)` **implicitly, as the A-selection cut**, and `2·ln k` explicitly in the
B correction. You never escape the logarithm — you only choose where to pay it. What survives on
the right-hand side is O(1) and N-independent: `+2·ln w ≈ +2.2` (window widening, ~+0.17σ) minus
`2·ln(z*·z_B·√2π/μ₁) ≈ +3.7` in favour of the split (selecting at z\* ≈ 3 is probabilistically
cheaper per look than claiming at ~6.5, because the Gaussian density is flatter there). Net: the
**median** two-stage reach is a wash with the single stage — actually ~0.3σ *better* — at every N.

**So why does the split still lose?** The whole real cost is the **two-coin penalty**: discovery
requires two *independent* half-datasets to succeed simultaneously,
`P = Φ(√f·Z − z*)·Φ(√(1−f)·Z − z_B)`. At 50% joint power each stage must run at ~71%, i.e. each
quadrature leg carries a buffer of ~Φ⁻¹(0.71) ≈ 0.55:
`μ₂² ≈ (z* + a)² + (z_B + b)²` with `Φ(a)·Φ(b) = ½`. The single-stage scan pays no such penalty
— its one measurement *is* its claim. This term is positive at every N, so **there is no
crossover**. Numerically (`ab_split_budget.py`, both f and Z_cut re-optimized at each N;
figure `results/plots/ab_split_crossover.png`):

| N (looks) | single-stage μ₁ | median 2-stage (w=1) | joint power (w=1) | joint power (w=3) | cost | R* |
|--:|--:|--:|--:|--:|--:|--:|
| 10 | 5.44 | 5.28 | 5.73 | 5.92 | +0.48 | 15 |
| 10³ | 6.23 | 5.93 | 6.45 | 6.62 | +0.39 | 12 |
| ~7.0·10³ (this program) | 6.54 | — | 6.86 | 7.02 | **+0.49\*** | **13** |
| 10⁶ | 7.25 | 6.93 | 7.46 | 7.61 | +0.36 | 14 |
| 10¹⁰ | 8.43 | 8.11 | 8.63 | 8.76 | +0.33 | 17 |

\*the program row uses the Z_cut-grid optimum of section 2; the decade rows use a finer joint
(f, Z_cut) optimization. The cost drifts down by only ~0.15σ over *nine decades* of N and stays
clearly positive — extrapolating, a crossover, if it exists at all, would sit beyond N ~ 10²⁰⁺,
i.e. nowhere physics will ever be.

**When the split IS the more sensitive procedure.** The comparison above assumes the single-stage
trials factor is *exactly known* (perfect Gross–Vitells). The split wins precisely when it isn't.
Equate `25 + 2·ln N_def = μ₂²`: the split out-performs the corrected scan iff the trials count
you would otherwise have to defend exceeds

    N_def  >  N_equiv = exp((μ₂² − 25)/2)  =  R*·N_true,    R* = exp((μ₂² − μ₁²)/2) ≈ 13

at this budget (R\* ≈ 12–17 across all N; **14 on the combinatorial scan**, and 5 if the
pre-registered windows are pinned rather than given ±2σ_M). So the operational criterion is: **if
honest accounting of analyst/ML freedom — NN trainings, selection variants, hyperparameter scans,
binning choices — would inflate the defendable effective trials factor by more than ~14× the
counted N, the split is not just more auditable but genuinely more sensitive.**

**And that inflation has now been measured, not imagined.** For a fixed, classical scan-window
program it is implausible and the corrected single stage wins by ~0.5σ. For a network-driven hunt it
is not a matter of opinion: the published BumpNet application mis-estimates **10⁻⁵–10⁻⁴ of its
looks**, i.e. **78–310×** the Gaussian tail — 6 to 22 times past the break-even, ~7× at the round
10² (`ESTIMATOR_DEFECTS.md`). For an estimator of that class the split therefore *wins on
sensitivity*, before any argument about auditability. Suppressing defects upstream is the way to
overturn that, and it has to buy back a factor ~7 to do so.

The split's own exposure to the same inflation is only logarithmic. If the mis-estimation persists
at the selection cut, the pre-registered list grows from 490 to 4.9·10⁴ entries and the claim bar
from 6.29 to 6.99 — 0.7σ for a hundredfold — and even that assumes the list length is *modelled*
rather than simply read off the frozen document.

## 2bb. "Why not swap the halves?" — the symmetrised two-fold scheme

The natural objection to §2b: half the luminosity sits idle at each stage, so run **both**
directions — explore A / confirm B, **and** explore B / confirm A — and claim if *either*
confirms. This is two-fold cross-validation of the *claim* (not of the training), and it is a
genuinely different procedure, not a relabelling.

**The trade.** Both directions pre-register, so the confirmation windows double, `k → 2k`, i.e.
`2·ln 2` on the B bar (Z_cut = 3: `k_eff` 9.9 → 18.9, `Z_B` 5.64 → **5.75**). What is bought is a
second chance at the same signal. Since `Z_B_req > Z_cut` the two directions overlap only where
*both* halves clear the claim bar, so the union power is

    P_sym = Φ_cut(A)·Φ_req(B) + Φ_cut(B)·Φ_req(A) − Φ_req(A)·Φ_req(B),
    Φ_x(A) = Φ(√f·Z_full − Z_x),   Φ_x(B) = Φ(√(1−f)·Z_full − Z_x)

a factor `2 − Φ_req/Φ_cut` over one direction (at most ×2). Toys validate it to 0.01σ: at 50/50,
Z_cut = 3, toys give **7.11** against **7.12** analytic (`ab_split_toys.py`, exact-window world).

**The catch: the gain needs f = ½, which is the split §2 rejects.** The second direction only pays
when both directions have comparable power. But confirmation is the luminosity-hungry stage, so
the one-way optimum is strongly asymmetric — and there the second direction is dead.

| design (N = 7,030) | one-way | swapped | swap gain |
|---|--:|--:|--:|
| 50/50, Z_cut = 3 | 7.99 (+1.46) | **7.40 (+0.86)** | **+0.60** |
| Z_cut = 2, f re-optimised (0.18) | 7.02 (+0.49) | 7.12 (+0.59) | −0.10 |
| Z_cut = 3, f re-optimised (0.29) | 7.09 (+0.56) | 7.15 (+0.63) | −0.07 |
| Z_cut = 4, f re-optimised | 7.36 (+0.83) | 7.06 (+0.54) | +0.30 |
| **best design**, Z_cut ∈ [2, 4.5] | **7.02 (+0.49)** @ f = 0.18 | 7.04 (+0.50) @ f = 0.50 | −0.02 |

So the swap does not beat the asymmetric scheme — it **rescues the naive 50/50 one** (+1.45σ →
+0.86σ), moving the optimum back to f = ½. Layered on top of the *recommended* working point
(f = 0.30, Z_cut = 3) it is actively harmful: direction 2 must confirm a 5σ-global result on 30%
of the data and essentially never fires, so one pays the `2·ln 2` for a dead direction. At
best-design level the two are a **wash** (+0.49 vs +0.51), and the gap stays within **0.11σ over
nine decades of N**, changing sign with the design box — there is no clean crossover in N, only a
swing in which (f, Z_cut) the optimiser prefers.

**Why no fold count can win — the geometric statement of §2b.** In the (z_A, z_B) plane the
Neyman–Pearson acceptance region for a mean shift is the **half-plane**
`√f·z_A + √(1−f)·z_B > √(25 + 2·ln N)` — which *is* the single-stage scan on the recombined
dataset. One-way splitting approximates it with an **L-shaped corner**; the swap with a **two-step
staircase**. More steps track the line more closely (hence the 50/50 rescue), but every step
multiplies the pre-registered windows. Conservation of the LEE is that staircase-vs-half-plane
inefficiency — a restatement of the `2·ln N = 2·ln(N/k) + 2·ln k` identity, and immune to adding
folds.

**The decisive objection is not statistical.** Of the three things the split buys (§6), the swap
keeps two and destroys the third:

- **countable trials** — kept, `2k` is as countable as `k`;
- **auditability** — kept, still the number of lines in a frozen document;
- **procedural blindness** — **gone**. After direction 1 the analyst has seen B, so direction 2's
  "confirmation" on A is made by someone who already knows both halves. Statistically valid *only*
  if the whole two-directional procedure is frozen and executed in one automated pass, and a
  review committee is entitled to discount it even then.

It also does nothing for the correlated-systematics mode of §3, which for a NN-driven hunt is the
leading one. **Not to be confused with cross-half *trainings*** (§3, §4.4: train on A → predict B
and vice versa) — symmetrising the *background model* is free and attacks the sculpting mode;
symmetrising the *claim* is what costs. **Verdict: not recommended** at this budget. Adopt it only
if a 50/50 split is forced for some external reason, where it is worth +0.6σ.

## 2c. Working point: "globally > 3σ in B" — the evidence trigger

Design requirement: **multiple regions are flagged in A** (everything above Z_cut, as in §1),
all of them are unblinded in B, and the B result is quoted as a **global significance over the
pre-registered windows**. The claim we want out of B is `Z_B,global ≥ 3`. The conversion is the
countable-trials relation of §1 with the 5 replaced by the global target:

    Z_B,global ≥ Z_g   ⇔   Z_B,local ≥ √(Z_g² + 2·ln(3·k_eff))

| Z_cut (A) | k_bkg | B local bar for **3σ global** | best f | reach (50% power) | single-stage 3σ-global |
|--:|--:|--:|--:|--:|--:|
| 2.5 | 41.1 | 4.32 | 0.32 | 5.71 | 5.16 |
| **3.0** | **8.9** | **3.97** | **0.40** | **5.73** | 5.16 |
| 3.5 | 1.5 | 3.61 | 0.49 | 5.80 | 5.16 |

Reading: with the standard A selection (Z_cut = 3, ≈9 flagged windows, ±2σ_M freedom → k_eff·3 ≈
30 counted looks in B), a **local Z_B ≈ 4.0 in any unblinded window is a 3σ-global result** —
the number to pre-register. A signal of full-dataset strength Z_full ≈ 5.7 crosses this with 50%
probability at the optimal split (vs 5.16 for a 3σ-global single-stage scan: the same ~+0.6σ
two-stage price as at the discovery bar; note the evidence bar prefers a slightly larger f ≈ 0.4
than the discovery bar, because the B requirement is cheaper). By construction the **false-evidence rate per B
opening is p₁(3) = 1.35·10⁻³ under background-only**, independent of how many windows A happened
to flag — that is the whole point of the countable correction.

The full ladder at Z_cut = 3, written down at pre-registration time:

| B outcome (local Z in any pre-registered window) | global meaning |
|---|---|
| Z_B < 3.97 | background-compatible; window closed |
| **Z_B ≥ 3.97** | **Z_global ≥ 3 — evidence**; triggers full-luminosity follow-up / next data |
| Z_B ≥ 5.64 | Z_global ≥ 5 — discovery (§2) |

## 2d. Flagging in A with Benjamini–Hochberg (FDR) instead of a fixed Z_cut

Alternative A rule: compute the N one-sided p-values of the full A scan, order them
`p_(1) ≤ … ≤ p_(N)`, find the largest j with `p_(j) ≤ q·j/N`, and flag the top j windows —
BH at FDR level q. What changes (all toy-validated, `ab_split_budget.py` BH section):

**1. Under background-only, B almost never opens.** Under the global null FDR = FWER = q, so
**P(BH flags anything at all) = q** (toys: 0.101 at q = 0.10, and when it does flag, k_obs ≈
1.3). This is the operationally biggest difference from the fixed cut: Z_cut = 3 *guarantees*
~9 background windows and therefore always spends the one-shot B unblinding; BH keeps B closed
with probability 1−q, preserving it for the growing dataset. The flag list also carries a
pre-registerable guarantee ("expected fraction of false flags ≤ q").

**2. For a single isolated signal, BH degenerates to Bonferroni** — with no other small
p-values, the signal is flagged iff `p ≤ q/N`, i.e. `z_A ≥ Φ⁻¹(1−q/N)` ≈ 4.0–4.7 — a *harder*
A bar than Z_cut = 3 (toys: P(flag) = 0.69 vs 0.68 Bonferroni-analytic at the reach point). The
B side gets cheaper (k_obs ~ 1 → evidence bar 3.4 instead of 4.0), but at the *evidence* working
point the exchange is not neutral — looser A cuts are the better end of the trade there:

| A rule | A bar (1 signal) | P(open B \| bkg) | B bar (3σ global) | reach (50%) | vs fixed cut |
|---|--:|--:|--:|--:|--:|
| fixed Z_cut = 3 | 3.00 | ≈ 1 (always, k≈9) | 3.97 | **5.73** | — |
| BH q = 0.25 | 3.95 | 0.25 | 3.41 | 5.99 | +0.26 |
| BH q = 0.10 | 4.16 | 0.10 | 3.37 | 6.12 | +0.39 |
| BH q = 0.05 | 4.32 | 0.05 | 3.36 | 6.23 | +0.50 |

**3. BH's adaptivity pays exactly where this program has structure.** One BSM model typically
lights up *several* of the 78 selection-level spectra at once (Z′ → ee *and* μμ; W′ → ℓν + tb;
VLQ multi-channel). With c channels above the bar, the step-up threshold relaxes to
`z(q·c/N)` — at q = 0.10: 4.16 (c=1) → 4.00 (c=2) → 3.91 (c=3) → 3.78 (c=5). BH flags
coherent multi-channel signals that a per-window Z_cut treats as independent 3σ-ish wiggles;
it is the natural rule when the alternative is "one model, many windows".

**Verdict.** BH is not a free upgrade — for a lone bump it costs ~0.3–0.5σ of evidence reach
because its single-signal bar is Bonferroni-like. It buys three real things: **B stays blind
under background** (false B-opening rate = q instead of certainty — the total false-evidence
rate becomes q·p₁(3) ≈ 10⁻⁴ at q = 0.10), an **FDR-certified flag list**, and **multi-channel
adaptivity**. A sensible hybrid keeps both: pre-register the union — the fixed-Z_cut = 3 list
(sensitivity for lone bumps, ~9 windows) *plus* whatever BH(q = 0.1) adds through multi-channel
step-up — since §2b guarantees that a few extra countable windows cost essentially nothing in
the B correction (`2·ln(3·k)` moves by ~0.1 in Z per doubling of k).

## 2e. The same design on the combinatorial scan

The derivations above are worked on N = 7,030 because that is where the toys were run. The scan a
network-driven hunt would face is the lensed combinatorial one, N = 362,815, and the design carries
over unchanged (`ab_split_scan.csv`):

| | model space (7,030) | scan + lenses (362,815) |
|---|--:|--:|
| single stage, exactly corrected | 6.54 | **7.11** |
| optimised split (f ≈ 0.2–0.3, Z_cut = 2–3) | 7.02 (+0.49) | **7.58–7.65 (+0.47)** |
| naive 50/50, Z_cut = 3 | 7.99 (+1.46) | 8.90 (+1.79) |
| windows pre-registered at Z_cut = 3 | 9 | **490** |
| claim bar √(25 + 2·ln(3k)) | 5.64 | **6.29** |
| break-even R\* (w = 3 / w = 1) | 13 / 4 | **14 / 5** |

Fifty times more looks moves the cost of the split by 0.02σ, which is the point: the split is priced
in `ln k`, so it hardly cares how wide the scan behind it is. The background-only check scales the
same way — in 2·10⁴ toys at Z_cut = 3 the best confirmation over the 490 windows reaches **5.30**
against the **6.29** bar, so **zero false claims**, the same verdict as on the narrower basis.

## 3. What the split does NOT protect against

- **Correlated systematics.** A background-model artefact — a turn-on, a trigger edge, and in
  the DDP context a *NN that has learnt to sculpt a bump* — appears identically in A and B.
  The confirmation stage validates the *statistical* nature of a fluctuation, nothing more.
  Mitigation, which marries naturally with the A/B machinery a data-directed scan already
  needs: **cross the halves and the trainings** — train the background
  model on A to predict B and vice versa. Then a bump that appears in both halves under
  *independent* trainings cannot be a training artefact of either. (A physics bump survives; a
  sculpting artefact would need both independent trainings to sculpt at the same mass — still
  possible if the *method* biases a particular mass, so keep the signal-injection closure tests.)
- **The one-shot property.** B can be opened once. If a tantalising 2.8σ sits just below Z_cut
  in A, there is no third dataset. Freeze Z_cut (and everything else) *before* looking at A, and
  resist re-tuning after A is seen — otherwise the countable trials factor silently grows by the
  number of choices you allowed yourself.
- **Split-induced pathologies.** Split at the *event* level, randomly (not by run period), so
  detector/pileup conditions are statistically identical in both halves; keep an event in exactly
  one half; verify with a closure test that A- and B-spectra agree in sidebands.

## 4. Recommended design (if adopted)

1. **Asymmetric split, f ≈ 0.25** (A) / 0.75 (B), event-level random, frozen before any data
   look. 50/50 only if A+B conditional combination is planned instead of B-only claims.
2. **Z_cut = 3.0** in A → expect ~9 background selections (scales with N: use the level of
   granularity actually scanned; at the inclusive 2,513-look level it is ~3.4).
3. Pre-register each selected window: spectrum, event selection, mass window ±2σ_M, test
   statistic, background model version. The B trials factor is 2·ln(3·k_eff) — write it down at
   pre-registration time, before unblinding B, together with the outcome ladder of §2c:
   local Z_B ≥ 3.97 ⇔ 3σ global (evidence), Z_B ≥ 5.64 ⇔ 5σ global (discovery).
4. **Cross-trained backgrounds**: NN trained on A predicts B and vice versa; a confirmed bump
   must appear in B under the A-trained model.
5. B-only p-value for discovery; A+B for measurement. Additionally unblind 2–3 random
   *control* windows in B (selected without reference to A) as a bias check of the procedure.
6. Report the stage-A excess population against the `EXCESS_COUNTING.md` expectation
   (k_obs vs k_bkg) — an anomalous *number* of selections is itself a (weak) signal, and a
   normal number is a powerful background-only validation.

## 5. Toy validation (`searchbudget/stages/ab_split_toys.py`)

The whole search program is modelled as its N = 7,030 effective independent looks (iid standard
normals under background-only), split per look into independent halves
`z_A = √f·μ + g_A`, `z_B = √(1−f)·μ + g_B`, so the full-dataset scan is recovered exactly as
`z_full = √f·z_A + √(1−f)·z_B`. 20,000 background-only toys and 1,500 signal toys per strength
point. Three figures in `results/plots/`:

- **`ab_toys_background.png`** — background-only, three panels:
  (a) the number of stage-A selections is Poisson with mean `N·p1(Z_cut)` — toys give
  **8.89 vs 8.93 predicted**; (b) the best confirmation significance among the unblinded
  windows never approaches the claim threshold — max Z_B = **4.69 in 20,000 toys** vs a ~5.4
  bar, i.e. **zero false claims**; (c) the corrected global p-value (`k·p1(Z_B_best)`) is
  uniform-to-conservative under background — the countable trials factor is *correct*, not just
  convenient (the mild bulk conservatism is ordinary Bonferroni; the small-p tail that matters
  for discovery hugs the diagonal).
- **`ab_toys_power.png`** — discovery power vs injected `Z_full` for the three procedures on
  the *same* toy data. The 50%-power crossings — single-stage **6.52**, optimized split
  **7.14**, 50/50 **7.72** — match the Φ-product analytics (6.54 / 7.10 / 7.72) to 0.04σ.
  This plot is also what killed the first draft's "median arithmetic" claim that the optimized
  split beats the single-stage scan. Note the sub-unity plateau of the two-stage curves at
  small f: even an arbitrarily strong signal is lost when its A-half fluctuates below Z_cut.
- **`ab_toys_spectrum.png`** — one toy falling mass spectrum, the procedure end-to-end: a
  background fluctuation at Z_A = 3.3 is pre-registered, B is opened *only there* and shows
  0.6σ — dead; an injected `Z_full = 7` signal at 1.2 TeV passes A at 3.0σ and confirms in B at
  **6.3σ**, above the √(25+2·ln k) ≈ 5.1 bar, while a second (background) window selected in
  the same toy dies quietly.

## 6. Bottom line

The two-stage scheme is worth adopting for its **auditability** (exact, countable trials; no
LEE modelling debate) and its **procedural blindness**, *not* for sensitivity — the toys and
the Φ-product analytics agree that at the optimal asymmetric split (f ≈ 0.2–0.3) it costs
**~0.5σ** of discovery reach, and at 50/50 **~1.4σ**. Running both directions and claiming on the
union (§2bb) does not buy that back — it repairs the 50/50 split but is a wash against the
asymmetric optimum, and it forfeits procedural blindness. That half-sigma is the price of turning
"trust our Gross–Vitells toys over 6,000 spectra" into "count the lines of a frozen document".
Its blind spot is correlated systematics, which for a NN-based bump hunt is exactly the
dangerous axis — so it should be deployed *together with* cross-half trainings, signal
injection, and the excess-counting sanity check, all of which exist in this project already.

**One qualification to "not for sensitivity".** That verdict holds where the single-stage trials
factor is honestly known. Where the estimator is a network it is not, and the measured defect rate
of a published one puts the required conservatism 6–22× past the break-even — so on the
combinatorial scan the split is *also* the more sensitive procedure, and the half-sigma is not a
price but a saving (`ESTIMATOR_DEFECTS.md`).

Sources: reach analytics `searchbudget/stages/ab_split_budget.py` (+ plot
`results/plots/ab_split_reach.png`); toys `searchbudget/stages/ab_split_toys.py`; trials budget
`SEARCH_BUDGET.md` / `SEARCH_BUDGET_SELECTIONS.md`; expected 3σ population
`EXCESS_COUNTING.md`; Gross & Vitells, arXiv:1005.1891 (single-stage LEE); the "selective
inference / data splitting" literature (e.g. Cox 1975) for why B-only inference after
A-selection is unbiased.
