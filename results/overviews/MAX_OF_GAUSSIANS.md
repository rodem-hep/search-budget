# The distribution of the maximum of $n$ Gaussian draws

**Setup.** Each pseudo-experiment draws $n = 30{,}000$ i.i.d. values $X_1,\dots,X_n \sim \mathcal{N}(0,1)$
and records only the largest, $M = \max_i X_i$. We repeat this $N = 10{,}000$ times and histogram the
resulting maxima. What is that histogram converging to?

**Answer in one line.** It converges to the *exact* density $f_M(x) = n\,\phi(x)\,\Phi(x)^{n-1}$, which
is asymptotically (in $n$, not in $N$) a **Gumbel** distribution with location $a_n = 4.005$ and scale
$b_n = 0.220$. The maxima pile up at $\approx 4.1$ with a spread of $\approx 0.28$, and the distribution
is strongly **right-skewed** — not Gaussian.

Throughout, $\phi(x) = \tfrac{1}{\sqrt{2\pi}}e^{-x^2/2}$ is the standard normal density and
$\Phi(x) = \int_{-\infty}^{x}\phi(t)\,dt$ its CDF.

**Contents.** *Part I* (§1–6) is background-only: the exact law of the maximum and its Gumbel limit.
*Part II* (§7–9) puts a signal in one bin and asks whether it can win the scan, and what a two-stage
A/B confirmation buys you. *Part III* (§10–12) replaces "take the maximum" with "take everything above
a threshold $t$", scans $t \in [3,6]$, and puts the two rules on a common ROC — power against average
number of false positives. The argmax turns out to be a **single dominated point** inside that curve:
ranking bins throws away information that their absolute size carries.

---

# Part I — Background only

## 1. Two limits, kept apart

This problem has two independent knobs, and almost all confusion comes from conflating them.

| | Limit | What it gives you | Exact? |
|---|---|---|---|
| **$N \to \infty$** | many pseudo-experiments | the histogram converges to $f_M(x) = n\phi\Phi^{n-1}$ | **exact for every $n$** |
| **$n \to \infty$** | many draws per pseudo-experiment | the rescaled $f_M$ converges to Gumbel | asymptotic, error $O(1/\ln n)$ |

The first limit is just the law of large numbers and involves no approximation at all. The second is
extreme value theory, and it is the one that is imperfect at $n = 30{,}000$.

**Consequence:** running more pseudo-experiments gives you a cleaner estimate of the blue curve below.
It will never turn the blue curve into the orange one. The gap between them is set by the finite $n$
*inside* each pseudo-experiment, and no amount of Monte Carlo statistics touches it.

---

## 2. The exact distribution (first principles)

### 2.1 The CDF, from independence

The maximum is below $x$ **if and only if every single draw is below $x$**. That is the whole
derivation. By independence the joint probability factorises:

$$F_M(x) \;=\; P(M \le x) \;=\; P(X_1 \le x,\; \dots,\; X_n \le x) \;=\; \prod_{i=1}^{n} P(X_i \le x) \;=\; \boxed{\;\Phi(x)^n\;}$$

### 2.2 The density, two ways

**By differentiation.** Apply the chain rule to $F_M$:

$$f_M(x) \;=\; \frac{d}{dx}\,\Phi(x)^n \;=\; n\,\Phi(x)^{n-1}\,\Phi'(x) \;=\; \boxed{\;n\,\phi(x)\,\Phi(x)^{n-1}\;}$$

**By direct counting** (no calculus needed, and more illuminating). Ask for the probability that the
maximum lands in $[x, x+dx]$. For that to happen, *exactly one* draw must sit at $x$ and *all the
others* must fall below it:

$$f_M(x)\,dx \;=\; \underbrace{n}_{\substack{\text{which of the } n \\ \text{draws is the max}}} \;\times\; \underbrace{\phi(x)\,dx}_{\substack{\text{that draw} \\ \text{lands at } x}} \;\times\; \underbrace{\Phi(x)^{n-1}}_{\substack{\text{the other } n-1 \\ \text{all land below } x}}$$

Both routes agree. This result is **exact for all $n$** — no asymptotics have been used.

### 2.3 Why it is skewed

The two factors fight each other. $\Phi(x)^{n-1}$ rises from 0 to 1 as $x$ grows, while $\phi(x)$ falls
off as $e^{-x^2/2}$. Their product is sharply peaked, but *asymmetrically*: on the left, $\Phi(x)^{n-1}$
collapses **extremely** fast (it takes all $30{,}000$ draws simultaneously coming out low, a
catastrophically unlikely conspiracy), while on the right the decay is merely the single-Gaussian tail
$\phi(x)$. A steep wall on the left, a long tail on the right — hence positive skew.

### 2.4 Exact quantiles (useful in practice)

Inverting $F_M(x) = \Phi(x)^n = q$ costs nothing:

$$\boxed{\;x_q \;=\; \Phi^{-1}\!\left(q^{1/n}\right)\;}$$

This is the formula to use for any real threshold or $p$-value. Example: the median of $M$ is
$\Phi^{-1}(2^{-1/30000}) = 4.086$. There is no reason to use an approximation for numbers like this.

---

## 3. The Gumbel limit (first principles)

Now the second limit: what *shape* does $\Phi(x)^n$ take as $n \to \infty$, once we recentre and
rescale so it doesn't just run off to infinity?

### 3.1 The tail-exponentiation trick

Write the upper tail as $\epsilon(x) \equiv 1 - \Phi(x)$, which is small out where the maximum lives.
Then

$$F_M(x) \;=\; \bigl(1 - \epsilon(x)\bigr)^{n} \;=\; \exp\Bigl[\,n\ln\bigl(1-\epsilon(x)\bigr)\Bigr] \;\approx\; \exp\bigl[-\,n\,\epsilon(x)\bigr]$$

using $\ln(1-\epsilon) \approx -\epsilon$. **Everything now hinges on the shape of the Gaussian tail
$\epsilon(x)$ — the rest is bookkeeping.**

### 3.2 Choosing the centring $a_n$

Define $a_n$ as the point where you expect **exactly one** draw to exceed it:

$$n\,\bigl(1 - \Phi(a_n)\bigr) \;=\; 1 \qquad\Longleftrightarrow\qquad \epsilon(a_n) = \tfrac{1}{n}$$

This is the natural centre of the distribution of the max: below it you expect more than one
upcrossing, above it fewer than one.

### 3.3 The tail is locally exponential

Probe near $a_n$ by setting $x = a_n + b_n z$ with $b_n = 1/a_n$ (the reason for this choice appears in
a moment). Use the standard Gaussian tail asymptotic — itself derivable by one integration by parts,
$\int_x^\infty \phi(t)dt = \phi(x)/x - \int_x^\infty \phi(t)/t^2\,dt$ — namely

$$1 - \Phi(x) \;\simeq\; \frac{\phi(x)}{x} \qquad (x \to \infty).$$

Expand the numerator. The exponent of $\phi$ at $x = a_n + z/a_n$ is

$$-\tfrac{1}{2}\left(a_n + \frac{z}{a_n}\right)^{2} \;=\; -\frac{a_n^2}{2} \;-\; z \;-\; \frac{z^2}{2a_n^2}$$

so that

$$\phi\!\left(a_n + \tfrac{z}{a_n}\right) \;=\; \phi(a_n)\, e^{-z}\, e^{-z^2/2a_n^2} \;\approx\; \phi(a_n)\,e^{-z}\,\bigl[1 + O(a_n^{-2})\bigr].$$

The denominator contributes only $1/x = \tfrac{1}{a_n}\bigl(1 + O(a_n^{-2})\bigr)$. Therefore

$$\epsilon\!\left(a_n + \tfrac{z}{a_n}\right) \;\approx\; \frac{\phi(a_n)}{a_n}\,e^{-z} \;\approx\; \epsilon(a_n)\, e^{-z} \;=\; \frac{e^{-z}}{n}.$$

**This is the crux:** measured in units of $b_n = 1/a_n$, the Gaussian tail decays *exponentially*. The
Gaussian's $e^{-x^2/2}$ is locally indistinguishable from a pure exponential once you zoom in around
$a_n$, because the quadratic term $z^2/2a_n^2$ is suppressed. And $b_n = 1/a_n$ was chosen for exactly
this reason: it is the scale on which the log-tail slope, $-d\ln\epsilon/dx \simeq x$, equals 1 at
$x = a_n$.

### 3.4 The limit law

Substitute back into §3.1:

$$P\!\left(\frac{M - a_n}{b_n} \le z\right) \;=\; \Phi(a_n + b_n z)^n \;\approx\; \exp\!\left[-\,n \cdot \frac{e^{-z}}{n}\right] \;=\; \boxed{\;\exp\bigl(-e^{-z}\bigr)\;}$$

the **standard Gumbel CDF**. The double exponential is the fingerprint of an exponentially decaying
tail. (This is one of only three possible limits — Gumbel, Fréchet, Weibull — by the Fisher–Tippett–
Gnedenko theorem; the Gaussian, with its exponential-family tail, flows to Gumbel.)

### 3.5 Solving for $a_n$ and $b_n$ explicitly

Take $n\,\phi(a_n)/a_n = 1$ and take logs:

$$\ln n \;-\; \tfrac{1}{2}\ln(2\pi) \;-\; \ln a_n \;-\; \frac{a_n^2}{2} \;=\; 0
\qquad\Longrightarrow\qquad
a_n^2 \;=\; 2\ln n \;-\; 2\ln a_n \;-\; \ln(2\pi)$$

Iterate. To leading order $a_n \approx \sqrt{2\ln n}$, so $\ln a_n \approx \tfrac12 \ln(2\ln n)$.
Feed that back in:

$$a_n^2 \;\approx\; 2\ln n \;-\; \ln(2\ln n) \;-\; \ln(2\pi) \;=\; 2\ln n \;-\; \ln\ln n \;-\; \ln(4\pi)$$

Take the square root and expand $\sqrt{1-u} \approx 1 - u/2$:

$$\boxed{\;a_n \;=\; \sqrt{2\ln n} \;-\; \frac{\ln\ln n \;+\; \ln 4\pi}{2\sqrt{2\ln n}}\;, \qquad b_n \;=\; \frac{1}{a_n} \;\simeq\; \frac{1}{\sqrt{2\ln n}}\;}$$

**For $n = 30{,}000$:** $\ln n = 10.309$, $\sqrt{2\ln n} = 4.5407$, so

$$a_n \;=\; 4.5407 \;-\; \frac{2.3331 + 2.5310}{9.0814} \;=\; 4.5407 - 0.5356 \;=\; \mathbf{4.005}, \qquad b_n \;=\; \mathbf{0.220}$$

Note the very slow growth: $a_n \sim \sqrt{2\ln n}$ means the maximum creeps upward only logarithmically,
and $b_n \sim 1/\sqrt{2\ln n}$ means the *spread shrinks* as $n$ grows. The max of many draws is a
remarkably stable quantity.

---

## 4. Moments

### 4.1 The Gumbel is minus the log of an exponential

The cleanest derivation. Let $E \sim \text{Exp}(1)$ and set $Z = -\ln E$. Then

$$P(Z \le z) \;=\; P(-\ln E \le z) \;=\; P\bigl(E \ge e^{-z}\bigr) \;=\; \exp\bigl(-e^{-z}\bigr)$$

— exactly standard Gumbel. So a Gumbel variate *is* a log-transformed exponential, and its moment
generating function follows from the Gamma integral:

$$\mathbb{E}\bigl[e^{tZ}\bigr] \;=\; \mathbb{E}\bigl[E^{-t}\bigr] \;=\; \int_0^\infty x^{-t}e^{-x}\,dx \;=\; \boxed{\;\Gamma(1-t)\;}$$

Differentiating at $t=0$ and using $\Gamma'(1) = -\gamma$, $\psi'(1) = \pi^2/6$:

$$\mathbb{E}[Z] = \gamma = 0.5772, \qquad \operatorname{Var}[Z] = \frac{\pi^2}{6}, \qquad \text{skew}[Z] = \frac{12\sqrt{6}\,\zeta(3)}{\pi^3} = 1.1395$$

### 4.2 Back to $M$

Since $M \approx a_n + b_n Z$:

$$\mathbb{E}[M] \;\approx\; a_n + \gamma\, b_n \;=\; 4.005 + 0.577(0.220) \;=\; \mathbf{4.13}$$
$$\operatorname{sd}[M] \;\approx\; \frac{\pi\, b_n}{\sqrt 6} \;=\; 1.2825 \times 0.220 \;=\; \mathbf{0.28}$$
$$\text{skew}[M] \;\approx\; \mathbf{1.14} \quad (\text{scale-invariant})$$

The skewness is the honest signature that this is not a Gaussian: it is a fixed positive constant,
independent of $n$.

**Important caveat.** The moments of the *exact* density have **no closed form** — $\int x\, n\phi(x)\Phi(x)^{n-1}dx$
is not elementary. That is precisely why the Gumbel approximation stays popular despite being imperfect:
it is the only route to a formula. If you want accurate numbers, integrate the exact density numerically;
if you want an expression on a slide, quote Gumbel.

---

## 5. Numerical validation

Simulation: $N = 10{,}000$ pseudo-experiments $\times$ $n = 30{,}000$ draws.

| | simulated | Gumbel asymptotic |
|---|---|---|
| mean | 4.110 | 4.132 |
| std | 0.284 | 0.283 |
| median | 4.072 | 4.086 |
| skewness | 0.89 | 1.14 |

Kolmogorov–Smirnov tests of the 10,000 simulated maxima against each candidate:

| Reference distribution | KS statistic | $p$-value |
|---|---|---|
| **exact** $\Phi(x)^{n}$ | 0.0077 | **0.59** ✓ |
| Gumbel$(a_n, b_n)$ | 0.041 | $2\times10^{-15}$ ✗ |
| best-fit Gaussian | 0.058 | $6\times10^{-30}$ ✗ |

The exact form is confirmed. **The asymptotic Gumbel is decisively rejected** — with 10,000
pseudo-experiments you have enough statistics to resolve its error. A best-fit Gaussian is rejected far
harder still, so never quote a symmetric $\pm$ on a maximum.

### Why Gumbel is off at $n = 30{,}000$

Look back at the neglected factor in §3.3: $e^{-z^2/2a_n^2}$. The relative error is $O(1/a_n^2) = O(1/2\ln n)$
— the convergence to Gumbel is only **logarithmically** slow. With $\ln(30000) \approx 10$ the correction is
of order $1/20 \approx 5\%$, which is exactly the discrepancy seen above (and visible in the skewness,
0.89 measured versus 1.14 asymptotic). Reaching, say, 1% would require $n \sim e^{50}$. In practice, for
any realistic $n$, **the Gaussian max is never fully Gumbel.**

---

## 6. Practical summary

- **For intuition / a formula:** Gumbel, location $a_n = 4.005$, scale $b_n = 0.220$, mean 4.13, sd 0.28,
  skew 1.14.
- **For any actual number** (thresholds, $p$-values, quantiles): use the exact
  $F_M(x) = \Phi(x)^{n}$ and $x_q = \Phi^{-1}(q^{1/n})$. It costs nothing and it is right.
- **Never** model the maximum as Gaussian.

### Connection to the look-elsewhere effect

If a bump hunt scans $N_{\rm eff}$ effectively independent bins, the largest observed local significance
$z$ has global $p$-value

$$p_{\rm global} \;=\; 1 - \Phi(z)^{N_{\rm eff}} \;\approx\; N_{\rm eff}\,\bigl(1-\Phi(z)\bigr) \;=\; N_{\rm eff} \times p_{\rm local}$$

the last step (the familiar **trials factor**) following from $(1-\epsilon)^N \approx 1 - N\epsilon$ when
$N\epsilon \ll 1$. Everything then reduces to estimating $N_{\rm eff}$ — the shape of the distribution
follows automatically from the derivation above. Note this is where the $\Phi(x)^n$ form earns its keep:
in the far tail, where global $p$-values actually live, the Gumbel approximation's few-percent error is
not something you want between you and a discovery claim.

---
---

# Part II — Can a real signal win the scan?

Everything so far was background-only. Now put a signal in.

**Setup.** One bin carries a signal of strength $\mu$: $S \sim \mathcal{N}(\mu, 1)$. The other
$n-1 = 29{,}999$ bins are pure background, $B_i \sim \mathcal{N}(0,1)$. All independent.

---

## 7. Probability the signal bin is the maximum

### 7.1 Derivation

The signal bin is the maximum iff it beats **every** background bin. Condition on the signal's realised
value $s$, then use independence exactly as in §2.1 — given $s$, the probability that all $n-1$
background bins fall below it is $\Phi(s)^{n-1}$. Integrating over $s$ against the signal's own density:

$$\boxed{\;p(\mu) \;\equiv\; P(\text{signal bin is the max}) \;=\; \int_{-\infty}^{\infty} \phi(s-\mu)\,\Phi(s)^{\,n-1}\,ds\;}$$

Note this is the same object as before, now smeared over the signal's fluctuation: the integrand is
"signal lands at $s$" $\times$ "the entire rest of the scan stays below $s$."

Sanity check: at $\mu = 0$ all $n$ bins are exchangeable, so $p(0)$ must equal $1/n = 3.3\times 10^{-5}$
— and it does.

### 7.2 Results ($n = 30{,}000$)

| $\mu$ | $p(\mu)$ = P(signal is the max) |
|---|---|
| 3σ | **14%** |
| 4σ | **46%** |
| 5σ | **80%** |
| 6σ | **96%** |

Crossing points: **50% at $\mu = 4.11\sigma$**, 90% at $5.45\sigma$, 95% at $5.83\sigma$, 99% at $6.55\sigma$.

*(Plot: [`signal_wins_the_max.png`](../plots/max_of_gaussians/signal_wins_the_max.png); verified against Monte Carlo to 4 decimal places.)*

### 7.3 Why a 4σ signal loses half the time — a width mismatch

This is the counter-intuitive part, and the reason is a **mismatch of widths**:

- The **background maximum** is a *tight* distribution: $4.09 \pm 0.28$. With 30,000 tries you are
  practically guaranteed to find something near 4 — recall from §3.5 that $b_n \sim 1/\sqrt{2\ln n}$
  actually *shrinks* with $n$.
- The **signal bin** is a *single draw* with the full unit width, $\sigma = 1$ — more than three times
  broader.

So the signal is the sloppy competitor in this race. A 4σ signal is centred slightly *below* the typical
background max (4.00 vs 4.09), and its own fluctuation is wide enough that it often lands at 3.5 while
some empty bin somewhere lands at 4.3. Hence $p(4\sigma) \approx 46\%$, and the 50% crossing sits
almost exactly at the background max's median, as it must.

**Consequence for analysis design.** If you take the global maximum and study only that, you discard a
real 5σ signal one time in five. This inefficiency is *distinct from* the trials factor — it is not
about false alarms at all — and it is the argument for examining the top-$k$ excesses rather than only
the winner.

---

## 8. Two-stage A/B confirmation

**Protocol.** Sample A and sample B are independent, with the same signal in the same bin.

1. **Stage 1** — find the argmax bin in sample A.
2. **Stage 2** — unblind *only that one bin* in sample B. Confirm if $z_B > t$, with $t = 3$.

### 8.1 The key structural fact: no trials factor in stage 2

The bin that gets unblinded is selected **using A alone**, and B is independent of A. Therefore, under
background-only, that bin's value in B is a plain $\mathcal{N}(0,1)$ draw. **You made one look, not
30,000.** The look-elsewhere effect is not "corrected for" — it is *structurally absent*, and the
stage-2 local $p$-value **is** the global $p$-value.

$$\boxed{\;P(\text{false confirmation}) \;=\; 1 - \Phi(t) \;=\; 1.35\times 10^{-3} \;=\; \mathbf{3.0\sigma\ \textbf{global}}\;}$$

For comparison: to reach that same global false-alarm rate in a *single-stage* scan of 30,000 bins you
would need a local threshold of $z^\star = \Phi^{-1}\!\bigl(1 - \tfrac{1-\Phi(3)}{n}\bigr) = 5.35$. The
A/B protocol buys the same global protection at a 3σ bar.

### 8.2 Confirmation probability

Split on who won stage 1 (the two cases are exclusive and exhaustive):

$$\boxed{\;P_{\rm confirm}(\mu) \;=\; \underbrace{p(\mu)\,\bigl[1-\Phi(t-\mu)\bigr]}_{\substack{\text{signal bin won A,} \\ \text{and repeats in B}}} \;+\; \underbrace{\bigl[1-p(\mu)\bigr]\,\bigl[1-\Phi(t)\bigr]}_{\substack{\text{a background bin won A,} \\ \text{and flukes past } t \text{ again in B}}}\;}$$

The first factor in each term is stage 1, the second is stage 2; they multiply because A and B are
independent. In the signal term the unblinded bin is $\mathcal{N}(\mu,1)$ in B; in the background term
it is $\mathcal{N}(0,1)$.

### 8.3 Results ($n = 30{,}000$, $t = 3$)

| $\mu$ | wins scan in A | passes $z_B > 3$ | **confirmed** | fake path |
|---|---|---|---|---|
| 3σ | 14% | 50.0% | **7%** | $1.2\times10^{-3}$ |
| 4σ | 46% | 84.1% | **39%** | $7.3\times10^{-4}$ |
| 5σ | 80% | 97.7% | **79%** | $2.7\times10^{-4}$ |
| 6σ | 96% | 99.9% | **96%** | $4.8\times10^{-5}$ |

Crossing points: 50% at $\mu = 4.26\sigma$, 90% at $5.48\sigma$, 95% at $5.85\sigma$.

*(Plot: [`ab_confirmation.png`](../plots/max_of_gaussians/ab_confirmation.png); verified against Monte Carlo, e.g. 0.3858 predicted vs 0.3861 simulated
at 4σ.)*

### 8.4 Stage 2 is nearly free

Compare columns 2 and 4: 46% → 39%, 80% → 79%, 96% → 96%. **The confirmation requirement costs almost
nothing.** A signal strong enough to win a 30,000-bin scan is, essentially by construction, strong
enough to clear a mere 3σ bar in a second sample of equal size — the bar sits at 3 while the signal is
centred at 4, 5 or 6. All of the inefficiency was already spent in stage 1.

**This is the central design lesson: the whole game is winning the scan.** You buy a 3σ *global*
false-alarm rate — the complete elimination of the look-elsewhere effect — for a few percent of signal
efficiency. That is an extraordinarily good trade.

(3σ is the exception that proves the rule: stage 2 halves it, 14% → 7%, because a 3σ signal in B is a
coin flip against a 3σ threshold. The protocol only looks free once $\mu$ comfortably exceeds $t$.)

---

## 9. Caveat: splitting one dataset in half is *not* free

§8 assumed two independent datasets each carrying the **same** signal strength $\mu$. If instead you
**split a single dataset in half**, each half has half the luminosity, so significance scales as
$\sqrt{L}$:

$$\mu_{\rm half} \;=\; \mu_{\rm full}/\sqrt{2}$$

A 5σ signal in the full dataset is only a 3.5σ signal in each half. Comparing at a **matched 3.0σ global
false-alarm rate** (so the single-stage analysis uses the local threshold $z^\star = 5.35$ from §8.1):

| $\mu_{\rm full}$ | $\mu_{\rm half}$ | A/B split, two halves | single-stage, full data |
|---|---|---|---|
| 4σ | 2.83 | 4.8% | **8.9%** |
| 5σ | 3.54 | 20.5% | **36.5%** |
| 6σ | 4.24 | 49.2% | **74.3%** |
| 7σ | 4.95 | 77.0% | **95.1%** |

50% power needs $\mu_{\rm full} = 6.02\sigma$ (A/B split) versus $5.35\sigma$ (single-stage); 90% power
needs $7.76\sigma$ versus $6.63\sigma$.

**At matched false-alarm rate, splitting strictly *loses* sensitivity** — by roughly $0.7\sigma$ to
$1.1\sigma$ of required signal strength. The $\sqrt{2}$ penalty applied twice (once in each half)
outweighs the trials-factor saving.

So the A/B protocol's value is **not** raw sensitivity. Its value is that the false-alarm rate is
*exact and assumption-free*: it does not require estimating $N_{\rm eff}$, does not depend on the
background-model fit surviving into the far tail, and cannot be gamed by post-hoc choices in stage 1.
Whether that robustness is worth ~1σ of reach is an analysis-design judgement, not a statistics one —
**but it is only a free lunch when A and B are genuinely separate datasets (e.g. two run periods, or a
new dataset arriving later), not two halves of the same one.**

---

# Part III — Threshold selection instead of the maximum

Everything so far selected **exactly one** bin per experiment: the argmax. Now change only the
selection rule — keep **every** bin whose $z$-score exceeds a threshold $t$, scanned over
$t \in [3, 6]$ — and leave the rest of the setup untouched ($n = 30{,}000$ bins, one of them carrying
signal $\mu$, and, where relevant, the same sample-B confirmation at $z_B > 3$).

This turns out not to be a small variation. **At its optimum the threshold rule strictly dominates the
argmax rule**, and the reason is structural.

---

## 10. Stage 1: what the cut selects

### 10.1 Background — a Poisson number of fakes

Each of the $n$ background bins independently exceeds $t$ with probability $1-\Phi(t)$, so the number
of surviving background bins is Binomial, and for these tiny probabilities, Poisson:

$$K \;\sim\; \mathrm{Binomial}\bigl(n,\,1-\Phi(t)\bigr) \;\simeq\; \mathrm{Poisson}\bigl(\lambda(t)\bigr),
\qquad \boxed{\;\lambda(t) \;=\; n\,[\,1-\Phi(t)\,]\;}$$

$\lambda(t)$ is the **expected number of fake bins** the cut hands downstream, and it is the single
number that controls everything:

| $t$ | $1-\Phi(t)$ | $\lambda(t) = n[1-\Phi(t)]$ | $P(\text{at least one fake})$ |
|---|---|---|---|
| 3.0 | $1.35\times10^{-3}$ | **40.5** | $\approx 1$ |
| 3.5 | $2.33\times10^{-4}$ | 6.98 | 0.999 |
| 4.0 | $3.17\times10^{-5}$ | 0.950 | 0.613 |
| 4.5 | $3.40\times10^{-6}$ | 0.102 | 0.097 |
| 5.0 | $2.87\times10^{-7}$ | $8.6\times10^{-3}$ | $8.6\times10^{-3}$ |
| 5.5 | $1.90\times10^{-8}$ | $5.7\times10^{-4}$ | $5.7\times10^{-4}$ |
| 6.0 | $9.87\times10^{-10}$ | $3.0\times10^{-5}$ | $3.0\times10^{-5}$ |

A cut at $t=3$ is useless on its own: it selects **forty** background bins. The whole range 3→6 spans
six orders of magnitude in fake rate, which is why the choice of $t$ is the entire game.

### 10.2 Signal efficiency

The signal bin is a single draw from $\mathcal{N}(\mu, 1)$, so it clears the bar with probability

$$\boxed{\;\varepsilon(\mu, t) \;=\; P(z_A > t) \;=\; 1 - \Phi(t-\mu) \;=\; \Phi(\mu - t)\;}$$

No integral, no competition — just a Gaussian tail. Note the immediate consequence: **each efficiency
curve crosses 50% exactly at $t = \mu$**, since $\Phi(0) = \tfrac12$.

| $t$ | $\mu=3$ | $\mu=4$ | $\mu=5$ | $\mu=6$ |
|---|---|---|---|---|
| 3.0 | 50.0% | 84.1% | 97.7% | 99.9% |
| 4.0 | 15.9% | 50.0% | 84.1% | 97.7% |
| 5.0 | 2.3% | 15.9% | 50.0% | 84.1% |
| 6.0 | 0.1% | 2.3% | 15.9% | 50.0% |

*See [`threshold_scan.png`](../plots/max_of_gaussians/threshold_scan.png).*

---

## 11. Two-stage A/B with a threshold scan

Protocol: every bin above $t$ in sample A is unblinded in sample B; any of them with $z_B > 3$ is
declared confirmed.

### 11.1 The trials factor comes back — but only linearly

This is the crucial difference from §8. There, stage 1 delivered **exactly one** bin, so stage 2 was a
single test and the false-alarm rate was exactly $1-\Phi(3)$ with *no* look-elsewhere penalty. Here
stage 1 delivers $K \sim \mathrm{Poisson}(\lambda(t))$ bins, and *each one* gets an independent shot at
faking $z_B > 3$. Thinning a Poisson process by an independent probability gives a Poisson process, so
the number of false confirmations is

$$K_{\rm fake} \;\sim\; \mathrm{Poisson}\Bigl(\lambda(t)\,[\,1-\Phi(3)\,]\Bigr)
\qquad\Longrightarrow\qquad
\boxed{\;P_{\rm fake}(t) \;=\; 1 - e^{-\lambda(t)\,[1-\Phi(3)]}\;}$$

The trials factor is back — but it is now just $\lambda(t)$, the number of bins that survived stage 1,
rather than the full $N_{\rm eff} \approx n$. Stage 1 is doing the look-elsewhere bookkeeping for you,
and $t$ is the dial that sets the price.

### 11.2 The matching threshold $t^\star$

Demand that the threshold rule have the *same* false-alarm rate as the argmax rule, which was exactly
$1-\Phi(3)$. Comparing the two expressions, this requires $\lambda(t) = 1$:

$$\boxed{\;\lambda(t^\star) = 1 \quad\Longleftrightarrow\quad t^\star \;=\; \Phi^{-1}\!\bigl(1 - \tfrac1n\bigr)
\;=\; 3.99 \quad (n = 30{,}000)\;}$$

**The matching threshold is the one that lets exactly one background bin through, on average.** That is
a satisfying result: the argmax rule *always* passes exactly one bin, so the fair comparison is the
threshold that passes one bin *in expectation*. Both then feed one candidate to sample B, and both
inherit the same $1-\Phi(3) = 1.35\times10^{-3}$ = **3.0σ global** false-alarm rate.

### 11.3 Power

Because the signal bin and the background bins are independent, and the signal bin's fate in A and in B
is independent too, the confirmation probability is just a product of two Gaussian tails:

$$\boxed{\;P_{\rm confirm}(\mu, t) \;=\; \underbrace{\Phi(\mu - t)}_{\text{clears }t\text{ in A}}\;\times\;
\underbrace{\Phi(\mu - 3)}_{z_B > 3}\;}$$

Contrast this with the argmax version from §8.2, which needed the integral $p(\mu)$ and carried an extra
term for the case where a *background* bin won the scan. Here there is no competition to lose.

### 11.4 Results — threshold at $t^\star$ vs argmax

Both columns sit at an identical **3.0σ global** false-alarm rate, so this is an apples-to-apples
comparison of power:

| $\mu$ | argmax rule (§8) | threshold at $t^\star = 3.99$ | gain |
|---|---|---|---|
| 3σ | 7.2% | **8.1%** | +0.9 pt |
| 4σ | 38.6% | **42.5%** | +3.9 pt |
| 5σ | 78.5% | **82.5%** | +4.0 pt |
| 6σ | 96.3% | **97.7%** | +1.3 pt |

**The threshold rule wins at every signal strength.** Equivalently, in terms of the signal strength
needed to reach a given power:

| power | argmax | threshold at $t^\star$ | saved |
|---|---|---|---|
| 50% | 4.26σ | **4.16σ** | 0.10σ |
| 90% | 5.48σ | **5.32σ** | 0.16σ |

The gain is real but modest — roughly $0.1$–$0.16\sigma$ of reach, peaking around $\mu \approx 4.5$σ
where the two rules disagree most. It is free: same data, same false-alarm rate, strictly better power.

### 11.5 Why the threshold wins

The argmax rule forces the signal to **beat** the largest background fluctuation. The threshold rule
only asks it to **clear a fixed bar**. Those are different demands, and the second is easier.

Concretely, at $t^\star$ the number of background bins passing stage 1 is $\mathrm{Poisson}(1)$:

$$P(K=0) = 37\%, \qquad P(K=1) = 37\%, \qquad P(K \ge 2) = 26\%$$

So the threshold rule adapts to the experiment it actually got. When the background happens to fluctuate
low ($K=0$, 37% of the time) the signal walks through unopposed. When it fluctuates high ($K \ge 2$) the
rule simply passes *several* candidates and lets sample B sort them out — the signal is not thrown away
merely because some background bin happened to be bigger. The argmax rule cannot do either: it always
passes exactly one bin, so a background upward fluctuation doesn't just add a fake, it **actively
destroys** the signal candidate.

That is the whole effect. The argmax is a *competition*; the threshold is an *exam*. Ranking throws away
information that the absolute scale carries.

### 11.6 The full trade-off in $t$

Sweeping $t$ over 3→6 moves you along a power/false-alarm curve ([`threshold_vs_argmax.png`](../plots/max_of_gaussians/threshold_vs_argmax.png)):

| $t$ | $E[\text{false conf.}]$ | $P_{\rm fake}$ | global $Z$ | $\mu{=}3$ | $\mu{=}4$ | $\mu{=}5$ | $\mu{=}6$ |
|---|---|---|---|---|---|---|---|
| 3.0 | $5.5\times10^{-2}$ | $5.3\times10^{-2}$ | 1.6σ | 25.0% | 70.8% | 95.5% | 99.7% |
| 3.5 | $9.4\times10^{-3}$ | $9.4\times10^{-3}$ | 2.4σ | 15.4% | 58.2% | 91.2% | 99.2% |
| **3.99** | $1.35\times10^{-3}$ | $1.35\times10^{-3}$ | **3.0σ** | 8.1% | 42.5% | 82.5% | 97.7% |
| 4.5 | $1.4\times10^{-4}$ | $1.4\times10^{-4}$ | 3.6σ | 3.3% | 26.0% | 67.6% | 93.2% |
| 5.0 | $1.2\times10^{-5}$ | $1.2\times10^{-5}$ | 4.2σ | 1.1% | 13.3% | 48.9% | 84.0% |
| 5.5 | $7.7\times10^{-7}$ | $7.7\times10^{-7}$ | 4.8σ | 0.3% | 5.6% | 30.2% | 69.1% |
| 6.0 | $4.0\times10^{-8}$ | $4.0\times10^{-8}$ | 5.4σ | 0.1% | 1.9% | 15.5% | 49.9% |

There is no free lunch *within* the scan — $t$ trades power against fake rate monotonically. The point
of $t^\star$ is only that it is where the comparison with the argmax rule is fair, and there the
threshold rule wins.

Two practical readings of this table:

- **$t$ also sets your unblinding cost.** At $t=3$ you must unblind ~40 bins in sample B; at $t^\star$
  you unblind about one. If unblinding is expensive or politically constrained, that matters as much as
  the statistics.
- **For a single-stage threshold scan** (no sample B at all), matching 3.0σ global requires
  $t = 5.35$ — the same number as §8.1, as it must be. Its power ($\mu=5$σ → 36.5%) is far below either
  two-stage rule, which is the look-elsewhere effect being paid in full.

---

## 12. The ROC: power versus false positives

§11 compared the two rules at one carefully chosen point ($t^\star$). The honest comparison is the whole
curve. Letting more than one bin through is fine — what we actually care about is **how much signal we
find per false positive we tolerate**, so put those two on the axes and let $t$ trace out a ROC.

### 12.1 The two axes

$$x(t) \;=\; E[K_{\rm fake}] \;=\; \underbrace{n\,[1-\Phi(t)]}_{\text{bins surviving A}}\times
\underbrace{[1-\Phi(3)]}_{\text{each fakes }z_B>3}, \qquad\qquad
y(\mu, t) \;=\; \Phi(\mu-t)\,\Phi(\mu-3)$$

$x$ is the **average number of false positives**, not a probability — it is allowed to exceed 1, and at
$t=3$ it does not, but the number of *candidates* certainly does (40 of them). Using the expected count
rather than $P(\ge 1)$ is the right choice precisely because we are no longer restricted to one bin.

### 12.2 The argmax is a single point, not a curve

The argmax rule has no dial. Under background only it *always* returns exactly one bin, and that bin is
confirmed with probability $1-\Phi(3)$ regardless of how large the background maximum happened to be
(sample B is independent). So

$$x_{\rm argmax} \;=\; 1-\Phi(3) \;=\; 1.35\times10^{-3} \qquad\text{(independent of } \mu\text{)},
\qquad y_{\rm argmax}(\mu) \;=\; P_{\rm confirm}(\mu) \text{ from §8.2}$$

It plots as one point per $\mu$, all at the same $x$. *See [`roc_threshold_vs_argmax.png`](../plots/max_of_gaussians/roc_threshold_vs_argmax.png).*

### 12.3 The argmax point lies strictly inside the ROC

Read the dominance in both directions:

**At equal false positives** ($x = 1.35\times10^{-3}$, which is exactly $t = t^\star$), the threshold
rule finds more signal:

| $\mu$ | argmax | threshold |
|---|---|---|
| 3σ | 7.2% | **8.1%** |
| 4σ | 38.6% | **42.5%** |
| 5σ | 78.5% | **82.5%** |
| 6σ | 96.3% | **97.7%** |

**At equal power**, the threshold rule needs fewer false positives:

| $\mu$ | power | $t$ matching argmax's power | $E[K_{\rm fake}]$ | vs argmax |
|---|---|---|---|---|
| 3σ | 7.2% | 4.06 | $9.9\times10^{-4}$ | **1.4× fewer** |
| 4σ | 38.6% | 4.10 | $8.2\times10^{-4}$ | **1.6× fewer** |
| 5σ | 78.5% | 4.15 | $6.9\times10^{-4}$ | **2.0× fewer** |
| 6σ | 96.3% | 4.19 | $5.6\times10^{-4}$ | **2.4× fewer** |

The argmax point sits **below and to the right** of the threshold curve for every $\mu$. It is dominated,
so no choice of cost trade-off can ever prefer it. Ranking bins destroys information that their absolute
size carries, and the ROC is where that shows up as a measurable loss.

### 12.4 Careful: the bare ratio $y/x$ has no interior maximum

The natural instinct is to maximise "signal found per false positive", $y/x$. **That objective is
degenerate** — it grows without bound as $t$ increases. Substituting the definitions and using the
Gaussian tail asymptotic $1-\Phi(t)\simeq\phi(t)/t$:

$$\frac{y}{x} \;=\; \frac{\Phi(\mu-t)\,\Phi(\mu-3)}{n\,[1-\Phi(t)]\,[1-\Phi(3)]}
\;\;\propto\;\; \frac{1-\Phi(t-\mu)}{1-\Phi(t)} \;\;\xrightarrow[t\to\infty]{}\;\; e^{\mu t - \mu^2/2}$$

The numerator is the *signal* tail and the denominator the *background* tail; their ratio is the
likelihood ratio, which diverges exponentially. Numerically, at $\mu = 4$:

| $t$ | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| $y$ | 0.708 | 0.421 | 0.133 | 0.019 | 0.0011 | 0.00003 |
| $y/x$ | 13 | 328 | $1.1\times10^4$ | $4.8\times10^5$ | $2.2\times10^7$ | $1.1\times10^9$ |

The ratio just says **"cut harder"** — it is maximised by finding nothing, very purely. It is a valid
*comparison* between two rules at the same operating point (and the threshold does beat the argmax on it
for $t \gtrsim 3.96$), but it is not an objective you can optimise.

### 12.5 The well-posed objective: the ROC tangent

The fix is standard: put a price on a false positive. Let $c$ be how much one false positive costs
relative to the value of one confirmed discovery, and maximise the linear utility

$$U(t) \;=\; y(\mu, t) \;-\; c\,x(t)$$

Setting $dU/dt = 0$ gives the classic condition **slope of the ROC $= c$**, and here the slope has a
closed form:

$$\frac{dy}{dx} \;=\; \frac{dy/dt}{dx/dt}
\;=\; \frac{\phi(\mu-t)\,\Phi(\mu-3)}{n\,\phi(t)\,[1-\Phi(3)]}
\;=\; \frac{\Phi(\mu-3)}{n\,[1-\Phi(3)]}\;e^{\mu t - \mu^2/2}$$

using $\phi(\mu-t)/\phi(t) = e^{\mu t - \mu^2/2}$. Setting this equal to $c$ and solving for $t$:

$$\boxed{\;t_{\rm opt}(\mu, c) \;=\; \frac{\mu}{2} \;+\; \frac{1}{\mu}\,
\ln\!\left[\frac{c\;n\,[1-\Phi(3)]}{\Phi(\mu-3)}\right]\;}$$

(Verified against a numerical derivative to five digits.) The slope increases monotonically with $t$,
which is exactly the statement that the ROC is **concave** — so this stationary point is the maximum, and
it is unique.

| $c$ (cost of a fake) | $\mu=3$ | $\mu=4$ | $\mu=5$ | $\mu=6$ |
|---|---|---|---|---|
| 1 | 2.96 | 2.97 | 3.24 | 3.62 |
| 10 | 3.73 | 3.54 | 3.71 | 4.00 |
| 100 | 4.50 | 4.12 | 4.17 | 4.38 |
| 1000 | 5.27 | 4.70 | 4.63 | 4.77 |
| $10^4$ | 6.03 | 5.27 | 5.09 | 5.15 |

Two things fall out of this table:

- **$t_{\rm opt}$ depends on $\mu$.** There is no universally optimal threshold — you have to commit to
  the signal strength you are trying to be sensitive to (or average over a prior on $\mu$). The
  dependence is weak, though: over the whole $\mu \in [4,6]$ range the optimum moves by only ~0.3, so a
  single cut is a defensible compromise.
- **$t^\star = 3.99$ was never optimal for anything in particular.** Reading its implied cost off the
  slope, $c = dy/dx|_{t^\star}$, gives 21 ($\mu{=}3$), 59 ($\mu{=}4$), 41 ($\mu{=}5$), 9 ($\mu{=}6$) —
  no single exchange rate, just whatever falls out of matching the argmax's fake budget. That is exactly
  the point: $t^\star$ is a *fair comparison* point, not a recommendation. If you want an actual
  operating point, state $c$ and read $t_{\rm opt}$ off the box above.

---

# Part IV — Benjamini–Hochberg: letting the data set the cut

Parts II and III used the two rules at the extremes of rigidity: the **argmax** always returns exactly
one bin, whatever the data look like; the **fixed threshold** always applies the same bar $t$, whatever
the data look like. The Benjamini–Hochberg step-up sits between them — it returns a *variable* number of
bins, at a bar the data themselves set, and it comes with an interpretable dial: the nominal false
discovery rate $q$.

Everything else is unchanged: $n = 30{,}000$ bins in A, one of them possibly signal at $\mu$; every bin
BH selects in A is unblinded in B and confirmed at $z_B > 3$.

## 13. The rule, and why it is a softened threshold

One-sided $p_i = 1-\Phi(z_i)$, sorted $p_{(1)} \le \dots \le p_{(n)}$. BH selects the $K$ smallest, with

$$K(q) \;=\; \max\Big\{k \;:\; p_{(k)} \;\le\; \frac{k\,q}{n}\Big\}$$

### 13.1 The rank-1 bar

Take $k=1$ alone. The condition $p_{(1)} \le q/n$ is *exactly a fixed threshold*, at

$$z_1(q) \;=\; \Phi^{-1}\!\left(1 - \frac{q}{n}\right)$$

and it is the Bonferroni cut at level $q$. So **$q$ enters BH the same way $\lambda(t) = n[1-\Phi(t)]$
enters the threshold rule**: it is a budget of expected fake bins, spent at the single-bin bar $z_1(q)$.
Everything BH adds on top of Bonferroni is the step-up: ranks $2, 3, \dots$ get the *looser* bars
$2q/n,\; 3q/n,\dots$, so a bin can be selected even if it is not the most significant one.

### 13.2 The null behaviour is exactly $q$ (Daniels' theorem)

Under the global null, for independent p-values,

$$P\big(\text{BH makes} \ge 1 \text{ rejection} \,\big|\, H_0\big) \;=\; q \qquad\text{exactly.}$$

Verified by MC ($3\times10^5$ experiments): $q = 0.05 \to 0.0505$, $q = 0.2 \to 0.1990$,
$q = 0.5 \to 0.5008$. Note the corollary: the *number* of rejections is not $q$ — see §13.3.

Also: $q \to 1$ is degenerate. $p_{(n)} \le 1$ always holds, so BH at $q=1$ rejects **every bin**. The
scan below stops at $q = 0.8$.

### 13.3 BH is bursty

The threshold rule fires on a Poisson($\lambda$) number of bins: at $\lambda = 1$ it fires on 63% of null
experiments, and when it fires it usually returns one bin. BH is not like that:

| $q$ | $E[K \mid H_0]$ | $P(K\ge1\mid H_0) = q$ | $E[K \mid K \ge 1]$ | $P(K \ge 2 \mid K \ge 1)$ |
|---|---|---|---|---|
| 0.01 | 0.010 | 0.010 | 1.01 | 0.014 |
| 0.05 | 0.056 | 0.050 | 1.11 | 0.101 |
| 0.10 | 0.124 | 0.100 | 1.24 | 0.188 |
| 0.20 | 0.312 | 0.200 | 1.57 | 0.348 |
| **0.381** | **1.00** | **0.382** | **2.60** | **0.579** |
| 0.50 | 2.00 | 0.500 | 4.00 | 0.698 |
| 0.80 | 20.0 | 0.800 | 24.9 | 0.911 |

BH rarely fires, but when it does it fires in a **clump** — the step-up is self-reinforcing, because one
small p-value licenses the next rank, which licenses the next. That clumping is what costs it against the
fixed threshold (§15): the rank-2, rank-3, … selections spend fake budget without helping the signal,
which is normally rank 1.

## 14. What $z$-cut is BH actually applying, per pseudo-experiment?

This is the question that makes BH concrete for a bump hunt. Unlike $t$, **BH's bar is a random
variable**: in a PE where it rejects $K$ bins, the weakest bin it accepted is $p_{(K)}$, so the realized
cut is $z_{\rm cut} = \Phi^{-1}(1 - p_{(K)})$. Constraints:

- if $K = 1$, then $p_{(1)} \le q/n$, so $z_{\rm cut} \ge z_1(q)$ — **never looser than the nominal bar**;
- if $K = k$, then $p_{(k)} \le kq/n$, so $z_{\rm cut} > \Phi^{-1}(1 - kq/n)$ — each extra rejection
  loosens the floor by one step.

Under $H_0$, $2\times10^5$ PEs, *see [`bh_zcut.png`](../plots/max_of_gaussians/bh_zcut.png)*:

| $q$ | nominal $z_1(q)$ | matched fixed $t$ | median $z_{\rm cut}$ | central 68% | loosest seen |
|---|---|---|---|---|---|
| 0.001 | 5.40 | 5.41 | 5.53 | [5.42, 5.70] | 5.40 |
| 0.01 | 4.97 | 4.97 | 5.11 | [5.01, 5.31] | 4.84 |
| 0.05 | 4.65 | 4.63 | 4.77 | [4.67, 4.99] | 4.37 |
| 0.10 | 4.50 | 4.46 | 4.61 | [4.51, 4.83] | 4.08 |
| 0.20 | 4.35 | 4.26 | 4.43 | [4.23, 4.66] | 3.74 |
| **0.381** | **4.21** | **3.99** | **4.18** | [3.91, 4.43] | 3.37 |
| 0.50 | 4.15 | 3.82 | 4.00 | [3.72, 4.31] | 3.11 |
| 0.80 | 4.04 | 3.21 | 3.40 | [3.03, 3.89] | 2.19 |

Read it in two pieces.

**Low $q$: BH *is* Bonferroni.** At $q \le 0.05$, $K \ge 2$ almost never happens, so essentially every
rejection is a lone rank-1 one and $z_{\rm cut} \ge z_1(q)$ by construction. The median sits ~0.12 *above*
$z_1$, which is not the rule being conservative but a selection effect: conditioned on $p_{(1)} \le q/n$,
the triggering p-value is roughly uniform on $(0, q/n)$, so its median is $\approx q/2n$, i.e. half a step
above the bar. In this regime BH and the matched fixed threshold are the same cut ($z_1 \simeq t$ to two
decimals) and, unsurprisingly, they have the same power.

**High $q$: the step-up takes over.** Past $q^\star$ the multiple-rejection PEs dominate, the realized
bar drops *below* $z_1(q)$, and the two columns separate: at $q = 0.8$ BH's nominal bar is still 4.04 but
the median cut it actually applies is 3.40, with a tail down to 2.19. Note that $z_1(q)$ is nearly flat
in this regime (4.35 → 4.04 across a 4× change in $q$) — all the extra fake budget is bought by the step-up,
not by lowering the nominal bar.

**With a signal present the bar loosens further**, because the signal's own small p-value licenses rank 2:

| $q$ | median $z_{\rm cut}$, $H_0$ | $\mu = 4\sigma$ | $\mu = 5\sigma$ |
|---|---|---|---|
| 0.01 | 5.11 | 5.35 | 5.64 |
| 0.05 | 4.77 | 5.00 | 5.35 |
| 0.20 | 4.43 | 4.52 | 4.77 |
| 0.381 | 4.18 | 4.19 | 4.24 |

(The medians go *up* with $\mu$ because the accepted set is now usually just the signal bin itself, which
sits at high $z$ — the *floor* has loosened, but the weakest accepted bin is the signal.)

## 15. Where BH lands: the three rules on one ROC

Same axes as §12.1, so all three rules are comparable:
$x = E[K_{\rm fake}]$, $y = P(\text{signal confirmed in B})$. *See [`bh_scan.png`](../plots/max_of_gaussians/bh_scan.png), [`bh_vs_argmax.png`](../plots/max_of_gaussians/bh_vs_argmax.png),
[`roc_bh_vs_threshold.png`](../plots/max_of_gaussians/roc_bh_vs_threshold.png).*

### 15.1 $q^\star$: the FDR that buys the argmax's fake budget

Solving $E[K \mid H_0] = 1$ (the argmax's budget, which is what defined $t^\star$) gives

$$\boxed{\;q^\star \;=\; 0.381\;}$$

and it matches the argmax on *both* false-alarm measures: $E[\text{false confirmations}] =
1.34\times10^{-3}$ (argmax: $1.35\times10^{-3}$) and $P(\ge1 \text{ false confirmation}) =
1.336\times10^{-3}$, i.e. **exactly $3.0\sigma$ global**. Its effective single-bin bar is
$z_1(q^\star) = 4.21$, *stricter* than $t^\star = 3.99$.

### 15.2 Power at a matched $3.0\sigma$ global false-alarm rate

| $\mu$ | argmax | **BH at $q^\star$** | threshold at $t^\star$ |
|---|---|---|---|
| 3σ | 7.2% | **7.5%** | **8.1%** |
| 4σ | 38.6% | **40.3%** | **42.5%** |
| 5σ | 78.5% | **80.6%** | **82.5%** |
| 6σ | 96.3% | **97.2%** | **97.7%** |

**BH lands strictly between the two.**

- *It beats the argmax* because, like the threshold, it never destroys the signal candidate merely because
  some background bin happened to be larger; and the step-up lets the signal ride along at rank 2+ even
  when it is not the most significant bin.
- *It loses to the fixed threshold* because of §13.3. To keep $E[K\mid H_0] = 1$ while spending part of
  that budget on rank-2+ clumps, it must raise the single-bin bar to 4.21 — and those extra rank-2+
  selections are pure waste: they are background bins, and the signal (normally rank 1) gains nothing from
  them.

### 15.3 The loss is a large-$q$ effect, not an indictment of BH

Match BH and the fixed threshold at the *same* fake budget all the way down the ROC, and the gap closes:

| $q$ | $E[K \mid H_0]$ | matched $t$ | BH, $\mu{=}4\sigma$ | threshold, $\mu{=}4\sigma$ |
|---|---|---|---|---|
| 0.001 | 0.0009 | 5.42 | **6.8%** | 6.6% |
| 0.008 | 0.008 | 5.02 | 13.0% | 13.0% |
| 0.05 | 0.060 | 4.61 | 22.5% | **22.7%** |
| 0.18 | 0.275 | 4.28 | 31.9% | **32.7%** |
| 0.38 | 1.00 | 3.99 | 40.3% | **42.5%** |
| 0.57 | 3.08 | 3.71 | 47.9% | **51.6%** |

Below $q \approx 0.01$ the two rules are indistinguishable — BH *is* the Bonferroni cut there (§14), and
if anything a hair ahead, because the rank-2 ride-along helps the signal marginally more often than the
clumping hurts. The gap only opens once $q$ is large enough for multiple rejections to be common. So the
honest statement is not "BH is worse than a threshold" — it is **"BH run loose is worse than a threshold
run loose, and BH run tight is a threshold."**

### 15.4 What the FDR guarantee costs in global significance

The dial is interpretable, which is BH's real selling point, but the conventional choice is far more
conservative than either of the other two rules' operating points:

| $q$ | $E[K\mid H_0]$ | $P(\ge1$ false confirmation$)$ | global $Z$ |
|---|---|---|---|
| 0.05 | 0.056 | $7.6\times10^{-5}$ | 3.79σ |
| 0.10 | 0.124 | $1.7\times10^{-4}$ | 3.59σ |
| 0.20 | 0.312 | $4.2\times10^{-4}$ | 3.34σ |
| 0.381 | 1.00 | $1.34\times10^{-3}$ | 3.00σ |
| 0.50 | 2.00 | $2.7\times10^{-3}$ | 2.78σ |

And note what the guarantee delivers *after* stage B: the realized FDR among **confirmed** candidates is
near zero across the whole range (dashed curves in [`roc_bh_vs_threshold.png`](../plots/max_of_gaussians/roc_bh_vs_threshold.png)), because the $z_B > 3$
confirmation kills fakes at $1.35\times10^{-3}$ each regardless of how many stage 1 let through. The FDR
control is doing its work at stage 1, where it is not really what we are buying — the false-alarm rate we
care about is set by stage B. That is the structural reason BH's headline guarantee is less useful here
than in a one-stage analysis.

---

## 16. Reproducing

Plots:

| File | Content |
|---|---|
| [`max_of_gaussians_light.png`](../plots/max_of_gaussians/max_of_gaussians_light.png) | left: one pseudo-experiment with its max marked. right: the distribution of the 10,000 maxima, with the exact density and the asymptotic Gumbel overlaid |
| [`signal_wins_the_max.png`](../plots/max_of_gaussians/signal_wins_the_max.png) | left: the race — signal densities (3–6σ) against the background-max distribution. right: $p(\mu)$, the probability the signal wins the scan |
| [`ab_confirmation.png`](../plots/max_of_gaussians/ab_confirmation.png) | left: the two-stage funnel, stage 1 vs confirmed. right: confirmation power vs $\mu$, with the 3σ-global false-alarm floor |
| [`threshold_scan.png`](../plots/max_of_gaussians/threshold_scan.png) | left: $\lambda(t)$, the expected background bins above $t$, with $t^\star$ marked. right: signal efficiency $\Phi(\mu-t)$, crossing 50% at $t=\mu$ |
| [`threshold_vs_argmax.png`](../plots/max_of_gaussians/threshold_vs_argmax.png) | left: confirmation power vs $t$, with the argmax rule as dashed reference. right: the false-alarm rate vs $t$, meeting the argmax rule at $t^\star$ |
| [`roc_threshold_vs_argmax.png`](../plots/max_of_gaussians/roc_threshold_vs_argmax.png) | left: the ROC — power vs average false positives, with the argmax as a single dominated point. right: the ratio $y/x$, showing it has no interior maximum |
| [`bh_scan.png`](../plots/max_of_gaussians/bh_scan.png) | left: $E[K]$ and $P(\text{any rejection}) = q$ vs the nominal FDR, with $q^\star$ marked. right: stage-1 signal efficiency, MC vs the rank-1 bar alone |
| [`bh_vs_argmax.png`](../plots/max_of_gaussians/bh_vs_argmax.png) | left: BH's confirmation power vs $q$, with the argmax rule as dashed reference. right: the global false-alarm rate a nominal $q$ buys, meeting the argmax at $q^\star$ |
| [`roc_bh_vs_threshold.png`](../plots/max_of_gaussians/roc_bh_vs_threshold.png) | left: all three rules on the §12 ROC. right: the realized FDR before and after the stage-B confirmation |
| [`bh_zcut.png`](../plots/max_of_gaussians/bh_zcut.png) | the realized $z$-cut per pseudo-experiment (median and central 68%) against the nominal bar $\Phi^{-1}(1-q/n)$ and the matched fixed threshold |

Every figure and every number quoted above is produced by three scripts (`make stats`), all seeded,
so they reproduce bit-for-bit:

| Script | What it does |
|---|---|
| `searchbudget/stages/max_of_gaussians.py` | Parts I–III: the six figures above and the §5 validation table (moments, KS tests) |
| `searchbudget/stages/bh_fdr_ab.py` | the $q$ scan: null fake budget, signal efficiency, power, the ROC. Writes `results/tables/bh_fdr_scan.csv` and the three BH figures. MC is cached in `results/tables/bh_fdr_mc.npz`; pass `--refit` to redo it |
| `searchbudget/stages/bh_zcut.py` | the per-PE realized cut of §14. Writes `results/tables/bh_zcut_per_pe.csv` and [`bh_zcut.png`](../plots/max_of_gaussians/bh_zcut.png) |

**Part IV is Monte Carlo** — $K(q)$ is data-adaptive, so there is no closed form. Both BH scripts use
one trick worth knowing: BH depends on the p-vector only through its smallest few order
statistics, and the $m$ smallest of $n$ iid uniforms have the exact representation
$U_{(k)} = \Gamma_k/\Gamma_{(n+1)}$ with $\Gamma_k$ a sum of $k$ iid Exp(1) — so one Gamma draw plus $m$
exponentials gives the exact joint law without ever generating $n$ values. And $K(q)$ vectorizes over
both experiments and the $q$-grid via a suffix-minimum: with $R_k = n\,p_{(k)}/k$ and
$C_k = \min_{j \ge k} R_j$ (non-decreasing in $k$), $K(q) = \#\{k : C_k \le q\}$.

The whole of Parts I–III condenses to this, if you want to check the numbers without the plotting:

```python
import numpy as np
from scipy.stats import norm, gumbel_r, kstest
from scipy.integrate import quad

n, nb, t = 30_000, 29_999, 3.0

# ---- Part I: the background maximum ----
rng = np.random.default_rng(0)
M   = rng.standard_normal((10_000, n)).max(axis=1)      # 10,000 pseudo-experiments

ln_n = np.log(n)
b_n  = 1/np.sqrt(2*ln_n)
a_n  = np.sqrt(2*ln_n) - (np.log(ln_n) + np.log(4*np.pi))/(2*np.sqrt(2*ln_n))

pdf_exact = lambda x: n*norm.pdf(x)*norm.cdf(x)**(n-1)   # exact density
quantile  = lambda q: norm.ppf(q**(1/n))                 # exact quantiles

print(kstest(M, lambda x: norm.cdf(x)**n))                       # p = 0.59  -> exact
print(kstest(M, lambda x: gumbel_r.cdf(x, loc=a_n, scale=b_n)))  # p = 2e-15 -> rejected

# ---- Part II: signal ----
p_win  = lambda mu: quad(lambda s: norm.pdf(s-mu)*norm.cdf(s)**nb, mu-12, mu+12, limit=400)[0]
p_conf = lambda mu: p_win(mu)*norm.sf(t-mu) + (1-p_win(mu))*norm.sf(t)

for mu in (3, 4, 5, 6):
    print(f"mu={mu}:  wins scan {p_win(mu):.3f}   confirmed {p_conf(mu):.3f}")
# mu=3:  wins scan 0.142   confirmed 0.072
# mu=4:  wins scan 0.458   confirmed 0.386
# mu=5:  wins scan 0.803   confirmed 0.785
# mu=6:  wins scan 0.964   confirmed 0.963

print(f"false-alarm rate = {norm.sf(t):.2e}  ({norm.isf(norm.sf(t)):.1f} sigma, GLOBAL)")

# ---- Part III: threshold selection instead of the argmax ----
lam    = lambda tt: nb*norm.sf(tt)                       # expected background bins above tt
p_fake = lambda tt: 1 - np.exp(-lam(tt)*norm.sf(t))      # P(>=1 false confirmation)
p_thr  = lambda mu, tt: norm.cdf(mu-tt)*norm.cdf(mu-t)   # P(signal confirmed)

tstar = norm.isf(1/nb)                                   # 3.9879 -> lambda = 1
print(f"t* = {tstar:.3f},  lambda = {lam(tstar):.3f},  P_fake = {p_fake(tstar):.2e}")

for mu in (3, 4, 5, 6):
    print(f"mu={mu}:  argmax {p_conf(mu):.3f}   threshold@t* {p_thr(mu, tstar):.3f}")
# mu=3:  argmax 0.072   threshold@t* 0.081
# mu=4:  argmax 0.386   threshold@t* 0.425
# mu=5:  argmax 0.785   threshold@t* 0.825
# mu=6:  argmax 0.963   threshold@t* 0.977

# ---- the ROC, and the operating point implied by a cost c ----
x_roc = lambda tt: nb*norm.sf(tt)*norm.sf(t)          # avg. number of false positives
y_roc = lambda mu, tt: norm.cdf(mu-tt)*norm.cdf(mu-t) # P(signal confirmed)
x_arg = norm.sf(t)                                    # argmax: one bin, always

slope = lambda mu, tt: norm.pdf(mu-tt)*norm.cdf(mu-t)/(nb*norm.pdf(tt)*norm.sf(t))
t_opt = lambda mu, c: mu/2 + np.log(c*nb*norm.sf(t)/norm.cdf(mu-t))/mu   # dy/dx = c

print(t_opt(4, 100), slope(4, t_opt(4, 100)))   # 4.12, 100.0  -> closed form is exact
```
