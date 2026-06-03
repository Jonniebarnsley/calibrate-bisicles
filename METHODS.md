# Methods

> Draft. §1–§2 are author-written (ensemble details); §3 onward drafted here in full.
> Citations are `[Author Year]` placeholders.

---

## 1. Overview

*[author to write]* — framing of the approach (single-model PPE, emulation, Bayesian
calibration against IMBIE, projection to 2300, variance decomposition); what we add
relative to forward-propagation studies; pipeline schematic; notation.

## 2. Perturbed-parameter ensemble

*[author to write]* — BISICLES set-up; the six perturbed parameters and their physical
roles; sampling design and ranges; the 4 GCM × 3 scenario forcing and how the
`thermal_forcing` (ocean) and `warming` (atmosphere) covariates are derived; output
conventions (cumulative SLC, m SLE, mass-loss positive, 2007 = 0).

Throughout we write a simulator input as $x = (\boldsymbol{\theta}, \mathbf{g}, \mathbf{f})$,
where $\boldsymbol{\theta} \in \mathbb{R}^{6}$ are the perturbed physical parameters,
$\mathbf{g}$ is the GCM (encoded by dummy variables), and $\mathbf{f}$ the two
forcing covariates. The simulator output is the cumulative sea-level contribution
$y(x, t)$ at year $t \in \{2007, \dots, 2300\}$, in mm SLE.

---

## 3. Observational constraint

We calibrate against the Antarctic contribution to sea level recorded by the IMBIE
assessment [IMBIE 2018], expressed as cumulative SLC relative to a 2007 baseline.
Rather than fitting the full annual series over the observational window — which is
strongly temporally autocorrelated and would require an explicit annual error
covariance — we reduce it to a single scalar: the cumulative SLC accumulated between
$t_0 = 2007$ and $t_1 = 2021$,

$$
Y \;=\; C_\mathrm{obs}(t_1) - C_\mathrm{obs}(t_0),
$$

where $C_\mathrm{obs}(t)$ is the IMBIE cumulative contribution. This matches the
quantity the ensemble reports (Section 2) and, for the IMBIE record, gives
$Y \approx 5.28$ mm SLE (mass loss taken positive).

The observational uncertainty must likewise be that of the *difference*. Writing
$\sigma_0, \sigma_1$ for the reported cumulative uncertainties at $t_0, t_1$, the
variance of $Y$ is

$$
\sigma_\mathrm{obs}^2
= \operatorname{Var}\!\big[C_\mathrm{obs}(t_1) - C_\mathrm{obs}(t_0)\big]
= \sigma_1^2 + \sigma_0^2 - 2\,\operatorname{Cov}\!\big[C_\mathrm{obs}(t_1), C_\mathrm{obs}(t_0)\big].
$$

IMBIE constructs the cumulative uncertainty as the root-sum-square of *independent*
monthly error contributions [IMBIE 2018], i.e. $C_\mathrm{obs}(t)$ is an accumulation
of independent increments. Then $C_\mathrm{obs}(t_1) = C_\mathrm{obs}(t_0) + I$ with the
increment $I$ over $(t_0, t_1]$ independent of $C_\mathrm{obs}(t_0)$, so
$\operatorname{Cov}[C_\mathrm{obs}(t_1), C_\mathrm{obs}(t_0)] = \sigma_0^2$ and

$$
\sigma_\mathrm{obs}^2 = \sigma_1^2 - \sigma_0^2,
\qquad
\sigma_\mathrm{obs} = \sqrt{1.47^2 - 0.96^2} \approx 1.11~\text{mm}.
$$

The two common shortcuts are special cases of the covariance term and are both
inappropriate here: taking $\sigma_\mathrm{obs} = \sigma_1 = 1.47$ mm assumes the 2007
baseline is error-free ($\sigma_0 = 0$), while the arithmetic difference
$\sigma_1 - \sigma_0 = 0.51$ mm assumes perfectly correlated endpoints
($\operatorname{Cov} = \sigma_0\sigma_1$, giving $(\sigma_1 - \sigma_0)^2$) — neither is
consistent with IMBIE's independent-increment construction.

---

## 4. Emulation

The ensemble is informative but too expensive to resample, so we build a Gaussian-
process (GP) emulator of the simulator [O'Hagan 2006; Gu et al. 2019]. We compare two
formulations — independent per-year emulators and a singular-value-decomposition (SVD)
basis-function emulator — and adopt the latter.

### 4.1 Gaussian-process emulator

Let a GP be fitted to $n$ training pairs $\{(x_j, y_j)\}_{j=1}^{n}$ with inputs
$x_j \in \mathbb{R}^{p}$ ($p = 11$: six parameters, three GCM dummies, two forcing
covariates). We use a universal-kriging GP with a linear-trend mean,

$$
y(x) = h(x)^{\!\top}\boldsymbol{\beta} + Z(x),
\qquad
h(x) = (1, x^{\!\top})^{\!\top},
$$

where $Z(\cdot)$ is a zero-mean GP with variance $\sigma^2$ and correlation function
$c(\cdot,\cdot)$. We use a separable Matérn-5/2 correlation,

$$
c(x, x') = \prod_{d=1}^{p}
\left(1 + \frac{\sqrt{5}\,|x_d - x'_d|}{\gamma_d}
+ \frac{5\,(x_d - x'_d)^2}{3\,\gamma_d^2}\right)
\exp\!\left(-\frac{\sqrt{5}\,|x_d - x'_d|}{\gamma_d}\right),
$$

with range parameters $\{\gamma_d\}$. Inputs are rescaled to $[0,1]$ (with aPhi, gamma,
UMV log-transformed first); the range parameters and a nugget $\eta$ are estimated by
RobustGaSP's robust marginal-posterior procedure [Gu et al. 2019], which guards against
the near-degenerate correlations that destabilise maximum likelihood.

Writing $\mathbf{H}$ for the $n \times 2$ trend matrix (rows $h(x_j)^{\!\top}$),
$\tilde{\mathbf{R}} = \mathbf{R} + \eta \mathbf{I}$ for the nugget-augmented correlation
matrix ($R_{jk} = c(x_j, x_k)$), and $r(x^\ast) = (c(x^\ast, x_j))_{j}$, the generalised
least-squares trend estimate and process variance are

$$
\hat{\boldsymbol{\beta}} = (\mathbf{H}^{\!\top}\tilde{\mathbf{R}}^{-1}\mathbf{H})^{-1}
\mathbf{H}^{\!\top}\tilde{\mathbf{R}}^{-1}\mathbf{y},
\qquad
\hat{\sigma}^2 = \tfrac{1}{n-2}
(\mathbf{y} - \mathbf{H}\hat{\boldsymbol{\beta}})^{\!\top}\tilde{\mathbf{R}}^{-1}
(\mathbf{y} - \mathbf{H}\hat{\boldsymbol{\beta}}).
$$

The predictive mean and variance at a new input $x^\ast$ are

$$
m(x^\ast) = h(x^\ast)^{\!\top}\hat{\boldsymbol{\beta}}
+ r(x^\ast)^{\!\top}\tilde{\mathbf{R}}^{-1}(\mathbf{y} - \mathbf{H}\hat{\boldsymbol{\beta}}),
$$

$$
v(x^\ast) = \hat{\sigma}^2\Big[\,1 + \eta
- r(x^\ast)^{\!\top}\tilde{\mathbf{R}}^{-1} r(x^\ast)
+ u(x^\ast)^{\!\top}(\mathbf{H}^{\!\top}\tilde{\mathbf{R}}^{-1}\mathbf{H})^{-1} u(x^\ast)\,\Big],
$$

with $u(x^\ast) = h(x^\ast) - \mathbf{H}^{\!\top}\tilde{\mathbf{R}}^{-1} r(x^\ast)$ the
trend-extrapolation correction.

**Forcing inputs.** The two forcing covariates are strongly collinear across the twelve
$(\text{GCM}, \text{scenario})$ cells (Pearson $r = 0.965$). A leave-one-out comparison
nonetheless shows that at least one is essential — with none, the emulator cannot
separate scenarios and its end-of-century skill collapses (normalised RMSE $\approx
0.87$) — and that retaining both is marginally superior to either alone at all lead
times; we therefore retain both.

### 4.2 Per-year emulation

The first formulation fits an independent GP $\hat{f}_t(x)$ to each year's cumulative-SLC
field, for $t = 2007, \dots, 2300$. Because the output is already cumulative from the
2007 baseline, $\hat{f}_t$ directly emulates the quantity of interest at each horizon
with no summation of annual increments. The cost scales with the number of years
($\sim\!300$ GPs), and, since the $\hat{f}_t$ are fitted independently, the emulated
trajectory $t \mapsto \hat{f}_t(x)$ of a fixed $x$ need not be smooth in time.

### 4.3 SVD basis-function emulation

The second formulation represents whole trajectories in a low-dimensional temporal
basis [Higdon et al. 2008]. Collect the ensemble in a matrix
$\mathbf{Y} \in \mathbb{R}^{n \times T}$ ($n$ runs, $T$ years), with column (per-year)
mean $\mu_t$ and standard deviation $s_t$. Standardise each year,

$$
\tilde{Y}_{jt} = \frac{Y_{jt} - \mu_t}{s_t},
$$

and take the singular value decomposition $\tilde{\mathbf{Y}} = \mathbf{U}\mathbf{D}\mathbf{V}^{\!\top}$.
The columns of $\mathbf{V}$ are orthonormal temporal basis functions $\phi_i(t) = V_{ti}$,
and the scores $a_{ji} = (\mathbf{U}\mathbf{D})_{ji}$ are the coordinates of run $j$ on
basis $i$; the scores are mutually orthogonal across components,
$\sum_{j} a_{ji} a_{ji'} = d_i^2\,\delta_{ii'}$.

We then emulate the leading $r$ coefficients as functions of the inputs — a GP
$g_i(\cdot)$ trained on $\{(x_j, a_{ji})\}_{j}$ for $i = 1, \dots, r$, each as in
Section 4.1 — giving predictive means $\hat{a}_i(x^\ast)$ and variances
$\hat{\sigma}_i^2(x^\ast)$. A trajectory at a new input is reconstructed as

$$
\hat{Y}(x^\ast, t) = \mu_t + s_t \sum_{i=1}^{r} \hat{a}_i(x^\ast)\,\phi_i(t).
$$

This replaces $\sim\!300$ annual GPs by $r$ coefficient GPs. **Per-year** (rather than
single, global) standardisation is important: cumulative Antarctic SLC variance grows
by roughly three orders of magnitude between the observational window and 2300, so a
single scale would let the leading components be dominated by the late-century signal
and reconstruct the early, calibration-relevant period poorly.

### 4.4 Predictive variance and the truncation term

The reconstruction variance has two parts. Treating the coefficient emulators as
mutually independent (the scores are uncorrelated by construction), the retained
components contribute $s_t^2 \sum_{i\le r} \phi_i(t)^2\,\hat{\sigma}_i^2(x^\ast)$. The
truncated components $i > r$ are absent from the reconstruction but contribute to the
true trajectory; their (standardised) variance at year $t$ is

$$
\tau_t = \operatorname{Var}\!\Big[\textstyle\sum_{i>r} a_{i}\phi_i(t)\Big]
= \sum_{i>r} \phi_i(t)^2 \operatorname{Var}[a_i]
= \frac{1}{n}\sum_{j}\Big(\tilde{Y}_{jt} - \textstyle\sum_{i\le r} a_{ji}\phi_i(t)\Big)^{2},
$$

where the second equality uses score orthogonality and the third evaluates $\tau_t$
directly as the mean-squared truncation residual. The total predictive variance is then

$$
\operatorname{Var}\hat{Y}(x^\ast, t)
= s_t^2\!\left[\sum_{i=1}^{r}\phi_i(t)^2\,\hat{\sigma}_i^2(x^\ast) \;+\; \tau_t\right].
$$

The truncation term $\tau_t$ is essential. Omitting it leaves the reconstruction blind
to the discarded modes and makes its predictive intervals over-confident: leave-one-out
$95\%$ coverage at the observational target falls to $\approx 0.83$, and the
standardised-residual (Mahalanobis) distance
$\mathrm{MD} = \big(\sum_j (\hat{Y}_{-j} - Y_j)^2 / v_j\big)^{1/2}$ sits well above its
expected value $\sqrt{n}$. Including $\tau_t$ restores coverage to $\approx 0.95$ and
$\mathrm{MD} \approx \sqrt{n}$. $\tau_t$ is the deterministic, per-year analogue of the
single i.i.d. residual precision estimated within the fully Bayesian formulation of
[Higdon et al. 2008].

### 4.5 The truncation–emulation trade-off and component count

Per-coefficient emulability is assessed by leave-one-out,

$$
R^2_i = 1 - \frac{\sum_{j}\big(\hat{a}_{i,-j} - a_{ji}\big)^2}{\sum_{j}\big(a_{ji} - \bar{a}_i\big)^2},
$$

where $\hat{a}_{i,-j}$ predicts $a_{ji}$ with run $j$ held out. The number of retained
components $r$ is then governed by a trade-off, not a fixed variance threshold. Adding a
component *reduces* the truncation variance $\tau_t$ (removing uncertainty), but
successive components carry progressively finer, less input-dependent structure that is
emulated less faithfully ($R^2_i$ declines with $i$), so each added component
*contributes* more emulation variance per unit of signal. The competing terms in
$\operatorname{Var}\hat{Y}$ thus balance: the end-to-end leave-one-out reconstruction
error and its calibration improve sharply over the first few components and then
plateau — for this ensemble by $r \approx 3$–$5$ — beyond which further components
remove little $\tau_t$ while injecting poorly-constrained modes. We select $r$ from
where the reconstruction error flattens; a cumulative-variance rule would over-count
components.

### 4.6 Comparison and choice

The two formulations agree closely: the calibrated effective sample sizes (Section 5)
and the projected SLC distributions (Section 6) from the SVD emulator reproduce those
from the per-year emulator to within sampling noise. Given this agreement, together
with the SVD emulator's far lower cost ($r \approx 5$ coefficient GPs versus $\sim\!300$
annual fits), its temporally coherent and smooth reconstructions, and its well-
calibrated predictive uncertainty once $\tau_t$ is included, we adopt the SVD basis-
function emulator with $r = 5$ for all subsequent calibration, projection and
sensitivity analyses. Hereafter $\hat{f}(x, t)$ and $v(x, t)$ denote its predictive
mean and variance.

---

## 5. Bayesian calibration (observational weighting)

We update a prior over the physical parameters $\boldsymbol{\theta}$ by their likelihood
given the observed scalar $Y$ (Section 3), conditional on the controllable inputs (GCM
and scenario). Let $M(x) = \hat{f}(x, t_1)$ be the emulated obs-period quantity and
$v(x, t_1)$ its emulator variance. The total residual variance combines observation,
model-discrepancy and emulator terms,

$$
\sigma_\mathrm{tot}^2(x) = \sigma_\mathrm{obs}^2 + \sigma_\mathrm{mod}^2 + v(x, t_1),
$$

with a Gaussian likelihood for the aggregated target,

$$
\log \mathcal{L}(x) = -\tfrac{1}{2}\log\sigma_\mathrm{tot}^2(x)
- \frac{\big(Y - M(x)\big)^2}{2\,\sigma_\mathrm{tot}^2(x)} + \text{const},
$$

the $\sigma$-dependent normalisation being retained because $\sigma_\mathrm{tot}$ varies
between samples (the emulator variance is heteroscedastic).

**Model discrepancy.** We set $\sigma_\mathrm{mod} = c\,\sigma_\mathrm{obs}$, anchoring
the structural model–reality discrepancy to the observation error rather than to the
ensemble spread (the latter being parametric *signal*, not noise); $c$ is the primary
sensitivity lever (Section 7).

**Prior and sampling.** The prior $\pi(\boldsymbol{\theta})$ is the explicit design prior
— independent uniform (n, m, uf) and log-uniform (aPhi, gamma, UMV) over the
perturbation ranges. We draw $K$ Monte-Carlo samples $\boldsymbol{\theta}_k \sim \pi$ and
evaluate $M, v$ through the emulator.

**Weighting.** Drawing from the prior and targeting the posterior
$p(\boldsymbol{\theta}\,|\,Y) \propto \pi(\boldsymbol{\theta})\,\mathcal{L}(\boldsymbol{\theta})$
makes this self-normalised importance sampling with weights proportional to the
likelihood. Computed *within each GCM stratum* and evaluated in log-space for stability,

$$
w_k = \frac{\exp\!\big(\log\mathcal{L}_k - \max_{k'}\log\mathcal{L}_{k'}\big)}
{\sum_{l}\exp\!\big(\log\mathcal{L}_l - \max_{k'}\log\mathcal{L}_{k'}\big)},
\qquad \sum_k w_k = 1.
$$

Weights are pooled across scenarios (the obs-period prediction is effectively
scenario-invariant, Section 4.1) and computed separately per GCM; the prior over GCMs is
held fixed (not data-updated). The weighted sample $\{(\boldsymbol{\theta}_k, w_k)\}$
represents the per-GCM posterior.

**Diagnostics.** Constraint strength is summarised by the Kish effective sample size,

$$
\mathrm{ESS} = \frac{\big(\sum_k w_k\big)^2}{\sum_k w_k^2} = \frac{1}{\sum_k w_k^2},
$$

reported as $\mathrm{ESS}/K$ per stratum, alongside a check that $Y$ lies within the
sampled range of $M$. This two-stage (fit-then-weight) procedure is a plug-in
approximation to the joint Bayesian treatment of [Kennedy & O'Hagan 2001; Higdon et al.
2008], in which the emulator hyperparameters, $\boldsymbol{\theta}$ and the discrepancy
are inferred together.

---

## 6. Projection

For each (GCM $g$, scenario $s$, future year $t$) we reconstruct the SLC of the
posterior-weighted parameter samples under the scenario-correct forcing,
$M_k = \hat{f}(\boldsymbol{\theta}_k, g, s, t)$, carrying their per-GCM weights
$w_k^{(g)}$ from Section 5 (the weights depend on GCM only and are reused across
scenarios). The calibrated distribution is summarised by weighted quantiles obtained
from the weighted empirical CDF: with the $M_k$ sorted ascending and Hazen plotting
positions

$$
W_{(k)} = \frac{\sum_{l \le k} w_{(l)} - \tfrac{1}{2} w_{(k)}}{\sum_l w_{(l)}},
$$

the quantile $Q(p)$ is obtained by linear interpolation of $W_{(k)} \mapsto M_{(k)}$ at
probability $p \in \{0.05, 0.17, 0.50, 0.83, 0.95\}$. Prior (unweighted) bands are
reported alongside for comparison. A marginal over GCMs is formed by pooling all
samples with an explicit, uniform GCM prior, $w_k^{(g)} / N_\mathrm{GCM}$ — a labelled
judgement, since the GCMs are not calibrated.

---

## 7. Sensitivity analysis

We test robustness of the calibrated projections to the main subjective choices, one at
a time about the baseline: the discrepancy multiplier $c$; the observation uncertainty
$\sigma_\mathrm{obs}$ (endpoint vs independent-increment, Section 3); and the likelihood
family. For the latter we replace the Gaussian by a scaled Student-$t$ with $\nu$
degrees of freedom,

$$
\log \mathcal{L}(x) = -\log\sigma_\mathrm{tot}(x)
+ \log t_\nu\!\left(\frac{Y - M(x)}{\sigma_\mathrm{tot}(x)}\right),
$$

which down-weights outlying runs. Each variant re-uses the fixed emulator predictions
and only recomputes the weights (Section 5), so the sweep is inexpensive; we report the
stability of the 2100 and 2300 posterior quantiles across settings.

---

## 8. Variance decomposition of ensemble uncertainty (methods)

To attribute the (prior) ensemble SLC spread to its sources and their evolution in time
(cf. [Seroussi et al. 2023; Coulon et al. 2025]) we use two complementary,
variance-based techniques at each year.

**Analysis of variance (discrete factors).** For the factorial GCM × scenario design we
model the per-run SLC at a fixed year as

$$
Y = \beta_0 + \alpha_g + \beta_s + (\alpha\beta)_{gs} + \varepsilon,
$$

with GCM main effect $\alpha_g$, scenario main effect $\beta_s$, their interaction
$(\alpha\beta)_{gs}$, and residual $\varepsilon$. The total sum of squares partitions as
$\mathrm{SS}_\mathrm{tot} = \mathrm{SS}_\mathrm{GCM} + \mathrm{SS}_\mathrm{scen} +
\mathrm{SS}_\mathrm{int} + \mathrm{SS}_\mathrm{res}$, and each factor's contribution to
the variance is the fraction $\mathrm{SS}_\bullet / \mathrm{SS}_\mathrm{tot}$. We use
**Type III** (marginal) sums of squares — each term assessed after all others — which
are appropriate when interactions are present and the design is only near-balanced; the
residual term measures the within-cell variance, i.e. that due to the perturbed
parameters. ANOVA operates on the raw ensemble and so is independent of the emulator.

**Sobol' variance-based sensitivity (continuous inputs).** For the parametric variation
we use Sobol' indices, which require no linearity assumption. For a square-integrable
$Y = f(X)$ with independent inputs $X = (X_1, \dots, X_d)$, the first-order and
total-order indices of input $i$ are

$$
S_i = \frac{\operatorname{Var}_{X_i}\!\big[\mathbb{E}_{X_{\sim i}}(Y \mid X_i)\big]}{\operatorname{Var}(Y)},
\qquad
S_{Ti} = \frac{\mathbb{E}_{X_{\sim i}}\!\big[\operatorname{Var}_{X_i}(Y \mid X_{\sim i})\big]}{\operatorname{Var}(Y)}
= 1 - \frac{\operatorname{Var}_{X_{\sim i}}\!\big[\mathbb{E}_{X_i}(Y \mid X_{\sim i})\big]}{\operatorname{Var}(Y)},
$$

where $X_{\sim i}$ denotes all inputs but $i$. $S_i$ is the variance attributable to
$X_i$ alone (its main effect) and $S_{Ti}$ that attributable to $X_i$ including all
interactions, so $\sum_i S_i \le 1 \le \sum_i S_{Ti}$ and $S_{Ti} - S_i$ measures
interaction. We estimate them with the Saltelli scheme [Saltelli et al. 2010]: two
independent input matrices $\mathbf{A}, \mathbf{B}$ ($N \times d$) and the hybrid
matrices $\mathbf{A}_B^{(i)}$ (matrix $\mathbf{A}$ with column $i$ replaced by that of
$\mathbf{B}$), with $N(d+2)$ model evaluations supplied by the emulator. The discrete
GCM dimension is represented continuously through the forcing-severity axis (within a
scenario the axis spans the GCMs' forcing values), and the emulator reconstruction —
including the truncation variance of Section 4.4 — provides the responses.

---

## 9. Implementation and reproducibility

The pipeline is implemented in R using RobustGaSP [Gu et al. 2019] for emulation,
`sensitivity` for the Sobol' indices, and `car` for Type III ANOVA, and is driven by a
single configuration file (all subjective choices in one place). *[Code and data
availability statement — author to complete.]*

---

### References to seed
Higdon et al. (2008); Kennedy & O'Hagan (2001); O'Hagan (2006); Gu et al. (2019,
RobustGaSP); Saltelli et al. (2010); Edwards et al. (2021); Seroussi et al. (2023);
Coulon et al. (2025); IMBIE Team (2018).
