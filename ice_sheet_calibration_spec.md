# Bayesian Calibration of Antarctic Ice Sheet Ensemble — Project Spec

This document specifies a Bayesian weighting workflow to calibrate a perturbed-parameter
ensemble of Antarctic ice sheet model runs against the IMBIE observational record
(2007–2021), using a per-year emulator, and to produce calibrated sea level
contribution (SLC) projections to 2300.

It is written to be handed to Claude Code as a working brief. Sections marked
**[DECISION NEEDED]** require a choice or a value before or during implementation.
Sections marked **[CHECK]** are diagnostics to run and inspect, not assume.

---

## 1. Scientific setup and key decisions already made

These decisions came out of a prior design discussion and should be treated as fixed
unless the diagnostics below contradict them.

1. **Aggregate the obs target over time.** Calibrate against the *cumulative*
   2007–2021 SLC (a scalar), not individual annual values. This sidesteps the
   temporal autocorrelation in the IMBIE timeseries. See §3.

2. **Do not calibrate the forcing GCM.** A GCM's fit to the 2007–2021 window reflects
   the realised internal variability of that 15-year period, not its multi-century
   predictive skill. Calibrating GCMs against this window would be unjustified.
   Instead, **stratify by GCM**: compute parameter weights independently within each
   GCM, and treat the prior over GCMs as fixed (not data-updated). See §4.

3. **Pool across scenarios when weighting.** The 2007–2021 SLC distributions are
   empirically indistinguishable across scenarios (KDEs overlap). The data cannot
   discriminate between scenarios over this window, so pooling introduces no bias and
   improves effective sample size. Weights therefore depend on **GCM only**. See §4.

4. **Calibrate only the 6 perturbed physical parameters.** These are the genuinely
   free, uncertain dimensions of the prior. All other emulator inputs are either
   encodings of the stratification variables or forcing-determined covariates, and
   must NOT be sampled freely. See §5.

5. **Project by (GCM, scenario).** Although weighting is scenario-agnostic, the
   projections diverge by scenario, so final outputs are presented per
   (GCM, scenario) cell, using the per-GCM parameter weights. See §7.

---

## 2. Inputs and data contract

Implementations should not assume column names; confirm against the actual files first.
The expected logical schema is:

**Ensemble SLC dataset** (long format)
- `run_id` — unique identifier per ensemble member
- `year` — integer, 2007–2300
- `slc` — sea level contribution. **[DECISION NEEDED]** confirm whether this is
  *cumulative* (from a 2007 baseline) or *annual increment*. The aggregation in §3
  differs accordingly.
- `scenario` — forcing scenario label
- `gcm` — forcing climate model label
- `p1 … p6` — the 6 perturbed parameter values for that run

**IMBIE observations**
- `year` — integer, 2007–2021
- `slc_obs` — observed SLC (same units, sign convention, and baseline as the ensemble)
- `slc_sd` — reported uncertainty (standard deviation)
- **[CHECK]** confirm whether IMBIE's reported uncertainty is for annual values or for
  the cumulative quantity, and whether it already accounts for systematic-error
  correlation. This determines how `sigma_obs_cum` is computed in §3.

**Emulators**
- One emulator per year, 2007–2300 (≈294 emulators).
- Inputs include: the 6 parameters, one or more **GCM one-hot / binary switches**, and
  a **2300 ocean thermal forcing anomaly** `T2300` (and possibly others — enumerate them).
- **[DECISION NEEDED]** record the exact input vector each emulator expects, in order,
  and which library produced them (e.g. `DiceKriging`, `RobustGaSP`, `mlegp`, a GP in
  `GPfit`, etc.), since the predict interface differs.

**Units convention.** Fix one sign convention (recommend: mass loss → positive SLC) and
one baseline (recommend: 2007 = 0). If converting IMBIE mass (Gt) to SLE, use
1 Gt ≈ 1/361.8 mm SLE. **[CHECK]** that ensemble and IMBIE agree on all three.

---

## 3. Aggregating the observational target

Collapse the 15-year window to a scalar.

If `slc` is cumulative:
```
Y   = slc_obs at 2021  (minus slc_obs at 2007 if baseline isn't already 2007)
M_i = emulated/observed cumulative slc at 2021 for run i
```

If `slc` is an annual increment:
```
Y   = sum_{t=2007}^{2021} slc_obs(t)
M_i = sum_{t=2007}^{2021} slc_i(t)
```

**Aggregate observational uncertainty.** Do NOT use sqrt(sum of annual variances) — that
assumes independence and will be overconfident given positive autocorrelation.
- Preferred: use IMBIE's own reported cumulative uncertainty if available.
- Otherwise: build the annual obs covariance with an AR(1) correlation and aggregate:
  `sigma_obs_cum^2 = 1^T Sigma_obs 1` where `Sigma_obs[t,t'] = sigma_t sigma_t' rho^|t-t'|`.
  **[DECISION NEEDED]** value of `rho` (try 0.5 / 0.7 / 0.9 in sensitivity).

**Total likelihood variance** (the scalar that controls weight sharpness):
```
sigma_total^2 = sigma_obs_cum^2 + sigma_discrep^2 + sigma_internal^2 + sigma_emulator^2
```
- `sigma_discrep` — structural model–reality discrepancy. **[DECISION NEEDED]**; start
  from a fraction (e.g. 30%) of the inter-run spread in `M_i`, and treat as the primary
  sensitivity lever.
- `sigma_internal` — internal climate variability (one GCM realisation vs. one Earth).
  **[DECISION NEEDED]**; estimate from a control run or inter-member spread if possible.
- `sigma_emulator` — emulator predictive sd at the evaluation point (see §6); this is
  per-prediction, so in practice it enters per-sample rather than as a single constant.

> The dominant methodological risk in this whole pipeline is setting `sigma_total` too
> small, which collapses all weight onto one or two runs. §6.3 ESS diagnostics catch this.

---

## 4. The weighting scheme (GCM-stratified, scenario-pooled)

For each GCM stratum `c`, over members/samples `i` in that stratum:

```
log L_i = -0.5 * (Y - M_i)^2 / sigma_total^2
w_i = softmax over the stratum: exp(log L_i - max_j log L_j) / sum_j exp(log L_j - max_j log L_j)
```

Weights are normalised **within each GCM stratum** and sum to 1 per stratum. Scenario is
NOT a grouping variable here — members sharing the same parameters but differing in
scenario receive (near-)identical weights, which is the intended behaviour.

Always work in log-space and subtract the per-stratum max before exponentiating.

**[CHECK] Balanced design.** Confirm the parameter sample is shared (or at least
overlapping) across scenarios within each GCM. If the Latin hypercube was drawn
separately per scenario with differing coverage, pooling is slightly biased — flag it.

---

## 5. Handling the "made-up" emulator inputs

These are emulator *inputs* but NOT free parameters. They must never be sampled
independently over their marginal ranges — doing so creates physically impossible input
combinations and pushes the emulator off the training manifold into confident nonsense.

**GCM binary/one-hot switches (Category A).** Pure encodings of GCM identity. Within a
GCM stratum, fix the switch to that GCM (set its one-hot to 1, others to 0). Never
perturb.

**`T2300` ocean thermal forcing anomaly (Category B).** A forcing-determined *response*,
not a knob. Established facts:
- It is **purely forcing-determined** (does NOT depend on the 6 parameters).
- It **does depend on (GCM, scenario)**.

Therefore it is an **exact lookup** keyed on (GCM, scenario):
```
T2300_lookup = distinct (gcm, scenario, T2300) rows from the ensemble design
```
No secondary emulator is needed (that would only be required if `T2300` varied with the
parameters). **[CHECK]** confirm `T2300` is single-valued per (gcm, scenario) cell.

**Crucial asymmetry between weighting and projection:**
- In **weighting** (an obs-period quantity), `T2300` is near-inert: because scenarios are
  indistinguishable in 2007–2021, the obs-period emulators are insensitive to it.
  **[CHECK]** inspect the obs-period emulators' length-scale / variable importance for
  `T2300` to confirm low sensitivity; if confirmed, any consistent value may be used and
  scenario need not enter the weight calculation at all.
- In **projection** (future years), `T2300` is highly influential and MUST be set to the
  scenario-correct lookup value for the cell being projected.

---

## 6. Emulator-based calibration loop (per GCM stratum)

Direct weighting of only the existing ensemble members often gives a small ESS. Using the
emulator to sample densely from the parameter prior produces a smoother parameter
posterior. Workflow per GCM stratum:

### 6.1 Dense prior sampling
- Draw `K` samples (e.g. K = 10,000) of `(p1…p6)` from the prior over the 6 parameters.
  **[DECISION NEEDED]** the prior form — match how the original ensemble was designed
  (e.g. uniform over the perturbation ranges, or a specified density). Use the SAME prior
  ranges the ensemble was built on; do not extrapolate beyond the emulator's training box.
- Build each sample's full emulator input vector: the 6 sampled params + this GCM's
  switch + a consistent `T2300` (any consistent value, since obs-period is insensitive —
  see §5).

### 6.2 Likelihood via the emulator
- Predict the (aggregated) 2007–2021 `M` for each sample using the obs-period emulator(s).
  - If emulating the cumulative-2021 value directly, one emulator call suffices.
  - If summing annual emulators, propagate predictive variance across years
    (with the same caution about temporal correlation — but since the target is a sum,
    the variance of the sum needs the cross-year covariance, which per-year-independent
    emulators do not give. **[DECISION NEEDED]** simplest robust route: train/evaluate a
    single emulator on the *aggregate* 2007–2021 quantity rather than summing 15 emulators).
- Use the emulator predictive mean as `M_i` and fold its predictive variance into
  `sigma_total^2` per sample.
- Compute weights per §4.

### 6.3 [CHECK] Diagnostics — run before trusting any result
- **Per-stratum ESS** (Kish): `ESS = 1 / sum(w_i^2)`. Report `ESS / K` per GCM. If a
  stratum's ESS is tiny (e.g. < a few % of K, or only 1–2 effective members in the raw
  ensemble), the conditional posterior is degenerate — say so explicitly and revisit
  `sigma_total` (most likely `sigma_discrep` is too small).
- **Coverage check.** Confirm the ensemble/sample `M` actually brackets the observed `Y`.
  If `Y` lies in the tail or outside the ensemble range, calibration is extrapolation and
  weights are meaningless — flag loudly.
- **Weight rank plot** per stratum (sorted weights) to eyeball collapse.
- **Posterior parameter marginals** vs. prior, per stratum, to see what IMBIE constrains.

---

## 7. Projection to 2300 (per GCM × scenario)

For each (GCM, scenario) cell and each future year `t`:
1. Take the parameter posterior for that GCM (the §6 weighted samples — weights depend on
   GCM only and are reused across scenarios).
2. Build inputs: posterior parameter samples + GCM switch + **scenario-correct `T2300`**
   from the lookup.
3. Predict `slc` at year `t` with that year's emulator.
4. Summarise with weighted statistics: weighted mean and weighted quantiles
   (e.g. 5/17/50/83/95%). Use a proper weighted-quantile routine (interpolate the
   weighted ECDF), not a naive sort.
5. Also compute the **prior** (unweighted) summary for comparison — the prior-vs-posterior
   band plot is the single most informative output.

Final deliverables:
- Calibrated SLC timeseries (median + bands) per (GCM, scenario).
- Optional marginal over GCMs using an **explicit, non-data-driven** GCM prior
  (e.g. uniform) — keep this separate and clearly labelled, since GCM weighting is a
  judgement call, not an inference.

---

## 8. Sensitivity analyses (required, not optional)

Re-run the pipeline varying:
1. `sigma_discrep` magnitude (e.g. 10 / 30 / 50% of ensemble spread) — dominant lever.
2. AR(1) `rho` for the obs covariance (0 / 0.5 / 0.9).
3. Likelihood form: Gaussian vs. Student-t (robustness to outlying runs).
4. Leave-one-year-out within 2007–2021 (only relevant if not fully aggregating).

If 2100 / 2300 posterior bands are stable across these, the result is robust. If they
swing, document which assumption drives it — that is itself a finding.

---

## 9. Suggested repository structure

```
.
├── R/
│   ├── 01_load_data.R         # read ensemble + IMBIE, enforce §2 data contract
│   ├── 02_aggregate_obs.R     # §3: build Y, M, sigma_total components
│   ├── 03_emulator_io.R       # input-vector builder, T2300 lookup, switch encoding (§5)
│   ├── 04_weights.R           # §4 + §6: dense sampling, likelihood, per-GCM weights
│   ├── 05_diagnostics.R       # §6.3: ESS, coverage, rank plots, posterior marginals
│   ├── 06_project.R           # §7: forward projection per (GCM, scenario)
│   └── 07_sensitivity.R       # §8
├── config.yml                 # all [DECISION NEEDED] values in one place
└── outputs/                   # figures + calibrated projection tables
```

Put every **[DECISION NEEDED]** value in `config.yml` so sensitivity runs only change
config, never code.

---

## 10. Open items to resolve with Claude Code (the [DECISION NEEDED] / [CHECK] list)

- [ ] Confirm `slc` is cumulative vs. annual increment, and unit/sign/baseline agreement.
- [ ] Confirm whether IMBIE uncertainty is annual or cumulative; set `sigma_obs_cum`.
- [ ] Enumerate exact emulator input vector(s) and the GP library / predict interface.
- [ ] Decide: aggregate-quantity emulator vs. summing 15 annual emulators (§6.2).
- [ ] Set `rho`, `sigma_discrep`, `sigma_internal` starting values.
- [ ] Confirm `T2300` is single-valued per (gcm, scenario) and obs-period-inert.
- [ ] Confirm balanced parameter design across scenarios.
- [ ] Specify the parameter prior (form + ranges) matching the ensemble design.
- [ ] Confirm `Y` lies inside the ensemble `M` range (coverage).
