# Bayesian Calibration of the Antarctic Ice-Sheet Ensemble — Summary

A working summary of the calibration pipeline: what it does, the decisions taken,
the key results, and the open questions. The original brief is in
`ice_sheet_calibration_spec.md`; this document records what was actually built and
learned.

---

## 1. Purpose

Calibrate a perturbed-parameter ensemble (PPE) of BISICLES Antarctic ice-sheet runs
against the IMBIE observational record (cumulative sea-level contribution,
2007–2021) using per-year Gaussian-process emulators, then produce calibrated
sea-level contribution (SLC) projections to 2300, stratified by GCM and scenario.

The approach is a **Bayesian observational calibration**: dense parameter samples
are weighted by their likelihood against the observed SLC, giving a posterior over
the physical parameters that is then propagated forward.

---

## 2. Data contract and key facts established

Resolved against the actual files (spec §10 checklist):

- **Ensemble** (`emulator_inputs_2007-2300_sle.csv`): wide format, one row per run,
  **381 usable runs** = 4 GCMs × 3 scenarios × ~32 LHS draws (3 dropped for NA).
  SLC is **cumulative, in metres, baselined to 2007 = 0, mass-loss positive**.
- **6 free parameters**: `aPhi, n, m, uf, gamma, UMV` (`aPhi, gamma, UMV`
  sampled/normalised on log10).
- **Forcing covariates** `thermal_forcing` and `warming` are **single-valued per
  (GCM, scenario)** — clean lookup tables, not free inputs.
- **Observation target**: IMBIE cumulative MB is already mass-loss-positive, so
  **Y = 7.42 − 2.14 = +5.28 mm** over 2007–2020.92 (the record ends 2020.92, not
  2021). (Uncertainties are stored with a spurious negative sign — magnitudes used.)
- **σ_obs = 1.11 mm**, the uncertainty on the 2007→2021 *difference*:
  √(σ₂₀₂₁² − σ₂₀₀₇²) = √(1.47² − 0.96²). IMBIE builds its cumulative error as a
  root-sum-square of independent monthly errors, so the cumulative variance is
  additive and the difference variance is σ₂₀₂₁² − σ₂₀₀₇². Using σ₂₀₂₁ = 1.47 alone
  would wrongly assume the 2007 baseline is exact; the arithmetic difference
  (1.47 − 0.96 = 0.51) would assume perfectly correlated errors, contradicting the
  RSS construction. No AR(1) construction needed.
- **Balanced design**: within each GCM the three scenarios share the identical
  parameter set, so pooling across scenarios for weighting is unbiased.
- **Coverage**: Y sits at the 84th percentile of the pooled ensemble. Per-GCM it is
  central for MIROC-ES2L, upper-tail for MRI-ESM2-0, and **above the entire range
  for the two low-forcing GCMs** (CESM2-WACCM, EC-Earth3-Veg). Dense emulator
  sampling widens coverage enough that Y falls (just) inside all four strata.

---

## 3. Pipeline architecture

All in R (RobustGaSP), `config.yml`-driven, run via `Rscript run_all.R` from the
`calibration/` directory.

| File | Role |
|---|---|
| `config.yml` | All decisions/levers; sensitivity runs change only this file |
| `R/00_setup.R` | Libraries, config, normalisation + weighted-quantile + likelihood helpers |
| `R/01_load_data.R` | Data contract: dummies, NA drop, baseline=0 check |
| `R/02_aggregate_obs.R` | Scalar target Y, σ_obs, (GCM,scenario) forcing lookup |
| `R/03_emulator_io.R` | Emulator dispatch (`predict_slc_mm`), per-year GP training (cached), input builders |
| `R/svd_emulator.R` | SVD basis-function emulator (used when `svd.enabled`): coefficient GPs + truncation variance + diagnostic CSV (§9) |
| `R/04_weights.R` | Dense prior sampling, likelihood, per-GCM weights |
| `R/05_diagnostics.R` | ESS, coverage, weight-rank + prior/posterior marginal plots |
| `R/06_project.R` | Projection per (GCM × scenario), weighted bands, GCM marginal |
| `R/07_sensitivity.R` | One-at-a-time sensitivity sweeps (spec §8) |
| `R/08_variance_decomposition.R` | Prior-ensemble variance attribution: ANOVA + Sobol' (driver `run_variance.R`) |

Outputs (`outputs/`): `weights.csv`, `projection.csv`, `sensitivity.csv`,
`svd_diagnostic.csv` (SVD mode), `variance_*` (from `run_variance.R`), and
diagnostic/band/sensitivity PDFs. Per-year mode caches annual GPs in
`outputs/emulators/` (first build ~22 min, ~10 min cached); SVD mode builds in ~1–2
min, so `run_all.R` ≈ 5 min and `run_variance.R` ≈ 3 min.

### Method details

- **Emulator** (two modes via `svd.enabled`; **SVD is currently on** — see §9): both
  share `rgasp`, Matérn-5/2, α = 1.9, estimated nugget, linear trend, and the same
  **11 inputs** — 6 parameters + 3 GCM dummies (EC-Earth3-Veg = reference) +
  `thermal_forcing` + `warming`, normalised to [0,1]. *Per-year mode* fits one GP per
  year on that year's cumulative SLC column (no summing of annuals); *SVD mode* (§9)
  emulates SVD coefficients instead. (A leave-one-out test confirmed both forcing
  inputs earn their place and are essential for the projection — see §8.)
- **Priors**: explicit **design priors** set in config (`priors:` block), drawn as
  true uniform (`n, m, uf`) and true log-uniform (`aPhi, gamma, UMV`) over raw-unit
  ranges, then mapped into the emulator's training-normalised space. Independent of
  the ensemble's realised min/max.
- **Weighting**: dense samples per GCM stratum → emulator predicts obs-period M and
  sd → log-likelihood → softmax **within each GCM** (scenario pooled). Uncertainty
  budget: `σ_total² = σ_obs² + σ_mod² + σ_emulator²`, with the per-sample emulator
  variance folded in and the σ-dependent normalisation term retained.
- **Projection**: reuse each GCM's posterior parameter weights, apply the
  **scenario-correct** forcing lookup, predict each future year, summarise with
  interpolated weighted-ECDF quantiles (5/17/50/83/95). Prior (unweighted) bands
  reported alongside; a uniform-GCM-prior marginal is provided, clearly labelled.

---

## 4. Key design decisions and their rationale

1. **Both forcing covariates as emulator inputs.** `thermal_forcing` and `warming`
   are included (not just GCM dummies) because they improve predictive skill and are
   pure (GCM, scenario) lookups. This also makes the emulator able to produce
   scenario divergence in projection.

2. **Single model-discrepancy term**, replacing separate discrepancy +
   internal-variability terms. Defined as a scalar multiple of the observation
   error: **σ_mod = `sigma_mod_mult` × σ_obs**.
   - Tried `sigma_mod_mult = 10` first → σ_mod dwarfs the few-mm spread in M →
     near-uniform weights (ESS ~100%), a near-no-op calibration.
   - Settled on **`sigma_mod_mult = 2`** with σ_obs = 1.11 mm → **σ_mod ≈ 2.22 mm**,
     a modest, live update (ESS 72–93%, no collapse). (An earlier run used mult = 1.5
     with σ_obs = 1.47, giving essentially the same σ_mod.)

3. **Discrepancy is anchored to the observation error, NOT the ensemble spread.**
   The earlier "30% of ensemble spread" heuristic is *circular*: the ensemble spread
   is parametric *signal* (sensitivity), not noise, so using it to set the noise
   floor conflates the two — and per-stratum it would erase the genuine
   identifiability differences between GCMs. Anchoring σ_mod to σ_obs (a noise-side,
   external quantity, constant across GCMs) lets the differing parametric signal play
   out correctly (MIROC, with the largest spread, updates most). Neither choice
   *measures* discrepancy — both are scale-setting judgements — but σ_obs is the
   right family of quantity to anchor to.

4. **Explicit, configurable priors** decoupled from the ensemble's realised range.
   The prior should be the design prior. Marginal plots draw the prior as an analytic
   rectangle (faithful to the uniform draw), not a KDE.

---

## 5. Key results

### Calibration (σ_mod_mult = 2, σ_obs = 1.11 mm)

| GCM | ESS / K | Y percentile in M | Y inside? |
|---|---|---|---|
| CESM2-WACCM | 78% | 100.0 | yes (edge) |
| EC-Earth3-Veg | 79% | 99.9 | yes (edge) |
| MIROC-ES2L | 72% | 32.3 | yes |
| MRI-ESM2-0 | 93% | 85.4 | yes |

Only **γ** is visibly constrained (for the GCMs where Y is informative); the other
parameters' posteriors stay close to their uniform priors.

### Projections (posterior median [5, 95] %, mm SLE)

Uniform-GCM-prior marginal:

| Year | SSP1-26 | SSP5-34-over | SSP5-85 |
|---|---|---|---|
| 2100 | 17 [−27, 92] | 30 [−29, 118] | 41 [−41, 158] |
| 2300 | 73 [−132, 483] | 126 [−142, 628] | 930 [129, 2017] |

Scenarios diverge correctly (overlap near-present, fan out later) — e.g. MIROC-ES2L
posterior median: 2050 = 21/21/20 → 2100 = 46/57/87 → 2300 = 74/127/1073 mm for
SSP1-26 / SSP5-34-over / SSP5-85.

Per-GCM, 2300 SSP5-85 (prior → posterior):

| GCM | prior | posterior |
|---|---|---|
| CESM2-WACCM | 1133 [505, 2336] | 1008 [471, 2085] |
| EC-Earth3-Veg | 1054 [381, 2334] | 1008 [366, 2262] |
| MIROC-ES2L | 1223 [622, 2412] | 1073 [600, 2009] |
| MRI-ESM2-0 | 453 [−77, 1477] | 468 [−76, 1508] |

### Sensitivity (spec §8) — the result is robust

One-at-a-time sweeps of the discrepancy multiplier, σ_obs, and Gaussian-vs-Student-t,
at 2100 and 2300:

- **Median essentially invariant**: 2100 marginal SSP5-85 median fixed at 41 mm
  across all settings; 2300 median moves ~12% (885 → 995 mm) across the full
  multiplier range 1→10.
- **Likelihood form and σ_obs choice negligible** (Student-t ≈ Gaussian → result not
  driven by a few outlying high-likelihood runs).
- The discrepancy lever's main effect is on the **upper tail** (q95 inflates as the
  calibration weakens toward the prior). ESS rises monotonically toward K as the
  lever loosens; MIROC is always the most-constrained; no collapse even at mult = 1.

---

## 6. Why the calibration shifts projections counterintuitively

Although the observation sits *above* most of the ensemble (so the weighting
up-weights high *near-term* members), several GCMs' future bands shift *down*. The
reason is that **near-term and committed (long-term) SLC are decoupled — for CESM2
even anti-correlated — in parameter space** (corr of SLC₂₀₂₁ vs SLC₂₃₀₀ under
SSP5-85):

| GCM | corr(2021, 2300) | γ vs 2021 | γ vs 2300 |
|---|---|---|---|
| CESM2-WACCM | −0.36 | −0.38 | +0.58 |
| EC-Earth3-Veg | −0.02 | +0.45 | +0.49 |
| MIROC-ES2L | +0.50 | +0.88 | +0.61 |
| MRI-ESM2-0 | +0.11 | +0.53 | +0.62 |

For CESM2, γ's effect *flips sign* across timescales: the highest-near-term-loss
members are the low-γ members, which have the *lowest* long-term loss — so
up-weighting them pulls the future down. The future-shift direction in general is
the product of (where Y sits in the near-term distribution) × (sign of the
near-term↔long-term coupling), which is GCM-dependent.

**Lesson**: a 15-year record constrains the *transient* response, while the
parameters that dominate the *committed* response (n, m, γ all ~+0.5 with 2300) are
only weakly — and sometimes oppositely — linked to it. The calibration therefore has
limited and occasionally backwards leverage on the projection. (Caveat: these are
n ≈ 32 correlations; the γ-2300 link is robust, the near-term signs noisier.)

---

## 7. Relationship to Edwards et al. (2021)

| | Edwards et al. 2021 | This pipeline |
|---|---|---|
| Conditioning on obs | **None** — forward propagation of prior parameter distributions | **Bayesian calibration** — likelihood-weighting against IMBIE |
| Ensemble | Multi-model (ISMIP6/GlacierMIP), model dummies + nugget absorb structural spread | Single-model BISICLES PPE, 6 perturbed parameters |
| Climate uncertainty | Continuous global-mean T from a 500-member FaIR ensemble | 4 discrete GCMs (fixed prior) + per-(GCM,scenario) forcing lookups |
| Parameter distributions | Process-derived priors (never updated by data) | LHS prior reweighted by the observation |
| Emulator | RobustGaSP, Matérn-5/2, nugget, per-year | Same library/kernel/structure |

The defining difference is **calibration vs. forward propagation**. Edwards
deliberately avoided calibrating ice-sheet parameters on short records — the same
transient-vs-committed tension surfaced in §6 is precisely why. This pipeline adds
that observational constraint, mitigated by a generous discrepancy term.

---

## 8. Variance decomposition of the prior ensemble

A complementary analysis (`R/08_variance_decomposition.R`, driver `run_variance.R`;
needs the `sensitivity` and `car` packages), **independent of the calibration
weighting** — it attributes the *prior* ensemble's SLC spread to its sources over
time, following Seroussi et al. (2023) and Coulon et al. (2025), reusing the cached
per-year emulators. Two complementary methods, every 5 years:

- **ANOVA (discrete factors)** on the raw ensemble: per-year Type III decomposition of
  `slr ~ gcm * scenario` into GCM / scenario / interaction / residual fractions
  (residual = within-cell parametric variance), plus a per-scenario `slr ~ gcm`. Uses
  sum-to-zero contrasts and drops the intercept row (variance = deviations from the
  mean) — both required for a correct Type III variance partition.
- **Sobol'/Saltelli (continuous)** via the emulator: first-order (**S**) and
  total-order (**T**) indices for the 6 parameters + a forcing-severity axis, run
  **per scenario**. N = 4000 base samples (converged; checked at N = 1000–8000).

### The forcing-input question
`thermal_forcing` and `warming` are correlated at **0.965** across the 12
(GCM, scenario) cells (both spike under SSP5-85). This splits into two answers:

- **For the emulator, keep both.** A leave-one-out test of
  {none / warming / thermal_forcing / both} input sets showed forcing inputs are
  *essential for projection* (2300 LOO normalised-RMSE ≈ 0.87 with GCM dummies only —
  useless — vs ≈ 0.11 with any forcing input); `both` is best at every horizon (a
  modest ~3–10% gain over the best single); and `thermal_forcing` is the stronger
  single input (clearly better at 2100, tied at 2300). In the obs period all configs
  are sub-mm (forcing near-inert), consistent with the weighting relying on GCM only.
- **For the Sobol', collapse to one axis.** Saltelli assumes *independent* inputs; at
  r = 0.965 they cannot be sampled independently without breaking that assumption and
  pushing the emulator off-manifold. So they are reduced to a single
  **forcing-severity axis** (warming, with `thermal_forcing` mapped onto the realised
  manifold).

### GCM represented through forcing intensity (per-scenario Sobol')
Within each scenario the forcing axis spans that scenario's GCMs' warming values (low-
to high-forcing GCM), with the GCM dummy **tracking the nearest GCM** so each GCM's
real behaviour is recovered at its forcing point. So per scenario the spread is
decomposed into {6 parameters + forcing (= GCM) + interactions}.

### S vs T, and the plots
- **S** (first-order) = variance from an input *alone*; Σ S ≤ 1.
- **T** (total-order) = variance involving an input *including all interactions*;
  Σ T ≥ 1; **T − S** = its interaction contribution.
- The stacked-area panels plot the first-order **S** per input, with a grey
  **interactions** slice = 1 − Σ S (not Σ(T − S), which would double-count). T is in
  the CSVs. `variance_decomposition.pdf` = page 1 ANOVA (combined stack + per-scenario
  GCM share), page 2 the three per-scenario Sobol' panels.
- The Sobol' panels apply a centered moving-average (config `sobol_smooth_window`,
  default 5 timesteps ≈ 25 yr; set to 1 to disable) to remove year-to-year jitter
  arising from the independent per-year GP fits. This is a visualisation choice only —
  the CSVs hold the unsmoothed indices.

### Findings
- **ANOVA**: scenario variance ≈ 0 near-term (0.002 at 2021); GCM + parametric
  dominate; **scenario takes over by 2300** (≈ 0.52). GCM's within-scenario share
  declines over time everywhere, most under SSP5-85.
- **Sobol' (per scenario)**:
  - **SSP1-26 and SSP5-34-over**: the **GCM/forcing axis dominates** the within-scenario
    spread (~0.5) — because CESM2 is a low-emissions warming outlier, *which GCM* you
    pick drives most of the spread — with gamma second.
  - **SSP5-85**: a clear shift — GCM/forcing dominates early, but by 2300 the
    **parameters take over** (gamma ≈ 0.32, then m and n) and forcing recedes. Under
    strong forcing the parametric response is amplified.
  - **gamma is the leading parameter throughout** — the same driver the calibration and
    the near-/long-term correlation analysis (§6) flagged.

A first `run_variance.R` takes ~15–20 min (57 years × 3 scenarios at N = 4000;
emulators load from cache, the Sobol' evaluations do not). Outputs:
`variance_anova_combined.csv`, `variance_anova_per_scenario.csv`,
`variance_sobol_by_scenario.csv`, `variance_decomposition.pdf`.

---

## 9. Optional emulator mode: SVD basis-function emulator

Selected by `svd.enabled` (currently **on**). Instead of one GP per year, this is a
Higdon et al. (2008) basis-function emulator: standardise each year, take the SVD
basis, emulate the leading `svd.n_components` (= 5) coefficients as GPs, and
reconstruct SLC(year) = mean + Σ coefficient·basis. All stages run in either mode via
a single `predict_slc_mm(Xnorm, year)` dispatch in `03`; module `R/svd_emulator.R`.

It *is* essentially the Higdon framework (SVD basis, a GP per coefficient, a residual
variance term). We depart in three ways, all suited to this problem: per-year
standardisation (below), a **scalar** discrepancy (§4.2) rather than Higdon's
functional δ(x), and two-stage importance weighting rather than full joint MCMC.

**Predictive variance carries a truncation-variance term**: the discarded components
enter as a per-year residual variance, so the reconstruction's uncertainty is honest
rather than over-confident.

Two problem-specific choices vs Higdon's standard recipe:
- **Per-year (per-output) standardisation.** Cumulative SLC variance spans ~3 orders of
  magnitude (mm at 2021 → m at 2300), so plain single-scale SVD lets the leading PCs
  ignore the early calibration period — its 2021 reconstruction error is ~2× σ_obs.
  Standardising each year fixes it.
- **Deterministic per-year truncation variance** (vs Higdon's single estimated iid
  residual precision λη). Without it the 2021 predictive intervals are over-confident
  (pass ~0.83); with it, ~0.95.

**Choosing the number of components.** Not the 99.9%-variance rule (which suggests 9)
but the end-to-end LOO reconstruction error, which plateaus by r ≈ 3–5 because
*emulation* error then dominates, not truncation. `svd.run_diagnostic` writes
`outputs/svd_diagnostic.csv` (r = 1..9): cumulative variance explained, per-coefficient
LOO R² (emulatability of each mode), and nRMSE / NED(÷sd²) / pass-rate / RMSE-mm at
2021 & 2300. At **r = 5**: 99.4% variance, 2021 reconstruction RMSE 0.63 mm (< σ_obs),
pass 0.95.

**Validation (r = 5).** Reproduces the per-year pipeline closely: ESS 75/80/75/95%
(per-year 78/79/72/93), 2300 SSP5-85 marginal 940 [162, 1975] mm (per-year
930 [129, 2017]), well-calibrated (pass ~0.95). The variance decomposition (§8) is also
qualitatively unchanged — the ANOVA is *identical* (it uses the raw ensemble, not the
emulator), and the Sobol' story holds, with slightly less weight on gamma and more on
n/m at SSP5-85.

**Speed.** ~5 coefficient GPs vs ~294 annual GPs. Both `06` and `08` predict each
(GCM, scenario)'s coefficients once and reconstruct all years; a year-by-year call
re-predicts them ~57× (the mistake that first made the projection hang for 30+ min).

**Trade-off.** SVD is slightly less accurate at the 2021 calibration target than the
direct per-year GP (0.63 vs 0.47 mm LOO) — still well under σ_obs — in exchange for
speed, temporally-coherent trajectories, and smooth time-resolved diagnostics. Revert
with `svd.enabled: false`.

---

## 10. Caveats and limitations

- **Short-record leverage**: see §6 — the constraint is modest and can move
  projections counterintuitively.
- **Coverage at the high edge** for CESM2 / EC-Earth3 (per the agreed stance, this is
  treated as variability, not GCM quality; handled by keeping σ_mod generous).
- **Projection bands carry parametric spread only** — the per-point emulator
  predictive variance is folded into the *weighting* but not into the projection
  bands (Edwards samples it into the PDF; adding it would be needed for a like-for-like
  uncertainty comparison).
- **σ_mod is a judgement, not a measurement**; the result's robustness to it (§5) is
  reassuring but the absolute band widths inherit that assumption.
- **GCM marginal** uses an explicit uniform prior over GCMs — a labelled judgement,
  not an inference (GCMs are not calibrated, by design).
- Marginal-plot KDEs roll off to zero at the prior edges — a kernel boundary
  artifact, not a real posterior feature.
- **Variance decomposition**: the Sobol' forcing axis assumes a *continuous* forcing
  distribution, but the real forcing is the discrete set of (GCM, scenario) values, so
  the forcing index's magnitude depends on that assumption; per scenario the GCM dummy
  tracks the nearest GCM, mildly extrapolating each GCM between forcing points. The
  ANOVA and Sobol' are complementary *views*, not a single merged stack (merging would
  double-count the forcing/scenario contribution).

---

## 11. Possible next steps

- Make σ_mod a relative term (∝ forcing or ∝ |M|) **if** a physical case is made that
  BISICLES is less trustworthy under stronger forcing — but not tied to ensemble
  spread (§4.3).
- Optionally add emulator predictive variance into the projection bands.
- Extend sensitivity to a full grid (currently one-at-a-time) if interactions are of
  interest.
- Optional ocean / ice-dynamics / GIA parameter grouping for the Sobol' indices
  (à la Coulon), once the physical categorisation of the 6 parameters is confirmed.
