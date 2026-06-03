# Variance Decomposition Analysis: Implementation Guidance

## Overview

This document describes the statistical analysis to be implemented in R. The goal is to decompose the sources of uncertainty in an Antarctic ice sheet ensemble into contributions from discrete factors (GCMs, emissions scenarios) and continuous parameters, following the approach of Seroussi et al. (2023) and Coulon et al. (2025). The output should show how the dominant drivers of ensemble spread evolve over time.

---

## Ensemble Structure

- **128 ensemble members per scenario**, sampled using a **Sobol' sequence** across 6 continuous parameters
- **4 GCMs**: each contributes 32 simulations per scenario, except MIROC which contributes 31 (treat as balanced — the difference is negligible)
- **3 scenarios**: all GCM–scenario combinations have been run
- **Total**: ~384 simulations per timestep across all scenarios
- **Emulator**: a single emulator trained jointly across all scenarios; can generate predictions at arbitrary parameter combinations at negligible cost

---

## Analysis Structure

Two complementary analyses are required:

1. **Combined analysis**: all scenarios together — quantifies how much variance is simply attributable to scenario choice versus within-scenario uncertainty
2. **Per-scenario analysis**: repeat for each scenario separately — quantifies what drives uncertainty conditional on a given emissions pathway

Both analyses should be repeated at each model timestep to produce temporal evolution of uncertainty sources.

---

## Part 1: Sobol' Sensitivity Indices (Continuous Parameters)

### Why Sobol' indices rather than ANOVA for continuous parameters

The 6 continuous parameters were sampled using a Sobol' sequence, which is the appropriate design for Sobol' sensitivity analysis. Unlike ANOVA regression terms, Sobol' indices make no assumption of linearity and capture the full variance contribution of each parameter including nonlinear effects and interactions.

Two types of index are computed:
- **First-order index** $S_i$: fraction of total variance explained by parameter $i$ alone
- **Total-order index** $S_{Ti}$: fraction explained by parameter $i$ including all its interactions with other parameters

The difference $S_{Ti} - S_i$ indicates how much of a parameter's influence comes through interactions with other parameters rather than its direct effect.

### Estimator

Use the **Saltelli (2010) estimator**, implemented in the R `sensitivity` package as `sobolSalt`. This estimator requires samples arranged in a specific paired design: two independent base matrices A and B (each of dimension $N \times k$, where $k = 6$), plus $k$ additional matrices where column $i$ of A is replaced by column $i$ of B.

Total emulator evaluations required: $N(k + 2) = N \times 8$

### Generating the paired design

The original ensemble was not structured as a Saltelli paired design, so new samples must be generated via the emulator. Steps:

1. Generate two independent $N \times 6$ Sobol' base matrices A and B (use `randtoolbox::sobol` or `sensitivity` built-ins with scrambling enabled)
2. Pass to `sobolSalt(model = NULL, X1 = A, X2 = B, scheme = "A")` to construct the full design matrix
3. Feed the design matrix through the emulator to get predictions
4. Call `tell(sobol_object, y = predictions)` to compute indices
5. Extract `$S` (first-order) and `$T` (total-order) indices

### Sample size

Use $N = 10{,}000$ as the base sample size (80,000 total emulator evaluations). Run a convergence check by computing indices at $N = 1000, 2000, 5000, 10000$ and confirming stabilisation. If indices have not stabilised by $N = 10{,}000$, increase to $N = 50{,}000$.

### Scenario handling for Sobol' analysis

Since the emulator was trained on all scenarios jointly, scenario should be treated as a discrete input variable when generating the Sobol' design. Compute:

- **Combined Sobol' indices**: pass scenario as an additional discrete column alongside the 6 continuous parameters (or run separately for each scenario and compare)
- **Per-scenario Sobol' indices**: fix the scenario input at each scenario value in turn, vary only the 6 continuous parameters, and compute indices separately for each scenario

The per-scenario indices are more interpretable and should be the primary output.

### Temporal evolution

Repeat the Sobol' index computation at each model timestep. This produces time series of $S_i(t)$ and $S_{Ti}(t)$ for each parameter, showing how parameter sensitivity evolves. This is computationally cheap since all evaluations go through the emulator.

---

## Part 2: ANOVA for Discrete Factors (GCM, Scenario)

### Purpose

ANOVA quantifies what fraction of ensemble variance is attributable to:
- GCM choice (main effect)
- Scenario choice (main effect, combined analysis only)
- GCM × Scenario interaction

### Why Type III sums of squares

Use **Type III (marginal) sums of squares** via `car::Anova(type = 3)`. This assesses each factor after accounting for all other terms in the model, which is the correct approach when interactions are present. Type I (sequential) sums of squares are order-dependent and less appropriate here.

The near-imbalance (MIROC has 31 members per scenario rather than 32) makes Type III preferable over Type I.

### Model specification

**Combined analysis** (all scenarios):

```r
model_combined <- lm(slr ~ gcm * scenario, data = df_year)
anova_combined <- car::Anova(model_combined, type = 3)
```

The interaction term `gcm * scenario` expands to main effects plus their interaction. This captures whether certain GCMs behave differently across scenarios — scientifically expected and should be included.

**Per-scenario analysis**:

```r
model_scenario <- lm(slr ~ gcm, data = df_year_scenario)
anova_scenario <- car::Anova(model_scenario, type = 3)
```

Only GCM is a factor within a single scenario.

### Extracting explained variance fractions

Divide each term's sum of squares by the total sum of squares (including residuals) to get the fraction of variance explained:

```r
ss <- anova_combined$`Sum Sq`
names(ss) <- rownames(anova_combined)
fractions <- ss / sum(ss)
```

The residual fraction represents variance not explained by GCM or scenario — this includes both continuous parameter variance and any unexplained noise.

### Temporal loop

Wrap the ANOVA in a loop over timesteps:

```r
timesteps <- unique(df$year)

anova_results <- map(timesteps, function(t) {
  df_t <- filter(df, year == t)
  model <- lm(slr ~ gcm * scenario, data = df_t)
  anova_table <- car::Anova(model, type = 3)
  ss <- anova_table$`Sum Sq`
  names(ss) <- rownames(anova_table)
  ss / sum(ss)
})

anova_df <- bind_rows(anova_results) %>%
  mutate(year = timesteps)
```

---

## Part 3: Combining the Two Analyses

The ANOVA and Sobol' analyses are complementary and should be combined into a single variance decomposition picture.

### Attribution of total variance

At each timestep, the total ensemble variance can be approximately partitioned as:

| Source | Method | Notes |
|---|---|---|
| Scenario | ANOVA main effect | Combined analysis only |
| GCM | ANOVA main effect | Both analyses |
| GCM × Scenario | ANOVA interaction | Combined analysis only |
| Continuous parameters (individual) | Sobol' first-order indices | Per-scenario analyses |
| Parameter interactions | Sobol' $S_{Ti} - S_i$ | Per-scenario analyses |
| Residual | Remainder | Unexplained variance |

Note that the Sobol' indices and ANOVA fractions operate on different slices of the ensemble, so they cannot simply be added together without care. The cleanest approach is:

1. Use ANOVA fractions to characterise discrete-factor contributions
2. Use per-scenario Sobol' indices to characterise within-scenario parameter contributions
3. Scale the Sobol' fractions by the within-scenario residual fraction from ANOVA to express them in terms of total ensemble variance

### Grouping parameters

Following Coulon et al. (2025), consider grouping the 6 continuous parameters into conceptual categories (e.g., ocean-related, atmosphere-related) and reporting group-level Sobol' indices alongside individual ones. Group-level first-order indices can be computed by summing individual first-order indices within a group (valid if parameters within a group are uncorrelated, which Sobol' sampling approximately ensures).

---

## Visualisation

The standard output for this type of analysis is a **stacked area chart** of explained variance fractions over time, one panel per scenario plus one combined panel. This directly mirrors the presentation in Seroussi et al. (2023, their Figures 6 and 8) and Coulon et al. (2025, their Figures 4 and 5).

### Suggested plot structure

- **x-axis**: time (model years)
- **y-axis**: fraction of total variance (0 to 1)
- **Stacked areas**: one colour per uncertainty source (GCM, scenario, each parameter or parameter group, residual)
- **Panels**: one for each scenario (per-scenario analysis) plus one combined

Use consistent colours across panels. The residual (unexplained) variance should be shown in a neutral colour at the top of the stack.

### Suggested R packages

- `ggplot2` for plotting
- `tidyr::pivot_longer` to reshape the results dataframe for `ggplot`
- `patchwork` or `facet_wrap` for multi-panel layout

---

## R Package Dependencies

| Package | Purpose |
|---|---|
| `sensitivity` | Sobol' index computation (`sobolSalt`, `tell`) |
| `randtoolbox` | Sobol' sequence generation |
| `car` | Type III ANOVA (`Anova`) |
| `tidyverse` | Data manipulation and plotting |
| `patchwork` | Multi-panel figure layout |

---

## Key References

- **Seroussi et al. (2023)**: ISMIP6 ensemble uncertainty decomposition using ANOVA across 13 ice sheet models, 9 GCMs, two melt parameterisation families. Uses a Gaussian process emulator to fill missing ensemble combinations. Decomposes variance into ice model, climate model, and ice–climate interaction terms at continental scale and for 11 individual glaciers. *The Cryosphere*, 17, 5197–5217.

- **Coulon et al. (2025)**: Two ice sheet models (Kori-ULB and PISM) with Latin Hypercube parameter sampling and Bayesian calibration against IMBIE observations. ANOVA with Type I sequential sums of squares (ice sheet model listed first). Separates ocean-related and atmosphere-related parametric uncertainty explicitly. Extends to year 2300 and shows shift in dominant uncertainty source from structural (ice model) to atmospheric parameters under SSP5-8.5. *Nature Communications*, 16, 10385.

- **Saltelli et al. (2010)**: Variance-based sensitivity analysis of model output — design and estimator for the total sensitivity index. *Computer Physics Communications*, 181, 259–270.
