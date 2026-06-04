# 04 — dense prior sampling, emulator likelihood, per-GCM weights (spec §4, §6).
# Weights are normalised within each GCM stratum (scenario is pooled / not a
# grouping variable). The dense prior is uniform on the normalised box [0,1]^6,
# which by construction equals the ensemble's LHS prior (log-uniform for log_cols,
# uniform otherwise).

if (!exists("predict_slc_mm")) source("R/03_emulator.R")

set.seed(cfg$weighting$seed)
K            <- cfg$weighting$n_samples
ref_scen     <- cfg$weighting$obs_ref_scenario
sigma_mod_mm <- cfg$weighting$sigma_mod_mult * sigma_obs_mm   # model discrepancy

priors <- cfg$priors

# Warn if any configured prior extends beyond the emulator's training box: samples
# there are extrapolation off the training manifold (spec §5/§6).
for (p in param_cols) {
  z <- norm_apply(as.numeric(priors[[p]]), norm_specs[[p]])
  if (min(z) < -0.02 || max(z) > 1.02)
    warning(sprintf("[PRIOR] %s prior maps to normalised [%.2f, %.2f] — outside emulator training box.",
                    p, min(z), max(z)))
}

# Draw K samples per parameter from its configured prior in raw units
# (log-uniform for log_cols, uniform otherwise). build_norm_inputs handles the
# normalisation into the emulator's training space.
sample_prior <- function(K) {
  out <- as.data.frame(lapply(param_cols, function(p) {
    rng <- as.numeric(priors[[p]])
    if (p %in% log_cols) 10^runif(K, log10(rng[1]), log10(rng[2]))
    else                 runif(K, rng[1], rng[2])
  }))
  names(out) <- param_cols
  out
}

weight_stratum <- function(gcm) {
  params <- sample_prior(K)
  X      <- build_norm_inputs(params, gcm, ref_scen)
  pr     <- predict_slc_mm(X, cfg$obs$target_year)

  # Per-sample total sigma: obs + model discrepancy + emulator predictive variance.
  # Heteroscedastic, so loglik_resid keeps the sigma-dependent normalisation term
  # (dropped in the spec's softmax shorthand but it matters when sigma varies).
  sigma <- sqrt(sigma_obs_mm^2 + sigma_mod_mm^2 + pr$sd^2)
  logL  <- loglik_resid(Y_mm - pr$mean, sigma,
                        cfg$weighting$likelihood, cfg$weighting$student_t_df)

  w <- exp(logL - max(logL))
  w <- w / sum(w)

  cbind(data.frame(gcm = gcm), params,
        data.frame(M_mm = pr$mean, emu_sd_mm = pr$sd, logL = logL, weight = w))
}

weights <- do.call(rbind, lapply(gcms, weight_stratum))
write.csv(weights, file.path(out_dir, "weights.csv"), row.names = FALSE)
message(sprintf("04: drew %d samples x %d GCM strata; wrote weights.csv", K, length(gcms)))
