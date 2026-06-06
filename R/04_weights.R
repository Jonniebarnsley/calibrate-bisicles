# 04 - Weights: Don't calibrate by gcm, since 2007-2021 period is so dominated
# interdecadal variability. Instead, calibrate each gcm one at a time and
# combine later using a uniform prior over gcm.

if (!exists("predict_slc_mm")) source("R/03_emulate.R")

set.seed(cfg$weighting$seed)
K        <- cfg$weighting$n_samples
ref_scen <- cfg$weighting$obs_ref_scenario

priors <- cfg$priors

# Warn if any configured prior extends beyond the emulator's training box:
# samples there are extrapolation off the training manifold (spec §5/§6).
for (p in param_cols) {
  z <- norm_apply(as.numeric(priors[[p]]), norm_specs[[p]])
  if (min(z) < -0.02 || max(z) > 1.02)
    warning(sprintf(
      "%s prior maps to normalised [%.2f, %.2f] — outside emulator training.",
      p, min(z), max(z)
    ))
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

# Per-basin emulator prediction + log-likelihood contribution for one set of
# normalised inputs X. Returns list(pr, logL). Per-sample total sigma is
# heteroscedastic (obs + model discrepancy + emulator uncertainty), so
# loglik_resid keeps the sigma-dependent normalisation (dropped in the spec's
# softmax shorthand but it matters when sigma varies).
basin_logL <- function(i, X) {
  pr        <- predict_slc_mm_list[[i]](X, cfg$obs$target_year)
  sigma_mod <- cfg$weighting$sigma_mod_mult * sigma_obs_mm_list[[i]]
  sigma     <- sqrt(sigma_obs_mm_list[[i]]^2 + sigma_mod^2 + pr$sd^2)
  logL      <- loglik_resid(Y_mm_list[[i]] - pr$mean, sigma,
                            cfg$weighting$likelihood, cfg$weighting$student_t_df)
  list(pr = pr, logL = logL)
}

weight_stratum <- function(gcm) {
  params <- sample_prior(K)
  X      <- build_norm_inputs(params, gcm, ref_scen)

  # Joint log-likelihood under independent basin errors: sum log-likelihoods.
  parts <- lapply(seq_len(n_basins), function(i) basin_logL(i, X))
  logL  <- Reduce("+", lapply(parts, `[[`, "logL"))

  # Per-basin emulator predictions stored alongside the joint logL.
  per_basin <- do.call(cbind, lapply(seq_len(n_basins), function(i) {
    setNames(data.frame(parts[[i]]$pr$mean, parts[[i]]$pr$sd),
             c(paste0("M_mm_",      basin_labels[i]),
               paste0("emu_sd_mm_", basin_labels[i])))
  }))

  w <- exp(logL - max(logL))
  w <- w / sum(w)

  cbind(data.frame(gcm = gcm), params, per_basin,
        data.frame(logL = logL, weight = w))
}

weights <- do.call(rbind, lapply(gcms, weight_stratum))
write.csv(weights, file.path(out_dir, "weights.csv"), row.names = FALSE)
message(sprintf(
  "04: drew %d samples x %d GCM strata; wrote weights.csv",
  K, length(gcms)
))
