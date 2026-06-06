# 00 — shared setup: libraries, config, output dir, normalisation helpers.
# Sourced by every stage; run the whole pipeline via ../run_all.R.

suppressPackageStartupMessages({
  library(RobustGaSP)
  library(dplyr)
  library(yaml)
})

# If a config path has not already been defined (e.g. by the
# run_diagnostics.R wrapper), use the default path.
if (!exists("config_path")) config_path <- "config.yml"
cfg <- yaml::read_yaml(config_path)

out_dir <- cfg$paths$outputs
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Helper functions ---

# Min-max normalisation to [0,1], optionally on a log10 scale. norm_fit()
# records the transform from the training data so norm_apply() can reapply it
# to fresh samples or to (gcm, scenario) forcing lookups.
norm_fit <- function(x, log = FALSE) {
  if (log) {
    stopifnot(all(x > 0))
    x <- log10(x)
  }
  list(lo = min(x), hi = max(x), log = log)
}
norm_apply <- function(x, spec) {
  if (spec$log) {
    stopifnot(all(x > 0))
    x <- log10(x)
  }
  (x - spec$lo) / (spec$hi - spec$lo)
}

# Weighted quantiles via linear interpolation of the weighted empirical
# Cumulative Density Function. Use for posterior bands.
weighted_quantile <- function(x, w, probs) {
  o  <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  stats::approx(cw, x, xout = probs, rule = 2)$y
}

# Per-sample log-likelihood of a residual (obs - model) under a Gaussian or
# scaled Student-t error model, with per-sample sigma. (sigma varies sample
# to sample once the emulator variance is folded in). df used only for t.
loglik_resid <- function(resid, sigma, likelihood = "gaussian", df = 4) {
  z <- resid / sigma
  if (identical(likelihood, "student_t"))
    -log(sigma) + stats::dt(z, df = df, log = TRUE)
  else
    -log(sigma) - 0.5 * z^2
}