# 00 — shared setup: libraries, config, output dir, normalisation helpers.
# Sourced by every stage; run the whole pipeline via ../run_all.R.

suppressPackageStartupMessages({
  library(RobustGaSP)
  library(dplyr)
  library(yaml)
})

if (!exists("CONFIG_PATH")) CONFIG_PATH <- "config.yml"
cfg <- yaml::read_yaml(CONFIG_PATH)

out_dir <- cfg$paths$outputs
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Min-max normalisation to [0,1], optionally on a log10 scale. norm_fit() records
# the transform from the training data so it can be reapplied (norm_apply) or
# inverted (norm_invert) for dense samples and (gcm, scenario) lookups.
norm_fit <- function(x, log = FALSE) {
  if (log) { stopifnot(all(x > 0)); x <- log10(x) }
  list(lo = min(x), hi = max(x), log = log)
}
norm_apply <- function(x, spec) {
  if (spec$log) { stopifnot(all(x > 0)); x <- log10(x) }
  (x - spec$lo) / (spec$hi - spec$lo)
}
norm_invert <- function(z, spec) {
  x <- z * (spec$hi - spec$lo) + spec$lo
  if (spec$log) x <- 10^x
  x
}

# Weighted quantiles via linear interpolation of the weighted ECDF (Hazen-style
# plotting positions). Use for posterior bands — never a naive sort (spec §7).
weighted_quantile <- function(x, w, probs) {
  o  <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  stats::approx(cw, x, xout = probs, rule = 2)$y
}

# Per-sample log-likelihood of a residual (obs - model) under a Gaussian or scaled
# Student-t error model, with per-sample sigma. Additive constants that cancel in a
# within-stratum softmax are dropped; sigma-dependent terms are kept (sigma varies
# sample to sample once the emulator variance is folded in). df used only for t.
loglik_resid <- function(resid, sigma, likelihood = "gaussian", df = 4) {
  z <- resid / sigma
  if (identical(likelihood, "student_t"))
    -log(sigma) + stats::dt(z, df = df, log = TRUE)
  else
    -log(sigma) - 0.5 * z^2
}
