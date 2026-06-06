# 07 — sensitivity of the calibrated projections to the weighting assumptions (§8).
# One-at-a-time sweeps around the baseline along two axes: the discrepancy
# multiplier (dominant) and the likelihood form. The obs-period predictions
# (M, sd) and the future per-sample predictions are FIXED across these settings
# — only the weights change — so they are computed once and re-weighted, making
# the sweep cheap. Leave-one-year-out (spec axis 4) is N/A: the target is a
# single aggregated 2021 scalar. sigma_obs is fixed at the derived value from
# IMBIE (02), no longer a sweep axis now that the difference-of-squares formula
# is settled.

if (!exists("weights", inherits = FALSE)) source("R/04_weights.R")

sens   <- cfg$sensitivity
ryears <- sens$report_years
probs  <- c(0.05, 0.17, 0.50, 0.83, 0.95)
df_t   <- if (is.null(cfg$weighting$student_t_df)) 4 else cfg$weighting$student_t_df

# --- Fixed per-GCM quantities (independent of the weighting setting) ---
# Per-basin obs-period M & sd per sample (from 04, used to recompute the joint
# logL under each setting) + future TOTAL SLC per (report-year, scenario).
fixed <- setNames(lapply(gcms, function(g) {
  wg  <- weights[weights$gcm == g, ]
  fut <- list()
  for (s in scenarios) {
    X   <- build_norm_inputs(wg[, param_cols], g, s)
    Msy <- Reduce("+", lapply(predict_slc_mm_list,
                              function(pred) pred(X, ryears)$mean))   # K x length(ryears)
    for (j in seq_along(ryears)) fut[[paste(ryears[j], s)]] <- Msy[, j]
  }
  per_basin <- lapply(basin_labels, function(bl) list(
    M  = wg[[paste0("M_mm_",      bl)]],
    sd = wg[[paste0("emu_sd_mm_", bl)]]
  ))
  list(per_basin = per_basin, fut = fut)
}), gcms)

# Recompute per-GCM weights for one (mult, likelihood) setting; sigma_obs is fixed.
# Joint log-likelihood under independent basin errors: sum across basins.
setting_weights <- function(g, mult, likelihood) {
  f    <- fixed[[g]]
  logL <- Reduce("+", lapply(seq_along(basin_labels), function(i) {
    pb        <- f$per_basin[[i]]
    sigma_obs <- sigma_obs_mm_list[[i]]
    sigma     <- sqrt(sigma_obs^2 + (mult * sigma_obs)^2 + pb$sd^2)
    loglik_resid(Y_mm_list[[i]] - pb$M, sigma, likelihood, df_t)
  }))
  w <- exp(logL - max(logL)); w / sum(w)
}

# Weighted bands at each (report-year, scenario) per GCM + uniform-GCM marginal.
summarise_setting <- function(axis, value, mult, likelihood) {
  W    <- setNames(lapply(gcms, setting_weights, mult = mult,
                          likelihood = likelihood), gcms)
  rows <- list()
  add  <- function(gcm, s, yr, ess, M, w) {
    q <- weighted_quantile(M, w, probs)
    rows[[length(rows) + 1L]] <<- data.frame(
      axis = axis, value = as.character(value), gcm = gcm, scenario = s, year = yr,
      ess = ess, q05 = q[1], q17 = q[2], q50 = q[3], q83 = q[4], q95 = q[5],
      mean = weighted.mean(M, w))
  }
  for (yr in ryears) for (s in scenarios) {
    Mpool <- numeric(0); wpool <- numeric(0)
    for (g in gcms) {
      M <- fixed[[g]]$fut[[paste(yr, s)]]; w <- W[[g]]
      add(g, s, yr, 1 / sum(w^2), M, w)
      Mpool <- c(Mpool, M); wpool <- c(wpool, w / length(gcms))
    }
    add("ALL (uniform GCM prior)", s, yr, NA, Mpool, wpool)
  }
  do.call(rbind, rows)
}

# --- One-at-a-time sweep (baseline values appear within each axis list) ---
base_mult <- cfg$weighting$sigma_mod_mult
base_lik  <- cfg$weighting$likelihood

res <- list()
for (v in sens$sigma_mod_mult)
  res[[length(res) + 1L]] <- summarise_setting("sigma_mod_mult", v, v, base_lik)
for (v in sens$likelihood)
  res[[length(res) + 1L]] <- summarise_setting("likelihood", v, base_mult, v)

sensitivity <- do.call(rbind, res)
write.csv(sensitivity, file.path(out_dir, "sensitivity.csv"), row.names = FALSE)

# --- Headline table: uniform-GCM marginal, ssp585, posterior median [5,95] ---
hl <- sensitivity[sensitivity$gcm == "ALL (uniform GCM prior)" &
                    sensitivity$scenario == "ssp585", ]
cat("\nSensitivity — uniform-GCM marginal, ssp585, posterior median [5,95] mm:\n")
for (yr in ryears) {
  cat(sprintf("  %d:\n", yr))
  d <- hl[hl$year == yr, ]
  for (i in seq_len(nrow(d)))
    cat(sprintf("    %-16s = %-9s  %5.0f [%5.0f, %5.0f]\n",
                d$axis[i], d$value[i], d$q50[i], d$q05[i], d$q95[i]))
}

# --- Plot: dominant lever (sigma_mod_mult) vs marginal band, per (year, scenario) ---
sweep <- sensitivity[sensitivity$axis == "sigma_mod_mult" &
                       sensitivity$gcm == "ALL (uniform GCM prior)", ]
sweep$value <- as.numeric(sweep$value)
pdf(file.path(out_dir, "sensitivity_mult.pdf"), width = 10, height = 6.5)
par(mfrow = c(length(ryears), length(scenarios)), mar = c(3.5, 3.5, 2, 1),
    mgp = c(2, 0.6, 0))
for (yr in ryears) for (s in scenarios) {
  d <- sweep[sweep$year == yr & sweep$scenario == s, ]; d <- d[order(d$value), ]
  plot(d$value, d$q50, type = "b", pch = 19, log = "x", ylim = range(d$q05, d$q95),
       xlab = "sigma_mod_mult", ylab = "SLC (mm)",
       main = sprintf("%d  %s  (marginal)", yr, s))
  lines(d$value, d$q05, lty = 2, col = "grey40")
  lines(d$value, d$q95, lty = 2, col = "grey40")
}
dev.off()

message("07: sensitivity sweep complete; wrote sensitivity.csv + sensitivity_mult.pdf")
