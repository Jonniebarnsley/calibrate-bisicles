# 03d — forcing-input tests: 9 emulator configurations differing only in their
# forcing inputs (combinations of GMSTa / tas / thetao / smb / thermal_forcing).
# For each (test, year in {2021, 2300}) trains a fresh per-year RobustGaSP in
# memory (no disk caching), computes LOO statistics, and builds main-effects
# panels with non-dummy inputs swept on [0,1] (GCM dummies fixed at 0.5).
#
# Outputs (under outputs/forcing_tests/):
#   forcing_tests.csv          one row per (test, year): n, nRMSE, pass_2sd_pct,
#                              rmse_mm, mean_sd_mm, NED
#   loo_<year>.pdf             3x3 grid, one panel per test
#   main_effects_<year>.pdf    9 pages, one per test
#
# Independent of the production emulator setup; does NOT touch outputs/emulators/.

if (!exists("ensemble")) source("R/01_load_data.R")

.fi_tests <- list(
  list(label = "tas",                   forcing = c("tas")),
  list(label = "thetao",                forcing = c("thetao")),
  list(label = "tas+thetao",            forcing = c("tas", "thetao")),
  list(label = "smb",                   forcing = c("smb")),
  list(label = "thermal_forcing",       forcing = c("thermal_forcing")),
  list(label = "smb+thermal_forcing",   forcing = c("smb", "thermal_forcing")),
  list(label = "GMSTa",                 forcing = c("GMSTa")),
  list(label = "GMSTa+thetao",          forcing = c("GMSTa", "thetao")),
  list(label = "GMSTa+thermal_forcing", forcing = c("GMSTa", "thermal_forcing"))
)
.fi_years <- c(2021, 2300)

.fi_colors <- c("orangered", "seagreen3", "royalblue1", "plum2",
                "slateblue1", "sienna1", "goldenrod1", "darkseagreen", "orchid")

.fi_train <- function(year, forcing_cols, include_dummies = TRUE) {
  cols <- if (include_dummies) c(param_cols, dummy_gcms, forcing_cols)
          else                 c(param_cols, forcing_cols)
  stopifnot(all(forcing_cols %in% names(ensemble)))
  specs <- list()
  for (col in c(param_cols, forcing_cols))
    specs[[col]] <- norm_fit(ensemble[[col]], log = col %in% log_cols)
  yname <- ycol(year)
  ok    <- stats::complete.cases(ensemble[, c(param_cols, forcing_cols, yname)])
  d     <- ensemble[ok, , drop = FALSE]
  X     <- matrix(0, nrow = nrow(d), ncol = length(cols),
                  dimnames = list(NULL, cols))
  for (col in param_cols)   X[, col] <- norm_apply(d[[col]], specs[[col]])
  if (include_dummies)
    for (g   in dummy_gcms)   X[, g]   <- as.integer(d[[gcm_col]] == g)
  for (col in forcing_cols) X[, col] <- norm_apply(d[[col]], specs[[col]])
  m <- NULL
  invisible(utils::capture.output(
    m <- RobustGaSP::rgasp(
      design = X, response = d[[yname]], trend = cbind(1, X),
      kernel_type = cfg$emulator$kernel_type,
      alpha       = cfg$emulator$alpha,
      nugget.est  = cfg$emulator$nugget_est)))
  list(model = m, specs = specs, emu_cols = cols, X = X,
       y_m = d[[yname]], forcing_cols = forcing_cols)
}

.fi_loo_df <- function(fit) {
  loo <- RobustGaSP::leave_one_out_rgasp(fit$model)
  data.frame(actual = fit$y_m * 1000,
             mean   = loo$mean * 1000,
             sd     = loo$sd   * 1000)
}

.fi_stats <- function(df) {
  df$UB <- df$mean + 2 * df$sd; df$LB <- df$mean - 2 * df$sd
  df$ok <- df$actual > df$LB & df$actual < df$UB
  rmse  <- sqrt(mean((df$actual - df$mean)^2))
  list(n            = nrow(df),
       rmse_mm      = rmse,
       nRMSE        = rmse / sd(df$actual),
       pass_2sd_pct = mean(df$ok) * 100,
       NED          = sqrt(sum((df$actual - df$mean)^2 / df$sd^2)),
       mean_sd_mm   = mean(df$sd))
}

.fi_main_effects <- function(fit, me_gcm = NULL, me_scenario = "ssp585") {
  cols            <- fit$emu_cols
  # Only sweep the 6 physical parameters. Sweeping a forcing input would push
  # it outside the (GCM, scenario)-specific value it occupies in training, which
  # is off-manifold and unphysical.
  sweep_cols      <- param_cols
  dummies_present <- intersect(dummy_gcms, cols)
  ntest           <- 100; grid <- seq(0, 1, length.out = ntest)

  # Forcing values from the (me_gcm, me_scenario) row of the ensemble, in
  # normalised space using the emulator's own specs. This keeps the (non-swept)
  # forcing coordinates on the training manifold for that GCM, avoiding the
  # MIROC-style off-manifold extrapolation that "all other inputs = 0.5" causes
  # for cold-end GCMs whose forcing range doesn't include 0.5.
  forcing_vals <- list()
  if (length(fit$forcing_cols) > 0 && !is.null(me_gcm)) {
    sel <- ensemble[[gcm_col]] == me_gcm & ensemble[[scen_col]] == me_scenario
    stopifnot(any(sel))
    row <- ensemble[which(sel)[1], ]
    for (fcol in fit$forcing_cols)
      forcing_vals[[fcol]] <- norm_apply(row[[fcol]], fit$specs[[fcol]])
  }

  out <- list()
  for (col in sweep_cols) {
    X <- matrix(0.5, nrow = ntest, ncol = length(cols),
                dimnames = list(NULL, cols))
    # GCM one-hot
    if (length(dummies_present)) {
      X[, dummies_present] <- 0
      if (!is.null(me_gcm) && me_gcm %in% dummies_present) X[, me_gcm] <- 1
    }
    # Fixed forcing at the (me_gcm, me_scenario) actual values
    for (fcol in names(forcing_vals)) X[, fcol] <- forcing_vals[[fcol]]
    # Sweep this column 0->1 (overrides any fixed value above)
    X[, col] <- grid
    pr <- predict(fit$model, X, testing_trend = cbind(1, X))
    out[[col]] <- data.frame(x = grid,
                             mean = pr$mean * 1000,
                             LB   = (pr$mean - 2 * pr$sd) * 1000,
                             UB   = (pr$mean + 2 * pr$sd) * 1000)
  }
  out
}

.fi_loo_panel <- function(df, label, stats) {
  df$UB <- df$mean + 2 * df$sd; df$LB <- df$mean - 2 * df$sd
  df$ok <- df$actual > df$LB & df$actual < df$UB
  good  <- df[df$ok, ]; bad <- df[!df$ok, ]
  ymax  <- max(df$UB); ymin <- min(df$LB)
  plot(good$actual, good$mean, xlim = c(ymin, ymax), ylim = c(ymin, ymax),
       xlab = "Actual SLC (mm)", ylab = "Predicted SLC (mm)",
       main = label, col = "cornflowerblue", pch = 1, cex = 0.5)
  points(bad$actual, bad$mean, col = "indianred", pch = 1, cex = 0.5)
  segments(good$actual, good$UB, good$actual, good$LB,
           col = rgb(0, 0, 1, alpha = 0.2), lwd = 2)
  segments(bad$actual, bad$UB, bad$actual, bad$LB,
           col = rgb(1, 0, 0, alpha = 0.2), lwd = 2)
  abline(0, 1, lwd = 1)
  legend("topleft", bty = "n", cex = 0.7, text.col = "grey25",
         legend = c(sprintf("RMSE = %.3f mm", stats$rmse_mm),
                    sprintf("pass = %.1f %%",  stats$pass_2sd_pct),
                    sprintf("NED = %.1f",      stats$NED)))
}

.fi_main_effects_page <- function(me, label, year, has_dummies,
                                  me_gcm = NULL, me_scenario = "ssp585") {
  N <- length(me)
  a <- floor(sqrt(N)); b <- ceiling(N / a)
  par(mfrow = c(a, b), mar = c(4, 4, 2, 1), oma = c(0, 0, 2.5, 0))
  ymax <- max(vapply(me, function(d) max(d$UB), numeric(1)))
  ymin <- min(vapply(me, function(d) min(d$LB), numeric(1)))
  for (i in seq_along(me)) {
    d   <- me[[i]]
    col <- .fi_colors[((i - 1L) %% length(.fi_colors)) + 1L]
    plot(d$x, d$mean, type = "n", xlim = c(0, 1), ylim = c(ymin, ymax),
         xlab = names(me)[i], ylab = "SLC (mm)")
    polygon(c(d$x, rev(d$x)), c(d$LB, rev(d$UB)),
            col = adjustcolor(col, 0.3), border = NA)
    lines(d$x, d$mean, col = col, lwd = 1.8)
  }
  shown_gcm <- if (is.null(me_gcm)) ref_gcm else me_gcm
  hdr <- if (has_dummies)
    sprintf("%s — parameter main effects at %d (other params at 0.5; forcing at %s, %s)",
            label, year, shown_gcm, me_scenario)
  else
    sprintf("%s — parameter main effects at %d (other params at 0.5)",
            label, year)
  mtext(hdr, outer = TRUE, line = 0.7, cex = 1.0, font = 2)
}

run_forcing_input_tests <- function(include_dummies = TRUE,
                                     me_gcm = ref_gcm,
                                     me_scenario = "ssp585",
                                     out_subdir = if (include_dummies) "forcing_tests"
                                                  else "forcing_tests_no_gcm") {
  if (include_dummies) {
    stopifnot(all(me_gcm %in% c(ref_gcm, dummy_gcms)))
    stopifnot(me_scenario %in% unique(ensemble[[scen_col]]))
  }
  fi_dir <- file.path(out_dir, out_subdir)
  dir.create(fi_dir, showWarnings = FALSE, recursive = TRUE)
  # display_gcms: list of GCMs to compute main effects for, sharing the same
  # trained emulators. NULL means "no dummies" — forcing stays at 0.5.
  display_gcms <- if (include_dummies) as.list(me_gcm) else list(NULL)
  message(sprintf("Forcing-input tests: include_dummies=%s, me_gcm=[%s], me_scenario=%s -> %s",
                  include_dummies, paste(me_gcm, collapse = ","), me_scenario, fi_dir))
  results <- vector("list", length(.fi_tests))
  for (i in seq_along(.fi_tests)) {
    t <- .fi_tests[[i]]
    results[[i]] <- list(label = t$label, forcing = t$forcing, per_year = list())
    for (yr in .fi_years) {
      t0    <- Sys.time()
      fit   <- .fi_train(yr, t$forcing, include_dummies = include_dummies)
      loo   <- .fi_loo_df(fit)
      stats <- .fi_stats(loo)
      me_by_gcm <- list()
      for (g in display_gcms) {
        key <- if (is.null(g)) "__no_gcm__" else g
        me_by_gcm[[key]] <- .fi_main_effects(fit, me_gcm = g,
                                             me_scenario = me_scenario)
      }
      dt    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      results[[i]]$per_year[[as.character(yr)]] <-
        list(loo = loo, stats = stats, me_by_gcm = me_by_gcm)
      message(sprintf("  [%-22s | %d]  n=%d  RMSE=%.3f mm  pass=%4.1f%%  NED=%5.1f  (%4.1fs)",
                      t$label, yr, stats$n, stats$rmse_mm,
                      stats$pass_2sd_pct, stats$NED, dt))
    }
  }

  # Summary CSV
  rows <- list()
  for (r in results) for (yr in .fi_years) {
    s <- r$per_year[[as.character(yr)]]$stats
    rows[[length(rows) + 1L]] <- data.frame(
      test = r$label, inputs = paste(r$forcing, collapse = "+"),
      year = yr, n = s$n, nRMSE = s$nRMSE,
      pass_2sd_pct = s$pass_2sd_pct, rmse_mm = s$rmse_mm,
      mean_sd_mm = s$mean_sd_mm, NED = s$NED)
  }
  summary_df <- do.call(rbind, rows)
  write.csv(summary_df, file.path(fi_dir, "forcing_tests.csv"), row.names = FALSE)

  # LOO PDFs — 3x3 grid per year
  for (yr in .fi_years) {
    pdf(file.path(fi_dir, sprintf("loo_%d.pdf", yr)),
        width = 13, height = 13)
    par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))
    for (r in results) {
      s <- r$per_year[[as.character(yr)]]
      .fi_loo_panel(s$loo, r$label, s$stats)
    }
    dev.off()
  }

  # Main effects PDFs — one PDF per (year, GCM). Filename includes the
  # (GCM, scenario) query point when dummies are present so different choices
  # don't collide.
  for (yr in .fi_years) for (g in display_gcms) {
    key        <- if (is.null(g)) "__no_gcm__" else g
    gcm_suffix <- if (include_dummies) paste0("_", g, "_", me_scenario) else ""
    pdf(file.path(fi_dir, sprintf("main_effects_%d%s.pdf", yr, gcm_suffix)),
        width = 13, height = 7)
    for (r in results) {
      s  <- r$per_year[[as.character(yr)]]
      me <- s$me_by_gcm[[key]]
      .fi_main_effects_page(me, r$label, yr,
                            has_dummies = include_dummies,
                            me_gcm = g,
                            me_scenario = me_scenario)
    }
    dev.off()
  }

  cat("\n=== Forcing-input test summary (RMSE in mm, pass in %, NED unitless) ===\n")
  for (yr in .fi_years) {
    cat(sprintf("\n--- %d ---\n", yr))
    d <- summary_df[summary_df$year == yr, ]
    d <- d[order(d$rmse_mm), ]
    for (i in seq_len(nrow(d)))
      cat(sprintf("  %-22s  n=%d  RMSE=%6.3f mm  pass=%5.1f%%  NED=%5.1f  nRMSE=%.3f\n",
                  d$test[i], d$n[i], d$rmse_mm[i], d$pass_2sd_pct[i],
                  d$NED[i], d$nRMSE[i]))
  }
  message("\n03d: wrote forcing_tests.csv + 4 PDFs to ", fi_dir)
}
