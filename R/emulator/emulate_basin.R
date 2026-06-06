# build_emulator(ens, label) — fit emulators for one basin, run configured
# diagnostics, and return a predict function (closure over the fitted state).
# Dispatches to the SVD or per-year builder based on cfg$svd$enabled.
# Called once per basin from 03_emulate.R.

build_emulator <- function(ens, label) {

  if (isFALSE(cfg$svd$enabled)) {

    cache_dir <- file.path(out_dir, "emulators", label)
    built     <- build_peryear_emulators(ens, cache_dir)

    # LOO PDFs at the configured diagnostic years. Per-year LOO is just
    # leave_one_out_rgasp on each year's trained emulator. write_loo_pdf and
    # diag_dir live in 03_emulate.R (sourced before this file).
    if (isTRUE(cfg$emulator$leave_one_out))
      for (year in cfg$emulator$diagnostic_years) {
        m   <- built$train(year)
        loo <- leave_one_out_rgasp(m)
        write_loo_pdf(data.frame(actual = ens[[ycol(year)]] * 1000,
                                 mean   = loo$mean * 1000,
                                 sd     = loo$sd   * 1000),
                      year, label)
        message(sprintf("emulate_by_year: wrote loo_%d_%s.pdf", year, label))
      }

    built$predict

  } else {

    state <- build_svd_state(ens)

    # Leave one out tests, including truncation error into variance
    if (isTRUE(cfg$emulator$leave_one_out)) {
      r             <- state$ncomp
      coef_mean_loo <- sapply(state$coef_loo[1:r], function(l) l$mean)
      coef_sd_loo   <- sapply(state$coef_loo[1:r], function(l) l$sd)
      for (year in cfg$emulator$diagnostic_years) {
        out    <- reconstruct_slc_mm(coef_mean_loo, coef_sd_loo, year, state)
        year_i <- match(year, state$years)
        write_loo_pdf(data.frame(actual = state$Y[, year_i],
                                 mean   = as.numeric(out$mean),
                                 sd     = as.numeric(out$sd)),
                      year, label)
        message(sprintf("emulate_by_svd: wrote loo_%d_%s.pdf", year, label))
      }
    }

    if (isTRUE(cfg$svd$run_diagnostic)) run_svd_diagnostic(state)

    function(X, year) svd_predict_mm(X, year, state)
  }
}
