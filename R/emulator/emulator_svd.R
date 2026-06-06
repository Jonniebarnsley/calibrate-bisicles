# SVD basis-function emulator (Higdon et al. 2008). Standardises each year,
# takes the SVD basis, emulates the leading n_components coefficients as GPs,
# and reconstructs SLC(year) with a truncation-variance term (the discarded
# components enter as a per-year residual variance).
#
# Provides:
#   build_emulator_svd(ens, label)  - constructor; returns an emulator_svd
#   predict.emulator_svd(...)        - S3 method
#   reconstruct_slc_mm(...)          - shared helper (also used by LOO)
#   run_svd_diagnostic(emu)          - optional rank-selection CSV

stopifnot(exists("X_train"))

build_emulator_svd <- function(ens, label) {
  ycols   <- grep("^[0-9]{4}$", names(ens), value = TRUE)
  Y       <- as.matrix(ens[, ycols]) * 1000   # mm

  # Standardise the slc curves to ensure that later years (with higher slc)
  # don't dominate the SVD.
  mu      <- colMeans(Y)
  Yc      <- sweep(Y, 2, mu)
  sds     <- apply(Yc, 2, sd)
  sds[sds == 0] <- 1  # 0 -> 1 to avoid dividing by zero in next step
  Ys      <- sweep(Yc, 2, sds, "/")

  svd_out <- svd(Ys)
  scores  <- svd_out$u %*% diag(svd_out$d)    # Coefficients
  varfrac <- svd_out$d^2 / sum(svd_out$d^2)   # Variance explained per component
  ncomp   <- cfg$svd$n_components

  # Per-year truncation uncertainty arising from dropping higher components.
  # Mean squared residual of first n components vs the full reconstructed slc.
  trunc_var <- colMeans((Ys - scores[, 1:ncomp, drop = FALSE] %*%
                           t(svd_out$v[, 1:ncomp, drop = FALSE]))^2)

  # Fit one GP per retained component and save its LOO results.
  coef_emu <- vector("list", ncomp)
  coef_loo <- vector("list", ncomp)
  for (i in seq_len(ncomp)) {
    invisible(utils::capture.output(
      coef_emu[[i]] <- rgasp(
        design      = X_train,
        response    = scores[, i],
        trend       = trend,
        kernel_type = cfg$emulator$kernel_type,
        alpha       = cfg$emulator$alpha,
        nugget.est  = cfg$emulator$nugget_est
      )
    ))
    coef_loo[[i]] <- as.data.frame(leave_one_out_rgasp(coef_emu[[i]]))
  }

  message(sprintf(
    "SVD [%s]: %d components (%.2f%% var).",
    label, ncomp, 100 * sum(varfrac[1:ncomp])
  ))

  emu <- list(
    label     = label,
    years     = as.integer(ycols),
    Y         = Y,           # Raw output matrix in mm
    mu        = mu,          # Per-column means of Y
    sds       = sds,         # Per-column sd of Y
    Ys        = Ys,          # Standardised output matrix
    V         = svd_out$v,   # Temporal basis functions
    scores    = scores,      # SVD coefficient values
    varfrac   = varfrac,     # Variance explained by each comp
    ncomp     = ncomp,       # Number of components retained
    coef_emu  = coef_emu,    # List of GP emulators for each SVD coeff
    coef_loo  = coef_loo,    # LOO results for each emulator
    trunc_var = trunc_var    # Per-year truncation uncertainty
  )
  class(emu) <- c("emulator_svd", "emulator")
  emu
}

# Reconstruct slc mean & sd (mm) at requested year(s) from per-sample
# coefficient mean/sd matrices (K x r) using the SVD basis, including the
# truncation variance. Shared by predict.emulator_svd and the SVD LOO method.
reconstruct_slc_mm <- function(coef_mean, coef_sd, year, emu,
                               trunc_var = emu$trunc_var) {
  year_i <- match(year, emu$years)
  stopifnot(!any(is.na(year_i)))                # Years must be in training data
  r      <- ncol(coef_mean)
  V_ir   <- emu$V[year_i, 1:r, drop = FALSE]

  # Reconstruct mean and total variance in standardised (z) space
  mean_z    <- coef_mean %*% t(V_ir)
  emu_var_z <- (coef_sd^2) %*% t(V_ir^2)
  tot_var_z <- sweep(emu_var_z, 2, trunc_var[year_i], "+")

  # Un-standardise to slc mean and sd in mm
  sds <- emu$sds[year_i]
  mu  <- emu$mu[year_i]
  mean_mm <- sweep(sweep(mean_z, 2, sds, "*"), 2, mu, "+")
  sd_mm   <- sweep(sqrt(tot_var_z), 2, sds, "*")
  list(mean = mean_mm, sd = sd_mm)
}

# Predict SLC mean & sd (mm) at year(s) for a normalised input matrix X.
predict.emulator_svd <- function(object, X, year, ...) {
  K <- nrow(X); r <- object$ncomp
  coef_mean <- matrix(0, K, r)   # dim = (nsamples, ncomp)
  coef_sd   <- matrix(0, K, r)

  # For each component, predict its coefficient and store its mean and sd.
  for (i in seq_len(r)) {
    pr <- predict(object$coef_emu[[i]], X, testing_trend = cbind(1, X))
    coef_mean[, i] <- pr$mean
    coef_sd[, i]   <- pr$sd
  }

  out <- reconstruct_slc_mm(coef_mean, coef_sd, year, object)
  if (length(year) == 1L)
    list(mean = as.numeric(out$mean), sd = as.numeric(out$sd))
  else
    out
}

# Diagnostics: for r = 1..9, compute variance explained, coefficient LOO R2,
# and the end-to-end LOO reconstruction metrics (RMSE in mm, nRMSE, NED,
# pass_rate) at 2021 and 2300.
run_svd_diagnostic <- function(emu) {

  N_DIAG <- 9     # rank-selection table covers r = 1..N_DIAG

  # Local aliases for readability of the matrix algebra below.
  Y       <- emu$Y
  Ys      <- emu$Ys
  V       <- emu$V
  scores  <- emu$scores
  varfrac <- emu$varfrac
  ncomp   <- emu$ncomp

  # Reuse the production coef_loo for the first ncomp components; fit any
  # additional components (ncomp+1..N_DIAG) locally on the fly.
  coef_loo <- emu$coef_loo
  if (ncomp < N_DIAG) for (i in (ncomp + 1L):N_DIAG) {
    invisible(utils::capture.output(
      m <- rgasp(
        design      = X_train,
        response    = scores[, i],
        trend       = trend,
        kernel_type = cfg$emulator$kernel_type,
        alpha       = cfg$emulator$alpha,
        nugget.est  = cfg$emulator$nugget_est
      )
    ))
    coef_loo[[i]] <- as.data.frame(leave_one_out_rgasp(m))
  }
  coef_mean_loo <- sapply(coef_loo, function(l) l$mean)
  coef_sd_loo   <- sapply(coef_loo, function(l) l$sd)

  # Years at which to report end-to-end LOO metrics, and their column indices.
  diagnostic_years    <- c(2021, 2300)
  diagnostic_year_i   <- match(diagnostic_years, emu$years)
  slc_year_sd         <- apply(Y, 2, sd)

  rows <- lapply(seq_len(N_DIAG), function(r) {

    # Truncate coef LOO arrays and scores/basis to the first r components.
    coef_mean_loo_r <- coef_mean_loo[, 1:r, drop = FALSE]
    coef_sd_loo_r   <- coef_sd_loo[, 1:r, drop = FALSE]
    scores_r        <- scores[, 1:r, drop = FALSE]
    V_r             <- V[, 1:r, drop = FALSE]

    # Per-year truncation variance at this r (in standardised space, all years).
    trunc_var_r <- colMeans((Ys - scores_r %*% t(V_r))^2)

    # Reconstruct slc mean and sd (mm) at the diagnostic years only, using the
    # r-specific truncation.
    recon <- reconstruct_slc_mm(coef_mean_loo_r, coef_sd_loo_r,
                                diagnostic_years, emu, trunc_var = trunc_var_r)

    # R-squared for the r'th SVD component's LOO test.
    loo_r2 <- 1 - mean((coef_loo[[r]]$mean - scores[, r])^2) / var(scores[, r])

    out <- data.frame(
      r       = r,
      var_pc  = 100 * varfrac[r],         # variance explained by r'th component
      var_cum = 100 * sum(varfrac[1:r]),  # cumulative variance explained
      loo_r2  = loo_r2
    )

    for (k in seq_along(diagnostic_years)) {    # 2021 and 2300
      yi      <- diagnostic_year_i[k]
      err     <- Y[, yi] - recon$mean[, k]      # emulator-BISICLES misfit
      pred_sd <- recon$sd[, k]                  # emulator uncertainty
      rmse    <- sqrt(mean(err^2))
      nrmse   <- rmse / slc_year_sd[yi]
      ned     <- sqrt(sum((err / pred_sd)^2))
      pass    <- mean(abs(err) <= 2 * pred_sd)
      out[[paste0("rmse_mm_", diagnostic_years[k])]] <- rmse
      out[[paste0("nRMSE_",   diagnostic_years[k])]] <- nrmse
      out[[paste0("NED_",     diagnostic_years[k])]] <- ned
      out[[paste0("pass_",    diagnostic_years[k])]] <- pass
    }
    out
  })

  diagnostic_csv_path <- file.path(out_dir,
                                   sprintf("svd_diagnostic_%s.csv", emu$label))
  write.csv(do.call(rbind, rows), diagnostic_csv_path, row.names = FALSE)
  message("SVD: diagnostic written to ", diagnostic_csv_path)
}
