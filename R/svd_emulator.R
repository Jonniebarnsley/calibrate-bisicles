# SVD basis-function emulator (Higdon et al. 2008 style) — used when svd.enabled.
# Standardises each year, takes the SVD basis, emulates the leading n_components
# coefficients as GPs, and reconstructs SLC(year) with a truncation-variance term
# (the discarded components enter as a per-year residual variance). Provides
# svd_predict_mm(Xnorm, year). Sourced by 03_emulator_io.R; uses ensemble, X_train,
# ycol, cfg, out_dir from there.

stopifnot(exists("X_train"))

svd_ycols <- grep("^X[0-9]{4}$", names(ensemble), value = TRUE)
svd_years <- as.integer(sub("^X", "", svd_ycols))

.Yall <- as.matrix(ensemble[, svd_ycols]) * 1000           # mm; rows align with X_train
.ok   <- stats::complete.cases(.Yall)
.Y    <- .Yall[.ok, , drop = FALSE]
.Xs   <- X_train[.ok, , drop = FALSE]
.mu   <- colMeans(.Y)
.Yc   <- sweep(.Y, 2, .mu)
.sds  <- apply(.Yc, 2, sd); .sds[.sds == 0] <- 1          # per-year standardisation
.Ys   <- sweep(.Yc, 2, .sds, "/")
.sv   <- svd(.Ys)
.V    <- .sv$v
.scores  <- .sv$u %*% diag(.sv$d)
.varfrac <- .sv$d^2 / sum(.sv$d^2)

ncomp <- cfg$svd$n_components
ndiag <- if (isTRUE(cfg$svd$run_diagnostic)) max(9L, ncomp) else ncomp
.tr   <- cbind(1, .Xs)

# Fit one GP per coefficient (for prediction); keep LOO for the diagnostic.
.coef_emu <- vector("list", ndiag)
.coef_loo <- vector("list", ndiag)
for (i in seq_len(ndiag)) {
  invisible(utils::capture.output(
    .coef_emu[[i]] <- rgasp(design = .Xs, response = .scores[, i], trend = .tr,
      kernel_type = cfg$emulator$kernel_type, alpha = cfg$emulator$alpha,
      nugget.est = cfg$emulator$nugget_est)))
  if (isTRUE(cfg$svd$run_diagnostic))
    .coef_loo[[i]] <- as.data.frame(leave_one_out_rgasp(.coef_emu[[i]]))
}

# Per-year truncation variance (standardised) for the chosen n_components.
.trunc_var <- colMeans((.Ys - .scores[, 1:ncomp, drop = FALSE] %*%
                          t(.V[, 1:ncomp, drop = FALSE]))^2)

message(sprintf("SVD: %d components (%.2f%% var); truncation-variance term included.",
                ncomp, 100 * sum(.varfrac[1:ncomp])))

# Predict SLC mean & sd (mm) at year(s) for a normalised input matrix Xnorm.
# Var(Yhat[t]) = sds[t]^2 * ( sum_i V[t,i]^2 sd_i^2  +  trunc_var[t] ).
svd_predict_mm <- function(Xnorm, year) {
  ti <- match(year, svd_years); stopifnot(!any(is.na(ti)))
  K  <- nrow(Xnorm); r <- ncomp
  ah <- as <- matrix(0, K, r)
  for (i in seq_len(r)) {
    pr <- predict(.coef_emu[[i]], Xnorm, testing_trend = cbind(1, Xnorm))
    ah[, i] <- pr$mean; as[, i] <- pr$sd
  }
  Vt       <- .V[ti, 1:r, drop = FALSE]                    # length(year) x r
  mean_std <- ah %*% t(Vt)                                 # K x length(year)
  var_std  <- sweep((as^2) %*% t(Vt^2), 2, .trunc_var[ti], "+")
  mean_mm  <- sweep(sweep(mean_std, 2, .sds[ti], "*"), 2, .mu[ti], "+")
  sd_mm    <- sweep(sqrt(var_std), 2, .sds[ti], "*")
  if (length(year) == 1L) list(mean = as.numeric(mean_mm), sd = as.numeric(sd_mm))
  else                    list(mean = mean_mm, sd = sd_mm)
}

# ---------------------------------------------------------------------------
# Diagnostic CSV: for r = 1..9, variance explained, coefficient LOO R2, and the
# end-to-end LOO reconstruction metrics (nRMSE, NED = sqrt(sum((err/sd)^2)),
# pass = within +/-2sd, RMSE in mm) at 2021 and 2300, with truncation variance.
# ---------------------------------------------------------------------------
if (isTRUE(cfg$svd$run_diagnostic)) {
  ahat  <- sapply(.coef_loo, function(l) l$mean)           # n_train x ndiag
  asd   <- sapply(.coef_loo, function(l) l$sd)
  dyrs  <- c(2021, 2300); dri <- match(dyrs, svd_years)
  sds_y <- apply(.Y, 2, sd)
  rows <- lapply(1:9, function(r) {
    Yhat <- sweep(sweep(ahat[, 1:r, drop = FALSE] %*% t(.V[, 1:r, drop = FALSE]),
                        2, .sds, "*"), 2, .mu, "+")
    tvar <- colMeans((.Ys - .scores[, 1:r, drop = FALSE] %*% t(.V[, 1:r, drop = FALSE]))^2)
    Ysd  <- sweep(sqrt(sweep((asd[, 1:r, drop = FALSE]^2) %*% t(.V[, 1:r, drop = FALSE]^2),
                             2, tvar, "+")), 2, .sds, "*")
    out <- data.frame(r = r, var_pc = 100 * .varfrac[r],
                      var_cum = 100 * sum(.varfrac[1:r]),
                      loo_r2 = 1 - mean((.coef_loo[[r]]$mean - .scores[, r])^2) / var(.scores[, r]))
    for (k in seq_along(dyrs)) {
      cc <- dri[k]; err <- .Y[, cc] - Yhat[, cc]; s <- Ysd[, cc]
      out[[paste0("nRMSE_", dyrs[k])]]   <- sqrt(mean(err^2)) / sds_y[cc]
      out[[paste0("NED_", dyrs[k])]]     <- sqrt(sum((err / s)^2))
      out[[paste0("pass_", dyrs[k])]]    <- mean(abs(err) <= 2 * s)
      out[[paste0("rmse_mm_", dyrs[k])]] <- sqrt(mean(err^2))
    }
    out
  })
  write.csv(do.call(rbind, rows), file.path(out_dir, "svd_diagnostic.csv"), row.names = FALSE)
  message("SVD: diagnostic written to ", file.path(out_dir, "svd_diagnostic.csv"))
}
