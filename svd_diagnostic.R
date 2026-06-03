# SVD / functional-PCA diagnostic for the ensemble SLC timeseries.
# Q: how many temporal components capture the ensemble, and (crucially) how well does
# a truncated reconstruction recover the EARLY calibration period (2021) vs late
# century? Plain (variance-weighted) SVD is compared with a column-standardised SVD
# that gives the low-variance early years fair weight.
# Run from calibration/:  Rscript svd_diagnostic.R

ens   <- read.csv("emulator_inputs_2007-2300_sle.csv")
ycols <- grep("^X[0-9]{4}$", names(ens), value = TRUE)
years <- as.integer(sub("^X", "", ycols))
Y     <- as.matrix(ens[, ycols]) * 1000                 # mm SLE
Y     <- Y[stats::complete.cases(Y), ]
ri    <- match(c(2021, 2100, 2300), years)
sds_y <- apply(Y, 2, sd)
cat(sprintf("Runs: %d   years: %d (%d-%d)\n", nrow(Y), ncol(Y), min(years), max(years)))

# Centre (+ optionally column-standardise), SVD, return pieces for reconstruction.
svd_diag <- function(Y, standardise = FALSE) {
  mu  <- colMeans(Y)
  Yc  <- sweep(Y, 2, mu)
  sds <- if (standardise) apply(Yc, 2, sd) else rep(1, ncol(Y)); sds[sds == 0] <- 1
  sv  <- svd(sweep(Yc, 2, sds, "/"))
  list(mu = mu, sds = sds, sv = sv, varexp = sv$d^2 / sum(sv$d^2))
}
# Reconstruct with r components, back in mm.
recon <- function(d, r) {
  s <- d$sv
  a <- s$u[, 1:r, drop = FALSE] %*% diag(s$d[1:r], r, r) %*% t(s$v[, 1:r, drop = FALSE])
  sweep(sweep(a, 2, d$sds, "*"), 2, d$mu, "+")
}

for (std in c(FALSE, TRUE)) {
  d  <- svd_diag(Y, std)
  cv <- cumsum(d$varexp)
  cat(sprintf("\n=== %s SVD ===\n", if (std) "Column-standardised" else "Plain (variance-weighted)"))
  cat(sprintf("components for 99%% / 99.9%% / 99.99%% variance: %d / %d / %d\n",
              which(cv >= .99)[1], which(cv >= .999)[1], which(cv >= .9999)[1]))
  cat("cumulative var (r = 1..6):", paste(sprintf("%.4f", cv[1:6]), collapse = " "), "\n\n")
  cat(sprintf("%-3s %12s %12s %12s\n", "r", "nRMSE@2021", "nRMSE@2100", "nRMSE@2300"))
  for (r in c(1, 2, 3, 4, 5, 8)) {
    rmse <- sqrt(colMeans((Y - recon(d, r))^2))[ri]
    cat(sprintf("%-3d %12.4f %12.4f %12.4f    (2021 abs: %6.3f mm)\n",
                r, (rmse / sds_y[ri])[1], (rmse / sds_y[ri])[2], (rmse / sds_y[ri])[3], rmse[1]))
  }
}
cat(sprintf("\nContext: SLC SD at 2021/2100/2300 = %.2f / %.1f / %.1f mm;  sigma_obs = 1.11 mm\n",
            sds_y[ri][1], sds_y[ri][2], sds_y[ri][3]))
cat("=> reconstruction RMSE at 2021 should be well below sigma_obs to leave the calibration unaffected.\n")
