# Leave-one-out diagnostics.
#
# leave_one_out(emu, year) returns a data.frame with columns (actual, mean,
# sd) in mm, classed as "loo_result" and carrying (label, year) as
# attributes. plot.loo_result draws the scatter (actual vs predicted with
# +/-2 sd error bars and headline metrics) to the active device.
#
# Sourced lazily from 03_emulate.R inside its `if (cfg$emulator$leave_one_out)`
# block. The two compute methods share a return shape so plot.loo_result
# works on either mode's output.

leave_one_out.emulator_peryear <- function(emu, year, ...) {
  m   <- fetch_peryear_gp(emu, year)
  loo <- leave_one_out_rgasp(m)
  df  <- data.frame(actual = emu$ens[[ycol(year)]] * 1000,
                    mean   = loo$mean * 1000,
                    sd     = loo$sd   * 1000)
  structure(df, class = c("loo_result", "data.frame"),
            label = emu$label, year = year)
}

leave_one_out.emulator_svd <- function(emu, year, ...) {
  r             <- emu$ncomp
  coef_mean_loo <- sapply(emu$coef_loo[1:r], function(l) l$mean)
  coef_sd_loo   <- sapply(emu$coef_loo[1:r], function(l) l$sd)
  out           <- reconstruct_slc_mm(coef_mean_loo, coef_sd_loo, year, emu)
  year_i        <- match(year, emu$years)
  df <- data.frame(actual = emu$Y[, year_i],
                   mean   = as.numeric(out$mean),
                   sd     = as.numeric(out$sd))
  structure(df, class = c("loo_result", "data.frame"),
            label = emu$label, year = year)
}

# LOO scatter PDF: actual vs predicted, +/-2sd error bars, headline metrics
# (RMSE, pass-rate, normalised error distance). Draws to the active device;
# caller manages pdf()/dev.off().
plot.loo_result <- function(x, ...) {
  label <- attr(x, "label")
  year <- attr(x, "year")

  x$UB <- x$mean + 2 * x$sd
  x$LB <- x$mean - 2 * x$sd
  x$ok <- x$actual > x$LB & x$actual < x$UB

  good <- x[x$ok, ]
  bad  <- x[!x$ok, ]

  rmse <- sqrt(mean((x$actual - x$mean)^2))
  ned  <- sqrt(sum((x$actual - x$mean)^2 / x$sd^2))
  pass <- mean(x$ok) * 100
  lim  <- range(x$LB, x$UB)

  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mar = c(4, 4, 2.5, 1))

  plot(good$actual, good$mean,
       xlim = lim, ylim = lim,
       xlab = "Actual SLC (mm)", ylab = "Predicted SLC (mm)",
       main = sprintf("LOO — %s, %d", label, year),
       col = "cornflowerblue", pch = 1, cex = 0.5)
  points(bad$actual, bad$mean, col = "indianred", pch = 1, cex = 0.5)

  segments(good$actual, good$UB, good$actual, good$LB,
           col = rgb(0, 0, 1, alpha = 0.2), lwd = 2)
  segments(bad$actual, bad$UB, bad$actual, bad$LB,
           col = rgb(1, 0, 0, alpha = 0.2), lwd = 2)

  abline(0, 1, lwd = 1)
  legend("topleft", bty = "n", cex = 0.75, text.col = "grey25",
         legend = c(sprintf("RMSE = %.3f mm", rmse),
                    sprintf("pass = %.1f %%",  pass),
                    sprintf("NED = %.1f",      ned)))
  invisible(x)
}
