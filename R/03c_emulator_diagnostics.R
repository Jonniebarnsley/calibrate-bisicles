# 03c — emulator diagnostics for the by_continent (ANT-only) pipeline. Produces
# per-year PDFs in outputs/diagnostics/:
#   loo_<year>.pdf          single LOO actual-vs-predicted panel; small-text
#                           RMSE / pass / NED in the upper-left corner
#   main_effects_<year>.pdf 11-panel grid (sweep each non-dummy input on [0,1]
#                           with all others — incl GCM dummies — at 0.5; circum-
#                           centre of the dummy points, equidistant from all 4 GCMs)
#
# Mode-agnostic: per-year GP -> rgasp leave_one_out; SVD mode -> reconstruct LOO
# from coefficient-LOO arrays + truncation variance (same recipe svd_emulator.R
# uses for its diagnostic CSV). Regeneration is gated on a meta sidecar; sigma_mod
# changes do NOT trigger regeneration (sigma_mod isn't an emulator property).

if (!exists("emu_cols")) source("R/03_emulator_io.R")

.diag_meta_keys <- function() list(
  ensemble     = cfg$paths$ensemble,
  forcing_cols = cfg$data$forcing_cols,
  emulator     = cfg$emulator,
  svd          = cfg$svd,
  years        = if (is.null(cfg$diagnostics$years)) c(2021, 2300)
                 else cfg$diagnostics$years,
  layout       = "by-year-v1"
)
.diag_dir       <- file.path(out_dir, "diagnostics")
.diag_meta_path <- file.path(.diag_dir, "meta.rds")

.diag_colors <- c("orangered", "seagreen3", "royalblue1", "plum2",
                  "slateblue1", "sienna1", "goldenrod1", "darkseagreen", "orchid",
                  "orangered", "seagreen3", "royalblue1", "plum2")

# Per-year-mode LOO at one year: rgasp leave_one_out -> data.frame(actual, mean, sd) in mm.
.diag_loo_peryear <- function(year) {
  m    <- train_emulator(year)
  resp <- ensemble[[ycol(year)]]
  ok   <- is.finite(resp)
  loo  <- leave_one_out_rgasp(m)
  data.frame(actual = resp[ok] * 1000,
             mean   = loo$mean  * 1000,
             sd     = loo$sd    * 1000)
}

# SVD-mode LOO at one year: reconstruct from coefficient-LOO + truncation variance.
# .V/.mu/.sds/.Y/.scores/.trunc_var/.coef_emu/svd_years live in R/svd_emulator.R.
.diag_loo_svd <- function(year) {
  r        <- cfg$svd$n_components
  coef_loo <- lapply(seq_len(r),
                     function(i) as.data.frame(leave_one_out_rgasp(.coef_emu[[i]])))
  ahat <- sapply(coef_loo, function(l) l$mean)        # n_train x r
  asd  <- sapply(coef_loo, function(l) l$sd)
  ti   <- match(year, svd_years); stopifnot(!is.na(ti))
  V_t  <- .V[ti, 1:r, drop = FALSE]                   # 1 x r
  Yhat_std <- as.numeric(ahat %*% t(V_t))
  Yvar_std <- as.numeric((asd^2) %*% t(V_t^2)) + .trunc_var[ti]
  data.frame(actual = .Y[, ti],
             mean   = Yhat_std * .sds[ti] + .mu[ti],
             sd     = sqrt(Yvar_std) * .sds[ti])
}

.diag_loo <- function(year) {
  if (isTRUE(cfg$svd$enabled)) .diag_loo_svd(year)
  else                         .diag_loo_peryear(year)
}

# Main effects: each non-dummy input swept on [0,1], all other inputs (incl 3 GCM
# dummies) at 0.5 — the circumcentre of the 4 GCM dummy vectors {(0,0,0), (1,0,0),
# (0,1,0), (0,0,1)}, equidistant from each so no GCM is favoured. Dummies aren't
# plotted (sweeping a fractional one-hot is meaningless). predict_slc_mm
# dispatches per-year / SVD internally.
.diag_main_effects <- function(year) {
  sweep_cols <- setdiff(emu_cols, dummy_gcms)
  ntest      <- 100
  grid       <- seq(0, 1, length.out = ntest)
  out <- list()
  for (col in sweep_cols) {
    X <- matrix(0.5, nrow = ntest, ncol = length(emu_cols),
                dimnames = list(NULL, emu_cols))
    X[, col] <- grid
    pr <- predict_slc_mm(X, year)
    out[[col]] <- data.frame(x = grid, mean = pr$mean,
                             LB = pr$mean - 2 * pr$sd,
                             UB = pr$mean + 2 * pr$sd)
  }
  out
}

.diag_loo_pdf <- function(year) {
  pdf(file.path(.diag_dir, sprintf("loo_%d.pdf", year)),
      width = 4.6, height = 4.6)
  par(mar = c(4, 4, 2.5, 1))
  df    <- .diag_loo(year)
  df$UB <- df$mean + 2 * df$sd; df$LB <- df$mean - 2 * df$sd
  df$ok <- df$actual > df$LB & df$actual < df$UB
  good  <- df[df$ok, ]; bad <- df[!df$ok, ]
  rmse  <- sqrt(mean((df$actual - df$mean)^2))
  ned   <- sqrt(sum((df$actual - df$mean)^2 / df$sd^2))
  pass  <- mean(df$ok) * 100
  ymax  <- max(df$UB);  ymin <- min(df$LB)
  plot(good$actual, good$mean,
       xlim = c(ymin, ymax), ylim = c(ymin, ymax),
       xlab = "Actual SLC (mm)", ylab = "Predicted SLC (mm)",
       main = sprintf("ANT — %d", year),
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
  dev.off()
}

.diag_main_effects_pdf <- function(year) {
  me   <- .diag_main_effects(year)
  N    <- length(me)
  a    <- floor(sqrt(N)); b <- ceiling(N / a)
  pdf(file.path(.diag_dir, sprintf("main_effects_%d.pdf", year)),
      width = 3.2 * b, height = 2.8 * a)
  ymax <- max(vapply(me, function(d) max(d$UB), numeric(1)))
  ymin <- min(vapply(me, function(d) min(d$LB), numeric(1)))
  par(mfrow = c(a, b), mar = c(4, 4, 2, 1), oma = c(0, 0, 2.5, 0))
  for (i in seq_along(me)) {
    d   <- me[[i]]
    col <- .diag_colors[((i - 1L) %% length(.diag_colors)) + 1L]
    plot(d$x, d$mean, type = "n",
         xlim = c(0, 1), ylim = c(ymin, ymax),
         xlab = names(me)[i], ylab = "SLC (mm)")
    polygon(c(d$x, rev(d$x)), c(d$LB, rev(d$UB)),
            col = adjustcolor(col, 0.3), border = NA)
    lines(d$x, d$mean, col = col, lwd = 1.8)
  }
  mtext(sprintf("ANT — main effects at %d (all other inputs fixed at 0.5)", year),
        outer = TRUE, line = 0.7, cex = 1.0, font = 2)
  dev.off()
}

run_emulator_diagnostics <- function(force = FALSE) {
  meta <- .diag_meta_keys()
  if (!force && file.exists(.diag_meta_path) &&
      identical(readRDS(.diag_meta_path), meta)) {
    message("03c: diagnostics up-to-date, skipping")
    return(invisible(NULL))
  }
  dir.create(.diag_dir, showWarnings = FALSE, recursive = TRUE)
  for (yr in meta$years) {
    .diag_loo_pdf(yr)
    .diag_main_effects_pdf(yr)
    message(sprintf("03c: wrote loo_%d.pdf + main_effects_%d.pdf", yr, yr))
  }
  saveRDS(meta, .diag_meta_path)
}

run_emulator_diagnostics()
