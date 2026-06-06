# Main-effects diagnostics.
#
# main_effects(emu, year, gcm, scenario) sweeps each of the 6 physical params
# on [0,1] with the others fixed at 0.5, GCM dummies one-hot for `gcm`, and
# forcing pinned to the (gcm, scenario) training row. Querying on-manifold
# avoids the off-data extrapolation that fixed dummies = 0.5 + forcing = 0.5
# would produce (no training sample sits at the dummy-space circumcentre,
# and forcing is GCM-determined, so the GP extrapolates wildly there).
#
# The compute method is mode-agnostic — it only calls predict() via S3
# dispatch — so one definition serves both emulator_svd and emulator_peryear.
# It returns a list of per-param data.frames classed as "me_result", with
# (label, year, gcm, scenario) as attributes; plot.me_result draws the
# panel grid to the active device.

diag_colors <- c("orangered", "seagreen3", "royalblue1", "plum2",
                 "slateblue1", "sienna1", "goldenrod1", "darkseagreen")

main_effects.emulator <- function(emu, year, gcm, scenario, ...) {
  stopifnot(gcm %in% gcms, scenario %in% scenarios)

  # Anchor forcing at the actual (gcm, scenario) training value, normalised.
  fl <- forcing_lookup[forcing_lookup[[gcm_col]] == gcm &
                         forcing_lookup[[scen_col]] == scenario, ]
  stopifnot(nrow(fl) == 1L)
  forcing_norm <- vapply(forcing_cols,
                         function(c) norm_apply(fl[[c]], norm_specs[[c]]),
                         numeric(1))

  ntest <- 100
  grid  <- seq(0, 1, length.out = ntest)

  panels <- list()
  for (col in param_cols) {
    X <- matrix(0.5, nrow = ntest, ncol = length(emu_cols),
                dimnames = list(NULL, emu_cols))
    for (g in gcm_input_cols)          X[, g] <- as.numeric(g == gcm)
    for (i in seq_along(forcing_cols)) X[, forcing_cols[i]] <- forcing_norm[i]
    X[, col] <- grid
    pr <- predict(emu, X, year)
    panels[[col]] <- data.frame(x = grid, mean = pr$mean,
                                LB = pr$mean - 2 * pr$sd,
                                UB = pr$mean + 2 * pr$sd)
  }

  structure(panels, class = c("me_result", "list"),
            label = emu$label, year = year, gcm = gcm, scenario = scenario)
}

# Panel grid of per-param sweeps, with +/-2sd ribbon and shared y-range.
# Draws to the active device; caller manages pdf()/dev.off().
plot.me_result <- function(x, ...) {
  label    <- attr(x, "label");    year     <- attr(x, "year")
  gcm      <- attr(x, "gcm");      scenario <- attr(x, "scenario")

  N    <- length(x)
  a    <- floor(sqrt(N)); b <- ceiling(N / a)
  ymax <- max(vapply(x, function(d) max(d$UB), numeric(1)))
  ymin <- min(vapply(x, function(d) min(d$LB), numeric(1)))

  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mfrow = c(a, b), mar = c(4, 4, 2, 1), oma = c(0, 0, 2.5, 0))

  for (i in seq_along(x)) {
    d   <- x[[i]]
    col <- diag_colors[((i - 1L) %% length(diag_colors)) + 1L]
    plot(d$x, d$mean, type = "n",
         xlim = c(0, 1), ylim = c(ymin, ymax),
         xlab = names(x)[i], ylab = "SLC (mm)")
    polygon(c(d$x, rev(d$x)), c(d$LB, rev(d$UB)),
            col = adjustcolor(col, 0.3), border = NA)
    lines(d$x, d$mean, col = col, lwd = 1.8)
  }

  mtext(sprintf("Main effects at %d — %s, %s, %s (other params at 0.5)",
                year, label, gcm, scenario),
        outer = TRUE, line = 0.7, cex = 1.0, font = 2)
  invisible(x)
}
