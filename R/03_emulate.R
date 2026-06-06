# 03 — emulator trunk: shared infrastructure + basin dispatch.
# For a single basin (cfg$paths$ensemble is a scalar string) this is identical
# to the previous single-basin behaviour. For multiple basins it calls
# build_emulator() once per basin and exposes predict_slc_mm_list.
# predict_slc_mm always points to the first (primary) basin's predict function
# so downstream stages that only need one emulator are unaffected.

if (!exists("Y_mm"))    source("R/02_load_obs.R")
if (!exists("ensemble")) source("R/01_load_ensemble.R")

emu_cols <- c(param_cols, gcm_input_cols, forcing_cols)
message("03: input order -> ", paste(emu_cols, collapse = ", "))

# Save normalisation specs for inverting later
norm_specs <- list()
for (col in emu_cols) {
  norm_specs[[col]] <- norm_fit(ensemble[[col]], log = col %in% log_cols)
}

# Raw data.frame -> normalised emulator-input matrix (using norm_specs)
normalise_df <- function(df) {
  sapply(emu_cols, function(col) norm_apply(df[[col]], norm_specs[[col]]))
}

# Helper function to build normalised emulator inputs for a given set of
# parameter samples, gcm, and scenario. Useful later when predicting large
# Monte Carlo sample.
build_norm_inputs <- function(params_raw, gcm, scenario) {
  fl <- forcing_lookup[forcing_lookup[[gcm_col]] == gcm &
                         forcing_lookup[[scen_col]] == scenario, ]
  stopifnot(nrow(fl) == 1)
  df <- params_raw
  for (col in gcm_input_cols) df[[col]] <- as.integer(col == gcm)
  for (col in forcing_cols)   df[[col]] <- fl[[col]]
  normalise_df(df)
}

# Emulator inputs and trend
X_train <- normalise_df(ensemble)
trend   <- cbind(1, X_train)

# Functions to make diagnostic plots investigating emulator accuracy:
#  - leave-one-out
#  - main-effects
# Outputs go to outputs/diagnostics/.
diag_dir <- file.path(out_dir, "diagnostics")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

# LOO scatter PDF: actual vs predicted, with +/-2sd error bars and headline
# metrics (RMSE, pass-rate, normalised error distance).
write_loo_pdf <- function(df, year, label) {
  df$UB <- df$mean + 2 * df$sd
  df$LB <- df$mean - 2 * df$sd
  df$ok <- df$actual > df$LB & df$actual < df$UB

  # Subset df into pass and fail samples
  good  <- df[df$ok, ]
  bad <- df[!df$ok, ]

  # Compute headline statistics
  rmse  <- sqrt(mean((df$actual - df$mean)^2))
  ned   <- sqrt(sum((df$actual - df$mean)^2 / df$sd^2))
  pass  <- mean(df$ok) * 100
  lim   <- range(df$LB, df$UB)

  # Make and save plot
  pdf(file.path(diag_dir, sprintf("loo_%d_%s.pdf", year, label)),
      width = 4.6, height = 4.6)
  par(mar = c(4, 4, 2.5, 1))
  plot(good$actual, good$mean,
       xlim = lim, ylim = lim,
       xlab = "Actual SLC (mm)", ylab = "Predicted SLC (mm)",
       main = sprintf("LOO — %d", year),
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

# Main-effects PDF: each non-dummy input is swept on [0,1] with all others
# (including GCM dummies) fixed at 0.5. predict_fn is passed explicitly so
# the correct basin's emulator is used.
diag_colors <- c("orangered", "seagreen3", "royalblue1", "plum2",
                 "slateblue1", "sienna1", "goldenrod1", "darkseagreen")

write_main_effects_pdf <- function(year, label, predict_fn) {
  sweep_cols <- setdiff(emu_cols, gcm_input_cols)
  ntest      <- 100
  grid       <- seq(0, 1, length.out = ntest)

  # Build a list of (x, mean, LB, UB) data frames, one per swept column.
  me <- list()
  for (col in sweep_cols) {
    X <- matrix(0.5, nrow = ntest, ncol = length(emu_cols),
                dimnames = list(NULL, emu_cols))
    X[, col] <- grid
    pr <- predict_fn(X, year) # predict function depends on svd.enabled
    me[[col]] <- data.frame(x = grid, mean = pr$mean,
                            LB = pr$mean - 2 * pr$sd,
                            UB = pr$mean + 2 * pr$sd)
  }

  # Panel grid + shared y-range across panels.
  N    <- length(me)
  a    <- floor(sqrt(N)); b <- ceiling(N / a)
  ymax <- max(vapply(me, function(d) max(d$UB), numeric(1)))
  ymin <- min(vapply(me, function(d) min(d$LB), numeric(1)))

  # Make and save plot
  pdf(file.path(diag_dir, sprintf("main_effects_%d_%s.pdf", year, label)),
      width = 3.2 * b, height = 2.8 * a)
  par(mfrow = c(a, b), mar = c(4, 4, 2, 1), oma = c(0, 0, 2.5, 0))
  for (i in seq_along(me)) {
    d   <- me[[i]]
    col <- diag_colors[((i - 1L) %% length(diag_colors)) + 1L]
    plot(d$x, d$mean, type = "n",
         xlim = c(0, 1), ylim = c(ymin, ymax),
         xlab = names(me)[i], ylab = "SLC (mm)")
    polygon(c(d$x, rev(d$x)), c(d$LB, rev(d$UB)),
            col = adjustcolor(col, 0.3), border = NA)
    lines(d$x, d$mean, col = col, lwd = 1.8)
  }
  mtext(sprintf("Main effects at %d (all other inputs fixed at 0.5)", year),
        outer = TRUE, line = 0.7, cex = 1.0, font = 2)
  dev.off()
  message(sprintf("03: wrote main_effects_%d_%s.pdf", year, label))
}

# Source the mode-specific builder and the basin dispatcher.
# emulate_by_svd/year expose build_svd_state / build_peryear_emulators.
# emulate_basin exposes build_emulator().
if (isFALSE(cfg$svd$enabled)) {
  source("R/emulator/emulate_by_year.R")
} else {
  source("R/emulator/emulate_by_svd.R")
}
source("R/emulator/emulate_basin.R")

# Build one emulator per basin. predict_slc_mm always points to the first
# (primary) basin for stages that only need one emulator (projection, variance
# decomposition); predict_slc_mm_list holds all basins for stages that need
# the joint likelihood (weighting).
predict_slc_mm_list <- setNames(
  lapply(seq_along(basin_labels),
         function(i) build_emulator(ensembles[[i]], basin_labels[i])),
  basin_labels
)
predict_slc_mm <- predict_slc_mm_list[[1]]

if (isTRUE(cfg$emulator$main_effects))
  for (year in cfg$emulator$diagnostic_years)
    for (i in seq_along(predict_slc_mm_list))
      write_main_effects_pdf(year, basin_labels[i], predict_slc_mm_list[[i]])
