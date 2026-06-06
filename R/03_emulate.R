# 03 — emulator control: builds per-basin emulators and writes optional
# diagnostic PDFs. Mode is chosen by cfg$svd$enabled.

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

# All basins share the same inputs. Sufficient to compute the inputs for
# just one emulator and use them for all basin emulators.
X_train <- normalise_df(ensemble)
trend   <- cbind(1, X_train)

# Source the emulator base (generics + dispatcher) and the relevant mode
# file (which defines the predict method for its subclass).
source("R/emulator/emulator.R")
if (isTRUE(cfg$svd$enabled)) {
  source("R/emulator/emulator_svd.R")
} else {
  source("R/emulator/emulator_peryear.R")
}

# Build one emulator per basin.
emulators <- setNames(
  lapply(seq_along(basin_labels),
         function(i) build_emulator(ensembles[[i]], basin_labels[i])),
  basin_labels
)

# Diagnostic PDFs go to outputs/diagnostics/.
diag_dir <- file.path(out_dir, "diagnostics")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

# Leave-one-out PDFs, if requested.
if (isTRUE(cfg$emulator$leave_one_out)) {
  source("R/emulator/leave_one_out.R")
  for (emu in emulators)
    for (year in cfg$emulator$diagnostic_years) {
      loo <- leave_one_out(emu, year)
      fname <- sprintf("loo_%d_%s.pdf", year, emu$label)
      pdf(file.path(diag_dir, fname), width = 4.6, height = 4.6)
      plot(loo)
      dev.off()
      message(sprintf("03: wrote %s", fname))
    }
}

# Main-effects PDFs, if requested. Anchored on-manifold at the configured
# (main_effects_gcm, main_effects_scenario) — see emulator/main_effects.R.
if (isTRUE(cfg$emulator$main_effects)) {
  source("R/emulator/main_effects.R")
  me_gcm      <- cfg$emulator$main_effects_gcm
  me_scenario <- cfg$emulator$main_effects_scenario
  np <- length(param_cols)
  a  <- floor(sqrt(np))
  b <- ceiling(np / a)
  for (emu in emulators)
    for (year in cfg$emulator$diagnostic_years) {
      me <- main_effects(emu, year, me_gcm, me_scenario)
      fname <- sprintf("main_effects_%d_%s_%s_%s.pdf",
                       year, emu$label, me_gcm, me_scenario)
      pdf(file.path(diag_dir, fname), width = 3.2 * b, height = 2.8 * a)
      plot(me)
      dev.off()
      message(sprintf("03: wrote %s", fname))
    }
}

# Optional SVD rank-selection diagnostic (writes one CSV per basin).
if (isTRUE(cfg$svd$enabled) && isTRUE(cfg$svd$run_diagnostic))
  for (emu in emulators) run_svd_diagnostic(emu)
