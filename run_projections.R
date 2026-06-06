# Driver — calibration + projection only (weights, calibration diagnostics,
# projection bands, sensitivity sweep). Stops before variance decomposition.
#
# Config overrides: turns off the emulator diagnostic PDFs (LOO, main-effects)
# and the SVD r=1..9 rank-selection diagnostic — those are emulator-quality
# checks, not part of producing projection outputs.
#
# Run from by_continent/:  Rscript run_projections.R

config_path <- "config.yml"
source("R/00_setup.R")

cfg$emulator$leave_one_out <- FALSE
cfg$emulator$main_effects  <- FALSE
cfg$svd$run_diagnostic     <- FALSE

source("R/01_load_ensemble.R")
source("R/02_load_obs.R")
source("R/03_emulate.R")
source("R/04_weights.R")
source("R/05_diagnostics.R")
source("R/06_project.R")
source("R/07_sensitivity.R")
message("Projections complete. Outputs in ", out_dir)
