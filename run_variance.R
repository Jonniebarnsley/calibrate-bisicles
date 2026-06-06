# Driver — variance decomposition of the prior ensemble (Seroussi/Coulon-style).
# Independent of the calibration pipeline; reuses any cached emulators.
#
# Config overrides: turns off the emulator diagnostic PDFs (LOO, main-effects)
# and the SVD r=1..9 rank-selection diagnostic — none of them affect the
# variance attribution outputs.
#
# Run from by_continent/:  Rscript run_variance.R

config_path <- "config.yml"
source("R/00_setup.R")

cfg$emulator$leave_one_out <- FALSE
cfg$emulator$main_effects  <- FALSE
cfg$svd$run_diagnostic     <- FALSE

source("R/01_load_ensemble.R")
source("R/02_load_obs.R")
source("R/03_emulate.R")
source("R/08_variance_decomposition.R")
message("Variance decomposition complete. Outputs in ", out_dir)
