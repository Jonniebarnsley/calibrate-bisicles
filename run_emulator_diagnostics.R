# Driver — emulator diagnostics only (LOO scatter + main-effects panels).
# Builds the emulator and stops; LOO + main-effects PDFs come out as side
# effects of sourcing 03_emulate.R. Useful when iterating on emulator inputs/
# kernel without waiting for calibration or projection.
#
# Config overrides: forces LOO + main-effects on (so the run is useful even if
# they're off in config); turns off the SVD r=1..9 rank-selection diagnostic
# (separate concern from emulator goodness-of-fit).
#
# Run from by_continent/:  Rscript run_diagnostics.R

config_path <- "config.yml"
source("R/00_setup.R")

cfg$emulator$leave_one_out <- TRUE
cfg$emulator$main_effects  <- TRUE
cfg$svd$run_diagnostic     <- FALSE

source("R/01_load_ensemble.R")
source("R/02_load_obs.R")
source("R/03_emulate.R")
message("Emulator diagnostics complete. PDFs in ", file.path(out_dir, "diagnostics"))
