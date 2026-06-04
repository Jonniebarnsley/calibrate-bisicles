# Driver — emulator diagnostics only (LOO scatter + main-effects panels).
# Builds the emulator and stops; LOO + main-effects PDFs come out as side
# effects of sourcing 03a/03b/03 when the cfg$emulator flags are set. Useful
# when iterating on emulator inputs/kernel without waiting for projection.
#
# Run from by_continent/:  Rscript run_diagnostics.R

config_path <- "config.yml"
source("R/00_setup.R")
source("R/01_load_ensemble.R")
source("R/02_load_obs.R")
source("R/03_emulator.R")
message("Emulator diagnostics complete. PDFs in ", file.path(out_dir, "diagnostics"))
