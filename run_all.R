# Driver — calibration + projection + sensitivity, end to end (spec §1-§8).
# Run from the calibration/ directory:  Rscript run_all.R

CONFIG_PATH <- "config.yml"
source("R/00_setup.R")
source("R/01_load_data.R")
source("R/02_aggregate_obs.R")
source("R/03_emulator_io.R")
source("R/03c_emulator_diagnostics.R")
source("R/04_weights.R")
source("R/05_diagnostics.R")
source("R/06_project.R")
source("R/07_sensitivity.R")
message("Pipeline complete. Outputs in ", out_dir)
