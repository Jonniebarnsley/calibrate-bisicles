# Driver — calibration + projection + sensitivity, end to end (spec §1-§8).
# Run from the calibration/ directory:  Rscript run_all.R

config_path <- "config.yml"
source("R/00_setup.R")
source("R/01_load_ensemble.R")
source("R/02_load_obs.R")
source("R/03_emulate.R")
source("R/04_weights.R")
source("R/05_diagnostics.R")
source("R/06_project.R")
source("R/07_sensitivity.R")
message("Pipeline complete. Outputs in ", out_dir)
