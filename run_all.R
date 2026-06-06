# Driver — full pipeline end to end.
# Run from pipline/:  Rscript run_all.R

config_path <- "config.yml"
source("R/00_setup.R")
source("R/01_load_ensemble.R")
source("R/02_load_obs.R")
source("R/03_emulate.R")
source("R/04_weights.R")
source("R/05_diagnostics.R")
source("R/06_project.R")
source("R/07_sensitivity.R")
source("R/08_variance_decomposition.R")
message("Pipeline complete. Outputs in ", out_dir)