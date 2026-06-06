# Driver — variance decomposition of the prior ensemble (Seroussi/Coulon-style).
# Run from the calibration/ directory:  Rscript run_variance.R
# Independent of the calibration pipeline (run_all.R); reuses the cached emulators.

config_path <- "config.yml"
source("R/00_setup.R")
source("R/03_emulate.R")
source("R/08_variance_decomposition.R")
message("Variance decomposition complete. Outputs in ", out_dir)
