# Driver — emulator diagnostics only (LOO scatter + main-effects panels).
# Sources up to R/03c_emulator_diagnostics.R and stops. In per-year mode this
# trains only the diagnostic-year emulators (cfg$diagnostics$years; 2021 and
# 2300 by default). Useful when iterating on emulator inputs/kernel without
# waiting for the full projection-grid build.
#
# Run from by_continent/:  Rscript run_diagnostics.R

CONFIG_PATH <- "config.yml"
source("R/00_setup.R")
source("R/01_load_data.R")
source("R/02_aggregate_obs.R")
source("R/03_emulator_io.R")
source("R/03c_emulator_diagnostics.R")
message("Emulator diagnostics complete. PDFs in ", file.path(out_dir, "diagnostics"))
