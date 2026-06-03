# Driver — forcing-input tests (9 emulator configurations).
# Trains fresh in-memory emulators for each forcing-col combination and produces
# LOO + main-effects diagnostics per (test, year). Independent of the production
# emulator setup; does NOT touch outputs/emulators/.
#
# Run from by_continent/:  Rscript run_forcing_tests.R

CONFIG_PATH <- "config.yml"
source("R/00_setup.R")
source("R/01_load_data.R")
source("R/02_aggregate_obs.R")
source("R/03d_forcing_input_tests.R")
# Optional command-line arguments:
#   [1] GCM         (default = ref_gcm; "all" = all 4 GCMs; comma-separated supported)
#   [2] scenario    (default = ssp585)
# When multiple GCMs are requested, the emulators are trained once and reused
# across all GCMs (one main_effects PDF per GCM). Examples:
#   Rscript run_forcing_tests.R                              -> EC-Earth3-Veg, ssp585
#   Rscript run_forcing_tests.R CESM2-WACCM                  -> CESM2, ssp585
#   Rscript run_forcing_tests.R all ssp585                   -> all 4 GCMs, ssp585
#   Rscript run_forcing_tests.R CESM2-WACCM,MIROC-ES2L       -> two GCMs, ssp585
.args      <- commandArgs(trailingOnly = TRUE)
.gcm_arg   <- if (length(.args) > 0) .args[1] else ref_gcm
me_scenario <- if (length(.args) > 1) .args[2] else "ssp585"
me_gcm     <- if (identical(.gcm_arg, "all")) {
                c(ref_gcm, dummy_gcms)
              } else {
                strsplit(.gcm_arg, ",")[[1]]
              }
run_forcing_input_tests(include_dummies = TRUE,
                        me_gcm = me_gcm,
                        me_scenario = me_scenario)
message("Forcing-input tests complete. Outputs in ", file.path(out_dir, "forcing_tests"))
