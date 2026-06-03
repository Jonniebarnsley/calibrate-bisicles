# Driver — forcing-input tests WITHOUT GCM dummy switches.
# Same 9 combinations as run_forcing_tests.R, but each emulator is trained on
# 6 params + forcing only (no 3 GCM dummies). This tests whether the dummies are
# actually adding information beyond what the forcing variables already encode.
# Outputs go to outputs/forcing_tests_no_gcm/ to keep them separate.
#
# Run from by_continent/:  Rscript run_forcing_tests_no_gcm.R

CONFIG_PATH <- "config.yml"
source("R/00_setup.R")
source("R/01_load_data.R")
source("R/02_aggregate_obs.R")
source("R/03d_forcing_input_tests.R")
run_forcing_input_tests(include_dummies = FALSE)
message("Forcing-input tests (no GCM dummies) complete. Outputs in ",
        file.path(out_dir, "forcing_tests_no_gcm"))
