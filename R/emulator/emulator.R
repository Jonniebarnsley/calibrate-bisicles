# Emulator base: S3 generics for diagnostics + a constructor that dispatches
# to the SVD or per-year subclass based on cfg$svd$enabled. Sourced from
# 03_emulate.R before the relevant mode file (emulator_svd.R or
# emulator_peryear.R), which defines the predict method for its subclass.
#
# Class layout:
#   emulator_svd      <- c("emulator_svd",     "emulator")
#   emulator_peryear  <- c("emulator_peryear", "emulator")
#
# Verbs:
#   predict(emu, X, year)                 # mode-specific (predict.emulator_<mode>)
#   leave_one_out(emu, year)              # mode-specific, returns loo_result
#   main_effects(emu, year, gcm, scen)    # mode-agnostic, returns me_result

leave_one_out <- function(emu, ...) UseMethod("leave_one_out")
main_effects  <- function(emu, ...) UseMethod("main_effects")

build_emulator <- function(ens, label) {
  if (cfg$svd$enabled) {
    build_emulator_svd(ens, label)
  } else {
    build_emulator_peryear(ens, label)
  }
}
