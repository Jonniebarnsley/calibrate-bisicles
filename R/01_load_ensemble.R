# 01 — load and clean the ensemble.

if (!exists("cfg")) source("R/00_setup.R")

ensemble_raw <- read.csv(cfg$paths$ensemble, check.names = FALSE)

param_cols     <- cfg$data$param_cols
log_cols       <- cfg$data$log_cols
forcing_cols   <- cfg$data$forcing_cols
gcm_col        <- cfg$data$gcm_col
scen_col       <- cfg$data$scenario_col
gcm_input_cols <- cfg$data$gcm_input_cols
if (is.null(gcm_input_cols)) gcm_input_cols <- character(0)

# Find all unique gcms and scenarios from the data
stopifnot(gcm_col %in% names(ensemble_raw),
          scen_col %in% names(ensemble_raw))
gcms      <- sort(unique(ensemble_raw[[gcm_col]]))
scenarios <- sort(unique(ensemble_raw[[scen_col]]))

# Year int -> string to match column headers
ycol         <- function(y) as.character(y)
target_col   <- ycol(cfg$obs$target_year)
baseline_col <- ycol(cfg$obs$baseline_year)

# Checks:
# - All data requested in config must exist in the ensemble
stopifnot(all(param_cols %in% names(ensemble_raw)),
          all(forcing_cols %in% names(ensemble_raw)),
          all(gcms %in% names(ensemble_raw)),
          target_col %in% names(ensemble_raw),
          baseline_col %in% names(ensemble_raw))

# - gcm_input_cols must be a subset of gcms
stopifnot(all(gcm_input_cols %in% gcms))

# - Number of gcm binary flags must be either 0 or k-1 for k total gcms
#   (k=4 in our case). Any other count makes the trend matrix
#   rank-deficient under cbind(1, X) and rgasp will crash.
if (!length(gcm_input_cols) %in% c(0L, length(gcms) - 1L))
  stop(sprintf(
    "gcm_input_cols must have 0 or %d entries (got %d).
    With k GCMs the trend cbind(1, X) requires k-1 dummies;
    any other count is rank-deficient.",
    length(gcms) - 1L, length(gcm_input_cols)
  ))

# Baseline year must be ~0 sle in all ensemble members.
stopifnot(max(abs(ensemble_raw[[baseline_col]]), na.rm = TRUE) < 1e-9)

# Drop runs that didn't complete
ensemble  <- na.omit(ensemble_raw)
ensemble$run_id <- seq_len(nrow(ensemble))

# Forcing covariates are single-valued per (gcm, scenario); precompute the
# lookup table for downstream consumers (build_norm_inputs in 03, 08).
forcing_lookup <- unique(ensemble[, c(gcm_col, scen_col, forcing_cols)])
stopifnot(nrow(forcing_lookup) ==
            nrow(unique(ensemble[, c(gcm_col, scen_col)])))

message(sprintf(
  "01: loaded %d runs (%d dropped for NA); %d GCMs x %d scenarios.",
  nrow(ensemble), nrow(ensemble_raw) - nrow(ensemble),
  length(unique(ensemble[[gcm_col]])), length(unique(ensemble[[scen_col]]))
))
