# 01 — load and clean the ensemble.

if (!exists("cfg")) source("R/00_setup.R")

param_cols     <- cfg$data$param_cols
log_cols       <- cfg$data$log_cols
forcing_cols   <- cfg$data$forcing_cols
gcm_col        <- cfg$data$gcm_col
scen_col       <- cfg$data$scenario_col
gcm_input_cols <- cfg$data$gcm_input_cols
if (is.null(gcm_input_cols)) gcm_input_cols <- character(0)

# Year int -> string to match column headers
ycol         <- function(y) as.character(y)
target_col   <- ycol(cfg$obs$target_year)
baseline_col <- ycol(cfg$obs$baseline_year)

ensemble_paths <- as.character(cfg$paths$ensemble)
n_basins       <- length(ensemble_paths)

# Derive basin labels from filenames (first underscore-delimited token).
# Single-basin (ANT_*.csv) gets the label "ANT" — there's no empty-label
# special case, so diagnostic filenames and weights columns are always
# basin-suffixed.
basin_labels <- sub("_.*$", "", basename(tools::file_path_sans_ext(ensemble_paths)))

load_one <- function(path) {
  raw <- read.csv(path, check.names = FALSE)

  # Checks:
  # - All data requested in config must exist in the ensemble
  stopifnot(all(param_cols %in% names(raw)),
            all(forcing_cols %in% names(raw)),
            gcm_col  %in% names(raw),
            scen_col %in% names(raw),
            target_col   %in% names(raw),
            baseline_col %in% names(raw))

  # - gcm_input_cols must be a subset of gcms found in the data
  gcms_in_data <- sort(unique(raw[[gcm_col]]))
  stopifnot(all(gcm_input_cols %in% gcms_in_data))

  # - Number of gcm binary flags must be either 0 or k-1 for k total gcms
  #   (k=4 in our case). Any other count makes the trend matrix
  #   rank-deficient under cbind(1, X) and rgasp will crash.
  if (!length(gcm_input_cols) %in% c(0L, length(gcms_in_data) - 1L))
    stop(sprintf(
      "gcm_input_cols must have 0 or %d entries (got %d).
    With k GCMs the trend cbind(1, X) requires k-1 dummies;
    any other count is rank-deficient.",
      length(gcms_in_data) - 1L, length(gcm_input_cols)
    ))

  # Baseline year must be ~0 sle in all ensemble members.
  stopifnot(max(abs(raw[[baseline_col]]), na.rm = TRUE) < 1e-9)

  # Drop runs that didn't complete
  out <- na.omit(raw)
  out$run_id <- seq_len(nrow(out))
  attr(out, "n_raw") <- nrow(raw)
  out
}

ensembles <- setNames(lapply(ensemble_paths, load_one), basin_labels)

ensemble <- ensembles[[1]]   # primary basin; shared input space
gcms      <- sort(unique(ensemble[[gcm_col]]))
scenarios <- sort(unique(ensemble[[scen_col]]))

# Forcing covariates are single-valued per (gcm, scenario); precompute the
# lookup table for downstream consumers (build_norm_inputs in 03, 08).
forcing_lookup <- unique(ensemble[, c(gcm_col, scen_col, forcing_cols)])
stopifnot(nrow(forcing_lookup) ==
            nrow(unique(ensemble[, c(gcm_col, scen_col)])))

n_dropped <- attr(ensemble, "n_raw") - nrow(ensemble)
message(sprintf(
  "01: loaded %d runs (%d dropped for NA) across %d basin(s).",
  nrow(ensemble), n_dropped, n_basins
))
