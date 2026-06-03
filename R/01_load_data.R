# 01 — load ensemble + IMBIE and enforce the data contract (spec §2).
# Confirmed: ensemble slc is CUMULATIVE, in metres, baselined to 2007 = 0,
# with mass loss POSITIVE (same convention as IMBIE's cumulative column).

if (!exists("cfg")) source("R/00_setup.R")

ensemble_raw <- read.csv(cfg$paths$ensemble)
imbie        <- read.csv(cfg$paths$imbie)

param_cols   <- cfg$data$param_cols
log_cols     <- cfg$data$log_cols
forcing_cols <- cfg$data$forcing_cols
gcm_col      <- cfg$data$gcm_col
scen_col     <- cfg$data$scenario_col
gcms         <- cfg$data$gcms
ref_gcm      <- cfg$data$gcm_reference

# Year columns are read as X#### (e.g. X2021).
ycol         <- function(y) paste0("X", y)
target_col   <- ycol(cfg$obs$target_year)
baseline_col <- ycol(cfg$obs$baseline_year)

stopifnot(all(param_cols %in% names(ensemble_raw)),
          all(forcing_cols %in% names(ensemble_raw)),
          target_col %in% names(ensemble_raw),
          baseline_col %in% names(ensemble_raw))

# Baseline year must be ~0 everywhere (cumulative slc referenced to it).
stopifnot(max(abs(ensemble_raw[[baseline_col]]), na.rm = TRUE) < 1e-9)

# GCM dummy switches; the reference level is dropped (dummy, not one-hot).
dummy_gcms <- setdiff(gcms, ref_gcm)
for (g in dummy_gcms) ensemble_raw[[g]] <- as.integer(ensemble_raw[[gcm_col]] == g)

# Drop runs with NA in any model input or the target year.
need     <- c(param_cols, forcing_cols, target_col)
ensemble <- ensemble_raw[stats::complete.cases(ensemble_raw[, need]), ]
ensemble$run_id <- seq_len(nrow(ensemble))

message(sprintf(
  "01: loaded %d runs (%d dropped for NA); %d GCMs x %d scenarios.",
  nrow(ensemble), nrow(ensemble_raw) - nrow(ensemble),
  length(unique(ensemble[[gcm_col]])), length(unique(ensemble[[scen_col]]))))
