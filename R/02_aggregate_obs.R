# 02 — aggregate the observational target to a scalar (spec §3).
# IMBIE's cumulative column is already mass-loss-positive (increasing), so the
# 2007->target change is taken directly (NOT negated). Reported uncertainties are
# stored with a spurious negative sign, so magnitudes are used.

if (!exists("ensemble")) source("R/01_load_data.R")

cum_col <- grep("Cumulative.*mass.*balance.*mm", names(imbie), value = TRUE)
cum_col <- cum_col[!grepl("uncertainty", cum_col, ignore.case = TRUE)][1]
unc_col <- grep("Cumulative.*uncertainty.*mm", names(imbie), value = TRUE)[1]
stopifnot(!is.na(cum_col), !is.na(unc_col))

nearest_row <- function(df, yr) df[which.min(abs(df$Year - yr)), ]
r0 <- nearest_row(imbie, cfg$obs$baseline_year)
r1 <- nearest_row(imbie, cfg$obs$target_year)

Y_mm <- r1[[cum_col]] - r0[[cum_col]]   # observed cumulative slc, mass-loss-positive

sigma_obs_mm <- cfg$obs$sigma_obs_cum_mm
if (is.null(sigma_obs_mm)) sigma_obs_mm <- abs(r1[[unc_col]])

message(sprintf(
  "02: Y = %.3f mm  (IMBIE %.2f -> %.2f over %.2f-%.2f);  sigma_obs = %.3f mm",
  Y_mm, r0[[cum_col]], r1[[cum_col]], r0$Year, r1$Year, sigma_obs_mm))

# Forcing lookup — thermal_forcing & warming are single-valued per (gcm, scenario).
forcing_lookup <- unique(ensemble[, c(gcm_col, scen_col, forcing_cols)])
stopifnot(nrow(forcing_lookup) ==
            nrow(unique(ensemble[, c(gcm_col, scen_col)])))

# Observed ensemble target (mm), for coverage diagnostics in stage 05.
ensemble$M_mm <- ensemble[[target_col]] * 1000
