# 02 — load IMBIE observations

if (!exists("cfg"))      source("R/00_setup.R")
if (!exists("n_basins")) source("R/01_load_ensemble.R")

cum_col <- "Cumulative mass balance (mm)"
unc_col <- "Cumulative mass balance uncertainty (mm)"

# Get nearest row in IMBIE. Searching for 2021 in Otosaka et al. (2023) data
# returns the 2020.92 row since the record ends there.
nearest_row <- function(df, yr) df[which.min(abs(df$Year - yr)), ]

load_one_obs <- function(path, label = "") {
  imbie <- read.csv(path, check.names = FALSE)
  r0    <- nearest_row(imbie, cfg$obs$baseline_year)
  r1    <- nearest_row(imbie, cfg$obs$target_year)

  Y_mm <- r1[[cum_col]] - r0[[cum_col]] # observed cumulative slc

  # Cumulative obs uncertainty. IMBIE constructs the cumulative error as a
  # root-sum-square of independent monthly errors,
  # so Var(C(t1) - C(t0)) = sigma_t1^2 - sigma_t0^2.
  sigma_obs_mm <- sqrt(r1[[unc_col]]^2 - r0[[unc_col]]^2)

  message(sprintf(
    "02%s: Y = %.3f mm  (IMBIE %.2f -> %.2f over %.2f-%.2f);  sigma_obs = %.3f mm",
    if (nchar(label) > 0) paste0(" [", label, "]") else "",
    Y_mm, r0[[cum_col]], r1[[cum_col]], r0$Year, r1$Year, sigma_obs_mm
  ))

  list(Y_mm = Y_mm, sigma_obs_mm = sigma_obs_mm)
}

imbie_paths <- as.character(cfg$paths$imbie)
stopifnot(length(imbie_paths) == n_basins)


obs_list          <- setNames(
  Map(load_one_obs, imbie_paths, basin_labels),
  basin_labels
)
Y_mm_list         <- lapply(obs_list, `[[`, "Y_mm")
sigma_obs_mm_list <- lapply(obs_list, `[[`, "sigma_obs_mm")

# Scalars from the primary basin for stages that only need a single obs
# (e.g. stage 05 coverage check).
Y_mm         <- Y_mm_list[[1]]
sigma_obs_mm <- sigma_obs_mm_list[[1]]
