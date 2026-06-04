# 02 — load IMBIE observations

if (!exists("cfg")) source("R/00_setup.R")

imbie <- read.csv(cfg$paths$imbie, check.names = FALSE)

cum_col <- "Cumulative mass balance (mm)"
unc_col <- "Cumulative mass balance uncertainty (mm)"

# Get nearest row in IMBIE. Searching for 2021 in Otosaka et al. (2023) data
# returns the 2020.92 row since the record ends there.
nearest_row <- function(df, yr) df[which.min(abs(df$Year - yr)), ]
r0 <- nearest_row(imbie, cfg$obs$baseline_year)
r1 <- nearest_row(imbie, cfg$obs$target_year)

Y_mm <- r1[[cum_col]] - r0[[cum_col]] # observed cumulative slc

# Cumulative obs uncertainty. IMBIE constructs the cumulative error as a
# root-sum-square of independent monthly errors,
# so Var(C(t1) - C(t0)) = sigma_t1^2 - sigma_t0^2.
sigma_obs_mm <- sqrt(r1[[unc_col]]^2 - r0[[unc_col]]^2)

message(sprintf(
  "02: Y = %.3f mm  (IMBIE %.2f -> %.2f over %.2f-%.2f);  sigma_obs = %.3f mm",
  Y_mm, r0[[cum_col]], r1[[cum_col]], r0$Year, r1$Year, sigma_obs_mm
))
