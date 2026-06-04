# 03 — obs-period emulator: train GP, build input vectors, predict (spec §5, §6).
# Emulator inputs (11, fixed order): 6 params + 3 GCM dummies + thermal_forcing +
# warming. Forcing covariates are included because they improve predictive skill
# and are pure (gcm, scenario) lookups. Trained on the target-year cumulative
# column (metres); the data is already cumulative-from-2007, so this directly
# emulates the aggregate quantity — no summing of annual emulators (spec §6.2).

if (!exists("Y_mm")) source("R/02_aggregate_obs.R")

emu_cols <- c(param_cols, dummy_gcms, forcing_cols)

# Normalisation fit on training data (dummies are 0/1 and left untouched).
norm_specs <- list()
for (col in c(param_cols, forcing_cols)) {
  norm_specs[[col]] <- norm_fit(ensemble[[col]], log = col %in% log_cols)
}

# Raw data.frame -> normalised emulator-input matrix (column order = emu_cols).
design_df <- function(df) {
  out <- matrix(0, nrow = nrow(df), ncol = length(emu_cols),
                dimnames = list(NULL, emu_cols))
  for (col in param_cols)   out[, col] <- norm_apply(df[[col]], norm_specs[[col]])
  for (col in dummy_gcms)   out[, col] <- as.numeric(df[[col]])
  for (col in forcing_cols) out[, col] <- norm_apply(df[[col]], norm_specs[[col]])
  out
}

X_train <- design_df(ensemble)             # input design — identical for every year
trend   <- cbind(1, X_train)

emu_cache_dir <- file.path(out_dir, "emulators")
if (!dir.exists(emu_cache_dir)) dir.create(emu_cache_dir, recursive = TRUE)

# Train (or load) the per-year GP for a year's cumulative SLC column. Disk cache
# (clear the dir if the `emulator:` block changes) + in-memory memoisation so
# repeated calls within a session are instant.
.emu_mem <- new.env(parent = emptyenv())
train_emulator <- function(year, cache = isTRUE(cfg$projection$cache_emulators)) {
  key <- as.character(year)
  if (!is.null(.emu_mem[[key]])) return(.emu_mem[[key]])
  f <- file.path(emu_cache_dir, sprintf("emulator_%d.rds", year))
  if (cache && file.exists(f)) { m <- readRDS(f); .emu_mem[[key]] <- m; return(m) }
  resp <- ensemble[[ycol(year)]]
  ok   <- is.finite(resp)
  m <- NULL
  invisible(utils::capture.output(            # silence the optimiser chatter
    m <- rgasp(
      design      = X_train[ok, , drop = FALSE],
      response    = resp[ok],
      trend       = trend[ok, , drop = FALSE],
      kernel_type = cfg$emulator$kernel_type,
      alpha       = cfg$emulator$alpha,
      nugget.est  = cfg$emulator$nugget_est
    )))
  if (cache) saveRDS(m, f)
  .emu_mem[[key]] <- m
  m
}

# Emulator mode: per-year GPs (default) or the SVD basis-function emulator.
if (isTRUE(cfg$svd$enabled)) {
  source("R/svd_emulator.R")                  # builds basis + coefficient emulators (+ diagnostic)
} else {
  emulator <- train_emulator(cfg$obs$target_year)   # obs-period GP
  # [CHECK §5] obs-period sensitivity to inputs (thermal_forcing/warming expected inert).
  message("03: input order -> ", paste(emu_cols, collapse = ", "))
  invisible(findInertInputs(emulator))
}

# Unified predictor: SLC mean & sd (mm) for a NORMALISED input matrix, at one year
# (returns vectors) or several (returns K x length(year) matrices). Dispatches on the
# emulator mode. In SVD mode a year vector reuses one coefficient prediction across
# all years, which is the whole point of the basis emulator — so callers that span
# many years should pass the full vector rather than looping year by year.
predict_slc_mm <- function(Xnorm, year) {
  if (isTRUE(cfg$svd$enabled)) return(svd_predict_mm(Xnorm, year))
  if (length(year) == 1L) {
    pr <- predict(train_emulator(year), Xnorm, testing_trend = cbind(1, Xnorm))
    return(list(mean = pr$mean * 1000, sd = pr$sd * 1000))
  }
  M <- S <- matrix(0, nrow(Xnorm), length(year))
  for (j in seq_along(year)) {
    pr <- predict(train_emulator(year[j]), Xnorm, testing_trend = cbind(1, Xnorm))
    M[, j] <- pr$mean * 1000; S[, j] <- pr$sd * 1000
  }
  list(mean = M, sd = S)
}

# Build a normalised input matrix from normalised parameter samples (K x 6, in
# [0,1], columns named param_cols) for one GCM stratum + scenario: sets the GCM
# dummies and the (gcm, scenario) forcing lookup, normalised consistently.
build_norm_inputs <- function(params_norm, gcm, scenario) {
  K <- nrow(params_norm)
  X <- matrix(0, nrow = K, ncol = length(emu_cols),
              dimnames = list(NULL, emu_cols))
  X[, param_cols] <- as.matrix(params_norm)[, param_cols]
  for (g in dummy_gcms) X[, g] <- as.numeric(g == gcm)
  fl <- forcing_lookup[forcing_lookup[[gcm_col]] == gcm &
                         forcing_lookup[[scen_col]] == scenario, ]
  stopifnot(nrow(fl) == 1)
  for (col in forcing_cols) X[, col] <- norm_apply(fl[[col]], norm_specs[[col]])
  X
}

# Raw parameter data.frame -> normalised parameter data.frame (training space),
# columns named param_cols. Used to push stored posterior samples back through the
# emulator (projection §7, sensitivity §8).
norm_params <- function(df) {
  out <- as.data.frame(lapply(param_cols,
                              function(c) norm_apply(df[[c]], norm_specs[[c]])))
  names(out) <- param_cols
  out
}
