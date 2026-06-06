# Per-year GP emulator: one rgasp model per output year, trained lazily on
# first request and cached in-memory + on disk. The on-disk cache lives at
# outputs/emulators/<label>/emulator_<year>.rds — delete the basin's
# subdirectory when anything that affects training (emulator: block, forcing
# cols, input layout) changes.
#
# Provides:
#   func build_emulator_peryear(ens, label) - constructor
#   predict.emulator_peryear(...)           - S3 method
#   fetch_peryear_gp(emu, year)             - shared helper (also used by LOO)

stopifnot(exists("X_train"))

build_emulator_peryear <- function(ens, label) {
  cache_dir <- file.path(out_dir, "emulators", label)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  emu <- list(
    label     = label,
    ens       = ens,
    cache_dir = cache_dir,
    cache     = new.env(parent = emptyenv())   # in-memory memoisation (reference)
  )
  class(emu) <- c("emulator_peryear", "emulator")
  emu
}

# Return the fitted rgasp model for `year`, fitting it if necessary. Uses
# the emulator's in-memory env first, then the on-disk cache (if enabled),
# and only fits from scratch as a last resort. The cache env has reference
# semantics, so writes persist across calls even though emu is passed by
# value.
fetch_peryear_gp <- function(emu, year) {
  key <- as.character(year)

  # In-memory hit
  if (!is.null(emu$cache[[key]])) return(emu$cache[[key]])

  # On-disk hit
  f <- file.path(emu$cache_dir, sprintf("emulator_%d.rds", year))
  if (isTRUE(cfg$projection$cache_emulators) && file.exists(f)) {
    m <- readRDS(f)
    emu$cache[[key]] <- m
    return(m)
  }

  # Fit from scratch
  invisible(utils::capture.output(            # silence the optimiser chatter
    m <- rgasp(
      design      = X_train,
      response    = emu$ens[[ycol(year)]],
      trend       = trend,
      kernel_type = cfg$emulator$kernel_type,
      alpha       = cfg$emulator$alpha,
      nugget.est  = cfg$emulator$nugget_est
    )
  ))

  if (isTRUE(cfg$projection$cache_emulators)) saveRDS(m, f)
  emu$cache[[key]] <- m
  m
}

# Predict SLC mean & sd (mm) at year(s) for a normalised input matrix X.
predict.emulator_peryear <- function(object, X, year, ...) {
  if (length(year) == 1L) {
    m  <- fetch_peryear_gp(object, year)
    pr <- predict(m, X, testing_trend = cbind(1, X))
    return(list(mean = pr$mean * 1000, sd = pr$sd * 1000))
  }

  M <- S <- matrix(0, nrow(X), length(year))
  for (j in seq_along(year)) {
    m  <- fetch_peryear_gp(object, year[j])
    pr <- predict(m, X, testing_trend = cbind(1, X))
    M[, j] <- pr$mean * 1000
    S[, j] <- pr$sd * 1000
  }
  list(mean = M, sd = S)
}
