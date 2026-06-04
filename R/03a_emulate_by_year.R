# 03a — Trains and caches one rgasp model per year.
# Sourced from 03_emulator.R only when cfg$svd$enabled is false.

stopifnot(exists("X_train"))

emu_cache_dir <- file.path(out_dir, "emulators")
if (!dir.exists(emu_cache_dir)) dir.create(emu_cache_dir, recursive = TRUE)

# Train (or load) the per-year GP for a year's cumulative SLC column. Disk cache
# (clear the dir if the `emulator:` block changes) + in-memory memoisation so
# repeated calls within a session are instant.
.emu_mem <- new.env(parent = emptyenv())
train_emulator <- function(year, cache = cfg$projection$cache_emulators) {
  key <- as.character(year)

  # If the emulator exists in memory already, load it
  if (!is.null(.emu_mem[[key]])) return(.emu_mem[[key]])

  # If the emulator is not in memory but is cached on disk, load it into memory
  f <- file.path(emu_cache_dir, sprintf("emulator_%d.rds", year))
  if (cache && file.exists(f)) {
    m <- readRDS(f)
    .emu_mem[[key]] <- m
    return(m)
  }

  # Otherwise, build an emulator using rgasp
  invisible(utils::capture.output(            # silence the optimiser chatter
    m <- rgasp(
      design      = X_train,
      response    = ensemble[[ycol(year)]],
      trend       = trend,
      kernel_type = cfg$emulator$kernel_type,
      alpha       = cfg$emulator$alpha,
      nugget.est  = cfg$emulator$nugget_est
    )
  ))

  # Save it to disk if cacheing is enabled in config
  if (cache) saveRDS(m, f)
  .emu_mem[[key]] <- m
  m
}

# Predict slc in mm for a given input set and year(s).
peryear_predict_mm <- function(X, year) {

  # If only predicting one year, return a list with mean and sd.
  if (length(year) == 1L) {
    pr <- predict(train_emulator(year), X, testing_trend = cbind(1, X))
    return(list(mean = pr$mean * 1000, sd = pr$sd * 1000))
  }

  # Otherwise (for multiple years) build and populate matrices for mean and sd.
  M <- S <- matrix(0, nrow(X), length(year))
  for (j in seq_along(year)) {
    pr <- predict(train_emulator(year[j]), X, testing_trend = cbind(1, X))
    M[, j] <- pr$mean * 1000
    S[, j] <- pr$sd * 1000
  }
  list(mean = M, sd = S)
}

# Runs when sourced
# Check for inputs to which the outputs at the obs target year are insensitive.
findInertInputs(train_emulator(cfg$obs$target_year))

# LOO PDFs at the configured diagnostic years. Per-year LOO is just
# leave_one_out_rgasp on each year's trained emulator. write_loo_pdf and
# diag_dir live in 03_emulator.R (sourced before this file).
if (isTRUE(cfg$emulator$leave_one_out)) for (year in cfg$emulator$diagnostic_years) {
  m   <- train_emulator(year)
  loo <- leave_one_out_rgasp(m)
  write_loo_pdf(data.frame(actual = ensemble[[ycol(year)]] * 1000,
                           mean   = loo$mean * 1000,
                           sd     = loo$sd   * 1000),
                year)
  message(sprintf("03a: wrote loo_%d.pdf", year))
}
