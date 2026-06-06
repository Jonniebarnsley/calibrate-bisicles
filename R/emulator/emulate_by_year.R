# Trains and caches one rgasp model per year.
# Sourced from emulate_basin.R; exposes build_peryear_emulators(ens, cache_dir)
# which returns list(train, predict). Side-effects (LOO PDFs) are triggered
# from emulate_basin.R, not here.

stopifnot(exists("X_train"))

# build_peryear_emulators: returns list(train = fn, predict = fn).
# train(year) lazily fits or loads the GP for that year's cumulative SLC column.
# Disk cache (clear cache_dir if the `emulator:` block changes) + in-memory
# memoisation so repeated calls within a session are instant.
build_peryear_emulators <- function(ens, cache_dir) {
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  .mem <- new.env(parent = emptyenv())

  train_fn <- function(year, cache = cfg$projection$cache_emulators) {
    key <- as.character(year)

    # If the emulator exists in memory already, load it
    if (!is.null(.mem[[key]])) return(.mem[[key]])

    # If the emulator is not in memory but is cached on disk, load it into memory
    f <- file.path(cache_dir, sprintf("emulator_%d.rds", year))
    if (cache && file.exists(f)) {
      m <- readRDS(f)
      .mem[[key]] <- m
      return(m)
    }

    # Otherwise, build an emulator using rgasp
    invisible(utils::capture.output(            # silence the optimiser chatter
      m <- rgasp(
        design      = X_train,
        response    = ens[[ycol(year)]],
        trend       = trend,
        kernel_type = cfg$emulator$kernel_type,
        alpha       = cfg$emulator$alpha,
        nugget.est  = cfg$emulator$nugget_est
      )
    ))

    # Save it to disk if cacheing is enabled in config
    if (cache) saveRDS(m, f)
    .mem[[key]] <- m
    m
  }

  # Predict slc in mm for a given input set and year(s).
  predict_fn <- function(X, year) {

    # If only predicting one year, return a list with mean and sd.
    if (length(year) == 1L) {
      pr <- predict(train_fn(year), X, testing_trend = cbind(1, X))
      return(list(mean = pr$mean * 1000, sd = pr$sd * 1000))
    }

    # Otherwise (for multiple years) build and populate matrices for mean and sd.
    M <- S <- matrix(0, nrow(X), length(year))
    for (j in seq_along(year)) {
      pr <- predict(train_fn(year[j]), X, testing_trend = cbind(1, X))
      M[, j] <- pr$mean * 1000
      S[, j] <- pr$sd * 1000
    }
    list(mean = M, sd = S)
  }

  list(train = train_fn, predict = predict_fn)
}
