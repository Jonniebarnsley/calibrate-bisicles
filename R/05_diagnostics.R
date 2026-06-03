# 05 — diagnostics: ESS, coverage, weight ranks, posterior marginals (spec §6.3).
# Run and INSPECT before trusting any result. Collapse/coverage issues are raised
# as warnings; the primary fix for collapse is raising sigma_discrep_mm in config.

if (!exists("weights", inherits = FALSE)) source("R/04_weights.R")

summ <- weights |>
  dplyr::group_by(gcm) |>
  dplyr::summarise(
    K        = dplyr::n(),
    ESS      = 1 / sum(weight^2),
    ESS_frac = (1 / sum(weight^2)) / dplyr::n(),
    M_min    = min(M_mm),
    M_med    = median(M_mm),
    M_max    = max(M_mm),
    Y_inside = (Y_mm >= min(M_mm) & Y_mm <= max(M_mm)),
    Y_pctile = mean(M_mm < Y_mm) * 100,
    .groups  = "drop"
  )
print(summ)
write.csv(summ, file.path(out_dir, "diagnostics_summary.csv"), row.names = FALSE)

for (i in seq_len(nrow(summ))) {
  s <- summ[i, ]
  if (!s$Y_inside)
    warning(sprintf("[COVERAGE] %s: Y=%.2f mm outside sampled M [%.2f, %.2f] -> extrapolation.",
                    s$gcm, Y_mm, s$M_min, s$M_max))
  if (s$ESS_frac < 0.02)
    warning(sprintf("[COLLAPSE] %s: ESS=%.0f (%.1f%% of K) -> raise sigma_mod_mult.",
                    s$gcm, s$ESS, 100 * s$ESS_frac))
}

# Sorted-weight (rank) plots per stratum — eyeball collapse onto a few samples.
pdf(file.path(out_dir, "weight_ranks.pdf"), width = 8, height = 6)
par(mfrow = c(2, 2))
for (g in gcms) {
  w <- sort(weights$weight[weights$gcm == g], decreasing = TRUE)
  plot(w, type = "h", log = "x", main = g, xlab = "rank", ylab = "weight")
}
dev.off()

# Posterior (weighted KDE) vs prior parameter marginals, per GCM. The prior is the
# true uniform/log-uniform sampled in 04, so it is drawn analytically as a rectangle
# from the configured bounds (on the log10 axis a log-uniform prior is rectangular
# too); both curves integrate to 1.
pdf(file.path(out_dir, "posterior_marginals.pdf"), width = 10, height = 7)
par(mfrow = c(2, 3))
for (g in gcms) {
  wg <- weights[weights$gcm == g, ]
  for (p in param_cols) {
    lg    <- p %in% log_cols
    xplot <- if (lg) log10(wg[[p]]) else wg[[p]]
    bnds  <- as.numeric(cfg$priors[[p]]); if (lg) bnds <- log10(bnds)
    lo <- bnds[1]; hi <- bnds[2]; prior_h <- 1 / (hi - lo)
    pos   <- density(xplot, bw = bw.nrd0(xplot), weights = wg$weight)
    plot(NA, xlim = range(bnds, pos$x), ylim = c(0, max(prior_h, pos$y)),
         main = paste(g, "-", p), ylab = "density",
         xlab = if (lg) paste0("log10(", p, ")") else p)
    lines(c(lo, lo, hi, hi), c(0, prior_h, prior_h, 0), col = "grey50", lwd = 2)
    lines(pos, col = "firebrick", lwd = 2)
  }
}
dev.off()

message("05: diagnostics complete; figures + summary in ", out_dir,
        " (grey = prior, red = posterior).")
