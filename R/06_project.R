# 06 — projection to 2300 per (GCM x scenario) with weighted bands (spec §7).
# Reuses the per-GCM parameter posterior from §6 (weights depend on GCM only) and
# applies the SCENARIO-CORRECT forcing lookup, so scenarios diverge in the future
# even though they were pooled for weighting. Reports prior (unweighted) and
# posterior (weighted) 5/17/50/83/95% bands at each grid year.

if (!exists("weights", inherits = FALSE)) source("R/05_diagnostics.R")

proj      <- cfg$projection
years     <- seq(proj$year_from, proj$year_to, by = proj$year_by)
probs     <- c(0.05, 0.17, 0.50, 0.83, 0.95)
scenarios <- cfg$data$scenarios

# Normalise the saved posterior parameters once per GCM (reused across years).
# norm_params() is defined in 03_emulator_io.R.
post_by_gcm <- setNames(lapply(gcms, function(g) {
  wg <- weights[weights$gcm == g, ]
  list(wg = wg, pn = norm_params(wg))
}), gcms)

row_qs <- function(gcm, scenario, year, dist, M, w) {
  q <- weighted_quantile(M, w, probs)
  data.frame(gcm = gcm, scenario = scenario, year = year, dist = dist,
             q05 = q[1], q17 = q[2], q50 = q[3], q83 = q[4], q95 = q[5],
             mean = weighted.mean(M, w))
}

# Predict each (GCM, scenario) once over ALL years -> K x n_years matrix (in SVD mode
# this reuses one coefficient prediction across years; in per-year mode it loops the
# cached annual GPs). Summaries then index the right column per year.
Mgs <- list()
for (g in gcms) for (s in scenarios) {
  X <- build_norm_inputs(post_by_gcm[[g]]$pn, g, s)
  Mgs[[paste(g, s)]] <- predict_slc_mm(X, years)$mean
}

records <- list()
for (yi in seq_along(years)) {
  yr <- years[yi]
  for (s in scenarios) {
    M_pool <- numeric(0); w_pool <- numeric(0)
    for (g in gcms) {
      pg <- post_by_gcm[[g]]
      M  <- Mgs[[paste(g, s)]][, yi]                                     # mm
      records[[length(records) + 1L]] <- rbind(
        row_qs(g, s, yr, "prior",     M, rep(1, length(M))),
        row_qs(g, s, yr, "posterior", M, pg$wg$weight))
      M_pool <- c(M_pool, M)
      w_pool <- c(w_pool, pg$wg$weight / length(gcms))   # uniform prior over GCMs
    }
    if (isTRUE(proj$gcm_marginal)) {
      gm <- "ALL (uniform GCM prior)"
      records[[length(records) + 1L]] <- rbind(
        row_qs(gm, s, yr, "prior",     M_pool, rep(1, length(M_pool))),
        row_qs(gm, s, yr, "posterior", M_pool, w_pool))
    }
  }
}
projection <- do.call(rbind, records)
write.csv(projection, file.path(out_dir, "projection.csv"), row.names = FALSE)

# Prior-vs-posterior band plots: rows = GCM (+ marginal), cols = scenario.
plot_gcms <- c(gcms, if (isTRUE(proj$gcm_marginal)) "ALL (uniform GCM prior)")
pdf(file.path(out_dir, "projection_bands.pdf"),
    width = 3.4 * length(scenarios), height = 2.6 * length(plot_gcms))
par(mfrow = c(length(plot_gcms), length(scenarios)),
    mar = c(3, 3.2, 2, 1), mgp = c(1.9, 0.6, 0))
for (g in plot_gcms) for (s in scenarios) {
  d   <- projection[projection$gcm == g & projection$scenario == s, ]
  pri <- d[d$dist == "prior", ]
  pos <- d[d$dist == "posterior", ]
  plot(NA, xlim = range(years), ylim = range(d$q05, d$q95),
       xlab = "year", ylab = "SLC (mm)", main = paste(g, s), cex.main = 0.9)
  polygon(c(pri$year, rev(pri$year)), c(pri$q05, rev(pri$q95)),
          col = adjustcolor("grey55", 0.40), border = NA)
  polygon(c(pos$year, rev(pos$year)), c(pos$q05, rev(pos$q95)),
          col = adjustcolor("firebrick", 0.35), border = NA)
  lines(pri$year, pri$q50, col = "grey25", lwd = 2)
  lines(pos$year, pos$q50, col = "firebrick", lwd = 2)
}
dev.off()

message(sprintf("06: projected %d years x %d GCMs x %d scenarios; wrote projection.csv + bands.",
                length(years), length(gcms), length(scenarios)))
