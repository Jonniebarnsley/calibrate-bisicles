# 06 — projection to 2300 per (GCM x scenario) with weighted bands (spec §7).
# Reuses the per-GCM parameter posterior from §6 (weights depend on GCM only) and
# applies the SCENARIO-CORRECT forcing lookup, so scenarios diverge in the future
# even though they were pooled for weighting. Reports prior (unweighted) and
# posterior (weighted) 5/17/50/83/95% bands at each grid year.

if (!exists("weights", inherits = FALSE)) source("R/05_diagnostics.R")

proj  <- cfg$projection
years <- seq(proj$year_from, proj$year_to, by = proj$year_by)
probs <- c(0.05, 0.17, 0.50, 0.83, 0.95)

# Per-GCM posterior subset (raw param values + their weights). build_norm_inputs
# accepts raw params directly, so no preprocessing needed.
post_by_gcm <- setNames(lapply(gcms,
                               function(g) weights[weights$gcm == g, ]),
                        gcms)

row_qs <- function(gcm, scenario, year, dist, M, w) {
  q <- weighted_quantile(M, w, probs)
  data.frame(gcm = gcm, scenario = scenario, year = year, dist = dist,
             q05 = q[1], q17 = q[2], q50 = q[3], q83 = q[4], q95 = q[5],
             mean = weighted.mean(M, w))
}

# Predict each (GCM, scenario) once over ALL years -> K x n_years TOTAL SLC matrix.
# Total SLC = sum of per-basin emulator means (linearity of expectation). For
# single-basin this reduces to the primary basin's predictions.
Mgs <- list()
for (g in gcms) for (s in scenarios) {
  X <- build_norm_inputs(post_by_gcm[[g]][, param_cols], g, s)
  Mgs[[paste(g, s)]] <- Reduce("+", lapply(predict_slc_mm_list,
                                           function(pred) pred(X, years)$mean))
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
        row_qs(g, s, yr, "posterior", M, pg$weight))
      M_pool <- c(M_pool, M)
      w_pool <- c(w_pool, pg$weight / length(gcms))   # uniform prior over GCMs
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

# ---------------------------------------------------------------------------
# Per-sample 2300 SLC for every (GCM, scenario) cell with its posterior weight
# (long format: model, scenario, weight, slc_mm). Weights are per-GCM, repeated
# unchanged across scenarios. Also a 5-panel weighted-KDE plot of (ssp126,
# ssp534-over) at 2300 — per-GCM + uniform-GCM marginal — using SIR resampling
# so MASS::kde2d (no native weight support) produces a weighted density.
# ---------------------------------------------------------------------------
i2300 <- match(2300L, years)
if (!is.na(i2300)) {
  K_res <- 10000L                                 # resample count per GCM / panel
  set.seed(cfg$weighting$seed)

  # SIR resample per GCM; pivot scenarios to columns so each row is one
  # posterior draw evaluated under all forcings (no weight column needed).
  sle_2300 <- do.call(rbind, lapply(gcms, function(g) {
    pg  <- post_by_gcm[[g]]
    idx <- sample.int(nrow(pg), K_res, replace = TRUE, prob = pg$weight)
    df  <- data.frame(model = rep(g, K_res))
    for (s in scenarios) df[[paste0(s, "_sle")]] <- Mgs[[paste(g, s)]][idx, i2300]
    df
  }))
  write.csv(sle_2300, file.path(out_dir, "sle_2300.csv"), row.names = FALSE)

  if (all(c("ssp126", "ssp534-over") %in% scenarios)) {
    suppressPackageStartupMessages(library(MASS))
    panels <- list()
    for (g in gcms) panels[[g]] <- list(
      x = Mgs[[paste(g, "ssp126")]][, i2300],
      y = Mgs[[paste(g, "ssp534-over")]][, i2300],
      w = post_by_gcm[[g]]$weight)
    panels[["ALL (uniform GCM prior)"]] <- list(
      x = unlist(lapply(panels[gcms], `[[`, "x"), use.names = FALSE),
      y = unlist(lapply(panels[gcms], `[[`, "y"), use.names = FALSE),
      w = unlist(lapply(panels[gcms],
                        function(p) p$w / length(gcms)), use.names = FALSE))

    lim <- range(unlist(lapply(panels, function(p) c(p$x, p$y))))

    pdf(file.path(out_dir, "sle_2300_ssp126_vs_ssp534over.pdf"),
        width = 11, height = 7.5)
    par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.7, 0))
    for (nm in names(panels)) {
      p   <- panels[[nm]]
      idx <- sample.int(length(p$w), K_res, replace = TRUE, prob = p$w)
      k   <- MASS::kde2d(p$x[idx], p$y[idx], n = 150, lims = c(lim, lim))
      image(k, col = hcl.colors(64, "viridis"), xlim = lim, ylim = lim,
            xlab = "SLC@2300 (mm), ssp126",
            ylab = "SLC@2300 (mm), ssp534-over",
            main = nm, cex.main = 0.95)
      contour(k, add = TRUE, drawlabels = FALSE,
              col = adjustcolor("white", 0.6), lwd = 0.5, nlevels = 6)
      abline(0, 1, lty = 2, col = "grey85")        # line of equality
    }
    dev.off()
    message(sprintf("06: wrote sle_2300.csv (%d rows) + sle_2300_ssp126_vs_ssp534over.pdf",
                    nrow(sle_2300)))
  }
}

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
