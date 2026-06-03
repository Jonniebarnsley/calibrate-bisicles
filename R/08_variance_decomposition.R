# 08 — variance decomposition of the PRIOR ensemble spread (uncertainty attribution).
# Independent of the calibration weighting. Two complementary analyses over time:
#   (A) ANOVA on the raw ensemble  -> variance from gcm / scenario / interaction
#   (B) Sobol' via the emulator, PER SCENARIO -> variance from the 6 params + a
#       forcing-severity axis that represents the GCM (within a scenario the axis
#       spans the GCMs' forcing values; the GCM dummy tracks the nearest GCM).
# Run via run_variance.R. Follows Seroussi et al. (2023) / Coulon et al. (2025).

if (!exists("train_emulator")) source("R/03_emulator_io.R")

suppressPackageStartupMessages({ library(car); library(sensitivity); library(randtoolbox) })

vcfg      <- cfg$variance
vyears    <- seq(vcfg$year_from, vcfg$year_to, by = vcfg$year_by)
scenarios <- cfg$data$scenarios

# ---------------------------------------------------------------------------
# (A) ANOVA — Type III SS fractions on the raw ensemble, per year.
# Type III requires sum-to-zero contrasts to give correct marginal SS. The
# (Intercept) row is excluded: variance is about deviations from the mean, so the
# decomposition is over {effects + Residuals}, where Residuals = within-cell
# (parametric) variance.
# ---------------------------------------------------------------------------
anova_fractions <- function(formula, df, factors) {
  contr <- setNames(rep(list(contr.sum), length(factors)), factors)
  m  <- lm(formula, data = df, contrasts = contr)
  a  <- car::Anova(m, type = 3)
  ss <- a$`Sum Sq`; names(ss) <- rownames(a)
  ss <- ss[!names(ss) %in% "(Intercept)"]
  ss / sum(ss)
}

anova_combined <- do.call(rbind, lapply(vyears, function(yr) {
  d <- data.frame(slr = ensemble[[ycol(yr)]],
                  gcm = factor(ensemble[[gcm_col]]),
                  scenario = factor(ensemble[[scen_col]]))
  d <- d[is.finite(d$slr), ]
  fr <- anova_fractions(slr ~ gcm * scenario, d, c("gcm", "scenario"))
  data.frame(year = yr, term = names(fr), fraction = as.numeric(fr))
}))

anova_per_scenario <- do.call(rbind, lapply(scenarios, function(s) {
  do.call(rbind, lapply(vyears, function(yr) {
    sel <- ensemble[[scen_col]] == s
    d <- data.frame(slr = ensemble[[ycol(yr)]][sel],
                    gcm = factor(ensemble[[gcm_col]][sel]))
    d <- d[is.finite(d$slr), ]
    fr <- anova_fractions(slr ~ gcm, d, "gcm")
    data.frame(scenario = s, year = yr, term = names(fr), fraction = as.numeric(fr))
  }))
}))

write.csv(anova_combined, file.path(out_dir, "variance_anova_combined.csv"), row.names = FALSE)
write.csv(anova_per_scenario, file.path(out_dir, "variance_anova_per_scenario.csv"), row.names = FALSE)
message("08A: ANOVA done.")

# ---------------------------------------------------------------------------
# (B) Sobol' — 6 params + 1 forcing-severity axis, PER SCENARIO, via the emulator.
# Within a scenario the forcing axis spans the GCMs' warming values: as severity
# rises we move from the low- to the high-forcing GCM. thermal_forcing is mapped from
# warming by a within-scenario linear fit, and the GCM dummy tracks the nearest GCM
# (by warming) so each GCM's real behaviour is recovered at its forcing point — i.e.
# GCM is "represented through forcing intensity".
# ---------------------------------------------------------------------------
sev_col   <- vcfg$forcing_severity              # "warming"
other_col <- setdiff(forcing_cols, sev_col)     # "thermal_forcing"
sobol_inputs <- c(param_cols, "forcing")
k <- length(sobol_inputs)                       # 7

forcing_map <- setNames(lapply(scenarios, function(s) {
  fl  <- forcing_lookup[forcing_lookup[[scen_col]] == s, ]
  fl  <- fl[order(fl[[sev_col]]), ]             # GCMs ordered by warming
  w   <- fl[[sev_col]]; fit <- lm(fl[[other_col]] ~ w)
  list(smin = min(w), smax = max(w), a = unname(coef(fit)[1]), b = unname(coef(fit)[2]),
       gcm_sorted = as.character(fl[[gcm_col]]), mids = (head(w, -1) + tail(w, -1)) / 2)
}), scenarios)

# A [0,1]^k Saltelli design block -> normalised emulator inputs for scenario s.
design_to_inputs <- function(D, s) {
  X <- matrix(0, nrow = nrow(D), ncol = length(emu_cols), dimnames = list(NULL, emu_cols))
  for (j in seq_along(param_cols)) {
    p <- param_cols[j]; rng <- as.numeric(cfg$priors[[p]])
    raw <- if (p %in% log_cols)
      10^(log10(rng[1]) + D[, j] * (log10(rng[2]) - log10(rng[1])))
    else rng[1] + D[, j] * (rng[2] - rng[1])
    X[, p] <- norm_apply(raw, norm_specs[[p]])
  }
  fm  <- forcing_map[[s]]
  sev <- fm$smin + D[, k] * (fm$smax - fm$smin)
  X[, sev_col]   <- norm_apply(sev, norm_specs[[sev_col]])
  X[, other_col] <- norm_apply(fm$a + fm$b * sev, norm_specs[[other_col]])
  nearest <- fm$gcm_sorted[findInterval(sev, fm$mids) + 1L]   # GCM dummy tracks forcing
  for (gg in dummy_gcms) X[, gg] <- as.numeric(nearest == gg)
  X
}

# Saltelli design ([0,1]^k) at base sample size N. Identical across scenarios/years
# (the scenario mapping happens in design_to_inputs), so it is built once below.
sobol_design <- function(N) {
  set.seed(1)
  S2 <- randtoolbox::sobol(n = N, dim = 2 * k)   # unscrambled; A/B = first/last k cols
  sensitivity::sobolSalt(model = NULL, X1 = as.data.frame(S2[, 1:k]),
                         X2 = as.data.frame(S2[, (k + 1):(2 * k)]), scheme = "A")
}
# Single (scenario, year) cell — used only by the convergence check.
sobol_one <- function(s, yr, N) {
  so <- sobol_design(N)
  so <- sensitivity::tell(so, predict_slc_mm(design_to_inputs(as.matrix(so$X), s), yr)$mean / 1000)
  data.frame(input = sobol_inputs, S = so$S$original, T = so$T$original)
}

# Main sweep: per scenario, predict ALL years in one pass (a single coefficient
# prediction in SVD mode; the cached annual GPs otherwise), then compute indices per
# year from the design's responses. tell() recomputes from the fixed design each call.
Nsob <- vcfg$sobol_n
so0  <- sobol_design(Nsob)
sobol_by_scenario <- do.call(rbind, lapply(scenarios, function(s) {
  Yall <- predict_slc_mm(design_to_inputs(as.matrix(so0$X), s), vyears)$mean / 1000
  do.call(rbind, lapply(seq_along(vyears), function(yi) {
    so <- sensitivity::tell(so0, Yall[, yi])
    data.frame(scenario = s, year = vyears[yi], input = sobol_inputs,
               S = so$S$original, T = so$T$original)
  }))
}))
write.csv(sobol_by_scenario, file.path(out_dir, "variance_sobol_by_scenario.csv"), row.names = FALSE)
message("08B: Sobol' done.")

# ---------------------------------------------------------------------------
# Convergence check at the final year + highest-forcing scenario (first-order S).
# ---------------------------------------------------------------------------
sconv <- scenarios[length(scenarios)]
cat("\nSobol' convergence check (", sconv, "@", max(vyears), "; first-order S):\n", sep = "")
conv <- lapply(vcfg$sobol_convergence, function(N)
  setNames(round(sobol_one(sconv, max(vyears), N)$S, 3), sobol_inputs))
conv <- do.call(rbind, conv); rownames(conv) <- paste0("N=", vcfg$sobol_convergence)
print(conv)

# ---------------------------------------------------------------------------
# Stacked-area plots: ANOVA (unchanged) + one Sobol' panel per scenario.
# ---------------------------------------------------------------------------
to_wide <- function(df, value, col) {
  m <- tapply(df[[value]], list(df$year, df[[col]]), function(x) x[1])
  m[is.na(m)] <- 0; m
}
stacked_area <- function(M, cols, main, show_legend = TRUE) {
  yrs <- as.numeric(rownames(M)); cum <- t(apply(M, 1, cumsum))
  plot(NA, xlim = range(yrs), ylim = c(0, 1), xlab = "year",
       ylab = "variance fraction", main = main)
  prev <- rep(0, length(yrs))
  for (j in seq_len(ncol(M))) {
    polygon(c(yrs, rev(yrs)), c(prev, rev(cum[, j])), col = cols[j], border = NA)
    prev <- cum[, j]
  }
  if (show_legend)
    legend("left", legend = colnames(M), fill = cols, bty = "n", cex = 0.7, inset = 0.01)
}
# Centered moving average with shrinking edges; k = 1 is a no-op.
smooth_ma <- function(x, k) {
  if (is.null(k) || k <= 1) return(x)
  n <- length(x); h <- (k - 1) %/% 2
  vapply(seq_len(n), function(i) mean(x[max(1, i - h):min(n, i + h)]), numeric(1))
}
sobol_wide <- function(s) {
  M <- to_wide(sobol_by_scenario[sobol_by_scenario$scenario == s, ], "S", "input")
  M[M < 0] <- 0; M <- M[, c(param_cols, "forcing")]
  kw <- vcfg$sobol_smooth_window           # smooth each input's S(t) before stacking
  for (j in seq_len(ncol(M))) M[, j] <- smooth_ma(M[, j], kw)
  M <- cbind(M, interactions = pmax(0, 1 - rowSums(M)))
  M / rowSums(M)
}

acol  <- c("#d73027", "#4575b4", "#984ea3", "#cccccc")
scol  <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#d73027", "#cccccc")
sccol <- setNames(c("#1a9850", "#fdae61", "#d73027"), scenarios)

pdf(file.path(out_dir, "variance_decomposition.pdf"), width = 13, height = 9)
# Page 1 — ANOVA (unchanged): combined stack + per-scenario GCM share.
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
Ma <- to_wide(anova_combined, "fraction", "term")[, c("scenario", "gcm", "gcm:scenario", "Residuals")]
stacked_area(Ma, acol, "ANOVA: discrete factors (all scenarios)")
plot(NA, xlim = range(anova_per_scenario$year), ylim = c(0, 1),
     xlab = "year", ylab = "GCM variance fraction", main = "Per-scenario ANOVA: GCM share")
for (s in scenarios) {
  d <- anova_per_scenario[anova_per_scenario$scenario == s & anova_per_scenario$term == "gcm", ]
  d <- d[order(d$year), ]; lines(d$year, d$fraction, col = sccol[s], lwd = 2)
}
legend("topright", legend = scenarios, col = sccol, lwd = 2, bty = "n")
# Page 2 — Sobol' per scenario (GCM represented through forcing intensity).
par(mfrow = c(1, length(scenarios)), mar = c(4, 4, 2.5, 1))
for (i in seq_along(scenarios))
  stacked_area(sobol_wide(scenarios[i]), scol, paste0("Sobol': ", scenarios[i]),
               show_legend = (i == 1))
dev.off()
message("08C: plots written to ", file.path(out_dir, "variance_decomposition.pdf"))
