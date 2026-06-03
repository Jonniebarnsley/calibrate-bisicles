library(RobustGaSP)
library(dplyr)
library(lhs)

setwd("/Users/jonniebarnsley/code/phd/local/R")
ensemble <- read.csv("emulator_inputs_2007-2300_sle.csv")
ensemble <- na.omit(ensemble)

# Add columns with binary switches for model and scenario
models <- c("CESM2-WACCM", "MRI-ESM2-0", "EC-Earth3-Veg", "MIROC-ES2L")
scenarios <- c("ssp126", "ssp534-over", "ssp585")
for (s in scenarios) {
  ensemble[[s]] <- as.integer(ensemble$scenario == s)
}
for (m in models) {
  ensemble[[m]] <- as.integer(ensemble$model == m)
}

# Select emulator inputs
inputs <- na.omit(
  dplyr::select(
    ensemble,
    aPhi,
    n,
    m,
    uf,
    gamma,
    UMV,
    # ssp126,
    # `ssp534-over`,
    # ssp585,
    `CESM2-WACCM`,
    `MRI-ESM2-0`,
    # `EC-Earth3-Veg`,
    `MIROC-ES2L`,
    # warming,
    # thermal_forcing
  )
)

labels <- c(
  expression(a[phi]), "n", "m", expression(u[f]), expression(gamma), "UMV"
)
plot_inputs <- dplyr::select(
  inputs,
  aPhi,
  n,
  m,
  uf,
  gamma,
  UMV
)
pairs(
  plot_inputs,
  log = c(1, 5, 6),
  lower.panel = NULL,
  cex = 0.5,
  labels = labels
)

normalise <- function(x, log = FALSE) {
  if (log) {
    stopifnot(all(x > 0))
    x <- log10(x)
  }
  (x - min(x)) / (max(x) - min(x))
}

log_cols <- c("aPhi", "gamma", "UMV")   # columns sampled on a log scale
normalised_inputs <- as.data.frame(
  lapply(names(inputs), function(col) {
    x <- inputs[col]
    normalise(
      x,
      log = col %in% log_cols
    )
  })
)
output <- na.omit(ensemble$X2021)
trend <- as.matrix(cbind(1, normalised_inputs))

model <- rgasp(
  design = normalised_inputs,
  response = output,
  trend = trend,
  kernel_type = "matern_5_2",
  alpha = 1.9,
  nugget.est = TRUE,
)

# check for inert inputs
p <- findInertInputs(model)
source("leaveoneout.R")
loo <- leave_one_out(model)
par(mfrow = c(1, 1))
plot(loo)
summary(loo)

source("main_effects.R")
me <- main_effects(model, normalised_inputs)
plot(me)

ensemble$success <- as.logical(loo$success)

# ---- Diagnose where leave-one-out struggles ----
ensemble <- ensemble |>
  mutate(
    failure = !success,
    success_group = factor(success, levels = c(FALSE, TRUE),
      labels = c("fail", "success")
    )
  )

overall_failure <- ensemble |>
  summarise(
    n = n(),
    failures = sum(failure),
    failure_rate = mean(failure)
  )
print(overall_failure)

failure_by_group <- ensemble |>
  group_by(scenario, model) |>
  summarise(
    n = n(),
    failures = sum(failure),
    failure_rate = mean(failure),
    .groups = "drop"
  ) |>
  arrange(desc(failure_rate))
print(failure_by_group)

param_cols <- c("aPhi", "n", "m", "uf", "gamma", "UMV", "thermal_forcing"
)

std_mean_diff <- function(x, fail) {
  x_fail <- x[fail]
  x_ok <- x[!fail]
  pooled_sd <- sqrt((var(x_fail) + var(x_ok)) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(0)
  }
  (mean(x_fail) - mean(x_ok)) / pooled_sd
}

smd <- sapply(param_cols, function(v) {
  std_mean_diff(ensemble[[v]], ensemble$failure)
})
smd_table <- data.frame(
  variable = names(smd),
  smd = as.numeric(smd),
  abs_smd = abs(as.numeric(smd))
) |>
  arrange(desc(abs_smd))
print(smd_table)

# Fit an interpretable model of failure propensity.
analysis_df <- ensemble |>
  mutate(
    log_aPhi = log10(aPhi),
    log_gamma = log10(gamma),
    log_UMV = log10(UMV)
  ) |>
  dplyr::select(
    failure,
    log_aPhi,
    n,
    m,
    uf,
    log_gamma,
    log_UMV,
    warming,
    thermal_forcing,
    scenario,
    model
  )

failure_glm <- glm(failure ~ ., data = analysis_df, family = binomial())
glm_coef <- summary(failure_glm)$coefficients
glm_table <- data.frame(
  term = rownames(glm_coef),
  estimate = glm_coef[, "Estimate"],
  std_error = glm_coef[, "Std. Error"],
  z_value = glm_coef[, "z value"],
  p_value = glm_coef[, "Pr(>|z|)"],
  odds_ratio = exp(glm_coef[, "Estimate"])
) |>
  arrange(p_value)
print(glm_table)

par(mfrow = c(2, 4))
for (v in param_cols) {
  boxplot(
    ensemble[[v]] ~ ensemble$success_group,
    main = v,
    xlab = "LOO outcome",
    ylab = v,
    col = c("indianred1", "grey80"),
    outline = FALSE
  )
}
