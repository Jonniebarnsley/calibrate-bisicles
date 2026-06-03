library(RobustGaSP)

leave_one_out <- function(model) {
  # Takes an rgasp model as an input and performs leave one out analysis.
  # Output is an object of the custom class 'leave_one_out', which can be
  # plotted and analysed.

  # extract information from the rgasp-class
  output <- slot(model, "output")

  # use RobustGaSP leave one out function
  loo <- leave_one_out_rgasp(model)

  # convert into dataframe and add some columns
  df <- as.data.frame(loo)
  df["actual"] <- as.numeric(output)
  df["UB"] <- df$mean + 2 * df$sd
  df["LB"] <- df$mean - 2 * df$sd
  df["success"] <- (df$actual > df$LB & df$actual < df$UB)

  class(df) <- "leave_one_out"
  return(df)
}

length(inputs)

plot.leave_one_out <- function(loo) {

  # Analyses the leave one out test by plotting the actual BISICLES (or whatever
  # numerical model) output against the emulator's predicted value, including
  # error bars. Highlights any cases where the actual value is outside of the
  # emulator's 95% confidence interval.

  # leave_one_out class can't be subsetted, so define an equivalent dataframe
  df <- data.frame(unclass(loo))

  # separate into points the emulator successfully predicted and points it didn't
  good <- subset(df, df$success == T)
  bad <- subset(df, df$success == F)

  # find bounds for axes
  ymax <- max(df$UB)
  ymin <- min(df$LB)

  #plot points and error bars for successes (in blue) and failures (in red)
  plot(
    good$actual,
    good$mean,
    xlim = c(ymin, ymax),
    ylim = c(ymin, ymax),
    xlab = "Actual sea level contribution (m)",
    ylab = "Predicted sea level contribution (m)",
    col = "cornflowerblue",
    pch = 1,
    cex = 0.5
  )
  points(
    bad$actual,
    bad$mean,
    col = "indianred",
    pch = 1,
    cex = 0.5
  )
  segments(
    good$actual,
    good$UB,
    good$actual,
    good$LB,
    col = rgb(0, 0, 1, alpha = 0.2),
    lwd = 2
  )
  segments(
    bad$actual,
    bad$UB,
    bad$actual,
    bad$LB,
    col = rgb(1, 0, 0, alpha = 0.2),
    lwd = 2
  )
  # line of equality
  lines(c(ymin - 2, ymax + 2), c(ymin - 2, ymax + 2))
}

summary.leave_one_out <- function(loo) {

  pass <- sum(loo$success)
  fail <- length(loo$success) - pass
  pass_rate <- pass / length(loo$success)

  fem <- loo$mean
  f   <- loo$actual
  std <- loo$sd
  n   <- length(loo$mean)

  ned   <- sqrt(sum((fem - f)^2 / std^2))
  rmse <- sqrt(sum((fem - f)^2 / n))

  cat("Pass:", pass, "\n")
  cat("Fail:", fail, "\n")
  cat("Pass rate:", pass_rate, "\n")
  cat("Normalized Euclidean Distance:", ned, "\n")
  cat("RMSE:", rmse)
}