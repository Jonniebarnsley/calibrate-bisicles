library(RobustGaSP)
library(stringr)

find_rectangle_dimensions <- function(N) {
  a <- floor(sqrt(N))
  b <- ceiling(N/a)
  return ( c(a, b) )
}

main_effects <- function(model, inputs) {
  
  N <- length(inputs)
  num_testpoints <- 100
  
  param_testpoints <- as.matrix(seq(0, 1, length.out=num_testpoints))
  df <- data.frame(x=param_testpoints)
  
  for (i in 1:N) {
    
    testing_inputs <- matrix(0.5, ncol=N, nrow=num_testpoints)
    testing_inputs[, i] <- param_testpoints

    testing_trend <- as.matrix(cbind(1, testing_inputs))
    
    model.predict <- predict(
      model,
      testing_inputs,
      testing_trend=testing_trend
    )
    
    name <- colnames(inputs)[i]
    data <- data.frame(
      LB=model.predict$lower95,
      mean=model.predict$mean,
      UB=model.predict$upper95
    )
    
    df[[name]] <- data
  }
  class(df) <- "main_effects"
  return(df)
}

plot.main_effects <- function(main_effects) {
  
  colors <- c('orangered', 'seagreen3', 'royalblue1', 'plum2',
  'slateblue1', 'sienna1', 'goldenrod1', 'darkseagreen', 'orchid',
  'orangered', 'seagreen3', 'royalblue1', 'plum2')
  
  # extract the names and number of inputs
  df <- data.frame(unclass(main_effects))
  ncols <- length(colnames(df))
  N <- (ncols - 1) / 3
  
  indices <- seq(2, 1+N*3, by=3)
  cols <- colnames(df)[indices]
  names <- sub("\\.LB$", "", cols)
  
  yvals <- dplyr::select(df, -x)
  ymin <- min(yvals)
  ymax <- max(yvals)
  
  x <- main_effects$x
  par(mfrow=find_rectangle_dimensions(N))
  for (i in 1:N) {
    var <- names[i]
    LB <- main_effects[[var]]$LB
    UB <- main_effects[[var]]$UB
    mean <- main_effects[[var]]$mean
    
    #ymax <- max(UB)
    #ymin <- min(LB)
    
    # empty plot
    plot(
      1,
      1,
      type='l',
      xlim=c(0, 1),
      ylim=c(ymin, ymax), #ylim=c(-6e11, 4e11), #ylim=c(-15, 20),
      xlab=var,
      ylab='sea level contribution (m)'
    )
    
    # shade 5-95% uncertainty region
    polygon(
      c(x, rev(x)),
      c(LB, rev(UB)),
      col=adjustcolor(colors[i], alpha.f=0.3),
      border=F
    )
    
    # plot emulator mean
    lines(x, mean, type='l', col=colors[i])
  }
}
