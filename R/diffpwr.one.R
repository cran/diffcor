diffpwr.one <- function(n.samples = 1000,
                        n,
                        emp.r,
                        hypo.r,
                        alpha = .05,
                        seed = 1234){
  # Monte Carlo simulation
  set.seed(seed)

  df <- data.frame(matrix(0, nrow = n.samples, ncol = 4))
  colnames(df) <- c("i", "point", "LL", "UL")
  df$i <- 1:n.samples

  for (j in 1:n.samples) {
    frame <- data.frame(mvrnorm(n = n,
                                mu = c(0, 0),
                                Sigma = matrix(c(1, emp.r,
                                                 emp.r, 1),
                                               2, 2)))
    df$point[j] <- cor(frame$X1, frame$X2)
  }

  df$LL <- atanh(df$point) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))
  df$UL <- atanh(df$point) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))

  df$check <- ifelse(atanh(hypo.r) > df$LL & atanh(hypo.r) < df$UL, 0, 1)

  pwr <- mean(df$check)

  res <- data.frame(n, pwr)

  # Visualization of the correlation difference

  LL <- tanh(atanh(emp.r) - (qnorm(1 - (alpha / 2)) / sqrt(n - 3)))
  UL <- tanh(atanh(emp.r) + (qnorm(1 - (alpha / 2)) / sqrt(n - 3)))

  plot(NA, ylim = c(0, 1),
       xlim = c(min(c(LL, hypo.r, UL)), max(c(LL, hypo.r, UL))),
       ylab = "", xlab = "Correlation",
       yaxt = "none")

  points(x = emp.r, y = .5, pch = 15)
  abline(v = hypo.r, lty = "dashed", col = "red")
  segments(x0 = LL,
           x1 = UL,
           y0 = .5)

  # Return output

  return(res)
}
