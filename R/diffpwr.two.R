diffpwr.two <- function(n.samples = 1000,
                        n1,
                        n2,
                        r1,
                        r2,
                        alpha = .05,
                        seed = 1234){

  # Monte Carlo simulation

  set.seed(seed)

  df1 <- data.frame(matrix(0, nrow = n.samples, ncol = 2))
  colnames(df1) <- c("i", "point")
  df1$i <- 1:n.samples

  for (j in 1:n.samples) {
    frame1 <- data.frame(mvrnorm(n = n1,
                                 mu = c(0, 0),
                                 Sigma = matrix(c(1, r1,
                                                  r1, 1),
                                                2, 2)))
    df1$point[j] <- cor(frame1$X1, frame1$X2)
  }

  df2 <- data.frame(matrix(0, nrow = n.samples, ncol = 2))
  colnames(df2) <- c("i", "point")
  df2$i <- 1:n.samples

  for (j in 1:n.samples) {
    frame2 <- data.frame(mvrnorm(n = n2,
                                 mu = c(0, 0),
                                 Sigma = matrix(c(1, r2,
                                                  r2, 1),
                                                2, 2)))
    df2$point[j] <- cor(frame2$X1, frame2$X2)
  }

  df_t <- data.frame(df1$point,
                     df2$point,
                     rep(1 / (n1 - 3), n.samples),
                     rep(1 / (n2 - 3), n.samples))

  colnames(df_t) <- c("point1", "point2", "SE1", "SE2")

  z <- (atanh(df_t$point1) - atanh(df_t$point2)) / sqrt(df_t$SE1 + df_t$SE2)

  check <- ifelse(abs(z) > abs(qnorm(1 - (alpha / 2))), 1, 0)

  pwr <- mean(check)

  out <- data.frame(r1, n1, r2, n2, pwr)

  # Visualization of the correlation difference

  LL1 <- tanh(atanh(r1) - (qnorm(1 - (alpha / 2)) / sqrt(n1 - 3)))
  UL1 <- tanh(atanh(r1) + (qnorm(1 - (alpha / 2)) / sqrt(n1 - 3)))

  LL2 <- tanh(atanh(r2) - (qnorm(1 - (alpha / 2)) / sqrt(n2 - 3)))
  UL2 <- tanh(atanh(r2) + (qnorm(1 - (alpha / 2)) / sqrt(n2 - 3)))


  plot(NA, ylim = c(.5, 2),
       xlim = c(min(LL1, LL2), max(UL1, UL2)),
       ylab = "", xlab = "Correlations",
       yaxt = "none")

  points(x = r1, y = 1, pch = 15)
  segments(x0 = LL1,
           x1 = UL1,
           y0 = 1)

  points(x = r2, y = 1.5, pch = 15)
  segments(x0 = LL2,
           x1 = UL2,
           y0 = 1.5)

  # Return output

  return(out)
}
