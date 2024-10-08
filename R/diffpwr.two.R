diffpwr.two <- function(n1,
                        n2,
                        rho1,
                        rho2,
                        alpha = .05,
                        n.samples = 1000,
                        seed = 1234){

  set.seed(seed)

  df1 <- data.frame(matrix(0, nrow = n.samples, ncol = 2))
  colnames(df1) <- c("i", "point")
  df1$i <- 1:n.samples

  for (i in 1:n.samples) {
    frame1 <- data.frame(mvrnorm(n = n1,
                                 mu = c(0, 0),
                                 Sigma = matrix(c(1, rho1,
                                                  rho1, 1),
                                                2, 2)))
    df1$point[i] <- cor(frame1$X1, frame1$X2)
  }

  set.seed(seed^2)

  df2 <- data.frame(matrix(0, nrow = n.samples, ncol = 2))
  colnames(df2) <- c("i", "point")
  df2$i <- 1:n.samples

  for (i in 1:n.samples) {
    frame2 <- data.frame(mvrnorm(n = n2,
                                 mu = c(0, 0),
                                 Sigma = matrix(c(1, rho2,
                                                  rho2, 1),
                                                2, 2)))
    df2$point[i] <- cor(frame2$X1, frame2$X2)
  }

  df_t <- data.frame(df1$point,
                     df2$point,
                     rep(1 / sqrt(n1 - 3), n.samples),
                     rep(1 / sqrt(n2 - 3), n.samples))

  colnames(df_t) <- c("point1", "point2", "SE1", "SE2")

  df_t$LL1 <- tanh(atanh(df_t$point1) - (qnorm(1 - (alpha / 2)) * df_t$SE1))
  df_t$UL1 <- tanh(atanh(df_t$point1) + (qnorm(1 - (alpha / 2)) * df_t$SE1))

  df_t$LL2 <- tanh(atanh(df_t$point2) - (qnorm(1 - (alpha / 2)) * df_t$SE2))
  df_t$UL2 <- tanh(atanh(df_t$point2) + (qnorm(1 - (alpha / 2)) * df_t$SE2))


  pwr <- round(ifelse(rho1 < rho2,
                      mean(df_t$UL1 < df_t$LL2),
                      mean(df_t$UL2 < df_t$LL1)),
               3)

  cov1 <- round(mean(ifelse(tanh(atanh(rho1)) > df_t$LL1
                            & tanh(atanh(rho1)) < df_t$UL1,
                            1, 0)),
                3)

  cov2 <- round(mean(ifelse(tanh(atanh(rho2)) > df_t$LL2
                            & tanh(atanh(rho2)) < df_t$UL2,
                            1, 0)),
                3)

  bias1_M <- round((mean(tanh(atanh(df_t$point1))) - tanh(atanh(rho1)))
                   / tanh(atanh(rho1)), 3)
  bias1_Md <- round((median(tanh(atanh(df_t$point1))) - tanh(atanh(rho1)))
                    / tanh(atanh(rho1)), 3)

  bias2_M <- round((mean(tanh(atanh(df_t$point2))) - tanh(atanh(rho2)))
                   / tanh(atanh(rho2)), 3)
  bias2_Md <- round((median(tanh(atanh(df_t$point2))) - tanh(atanh(rho2)))
                    / tanh(atanh(rho2)), 3)


  res <- data.frame(
    Parameters = c("rho1", "n1", "cov1", "bias1_M", "bias1_Md",
                   "rho2", "n2", "cov2", "bias2_M", "bias2_Md",
                   "pwr"),

    Estimates = c(rho1, n1, cov1, bias1_M, bias1_Md,
                  rho2, n2, cov2, bias2_M, bias2_Md,
                  pwr))

  return(res)
}
