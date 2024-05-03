diffpwr.one <- function(n,
                        r,
                        rho,
                        alpha = .05,
                        n.samples = 1000,
                        seed = 1234){

  set.seed(seed)

  df <- data.frame(matrix(0, nrow = n.samples, ncol = 4))
  colnames(df) <- c("i", "point", "LL", "UL")
  df$i <- 1:n.samples

  for (j in 1:n.samples) {
    frame <- data.frame(mvrnorm(n = n,
                                mu = c(0, 0),
                                Sigma = matrix(c(1, r,
                                                 r, 1),
                                               2, 2)))
    df$point[j] <- cor(frame$X1, frame$X2)
  }

  df$LL <- tanh(atanh(df$point) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))
  df$UL <- tanh(atanh(df$point) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))

  pwr <- round(mean(ifelse(tanh(atanh(rho)) > df$LL &
                             tanh(atanh(rho)) < df$UL,
                           0, 1)),
               3)

  cov <- round(mean(ifelse(tanh(atanh(r)) > df$LL
                           & tanh(atanh(r)) < df$UL,
                           1, 0)),
               3)

  bias_M <- round((mean(tanh(atanh(df$point))) - tanh(atanh(r))) /
                    tanh(atanh(r)), 3)

  bias_Md <- round((median(tanh(atanh(df$point))) - tanh(atanh(r))) /
                     tanh(atanh(r)), 3)

  res <- data.frame(
    Parameters = c("r", "rho", "n", "cov", "bias_M", "bias_Md", "pwr"),

    Estimates = c(r, rho, n, cov, bias_M, bias_Md, pwr))

  return(res)
}
