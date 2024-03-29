diffpwr.one <- function(n,
                        emp.r,
                        hypo.r,
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
                                Sigma = matrix(c(1, emp.r,
                                                 emp.r, 1),
                                               2, 2)))
    df$point[j] <- cor(frame$X1, frame$X2)
  }

  df$LL <- tanh(atanh(df$point) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))
  df$UL <- tanh(atanh(df$point) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))

  pwr <- round(mean(ifelse(tanh(atanh(hypo.r)) > df$LL &
                             tanh(atanh(hypo.r)) < df$UL,
                           0, 1)),
               3)

  cov <- round(mean(ifelse(tanh(atanh(emp.r)) > df$LL
                           & tanh(atanh(emp.r)) < df$UL,
                           1, 0)),
               3)

  bias <- round((mean(tanh(atanh(df$point))) - tanh(atanh(emp.r))) /
                  tanh(atanh(emp.r)), 3)

  res <- data.frame(emp.r, hypo.r, n, cov, bias, pwr)

  return(res)
}
