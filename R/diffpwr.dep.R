diffpwr.dep <- function(n,
                        rho12,
                        rho13,
                        rho23,
                        alpha = .05,
                        n.samples = 1000,
                        seed = 1234){

  set.seed(seed)
  df <- data.frame(matrix(0, nrow = n.samples, ncol = 4))
  colnames(df) <- c("i", "rho12", "rho13", "rho23")
  df$i <- 1:n.samples
  for (i in 1:n.samples){
    frame <- data.frame(mvrnorm(n = n, mu = c(0, 0, 0),
                                Sigma = matrix(c(1, rho12, rho13,
                                                 rho12, 1, rho23,
                                                 rho13, rho23, 1),
                                               3, 3)))
    df$rho12[i] <- cor(frame$X1, frame$X2)
    df$rho13[i] <- cor(frame$X1, frame$X3)
    df$rho23[i] <- cor(frame$X2, frame$X3)

    df$z12[i] <- atanh(df$rho12[i])
    df$z13[i] <- atanh(df$rho13[i])
    df$r1[i] <- (df$rho12[i] + df$rho13[i]) / 2
    df$Cov.dep.a[i] <- 1 / ((1 - (df$r1[i]^2))^2)
    df$Cov.dep.b[i] <- df$rho23[i] * (1 - (2 * (df$r1[i]^2)))
    df$Cov.dep.c[i] <- 0.5 * (df$r1[i]^2)
    df$Cov.dep.d[i] <- 1 - (2 * (df$r1[i]^2)) - (df$rho23[i]^2)
    df$Cov.dep[i] <- (df$Cov.dep.a[i] * df$Cov.dep.b[i]) -
      (df$Cov.dep.c[i] * df$Cov.dep.d[i])
    df$SE.dep[i] <- sqrt((2 - (2 * df$Cov.dep[i]))/(n - 3))
    df$diff.z.dep[i] <- (df$z12[i] - df$z13[i])/df$SE.dep[i]
  }

  pwr <- round(mean(ifelse(abs(df$diff.z.dep) > qnorm(1 - (alpha / 2)),
                           1, 0)),
               3)

  df$LL12 <- tanh(atanh(df$rho12) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))
  df$UL12 <- tanh(atanh(df$rho12) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))

  df$LL13 <- tanh(atanh(df$rho13) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))
  df$UL13 <- tanh(atanh(df$rho13) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))

  df$LL23 <- tanh(atanh(df$rho23) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))
  df$UL23 <- tanh(atanh(df$rho23) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3))))

  cov12 <- round(mean(ifelse(tanh(atanh(rho12)) > df$LL12
                             & tanh(atanh(rho12)) < df$UL12,
                             1, 0)),
                 3)

  cov13 <- round(mean(ifelse(tanh(atanh(rho13)) > df$LL13
                             & tanh(atanh(rho13)) < df$UL13,
                             1, 0)),
                 3)

  cov23 <- round(mean(ifelse(tanh(atanh(rho23)) > df$LL23
                             & tanh(atanh(rho23)) < df$UL23,
                             1, 0)),
                 3)

  bias12_M <- round((mean(tanh(atanh(df$rho12))) - tanh(atanh(rho12)))
                  / tanh(atanh(rho12)), 3)
  bias12_Md <- round((median(tanh(atanh(df$rho12))) - tanh(atanh(rho12)))
                     / tanh(atanh(rho12)), 3)

  bias13_M <- round((mean(tanh(atanh(df$rho13))) - tanh(atanh(rho13)))
                  / tanh(atanh(rho13)), 3)
  bias13_Md <- round((median(tanh(atanh(df$rho13))) - tanh(atanh(rho13)))
                     / tanh(atanh(rho13)), 3)

  bias23_M <- round((mean(tanh(atanh(df$rho23))) - tanh(atanh(rho23)))
                    / tanh(atanh(rho23)), 3)
  bias23_Md <- round((median(tanh(atanh(df$rho23))) - tanh(atanh(rho23)))
                     / tanh(atanh(rho23)), 3)

  res <- data.frame(
    Parameters = c("rho12", "cov12", "bias12_M", "bias12_Md",
                   "rho13", "cov13", "bias13_M", "bias13_Md",
                   "rho23", "cov23", "bias23_M", "bias23_Md",
                   "n", "pwr"),

    Estimates = c(rho12, cov12, bias12_M, bias12_Md,
                  rho13, cov13, bias13_M, bias13_Md,
                  rho23, cov23, bias23_M, bias23_Md,
                  n, pwr))

  return(res)
}
