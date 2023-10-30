diffpwr.dep <- function(n,
                        r12,
                        r13,
                        r23,
                        alpha = 0.05,
                        n.samples = 1000,
                        seed = 1234){

  set.seed(seed)
  df <- data.frame(matrix(0, nrow = n.samples, ncol = 4))
  colnames(df) <- c("i", "r12", "r13", "r23")
  df$i <- 1:n.samples
  for (i in 1:n.samples){
    frame <- data.frame(mvrnorm(n = n, mu = c(0, 0, 0),
                                Sigma = matrix(c(1, r12, r13,
                                                 r12, 1, r23,
                                                 r13, r23, 1),
                                               3, 3)))
    df$r12[i] <- cor(frame$X1, frame$X2)
    df$r13[i] <- cor(frame$X1, frame$X3)
    df$r23[i] <- cor(frame$X2, frame$X3)

    df$z12[i] <- atanh(df$r12[i])
    df$z13[i] <- atanh(df$r13[i])
    df$r1[i] <- (df$r12[i] + df$r13[i]) / 2
    df$Cov.dep.a[i] <- 1 / ((1 - (df$r1[i]^2))^2)
    df$Cov.dep.b[i] <- df$r23[i] * (1 - (2 * (df$r1[i]^2)))
    df$Cov.dep.c[i] <- 0.5 * (df$r1[i]^2)
    df$Cov.dep.d[i] <- 1 - (2 * (df$r1[i]^2)) - (df$r23[i]^2)
    df$Cov.dep[i] <- (df$Cov.dep.a[i] * df$Cov.dep.b[i]) -
      (df$Cov.dep.c[i] * df$Cov.dep.d[i])
    df$SE.dep[i] <- sqrt((2 - (2 * df$Cov.dep[i]))/(n - 3))
    df$diff.z.dep[i] <- (df$z12[i] - df$z13[i])/df$SE.dep[i]
  }

  pwr <- round(mean(ifelse(abs(df$diff.z.dep) > qnorm(1 - (alpha / 2)),
                           1, 0)),
               3)

  df$LL12 <- atanh(df$r12) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))
  df$UL12 <- atanh(df$r12) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))

  df$LL13 <- atanh(df$r13) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))
  df$UL13 <- atanh(df$r13) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))

  df$LL23 <- atanh(df$r23) - (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))
  df$UL23 <- atanh(df$r23) + (qnorm(1 - (alpha / 2)) * (1 / sqrt(n - 3)))

  cov12 <- round(mean(ifelse(atanh(r12) > df$LL12 & atanh(r12) < df$UL12,
                             1, 0)),
                 3)

  cov13 <- round(mean(ifelse(atanh(r13) > df$LL13 & atanh(r13) < df$UL13,
                             1, 0)),
                 3)

  cov23 <- round(mean(ifelse(atanh(r23) > df$LL23 & atanh(r23) < df$UL23,
                             1, 0)),
                 3)

  bias12 <- round((mean(atanh(df$r12)) - atanh(r12)) / atanh(r12), 3)
  bias13 <- round((mean(atanh(df$r13)) - atanh(r13)) / atanh(r13), 3)
  bias23 <- round((mean(atanh(df$r23)) - atanh(r23)) / atanh(r23), 3)

  res <- data.frame(r12, cov12, bias12,
                    r13, cov13, bias13,
                    r23, cov23, bias23,
                    n, pwr)

  return(res)
}
