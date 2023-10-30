visual_mc <- function(emp.r,
                      n,
                      alpha = .05,
                      n.intervals = 100,
                      seed = 1234){

  set.seed(seed)

  df <- data.frame(matrix(0, nrow = n.intervals, ncol = 4))
  colnames(df) <- c("i", "point", "LL", "UL")

  df$i <- 1:n.intervals
  for (j in 1:n.intervals) {
    frame <- data.frame(mvrnorm(n = n, mu = c(0, 0), Sigma = matrix(c(1, emp.r,
                                                                      emp.r, 1),
                                                                    2, 2)))

    df$point[j] <- cor(frame$X1, frame$X2)
  }

  df$LL <- atanh(df$point) - (qnorm(1 - (alpha/2)) * (1/sqrt(n - 3)))
  df$UL <- atanh(df$point) + (qnorm(1 - (alpha/2)) * (1/sqrt(n - 3)))

  plot(NA, ylim = c(1, n.intervals), xlim = c(min(c(min(df$LL - .05),
                                                    max(df$UL + .05))),
                                              max(c(min(df$LL - .05),
                                                    max(df$UL + .05)))),
       ylab = "", xlab = "r", yaxt = "none", frame.plot = FALSE)

  abline(v = emp.r, lty = "dashed", col = "red")

  for(i in 1:n.intervals){
    points(x = df$point[i], y = i, pch = 15)
    segments(x0 = df$LL[i],
             x1 = df$UL[i],
             y0 = i)
  }
}
