bootcor.two <- function(x1,
                        y1,
                        x2, 
                        y2,
                        k = 5000,
                        alpha = .05,
                        digit = 3,
                        seed = 1234){
  set.seed(seed)
  
  r1 <- cor(x1, y1)
  r2 <- cor(x2, y2)
  
  diff <- NULL
  
  for(i in 1:k){
    order1 <- sample(1:length(x1), 
                     size = length(x1),
                     replace = TRUE)
    order2 <- sample(1:length(x2), 
                     size = length(x2),
                     replace = TRUE)

    z1 <- atanh(cor(x1[order1], y1[order1]))
    z2 <- atanh(cor(x2[order2], y2[order2]))
    
    diff[i] <- tanh(z1 - z2)
  }
  
  out <- data.frame(round(r1, digit),
                    round(r2, digit),
                    round(mean(diff), digit), 
                    round(quantile(diff, alpha / 2), digit), 
                    round(quantile(diff, 1 - (alpha / 2)), digit))
  
  colnames(out) <- c("r1",
                     "r2",
                     "M", 
                     "LL", 
                     "UL")
  rownames(out) <- ""
  
  return(out)
}
