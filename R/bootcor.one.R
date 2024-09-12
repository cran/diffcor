bootcor.one <- function(x, 
                        y,
                        r_target,
                        k = 5000,
                        alpha = .05,
                        digit = 3,
                        seed = 1234){
  set.seed(seed)
  
  r_emp <- cor(x, y)
  
  diff <- NULL
  
  for(i in 1:k){
    order <- sample(1:length(x), 
                    size = length(x),
                    replace = TRUE)
    
    z_emp <- atanh(cor(x[order], y[order]))
    z_exp <- atanh(r_target)
    
    diff[i] <- tanh(z_emp - z_exp)
  }
  
  out <- data.frame(round(r_emp, digit),
                    r_target,
                    round(mean(diff), digit), 
                    round(quantile(diff, alpha / 2), digit), 
                    round(quantile(diff, 1 - (alpha / 2)), digit))
  
  colnames(out) <- c("r_emp",
                     "r_target",
                     "M", 
                     "LL", 
                     "UL")
  rownames(out) <- ""
  
  return(out)
}
