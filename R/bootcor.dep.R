bootcor.dep <- function(target,
                        x1, 
                        x2, 
                        k = 5000,
                        alpha = .05,
                        digit = 3,
                        seed = 1234){
  
  set.seed(seed)
  
  r_target_1 = cor(target, x1)
  r_target_2 = cor(target, x2)
  r_1_2 = cor(x1, x2)
  
  diff <- NULL
  
  for(i in 1:k){
    order <- sample(1:length(x1), size = length(x1), replace = TRUE)
    
    z_target_1 <- atanh(cor(target[order], x1[order]))
    z_target_2 <- atanh(cor(target[order], x2[order]))
    
    diff[i] <- tanh(z_target_1 - z_target_2)
    
  }
  
  out <- data.frame(round(r_target_1, digit),
                    round(r_target_2, digit),
                    round(r_1_2, digit),
                    round(mean(diff), digit), 
                    round(quantile(diff, alpha / 2), digit), 
                    round(quantile(diff, 1 - (alpha / 2)), digit))
  
  colnames(out) <- c("r_target_1",
                     "r_target_2",
                     "r_1_2",
                     "M", 
                     "LL", 
                     "UL")
  rownames(out) <- ""
  
  return(out)
}
