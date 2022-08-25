#' Test of differences between two dependent correlations
#'
#' @param r12 empirically observed correlation between first and second construct
#' @param r13 empirically observed correlation between first and third construct
#' @param r23 empirically observed correlation between second and third construct
#' @param n sample size the correlations are based on
#'
#' @return Fisher's z-value and corresponding p-values
#' @export

diffcor.dep <- function(r12, r13, r23, n, cor.names = NULL,
                        alternative = c("one.sided", "two.sided"), digit = 3){
  
  z12 <- atanh(r12)
  z13 <- atanh(r13)
  r1 <- (r12 + r13) / 2
  
  Cov.dep.a <- 1 / (1 - r1 ^ 2) ^ 2
  Cov.dep.b <- r23 * (1 - 2 * r1 ^ 2)
  Cov.dep.c <- -.50 * (r1 ^ 2)
  Cov.dep.d <- 1 - (2 * (r1 ^ 2)) - (r23 ^ 2)
  Cov.dep <- (Cov.dep.a * Cov.dep.b) + (Cov.dep.c * Cov.dep.d)
  SE.dep <- sqrt((2 - (2 * Cov.dep)) / (n - 3))
  diff.z.dep <- round(((z12 - z13) / SE.dep), digit)
  
  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  
  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.dep)), digit), scientific = F)}
  
  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.dep)), digit), scientific = F)}
  
  res.dep <- data.frame(round(r12, digit), round(r13, digit), round(r23, digit),
                        diff.z.dep, p)
  rownames(res.dep) <- cor.names
  colnames(res.dep) <- c("r12", "r13", "r23", "z", "p")
  
  return(res.dep)
}