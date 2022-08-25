#' Test of differences between two independent correlations
#'
#' @param r1 empirically observed correlation in first group
#' @param r2 empirically observed correlation in second group
#' @param n1 sample size the first correlation is based on
#' @param n2 sample size the second correlation is based on
#'
#' @return Fisher's z-value, corresponding p-values, and Cohen's q
#' @export

diffcor.two <- function(r1, r2, n1, n2, alpha = .05, cor.names = NULL,
                        alternative = c("one.sided", "two.sided"), digit = 3){
  
  fisher1 <- atanh(r1)
  fisher2 <- atanh(r2)
  Var1 <- 1 / (n1 - 3)
  Var2 <- 1 / (n2 - 3)
  SE <- sqrt((Var1 + Var2))
  
  diff.z.two <- round(((fisher1 - fisher2) / SE), digit)
  Cohen.q.two <- round((fisher1 - fisher2), digit)
  
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  
  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.two)), digit), scientific = F)}
  
  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.two)), digit), scientific = F)}
  
  LL.z.1 <- fisher1 - abs((qnorm(alpha / 2)) * sqrt(Var1))
  LL1 <- round(tanh(LL.z.1), digit)
  UL.z.1 <- fisher1 + abs((qnorm(alpha / 2)) * sqrt(Var1))
  UL1 <- round(tanh(UL.z.1), digit)
  
  LL.z.2 <- fisher2 - abs((qnorm(alpha / 2)) * sqrt(Var2))
  LL2 <- round(tanh(LL.z.2), digit)
  UL.z.2 <- fisher2 + abs((qnorm(alpha / 2)) * sqrt(Var2))
  UL2 <- round(tanh(UL.z.2), digit)
  
  res.two <- data.frame(round(r1, digit), round(r2, digit), LL1, UL1, LL2, UL2,
                        diff.z.two, p, Cohen.q.two)
  rownames(res.two) <- cor.names
  colnames(res.two) <- c("r1", "r2", "LL1", "UL1", "LL2", "UL2", "z", "p",
                         "Cohen_q")
  return(res.two)
}