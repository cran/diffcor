#' Fisher's z-test of difference between an observed correlation from a given value
#'
#' @param emp.r empirically observed correlation
#' @param hypo.r expected correlation
#' @param n sample size the observed correlation is based on
#' @return Fisher's z-values, confidence intervals of the empirical correlations,
#' corresponding p-values, and Cohen's q
#' @export

diffcor.one <- function(emp.r, hypo.r, n, alpha = .05, cor.names = NULL,
                        alternative = c("one.sided", "two.sided"), digit = 3){
  
  fisher.emp <- atanh(emp.r)
  exp.r <- atanh(hypo.r) + (hypo.r / ((2 * n) - 1))
  Var.z <- (n - 3)
  
  diff.z.one <- round((fisher.emp - exp.r) * sqrt(Var.z), digit);
  Cohen.q.one <- round((fisher.emp - exp.r), digit)
  
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  
  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.one)), digit), scientific = F)}
  
  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.one)), digit), scientific = F)}
  
  LL.z <- fisher.emp - abs((qnorm(alpha / 2)) / sqrt(Var.z))
  LL <- round(tanh(LL.z), digit)
  UL.z <- fisher.emp + abs((qnorm(alpha / 2)) / sqrt(Var.z))
  UL <- round(tanh(UL.z), digit);
  
  res.one <- data.frame(round(hypo.r, digit), round(emp.r, digit), round(LL, digit),
                        round(UL, digit), diff.z.one, p, Cohen.q.one)
  rownames(res.one) <- cor.names
  colnames(res.one) <- c("r_exp", "r_obs", "LL", "UL", "z", "p", "Cohen_q")
  
  return(res.one)}