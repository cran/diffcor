#' Fisher's z-tests concerning difference of an observed correlation from given value or difference between independent correlations
#'
#' Test concerning differences between hypothesis and observation
#'
#' @param emp.r empirically observed correlation
#' @param hypo.r expected correlation
#' @param n sample size the observed correlation is based on
#' @return Fisher's z-values, confidence intervals of the empirical correlations,
#' corresponding p-values, and Cohen's q
#' @export

diffcor.one <- function(emp.r, hypo.r, n, alpha = .05, cor.names = NULL,
                        alternative = c("one.sided", "two.sided"), digit = 3){

  fisher.emp <- atanh(emp.r);
  exp.r <- atanh(hypo.r) + (hypo.r / ((2 * n) - 1));
  Var.z <- (n - 3);

  diff.z.one <- round((fisher.emp - exp.r) * sqrt(Var.z), 2);
  Cohen.q.one <- round((fisher.emp - exp.r), 2);

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)

  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.one)), digit), scientific = F)}

  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.one)), digit), scientific = F)};

  LL.z <- fisher.emp - (qnorm(alpha/2) / sqrt(Var.z)); LL <- round(tanh(LL.z), digit)
  UL.z <- fisher.emp + (qnorm(alpha/2) / sqrt(Var.z)); UL <- round(tanh(UL.z), digit);

  res.one <- data.frame(round(hypo.r, digit), round(emp.r, digit), round(LL, digit),
                        round(UL, digit), diff.z.one, p, Cohen.q.one);
  rownames(res.one) <- cor.names;
  colnames(res.one) <- c("r_exp", "r_obs", "LL", "UL", "z", "p", "Cohen_q");

  return(res.one)}

#' Test concerning differences between two independent correlations
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

  fisher1 <- atanh(r1); fisher2 <- atanh(r2);
  Var1 <- 1/(n1 - 3); Var2 <- 1/(n2 - 3); SE <- sqrt((Var1 + Var2));

  diff.z.two <- round(((fisher1 - fisher2) / SE), 2);
  Cohen.q.two <- round((fisher1 - fisher2), 2)

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.two)), digit), scientific = F)};

  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.two)), digit), scientific = F)};

  LL.z.1 <- fisher1 - (qnorm(alpha/2) * Var1); LL1 <- round(tanh(LL.z.1), digit);
  UL.z.1 <- fisher1 + (qnorm(alpha/2) * Var1); UL1 <- round(tanh(UL.z.1), digit);

  LL.z.2 <- fisher2 - (qnorm(alpha/2) * Var2); LL2 <- round(tanh(LL.z.2), digit);
  UL.z.2 <- fisher2 + (qnorm(alpha/2) * Var2); UL2 <- round(tanh(UL.z.2), digit);

  res.two <- data.frame(round(r1, digit), round(r2, digit), LL1, UL1, LL2, UL2,
                        diff.z.two, p, Cohen.q.two)
  rownames(res.two) <- cor.names;
  colnames(res.two) <- c("r1", "r2", "LL1", "UL1", "LL2", "UL2", "z", "p",
                         "Cohen_q");
  return(res.two)
  }

#' Test concerning differences between two dependent correlations
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

  z12 <- atanh(r12); z13 <- atanh(r13); r1 <- (r12 + r13) / 2;
  Cov.dep.a <- 1 / (1 - r1^2)^2; Cov.dep.b <- r23 * (1 - 2 * r1^2)
  Cov.dep.c <- -.50 * (r1^2); Cov.dep.d <- 1 - (2 * (r1^2)) - r23^2;
  Cov.dep <- (Cov.dep.a * Cov.dep.b) + (Cov.dep.c * Cov.dep.d);
  SE.dep <- sqrt((2 - (2 * Cov.dep)) / (n - 3));
  diff.z.dep <- round(((z12 - z13) / SE.dep), 2);

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.dep)), digit), scientific = F)};
  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.dep)), digit), scientific = F)};

  res.dep <- data.frame(round(r12, digit), round(r13, digit), round(r23, digit),
                        diff.z.dep, p);
  rownames(res.dep) <- cor.names;
  colnames(res.dep) <- c("r12", "r13", "r23", "z", "p");

  return(res.dep)
}
