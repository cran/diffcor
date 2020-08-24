#' Fisher's z-tests concerning difference of an observed correlation from given value or difference between independent correlations
#'
#' Test concerning differences between hypothesis and observation
#'
#' @param emp.r empirically observed correlation
#' @param hypo.r expected correlation
#' @param n sample size the observed correlation is based on
#' @return Fisher's z-values, corresponding p-values, and Cohen's q
#' @examples cor.1 <- diffcor.one(.76, .70, 271)
#' @export

diffcor.one <- function(emp.r, hypo.r, n, cor.names = NULL, alternative = c("one.sided", "two.sided")){
  fisher.emp <- atanh(emp.r);
  exp.r <- atanh(hypo.r)+(hypo.r/((2*n)-1));
  Var.z <- (n-3);

  diff.z.one <- (fisher.emp-exp.r)*sqrt(Var.z);
  Cohen.q.one <- fisher.emp-exp.r;

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 1)
    p <- 1-pnorm(abs(diff.z.one));
  if (tside == 2)
    p <- 2*pnorm(-abs(diff.z.one));

  res.one <- data.frame(diff.z.one, p, Cohen.q.one);
  colnames(res.one) <- c("Fisher z", "p value", "Cohen's q");
  row.names(res.one) <- cor.names;

  return(res.one)}

#' Test concerning differences between two independent correlations
#'
#' @param r1 empirically observed correlation in first group
#' @param r2 empirically observed correlation in second group
#' @param n1 sample size the first correlation is based on
#' @param n2 sample size the second correlation is based on
#'
#' @return Fisher's z-value, corresponding p-values, and Cohen's q
#' @examples diffcor.two(r1 = .76, r2 = .70, n1 = 271, n2 = 323)
#' @export

diffcor.two <- function(r1, r2, n1, n2, cor.names = NULL, alternative = c("one.sided", "two.sided")){
  fisher1 <- atanh(r1); fisher2 <- atanh(r2);
  Var1 <- 1/(n1-3); Var2 <- 1/(n2-3); SE <- sqrt((Var1+Var2));

  diff.z.two <- (fisher1-fisher2)/SE;
  Cohen.q.two <- (fisher1-fisher2)

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 1)
    p <- 1-pnorm(abs(diff.z.two));
  if (tside == 2)
    p <- 2*pnorm(-abs(diff.z.two));


  res.two <- data.frame(diff.z.two, p, Cohen.q.two);
  colnames(res.two) <- c("Fisher z", "p value", "Cohen's q");
  row.names(res.two) <- cor.names

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
#' @examples diffcor.dep(r12 = .76, r13 = .70, r23 = .50, n = 271)
#' @export

diffcor.dep <- function(r12, r13, r23, n, cor.names = NULL, alternative = c("one.sided", "two.sided")){
  z12 <- atanh(r12); z13 <- atanh(r13); r1 <- (r12+r13)/2;
  Cov.dep.a <- 1/(1-r1^2)^2; Cov.dep.b <- r23*(1-2*r1^2); Cov.dep.c <- -.50*(r1^2); Cov.dep.d <- 1-(2*(r1^2))-r23^2;
  Cov.dep <- (Cov.dep.a*Cov.dep.b)+(Cov.dep.c*Cov.dep.d);
  SE.dep <- sqrt((2-(2*Cov.dep))/(n-3));
  diff.z.dep <- (z12-z13)/SE.dep;

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 1)
    p <- 1-pnorm(abs(diff.z.dep));
  if (tside == 2)
    p <- 2*pnorm(-abs(diff.z.dep));


  res.dep <- data.frame(diff.z.dep, p);
  colnames(res.dep) <- c("Fisher z", "p value");
  row.names(res.dep) <- cor.names
  return(res.dep)
}
