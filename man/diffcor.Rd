\name{diffcor}
\alias{diffcor.one}
\alias{diffcor.two}
\alias{diffcor.dep}
\title{Fisher's z-Tests Concerning Difference of Correlations}
\usage{
diffcor.one(emp.r, hypo.r, n, cor.names = NULL,
alternative = c("one.sided", "two.sided"), bonferroni = FALSE, digit = 3, alpha = .05)

diffcor.two(r1, r2, n1, n2, cor.names = NULL,
alternative = c("one.sided", "two.sided"), bonferroni = FALSE, digit = 3, alpha = .05)

diffcor.dep(r12, r13, r23, n, cor.names = NULL,
alternative = c("one.sided", "two.sided"), bonferroni = FALSE, digit = 3, alpha = .05)
}

\arguments{
\item{emp.r}{empirically observed correlation}
\item{hypo.r}{hypothesized correlation which shall be tested}
\item{n}{sample size in which the observed effect was found}
\item{cor.names}{Label for the respective correlation (e.g., "IQ-performance"). This is especially advantageous if you test a series of correlations simultaneously, see examples. Per default, cor.names is NULL}
\item{bonferroni}{a character string to specify of you wish to correct for multiple testing, using the Bonferroni correction}
\item{digit}{number of decimal places for p-values}
\item{alternative}{a character string to specify if you wish to test one-sided or two-sided differences}
\item{alpha}{level of significance. Per default, alpha is set .05}
\item{r1}{first correlation coefficient}
\item{r2}{second correlation coefficient}
\item{n1}{sample size the first correlation coeffient is based on}
\item{n2}{sample size the second correlation coeffient is based on}
\item{r12}{correlation between first and second construct}
\item{r13}{correlation between first and third construct}
\item{r23}{correlation between second and third construct}
}
\description{
Computations of Fisher's z-tests concerning differences between correlations.

diffcor.one() could be used to test for differences regarding an expected value, e.g., in construct validation.

diffcor.two() may be useful in replication studies, to test if the original study and the replication study differed in terms of effects.

diffcor.dep() can be applied to check if the correlation between one construct with another one (r12) is significantly different/higher/smaller than the correlation of one of the constructs with a third construct (r13), given the correlation of the constructs that are compared (r23).

The outputs for all the three functions provide the test statistic in z-units, p-values, and alpha levels that were corrected in terms of multiple testing according to Bonferroni. To help interpret the output, the procedure prompts if a single p value is smaller than the corrected alpha. For diffcor.one() and diffcor.two(), the effect size Cohens q is additionally provided. According to Cohen (1988), q = |.10|, |.30| and |.50| are considered small, moderate, and large differences, respectively.
}

\value{
\item{Fisher_z}{test statistic for correlation difference in units of the standard normal distribution}
\item{p}{p-value for one- or two-sided testing, depending on alternative = c("one.sided", "two.sided)}
\item{alpha_Bonferroni}{alpha level which is corrected for multiple testing; not displayed if bonferroni = FALSE; per default, bonferroni is set TRUE}
\item{Bonferroni_significant}{prompts if a single correlation is significantly different, considering the chosen alpha level and the number of correlations tested}
\item{Cohen_q}{effect size measure for differences of independent correlations}
}

\examples{
# empirically observed r = .76, expected r = .70,
# sample size n = 271, alpha level = .05 (default)
diffcor.one(.76, .70, 271, cor.names = NULL, alternative = "two.sided")

# works also with a dataframe of correlations.

diffcor.one(c(.76, .53, -.32), c(.70, .35, -.40),
  c(225, 250, 210),
  cor.names = c("a-b", "c-d", "e-f"), digit = 2, alternative = "one.sided")

# same applies to the functions diffcor.two() and diffcor.dep()

diffcor.two(r1 = .76, r2 = .70, n1 = 271, n2 = 323,
cor.names = NULL, alternative = "one.sided", alpha = .10)

diffcor.two(r1 = c(.39, .52, .22),
  r2 = c(.29, .44, .12),
  n1 = c(66, 66, 66), n2 = c(96, 96, 96),
  cor.names = c("a-b", "c-d", "e-f"), alternative = "one.sided")

diffcor.two(r1 = c(.39, .52, .22),
  r2 = c(.29, .44, .12),
  n1 = c(66, 66, 66),
  n2 = c(96, 96, 96),
  cor.names = NULL, bonferroni = FALSE, alternative = "one.sided")

diffcor.dep(r12 = .76, r13 = .70, r23 = .50, n = 271, digit = 4,
cor.names = NULL, alternative = "two.sided")

diffcor.dep(r12 = c(.76, .54, .22),
  r13 = c(.24, .01, .07), r23 = c(.25, .65, .90),
  n = c(500, 392, 92), cor.names = c("a-b", "c-d", "e-f"),
  alternative = "one.sided", alpha = .01)

diffcor.dep(r12 = c(.76, .54, .22),
  r13 = c(.24, .01, .07),
  r23 = c(.25, .65, .90),
  n = c(500, 392, 92), cor.names = NULL,
  alternative = "one.sided")}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}
\references{Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Eid, M., Gollwitzer, M., & Schmitt, M. (2015). Statistik und Forschungsmethoden (4.Auflage) [Statistics and research methods (4th ed.)]. Weinheim Basel: Beltz.

Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87, 245-251.}
