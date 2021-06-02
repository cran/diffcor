\name{diffcor}
\alias{diffcor.one}
\alias{diffcor.two}
\alias{diffcor.dep}
\title{Fisher's z-Tests Concerning Difference of Correlations}
\usage{
diffcor.one(emp.r, hypo.r, n, alpha = .05, cor.names = NULL,
alternative = c("one.sided", "two.sided"), digit = 3)

diffcor.two(r1, r2, n1, n2, alpha = .05, cor.names = NULL,
alternative = c("one.sided", "two.sided"), digit = 3)

diffcor.dep(r12, r13, r23, n, cor.names = NULL,
alternative = c("one.sided", "two.sided"), digit = 3)
}

\arguments{
\item{emp.r}{empirically observed correlation}
\item{hypo.r}{hypothesized correlation which shall be tested}
\item{n}{sample size in which the observed effect was found}
\item{alpha}{Likelihood of Type I error, DEFAULT = .05}
\item{cor.names}{Label for the respective correlation (e.g., "IQ-performance").
This is especially advantageous if you test a series of correlations
simultaneously, see examples. Per default, cor.names is NULL}
\item{digit}{number of decimal places for p-values, DEFAULT = 3}
\item{alternative}{a character string to specify if you wish to test one-sided
or two-sided differences}
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

diffcor.one() could be used to test for differences regarding an expected value,
e.g., in construct validation.

diffcor.two() may be useful in replication studies, to test if the original study
and the replication study differed in terms of effects.

diffcor.dep() can be applied to check if the correlation between one construct
with another one (r12) is significantly different/higher/smaller than the
correlation of one of the constructs with a third construct (r13), given the
correlation of the constructs that are compared (r23).

The outputs for all the three functions provide the compared correlations, test
statistic in z-units, and p-values. For diffcor.one() and diffcor.two(), the
output further provides confidence intervals of the empirical correlations and
the effect size Cohens q. According to Cohen (1988), q = |.10|, |.30| and |.50|
are considered small, moderate, and large differences, respectively.}

\value{
\item{r_exp}{Vector of the expected correlations in diffcor.one}
\item{r_obs}{Vector of the empirically observed correlations in diffcor.one}
\item{r1}{Vector of the empirically observed correlations in the first sample}
\item{r2}{Vector of the empirically observed correlations in the second sample}
\item{r12}{Vector of the empirically observed correlations between the first and
the second construct in diffcor.dep}
\item{r13}{Vector of the empirically observed correlations between the first and
the third construct in diffcor.dep}
\item{r23}{Vector of the empirically observed correlations between the second and
the third construct in diffcor.dep}
\item{LL}{Lower limit of the confidence interval of the empirical correlation in
diffcor.one, given the specified alpha level, DEFAULT = 95 percent}
\item{UL}{Upper limit of the confidence interval of the empirical correlation in
diffcor.one, given the specified alpha level, DEFAULT = 95 percent}
\item{LL1}{Lower limit of the confidence interval of the first empirical
correlation in diffcor.two, given the specified alpha level, DEFAULT = 95 percent}
\item{UL1}{Upper limit of the confidence interval of the first empirical
correlation in diffcor.two, given the specified alpha level, DEFAULT = 95 percent}
\item{LL2}{Lower limit of the confidence interval of the second empirical
correlation in diffcor.two, given the specified alpha level, DEFAULT = 95 percent}
\item{UL2}{Upper limit of the confidence interval of the second empirical
correlation in diffcor.two, given the specified alpha level, DEFAULT = 95 percent}
\item{z}{Test statistic for correlation difference in units of z distribution}
\item{p}{p value for one- or two-sided testing, depending on alternative =
c("one.sided", "two.sided)}
\item{Cohen_q}{Effect size measure for differences of independent correlations}
}

\examples{
diffcor.one(c(.76, .53, -.32), c(.70, .35, -.40),
  c(225, 250, 210),
  cor.names = c("a-b", "c-d", "e-f"), digit = 2, alternative = "one.sided")

diffcor.two(r1 = c(.39, .52, .22),
  r2 = c(.29, .44, .12),
  n1 = c(66, 66, 66), n2 = c(96, 96, 96), alpha = .01,
  cor.names = c("a-b", "c-d", "e-f"), alternative = "one.sided")

diffcor.dep(r12 = .76, r13 = .70, r23 = .50, n = 271, digit = 4,
cor.names = NULL, alternative = "two.sided")
}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}
\references{Cohen, J. (1988). Statistical power analysis for the behavioral
sciences (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Eid, M., Gollwitzer, M., & Schmitt, M. (2015). Statistik und Forschungsmethoden
(4.Auflage) [Statistics and research methods (4th ed.)]. Weinheim Basel: Beltz.

Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix.
Psychological Bulletin, 87, 245-251.}
