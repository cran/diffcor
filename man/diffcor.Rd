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
\item{emp.r}{Empirically observed correlation}
\item{hypo.r}{Hypothesized correlation which shall be tested}
\item{n}{Sample size in which the observed effect was found}
\item{alpha}{Likelihood of Type I error, DEFAULT = .05}
\item{cor.names}{Optional, label for the correlation (e.g., "IQ-performance").
Per default, cor.names is NULL}
\item{digit}{Number of digits in the output for all parameters, DEFAULT = 3}
\item{alternative}{A character string specifying if you wish to test one-sided
or two-sided differences}
\item{r1}{First correlation coefficient}
\item{r2}{Second correlation coefficient}
\item{n1}{Sample size in which the first correlation coeffient was observed}
\item{n2}{Sample size in which the second correlation coeffient was observed}
\item{r12}{Correlation between the first and the second construct}
\item{r13}{Correlation between the first and the third construct}
\item{r23}{Correlation between the second and the third construct}
}

\description{
Computations of Fisher's z-tests concerning differences between correlations.

diffcor.one() can be used to test if an observed correlation differs from an
expected one, for example, in construct validation. diffcor.two() can be used
to test if the correlation between two variables differs across two independent
studies. diffcor.dep() can be applied to check if the correlation between two
variables (r12) differs from the correlation between the first and a third one
(r13), given the intercorrelation of the compared constructs (r23). All outputs
provide the compared correlations, test statistic as z-score, and p-values. For
diffcor.one() and diffcor.two(), the output further provides confidence intervals
of the empirical correlations and the effect size Cohens q. According to Cohen
(1988), q = |.10|, |.30| and |.50| are considered small, moderate, and large
differences, respectively.}

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
sciences (2nd ed.). Lawrence Erlbaum.

Eid, M., Gollwitzer, M., & Schmitt, M. (2015). Statistik und Forschungsmethoden
(4.Auflage) [Statistics and research methods (4th ed.)]. Beltz.

Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix.
Psychological Bulletin, 87, 245-251.}
