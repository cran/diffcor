\name{diffcor.one}
\alias{diffcor.one}
\title{Fisher's z-test of difference between an empirical and a hypothesized correlation}
\usage{
diffcor.one(emp.r, hypo.r, n, alpha = .05, cor.names = NULL,
alternative = c("one.sided", "two.sided"), digit = 3)}

\arguments{
\item{emp.r}{Empirically observed correlation}
\item{hypo.r}{Hypothesized correlation which shall be tested}
\item{n}{Sample size in which the observed effect was found}
\item{alpha}{Likelihood of Type I error, DEFAULT = .05}
\item{cor.names}{OPTIONAL, label for the correlation (e.g., "IQ-performance").
DEFAULT is NULL}
\item{digit}{Number of digits in the output for all parameters, DEFAULT = 3}
\item{alternative}{A character string specifying if you wish to test one-sided
or two-sided differences}}

\description{
The function tests whether an observed correlation differs from an expected one,
for example, in construct validation. The output provides the compared
correlations, a z-score, a p-value, a confidence interval, and the effect size
Cohens q. According to Cohen (1988), q = |.10|, |.30| and |.50| are considered
small, moderate, and large differences, respectively.}

\value{
\item{r_exp}{Vector of the expected correlations}
\item{r_obs}{Vector of the empirically observed correlations}
\item{LL}{Lower limit of the confidence interval of the empirical correlation,
given the specified alpha level, DEFAULT = 95 percent}
\item{UL}{Upper limit of the confidence interval of the empirical correlation,
given the specified alpha level, DEFAULT = 95 percent}
\item{z}{Test statistic for correlation difference in units of z distribution}
\item{p}{p value for one- or two-sided testing, depending on alternative =
c("one.sided", "two.sided)}
\item{Cohen_q}{Effect size measure for differences of independent correlations}
}

\examples{
diffcor.one(c(.76, .53, -.32), c(.70, .35, -.40),
  c(225, 250, 210),
  cor.names = c("a-b", "c-d", "e-f"), digit = 2, alternative = "one.sided")
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
