\name{diffcor.two}
\alias{diffcor.two}
\title{Fisher's z-Tests for differences of correlations in two independent samples}
\usage{
diffcor.two(r1, r2, n1, n2, alpha = .05, cor.names = NULL,
alternative = c("one.sided", "two.sided"), digit = 3)
}

\arguments{
\item{r1}{Correlation coefficient in first sample}
\item{r2}{Correlation coefficient in second sample}
\item{n1}{First sample size}
\item{n2}{Second sample size}
\item{alpha}{Likelihood of Type I error, DEFAULT = .05}
\item{cor.names}{OPTIONAL, label for the correlation (e.g., "IQ-performance").
DEFAULT is NULL}
\item{digit}{Number of digits in the output for all parameters, DEFAULT = 3}
\item{alternative}{A character string specifying if you wish to test one-sided
or two-sided differences}
}

\description{
The function tests whether the correlation between two variables differs across
two independent studies/samples. The output provides the compared correlations,
test statistic as z-score, p-values, confidence intervals of the empirical
correlations, and the effect size Cohens q. According to Cohen (1988), q = |.10|,
|.30| and |.50| are considered small, moderate, and large differences,
respectively.}

\value{
\item{r1}{Vector of the empirically observed correlations in the first sample}
\item{r2}{Vector of the empirically observed correlations in the second sample}
\item{LL1}{Lower limit of the confidence interval of the first empirical
correlation, given the specified alpha level, DEFAULT = 95 percent}
\item{UL1}{Upper limit of the confidence interval of the first empirical
correlation, given the specified alpha level, DEFAULT = 95 percent}
\item{LL2}{Lower limit of the confidence interval of the second empirical
correlation, given the specified alpha level, DEFAULT = 95 percent}
\item{UL2}{Upper limit of the confidence interval of the second empirical
correlation, given the specified alpha level, DEFAULT = 95 percent}
\item{z}{Test statistic for correlation difference in units of z distribution}
\item{p}{p value for one- or two-sided testing, depending on alternative =
c("one.sided", "two.sided)}
\item{Cohen_q}{Effect size measure for differences of independent correlations}
}

\examples{
diffcor.two(r1 = c(.39, .52, .22),
  r2 = c(.29, .44, .12),
  n1 = c(66, 66, 66), n2 = c(96, 96, 96), alpha = .01,
  cor.names = c("a-b", "c-d", "e-f"), alternative = "one.sided")
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
