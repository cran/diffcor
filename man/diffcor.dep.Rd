\name{diffcor.dep}
\alias{diffcor.dep}
\title{Fisher's z-Tests of dependent correlations}
\usage{
diffcor.dep(r12, r13, r23, n, cor.names = NULL,
alternative = c("one.sided", "two.sided"), digit = 3)
}

\arguments{
\item{r12}{Correlation between the first and the second construct}
\item{r13}{Correlation between the first and the third construct}
\item{r23}{Correlation between the second and the third construct}
\item{n}{Sample size in which the observed effect was found}
\item{cor.names}{OPTIONAL, label for the correlation. DEFAULT is NULL}
\item{alternative}{A character string specifying if you wish to test one-sided
or two-sided differences}
\item{digit}{Number of digits in the output for all parameters, DEFAULT = 3}
}

\description{
The function checks if the correlation between two variables (r12) differs from
the correlation between the first and a third one (r13), given the
intercorrelation of the compared constructs (r23). The output provides the
compared correlations, test statistic as z-score, and p-values.}

\value{
\item{r12}{Vector of the empirically observed correlations between the first and
the second construct}
\item{r13}{Vector of the empirically observed correlations between the first and
the third construct}
\item{r23}{Vector of the empirically observed correlations between the second and
the third construct}
\item{z}{Test statistic for correlation difference in units of z distribution}
\item{p}{p value for one- or two-sided testing, depending on alternative =
c("one.sided", "two.sided)}
}

\examples{
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
