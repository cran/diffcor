\name{diffpwr.two}
\alias{diffpwr.two}
\title{Monte Carlo Simulation for the correlation difference between two
correlations that were observed in two independent samples}
\description{Computation of a Monte Carlo simulation to estimate the statistical
power the correlation difference between the correlation coefficients detected
in two indepdenent samples (e.g., original study and replication study).}
\usage{
diffpwr.two(n.samples = 1000, n1, n2, r1, r2, alpha = .05, seed = 1234)}

\arguments{
  \item{n.samples}{Number of samples generated in the Monte Carlo simulation.
  The recommended default is 1000 iterations.}
  \item{n1}{Sample size to be tested in the Monte Carlo simulation for the
  first sample.}
  \item{n2}{Sample size to be tested in the Monte Carlo simulation for the
  second sample.}
  \item{r1}{Correlarion observed in the first sample.}
  \item{r2}{Correlarion observed in the second sample.}
  \item{alpha}{Type I error. Default is .05.}
  \item{seed}{To make the results reproducible, a random seed is specified.}
}
\details{Depending on the number of generated samples (`n.samples`), correlation
coefficients of the sizes `r1` and `r2` are simulated. For each simulated pair
of coefficients, it is then checked whether the confidence intervals (with
given alpha level) of the correlations overlap. The ratio of simulated non-
overlapping confidence intervals equals the statistical power, given the
actual sample sizes (`n1` and `n2`; see Robert & Casella, 2010
<doi:10.1007/978-1-4419-1576-4>, for an overview of the Monte Carlo method).}

\value{
\item{r1}{The correlation coefficient from the first sample.}
\item{n1}{The sample size of the first sample.}
\item{r2}{The correlation coefficient from the second sample.}
\item{n2}{The sample size of the second sample.}
\item{pwr}{Statistical power as the ratio of simulated non-verlapping confidence
intervals.}
In addition, the correlation difference is visualized as a plot.
}

\references{
Robert, C., & Casella, G. (2010). Introducing Monte Carlo methods with R.
Springer. https://doi.org/10.1007/978-1-4419-1576-4}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{
diffpwr.two(n.samples = 1000, n1 = 1000, n2 = 594, r1 = .45, r2 = .39,
            alpha = .05, seed = 1234)}
