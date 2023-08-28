\name{diffpwr.one}
\alias{diffpwr.one}
\title{Monte Carlo Simulation for the correlation difference between an expected
and an observed correlation}
\description{Computation of a Monte Carlo simulation to estimate the statistical
power the correlation difference between an observed correlation coefficient and
an a fixed value against which the correlation should be tested.}
\usage{
diffpwr.one(n.samples = 1000, n, emp.r, hypo.r, alpha = .05, seed = 1234)}

\arguments{
  \item{n.samples}{Number of samples generated in the Monte Carlo simulation.
  The recommended default is 1000 iterations.}
  \item{n}{Sample size to be tested in the Monte Carlo simulation.}
  \item{emp.r}{Assumed observed correlation.}
  \item{hypo.r}{Correlation coefficient against which to test.}
  \item{alpha}{Type I error. Default is .05.}
  \item{seed}{To make the results reproducible, it is recommended to set a
  random seed.}
}
\details{Depending on the number of generated samples (`n.samples`), correlation
coefficients of the size `emp.r` are simulated. Confidence intervals are built
around the simulated correlation coefficients. For each simulated coefficient,
it is then checked whether the hypothesized correlation cofficient (`hypo.r`)
falls within this interval. The ratio of simulated confidence intervals
excluding the hypothesized coefficient equals the statistical power, given the
actual sample size (`n`; see Robert & Casella, 2010
<doi:10.1007/978-1-4419-1576-4>, for an overview of the Monte Carlo method).}

\value{
\item{n}{The sample size entered in the function.}
\item{pwr}{Statistical power as the ratio of simulated confidence intervals
excluding the hypothesized correlation.}
In addition, the correlation difference is visualized as a plot.
}

\references{
Robert, C., & Casella, G. (2010). Introducing Monte Carlo methods with R.
Springer. https://doi.org/10.1007/978-1-4419-1576-4}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{
diffpwr.one(n.samples = 1000, n = 500, emp.r = .30, hypo.r = .40, alpha = .05,
            seed = 1234)
}
