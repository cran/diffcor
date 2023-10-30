\name{diffpwr.one}
\alias{diffpwr.one}
\title{Monte Carlo Simulation for the correlation difference between an expected
and an observed correlation}
\description{Computation of a Monte Carlo simulation to estimate the statistical
power the correlation difference between an observed correlation coefficient and
an a fixed value against which the correlation should be tested.}
\usage{
diffpwr.one(n, emp.r, hypo.r, alpha = .05, n.samples = 1000, seed = 1234)}

\arguments{
\item{n}{Sample size to be tested in the Monte Carlo simulation.}
\item{emp.r}{Assumed observed correlation.}
\item{hypo.r}{Correlation coefficient against which to test.}
\item{alpha}{Type I error. Default is .05.}
\item{n.samples}{Number of samples generated in the Monte Carlo simulation. The
recommended minimum is 1000 iterations, which is also the default.}
\item{seed}{To make the results reproducible, it is recommended to set a random
seed.}
}

\details{Depending on the number of generated samples (`n.samples`), correlation
coefficients of the size `emp.r` are simulated. Confidence intervals are built
around the simulated correlation coefficients. For each simulated coefficient,
it is then checked whether the hypothesized correlation cofficient (`hypo.r`)
falls within this interval. All correlations are automatically transformed with
the Fisher z-transformation prior to computations. The ratio of simulated
confidence intervals excluding the hypothesized coefficient equals the
statistical power, given the actual sample size (`n`; see Robert & Casella, 2010
<doi:10.1007/978-1-4419-1576-4>, for an overview of the Monte Carlo method).}

\value{
\item{emp.r}{Empirically observed correlation.}
\item{hypo.r}{Correlation against which `emp.r` should be tested.}
\item{n}{The sample size entered in the function.}
\item{cov}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed correlation `emp.r`. Should be between .91 and .98 (Muthén
& Muthén, 2002 <doi:10.1207/S15328007SEM0904_8>).}
\item{bias}{Average relative difference between the assumed `emp.r` and the
simulated correlations.}
\item{pwr}{Statistical power as the ratio of simulated confidence intervals
excluding the hypothesized correlation.}
}

\references{
Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to decide
on sample size and determine power. Structural Equation Modeling: A
Multidisciplinary Journal, 9(4), 599–620.
https://doi.org/10.1207/S15328007SEM0904_8

Robert, C., & Casella, G. (2010). Introducing Monte Carlo methods with R.
Springer. https://doi.org/10.1007/978-1-4419-1576-4
}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{
diffpwr.one(n.samples = 1000, n = 500, emp.r = .30, hypo.r = .40, alpha = .05,
            seed = 1234)
}
