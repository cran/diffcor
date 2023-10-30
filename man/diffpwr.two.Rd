\name{diffpwr.two}
\alias{diffpwr.two}
\title{Monte Carlo Simulation for the correlation difference between two
correlations that were observed in two independent samples}
\description{Computation of a Monte Carlo simulation to estimate the statistical
power the correlation difference between the correlation coefficients detected
in two indepdenent samples (e.g., original study and replication study).}
\usage{
diffpwr.two(n1, n2, r1, r2, alpha = .05, n.samples = 1000, seed = 1234)}

\arguments{
\item{n1}{Sample size to be tested in the Monte Carlo simulation for the first
sample.}
\item{n2}{Sample size to be tested in the Monte Carlo simulation for the
second sample.}
\item{r1}{Correlarion observed in the first sample.}
\item{r2}{Correlarion observed in the second sample.}
\item{alpha}{Type I error. Default is .05.}
\item{n.samples}{Number of samples generated in the Monte Carlo simulation. The
recommended minimum is 1000 iterations, which is also the default.}
\item{seed}{To make the results reproducible, a random seed is specified.}
}

\details{Depending on the number of generated samples (`n.samples`), correlation
coefficients of the sizes `r1` and `r2` are simulated. For each simulated pair
of coefficients, it is then checked whether the confidence intervals (with
given alpha level) of the correlations overlap. All correlations are
automatically transformed with the Fisher z-transformation prior to
computations. The ratio of simulated non-overlapping confidence intervals equals
the statistical power, given the actual sample sizes (`n1` and `n2`; see Robert
& Casella, 2010 <doi:10.1007/978-1-4419-1576-4>, for an overview of the Monte
Carlo method).}

\value{
\item{r1}{Correlation observed in sample 1.}
\item{n1}{The sample size of the first sample.}
\item{cov1}{Coverage. Ratio of simulated confidence intervals including r1.}
\item{bias1}{Average relative difference between r1 and simulated correlations.}
\item{r2}{Correlation observed in sample 2.}
\item{n2}{The sample size of the second sample.}
\item{cov2}{Coverage. Ratio of simulated confidence intervals including r2.}
\item{bias2}{Average relative difference between r2 and simulated correlations.}
\item{pwr}{Statistical power as the ratio of simulated non-verlapping confidence
intervals.}

Biases should be as close to zero as possible and coverage should be ideally
between .91 and .98 (Muthén & Muthén, 2002 <doi:10.1207/S15328007SEM0904_8>).}

\references{
Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to decide
on sample size and determine power. Structural Equation Modeling: A
Multidisciplinary Journal, 9(4), 599–620.
https://doi.org/10.1207/S15328007SEM0904_8

Robert, C., & Casella, G. (2010). Introducing Monte Carlo methods with R.
Springer. https://doi.org/10.1007/978-1-4419-1576-4}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{
diffpwr.two(n.samples = 1000, n1 = 1000, n2 = 594, r1 = .45, r2 = .39,
            alpha = .05, seed = 1234)}
