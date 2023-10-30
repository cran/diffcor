\name{diffpwr.dep}
\alias{diffpwr.dep}
\title{Monte Carlo Simulation for the correlation difference between dependent
correlations}
\description{
Computation of a Monte Carlo simulation to estimate the statistical power of the
comparison between the correlations of a variable with two competing variables
that are also correlated with each other.}

\usage{
diffpwr.dep(n, r12, r13, r23, alpha = 0.05, n.samples = 1000, seed = 1234)}

\arguments{
\item{n}{Sample size to be tested in the Monte Carlo simulation.}
\item{r12}{Correlation between the criterion with which both competing variables
are correlated and the first of the two competing variables.}
\item{r13}{Correlation between the criterion with which both competing variables
are correlated and the second of the two competing variables.}
\item{r23}{Intercorrelation between the two competing variables.}
\item{alpha}{Type I error. Default is .05.}
\item{n.samples}{Number of samples generated in the Monte Carlo simulation. The
recommended minimum is 1000 iterations, which is also the default.}
\item{seed}{To make the results reproducible, it is recommended to set a random
seed.}
}

\details{Depending on the number of generated samples (`n.samples`), correlation
coefficients of the sizes `r12`, `r13`, and `r23` are simulated. For each
simulated sample, it is checked whether the correlations r12 and r13 differ,
given the correlation `r23`. The ratio of simulated z-tests of the correlation
difference tests exceeding the critical z-value, given the intended alpha-level,
equals the achieved statistical power (`n`; see Muthén & Muthén, 2002
<doi:10.1207/S15328007SEM0904_8>; Robert & Casella, 2010
<doi:10.1007/978-1-4419-1576-4>, for overviews of the Monte Carlo method).}

\value{
\item{r12}{Correlation between the criterion with which both competing variables
are correlated and the first of the two competing variables.}
\item{cov12}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed effect size r12.}
\item{bias12}{Average relative deviation of the simulated correlations r12 from
the intended one.}
\item{r13}{Correlation between the criterion with which both competing variables
are correlated and the second of the two competing variables.}
\item{cov13}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed effect size r13.}
\item{bias13}{Average relative deviation of the simulated correlations r13 from
the intended one.}
\item{r23}{Intercorrelation between the two competing variables.}
\item{cov23}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed effect size r23.}
\item{bias23}{Average relative deviation of the simulated correlations r23 from
the intended one.}
\item{n}{Sample size to be tested in the Monte Carlo simulation.}
\item{pwr}{Statistical power as the ratio of simulated difference tests that
yielded significance.}

Biases should be as close to zero as possible and coverage should be ideally
between .91 and .98 (Muthén & Muthén, 2002 <doi:10.1207/S15328007SEM0904_8>).}

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

\examples{diffpwr.dep(n.samples = 1000, n = 250, r12 = .30, r13 = .45,
                      r23 = .50, alpha = .05, seed = 1234)}
