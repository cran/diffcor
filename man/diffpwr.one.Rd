\name{diffpwr.one}
\alias{diffpwr.one}
\title{Monte Carlo Simulation for the correlation difference between an expected
and an observed correlation}
\description{Computation of a Monte Carlo simulation to estimate the statistical
power the correlation difference between an observed correlation coefficient and
an a fixed value against which the correlation should be tested.}
\usage{
diffpwr.one(n,
            r,
            rho,
            alpha = .05,
            n.samples = 1000,
            seed = 1234)}

\arguments{
\item{n}{Sample size to be tested in the Monte Carlo simulation.}
\item{r}{Assumed observed correlation.}
\item{rho}{Correlation coefficient against which to test (reflects the null
hypothesis).}
\item{alpha}{Type I error. Default is .05.}
\item{n.samples}{Number of samples generated in the Monte Carlo simulation. The
recommended minimum is 1,000 iterations, which is also the default.}
\item{seed}{To make the results reproducible, it is recommended to set a random
seed.}
}

\details{Depending on the number of generated samples (n.samples), correlation
coefficients of size r are simulated. Confidence intervals are constructed
around the simulated correlation coefficients. For each simulated coefficient,
it is then checked whether the hypothesized correlation cofficient (rho) falls
within this interval. All correlations are automatically transformed with the
Fisher z-transformation prior to computations. The ratio of simulated confidence
intervals excluding the hypothesized coefficient equals the statistical power,
given the intended alpha-level and sample size (see Robert & Casella, 2010
<doi:10.1007/978-1-4419-1576-4>, for an overview of the Monte Carlo method).

It should be noted that the Pearson correlation coefficient is sensitive to
linear association, but also to a host of statistical issues such as univariate
and bivariate outliers, range restrictions, and heteroscedasticity (e.g., Duncan
& Layard, 1973 <doi:10.1093/BIOMET/60.3.551>; Wilcox, 2013
<doi:10.1016/C2010-0-67044-1>). Thus, every power analysis requires that specific
statistical prerequisites are fulfilled and can be invalid with regard to the
actual data if the prerequisites do not hold, potentially biasing Type I error
rates.}

\value{
As dataframe with the following parameters
\item{r}{Empirically observed correlation.}
\item{rho}{Correlation against which r should be tested.}
\item{n}{The sample size entered in the function.}
\item{cov}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed correlation r. Should be between .91 and .98 (Muthén
& Muthén, 2002 <doi:10.1207/S15328007SEM0904_8>).}
\item{bias_M}{Difference between the mean of the distribution of the simulated
correlations and rho, divided by rho.}
\item{bias_Md}{Difference between the median of the distribution of the
simulated correlations and rho, divided by rho.}
\item{pwr}{Statistical power as the ratio of simulated confidence intervals
excluding rho.}
}

\references{
Duncan, G. T., & Layard, M. W. (1973). A Monte-Carlo study of asymptotically
robust tests for correlation coefficients. Biometrika, 60, 551–558.
https://doi.org/10.1093/BIOMET/60.3.551

Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to decide
on sample size and determine power. Structural Equation Modeling: A
Multidisciplinary Journal, 9(4), 599–620.
https://doi.org/10.1207/S15328007SEM0904_8

Robert, C., & Casella, G. (2010). Introducing Monte Carlo methods with R.
Springer. https://doi.org/10.1007/978-1-4419-1576-4

Wilcox, R. (2013). Introduction to robust estimation and hypothesis testing.
Elsevier. https://doi.org/10.1016/C2010-0-67044-1}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{
diffpwr.one(n = 500,
            r = .30,
            rho = .40,
            alpha = .05,
            n.samples = 1000,
            seed = 1234)
}
