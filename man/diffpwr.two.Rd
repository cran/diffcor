\name{diffpwr.two}
\alias{diffpwr.two}
\title{Monte Carlo Simulation for the correlation difference between two
correlations that were observed in two independent samples}
\description{Computation of a Monte Carlo simulation to estimate the statistical
power the correlation difference between the correlation coefficients detected
in two indepdenent samples (e.g., original study and replication study).}
\usage{
diffpwr.two(n1,
            n2,
            rho1,
            rho2,
            alpha = .05,
            n.samples = 1000,
            seed = 1234)}

\arguments{
\item{n1}{Sample size to be tested in the Monte Carlo simulation for the first
sample.}
\item{n2}{Sample size to be tested in the Monte Carlo simulation for the
second sample.}
\item{rho1}{Assumed population correlation to be observed in the first sample.}
\item{rho2}{Assumed population correlation to be observed in the second sample.}
\item{alpha}{Type I error. Default is .05.}
\item{n.samples}{Number of samples generated in the Monte Carlo simulation. The
recommended minimum is 1,000 iterations, which is also the default.}
\item{seed}{To make the results reproducible, a random seed is specified.}
}

\details{Depending on the number of generated samples (n.samples), correlation
coefficients are simulated. For each simulated pair of coefficients, it is then
checked whether the confidence intervals (with given alpha level) of the
correlations overlap. All correlations are automatically transformed with the
Fisher z-transformation prior to computations. The ratio of simulated non-
overlapping confidence intervals equals the statistical power, given the
alpha-level and sample sizes (see Robert & Casella, 2010
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
\item{rho1}{Assumed population correlation to be observed in the first sample.}
\item{n1}{Sample size of the first sample.}
\item{cov1}{Coverage. Ratio of simulated confidence intervals including rho1.}
\item{bias1_M}{Difference between the mean of the distribution of the simulated
correlations and rho1, divided by rho1.}
\item{bias1_Md}{Difference between the median of the distribution of the
simulated correlations and rho1, divided by rho1.}
\item{rho2}{Assumed population correlation to be observed in the second sample.}
\item{n2}{The sample size of the second sample.}
\item{cov2}{Coverage. Ratio of simulated confidence intervals including rho2.}
\item{bias2_M}{Difference between the mean of the distribution of the simulated
correlations and rho2, divided by rho2.}
\item{bias2_Md}{Difference between the median of the distribution of the
simulated correlations and rho2, divided by rho2.}
\item{pwr}{Statistical power as the ratio of simulated non-verlapping confidence
intervals.}

Biases should be as close to zero as possible and coverage should be ideally
between .91 and .98 (Muthén & Muthén, 2002 <doi:10.1207/S15328007SEM0904_8>).}

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
diffpwr.two(n1 = 1000,
            n2 = 594,
            rho1 = .45,
            rho2 = .39,
            alpha = .05,
            n.samples = 1000,
            seed = 1234)}
