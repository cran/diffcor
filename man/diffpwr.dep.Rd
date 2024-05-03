\name{diffpwr.dep}
\alias{diffpwr.dep}
\title{Monte Carlo Simulation for the correlation difference between dependent
correlations}
\description{
Computation of a Monte Carlo simulation to estimate the statistical power of the
comparison between the correlations of a variable with two competing variables
that are also correlated with each other.}

\usage{
diffpwr.dep(n,
            rho12,
            rho13,
            rho23,
            alpha = 0.05,
            n.samples = 1000,
            seed = 1234)}

\arguments{
\item{n}{Sample size to be tested in the Monte Carlo simulation.}
\item{rho12}{Assumed population correlation between the criterion with which
both competing variables are correlated and the first of the two competing
variables.}
\item{rho13}{Assumed population correlation between the criterion with which
both competing variables are correlated and the second of the two competing
variables.}
\item{rho23}{Assumed population correlation between the two competing variables.}
\item{alpha}{Type I error. Default is .05.}
\item{n.samples}{Number of samples generated in the Monte Carlo simulation. The
recommended minimum is 1,000 iterations, which is also the default.}
\item{seed}{To make the results reproducible, it is recommended to set a random
seed.}
}

\details{Depending on the number of generated samples (n.samples), correlation
coefficients simulated. For each simulated sample, it is checked whether the
correlations r12 and r13 differ, given the correlation r23. The ratio of
simulated z-tests of the correlation difference tests exceeding the critical
z-value, given the intended alpha-level and sample size, equals the achieved
statistical power(see Muthén & Muthén, 2002 <doi:10.1207/S15328007SEM0904_8>;
Robert & Casella, 2010 <doi:10.1007/978-1-4419-1576-4>, for overviews of the
Monte Carlo method).

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
\item{rho12}{Assumed population correlation between the criterion with which
both competing variables are correlated and the first of the two competing
variables.}
\item{cov12}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed effect size rho12.}
\item{bias12_M}{Difference between the mean of the distribution of the simulated
correlations and rho12, divided by rho12.}
\item{bias12_Md}{Difference between the median of the distribution of the
simulated correlations and rho12, divided by rho12.}
\item{rho13}{Assumed population correlation between the criterion with which
both competing variables are correlated and the second of the two competing
variables.}
\item{cov13}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed effect size rho13.}
\item{bias13_M}{Difference between the mean of the distribution of the simulated
correlations and rho13, divided by rho13.}
\item{bias13_Md}{Difference between the median of the distribution of the
simulated correlations and rho13, divided by rho13.}
\item{rho23}{Assumed population correlation between the two competing variables.
}
\item{cov23}{Coverage. Indicates the ratio of simulated confidence intervals
including the assumed effect size rho23.}
\item{bias23_M}{Difference between the mean of the distribution of the simulated
correlations and rho23, divided by rho23.}
\item{bias23_Md}{Difference between the median of the distribution of the
simulated correlations and rho23, divided by rho23.}
\item{n}{Sample size to be tested in the Monte Carlo simulation.}
\item{pwr}{Statistical power as the ratio of simulated difference tests that
yielded statistical significance.}

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

\examples{diffpwr.dep(n.samples = 1000,
                      n = 250,
                      rho12 = .30,
                      rho13 = .45,
                      rho23 = .50,
                      alpha = .05,
                      seed = 1234)}
