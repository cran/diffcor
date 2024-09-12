\name{bootcor.one}
\alias{bootcor.one}
\title{Bootstrapped Correlation Difference Test between an Empirical and an
Expected Correlation}
\description{Derivation of bootstrap confidence intervals for the calculation of
correlation differences between the empirically observed correlation coefficient
and a threshold against which this coefficient is tested.}

\usage{
bootcor.one(x,
            y,
            r_target,
            k = 5000,
            alpha = .05,
            digit = 3,
            seed = 1234)}

\arguments{
  \item{x}{A vector containing the values of the first variable being involved
  in the correlation.}
  \item{y}{A vector containing the values of the second variable being involved
  in the correlation.}
  \item{r_target}{A single value against which the correlation between x and y
  is tested.}
  \item{k}{The number of bootstrap samples to be drawn. The default is 5000.}
  \item{alpha}{Likelihood of Type I error. The default is .05.}
  \item{digit}{Number of digits in the output. The default is 3.}
  \item{seed}{A random seed to make the results reproducible.}
}

\details{
Bivariate correlation analyses as well as correlation difference tests possess
very strict statistical requirements that are not necessarily fulfilled when
using the basic \code{diffcor.one()} function from this package (Wilcox, 2013
<doi:10.1016/C2010-0-67044-1>). For instance, if the assumption of a normal
distribution does not hold, the significance test can lead to false positive or
false negative conclusions. To address potential deviations from normal
distribution, the present function applies bootstrapping to the data. The output
provides a confidence interval for the difference between the empirically
observed correlation coefficient and the threshold against which this
coefficient should be tested, whereby the interval is derived from bootstrapping
samples.}

\value{
\item{r_emp}{The empircally observed correlation between x and y.}
\item{r_target}{The threshold against which r_emp is tested.}
\item{M}{Mean of the confidence interval of the correlation difference between
\code{r_emp} and \code{r_target}.}
\item{LL}{Lower limit of the confidence interval of the correlation difference
between \code{r_emp} and \code{r_target}, given the entered Type I-level.}
\item{UL}{Upper limit of the confidence interval of the correlation difference
between \code{r_emp} and \code{r_target}, given the entered Type I-level.}
}

\references{
Wilcox, R. (2013). Introduction to robust estimation and hypothesis testing.
Elsevier. https://doi.org/10.1016/C2010-0-67044-1
}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{

df <- data.frame(a = rnorm(1000),
                 b = rnorm(1000))

bootcor.one(x = df$a,
            y = df$b,
            r_target = .10,
            k = 5000,
            alpha = .05,
            digit = 3,
            seed = 1234)
}
