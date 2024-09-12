\name{bootcor.two}
\alias{bootcor.two}
\title{Bootstrapped Correlation Difference Test between Correlations from Two
Independent Samples}
\description{Derivation of bootstrap confidence intervals for the calculation of
correlation differences between the empirically observed correlations obtained
from two independent samples.}

\usage{
bootcor.two(x1,
            y1,
            x2,
            y2,
            k = 5000,
            alpha = .05,
            digit = 3,
            seed = 1234)}

\arguments{
  \item{x1}{A vector containing the values of the first variable being involved
  in the correlation in Sample 1.}
  \item{y1}{A vector containing the values of the second variable being involved
  in the correlation in Sample 1.}
  \item{x2}{A vector containing the values of the first variable being involved
  in the correlation in Sample 2.}
  \item{y2}{A vector containing the values of the second variable being involved
  in the correlation in Sample 2.}
  \item{k}{The number of bootstrap samples that should be drawn. The default is
  5000.}
  \item{alpha}{Likelihood of Type I error. The default is .05.}
  \item{digit}{Number of digits in the output. The default is 3.}
  \item{seed}{A random seed to make the results reproducible.}
}

\details{
Bivariate correlation analyses as well as correlation difference tests possess
very strict statistical requirements that are not necessarily fulfilled when
using the basic \code{diffcor.two()} function from this package (Wilcox, 2013
<doi:10.1016/C2010-0-67044-1>). For instance, if the assumption of a normal
distribution does not hold, the significance test can lead to false positive or
false negative conclusions. To address potential deviations from normal
distribution, the present function applies bootstrapping to the data. The output
provides a confidence interval for the difference between the empirically
observed correlation coefficients obtained from two independent samples, whereby
the interval is derived from bootstrapping.}

\value{
\item{r1}{The empircally observed correlation between x and y in Sample 1.}
\item{r2}{The empircally observed correlation between x and y in Sample 2.}
\item{M}{Mean of the confidence interval of the correlation difference between
the correlations from the two samples.}
\item{LL}{Lower limit of the confidence interval of the correlation difference
between the correlations from the two samples, given the entered Type I-level.}
\item{UL}{Upper limit of the confidence interval of the correlation difference
between the correlations from the two samples, given the entered Type I-level.}
}

\references{
Wilcox, R. (2013). Introduction to robust estimation and hypothesis testing.
Elsevier. https://doi.org/10.1016/C2010-0-67044-1
}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{
df1 <- data.frame(a = rnorm(1000),
                  b = rnorm(1000))

df2 <- data.frame(x = rnorm(600),
                  y = rnorm(600))

bootcor.two(x1 = df1$a,
            y1 = df1$b,
            x2 = df2$x,
            y2 = df2$y,
            k = 5000,
            alpha = .05,
            digit = 3,
            seed = 1234)}
