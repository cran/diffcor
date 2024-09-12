\name{bootcor.dep}
\alias{bootcor.dep}
\title{Bootstrapped Correlation Difference Test for Dependent Correlations}
\description{Derivation of bootstrap confidence intervals for the calculation of
correlation differences for dependent correlations.}

\usage{
bootcor.dep(target,
            x1,
            x2,
            k = 5000,
            alpha = .05,
            digit = 3,
            seed = 1234)}

\arguments{
  \item{target}{A vector containing the values for the target variable for
  which the correlations of the two competing variables x1 and x2 should be
  compared.}
  \item{x1}{A vector containing the values of the first variable being
  correlated with the target variable.}
  \item{x2}{A vector containing the values of the second variable being
  correlated with the target variable.}
  \item{k}{The number of bootstrap samples that should be drawn. The default is
  5000.}
  \item{alpha}{Likelihood of Type I error. The default is .05.}
  \item{digit}{Number of digits in the output. The default is 3.}
  \item{seed}{A random seed to make the results reproducible.}
}

\details{
Bivariate correlation analyses as well as correlation difference tests possess
very strict statistical requirements that are not necessarily fulfilled when
using the basic \code{diffcor.dep()} function from this package (Wilcox, 2013
<doi:10.1016/C2010-0-67044-1>). For instance, if the assumption of a normal
distribution does not hold, the significance test can lead to false positive or
false negative conclusions. To address potential deviations from normal
distribution, the present function applies bootstrapping to the data. The output
provides a confidence interval for the difference between the empirically
observed correlations of two competing variables with a target variable,
whereby the interval is derived from bootstrapping..}

\value{
\item{r_target_1}{The empircally observed correlation between the first
variable and the target variable.}
\item{r_target_2}{The empircally observed correlation between the second
variable and the target variable.}
\item{M}{Mean of the confidence interval of the correlation difference between
\code{r_target_1} and \code{r_target_2}.}
\item{LL}{Lower limit of the confidence interval of the correlation difference
between \code{r_target_1} and \code{r_target_2}, given the entered Type I-level.
}
\item{UL}{Upper limit of the confidence interval of the correlation difference
between \code{r_target_1} and \code{r_target_2}, given the entered Type I-level.
}
}

\references{
Wilcox, R. (2013). Introduction to robust estimation and hypothesis testing.
Elsevier. https://doi.org/10.1016/C2010-0-67044-1
}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{&htest}

\examples{

df <- data.frame(target = rnorm(1000),
                 var1 = rnorm(1000),
                 var2 = rnorm(1000))

bootcor.dep(target = df$target,
            x1 = df$var1,
            x2 = df$var2,
            k = 5000,
            alpha = .05,
            digit = 3,
            seed = 1234)}
