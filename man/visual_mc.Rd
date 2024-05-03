\name{visual_mc}
\alias{visual_mc}
\title{Visualization of the simulated parameters}
\description{To evaluate the quality of the Monte Carlo simulation beyond bias
and coverage parameters (Muthén & Muthén, 2002), it can be helpful to also
inspect the simulated parameters visually. To this end, visual_mc() can be used
to visualize the simulated parameters (including corresponding confidence
intervals) in relation to the targeted parameter.}
\usage{
visual_mc(rho,
          n,
          alpha = .05,
          n.intervals = 100,
          seed = 1234)}

\arguments{
  \item{rho}{Targeted correlation coefficient of the simulation.}
  \item{n}{An integer reflecting the sample size.}
  \item{alpha}{Type I error. Default is .05.}
  \item{n.intervals}{An integer reflecting the number of simulated parameters
  that should be visualized in the graphic. Default is 100.}
  \item{seed}{To make the results reproducible, a random seed is specified.}
}

\value{
A plot in which the targeted correlation coefficient is visualized with a dashed
red line and the simulated correlation coefficients are visualized by black
squares and confidence intervals (level depending on the specification made in
the argument alpha).}

\references{
Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to decide
on sample size and determine power. Structural Equation Modeling: A
Multidisciplinary Journal, 9(4), 599–620.
https://doi.org/10.1207/S15328007SEM0904_8}

\author{Christian Blötner
\email{c.bloetner@gmail.com}}
\keyword{utilities}
\keyword{graphs}

\examples{visual_mc(rho = .25,
                    n = 300,
                    alpha = .05,
                    n.intervals = 100,
                    seed = 1234)}
