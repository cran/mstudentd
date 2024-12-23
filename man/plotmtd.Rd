% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotmtd.R
\name{plotmtd}
\alias{plotmtd}
\title{Plot of the Bivariate \eqn{t} Density}
\usage{
plotmtd(nu, mu, Sigma, xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
                xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
                zlab = "f(x,y)", col = "gray", tol = 1e-6, ...)
}
\arguments{
\item{nu}{numeric. The degrees of freedom.}

\item{mu}{length 2 numeric vector. The mean vector.}

\item{Sigma}{symmetric, positive-definite square matrix of order 2. The correlation matrix.}

\item{xlim, ylim}{x-and y- limits.}

\item{n}{A one or two element vector giving the number of steps in the x and y grid, passed to \code{\link[rgl]{plot3d.function}}.}

\item{xvals, yvals}{The values at which to evaluate \code{x} and \code{y}. If used, \code{xlim} and/or \code{ylim} are ignored.}

\item{xlab, ylab, zlab}{The axis labels.}

\item{col}{The color to use for the plot. See \code{\link[rgl]{plot3d.function}}.}

\item{tol}{tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma, for the estimation of the density. see \code{\link{dmtd}}.}

\item{...}{Additional arguments to pass to \code{\link[rgl]{plot3d.function}}.}
}
\value{
Returns invisibly the probability density function.
}
\description{
Plots the probability density of the multivariate \eqn{t} distribution with 2 variables
with location parameter \code{mu} and scatter matrix \code{Sigma}.
}
\examples{
nu <- 1
mu <- c(1, 4)
Sigma <- matrix(c(0.8, 0.2, 0.2, 0.2), nrow = 2)
plotmtd(nu, mu, Sigma)

}
\references{
S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
}
\seealso{
\code{\link{dmtd}}: probability density of a multivariate \eqn{t} density

\code{\link{contourmtd}}: contour plot of a bivariate \eqn{t} density.

\code{\link[rgl]{plot3d.function}}: plot a function of two variables.
}
\author{
Pierre Santagostini, Nizar Bouhlel
}
