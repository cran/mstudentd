% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mstudentd-package.R
\docType{package}
\name{mstudentd-package}
\alias{mstudentd-package}
\alias{mstudentd}
\title{Tools for Multivariate \eqn{t} Distributions}
\description{
This package provides tools for multivariate \eqn{t} distributions (MTD):
\itemize{
\item Calculation of distances/divergences between MTD:
\itemize{
\item Renyi divergence, Bhattacharyya distance, Hellinger distance: \code{\link{diststudent}}
\item Kullback-Leibler divergence: \code{\link{kldstudent}}
}
\item Tools for MTD:
\itemize{
\item Probability density: \code{\link{dmtd}}
\item Simulation from a MTD: \code{\link{rmtd}}
\item Plot of the density of a MTD with 2 variables: \code{\link{plotmtd}}, \code{\link{contourmtd}}
}
}
}
\references{
S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.

N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
\doi{10.1109/LSP.2023.3324594}
#' @keywords internal
}
\seealso{
Useful links:
\itemize{
  \item \url{https://forgemia.inra.fr/imhorphen/mstudentd}
  \item Report bugs at \url{https://forgemia.inra.fr/imhorphen/mstudentd/-/issues}
}

}
\author{
Pierre Santagostini \href{mailto:pierre.santagostini@agrocampus-ouest.fr}{pierre.santagostini@agrocampus-ouest.fr},
Nizar Bouhlel \href{mailto:nizar.bouhlel@agrocampus-ouest.fr}{nizar.bouhlel@agrocampus-ouest.fr}
}
