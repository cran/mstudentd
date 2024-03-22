#' Tools for Multivariate \eqn{t} Distributions
#'
#' This package provides tools for multivariate \eqn{t} distributions (MTD):
#' \itemize{
#' \item Calculation of distances/divergences between MTD:
#' \itemize{
#' \item Renyi divergence, Bhattacharyya distance, Hellinger distance: \code{\link{diststudent}}
#' \item Kullback-Leibler divergence: \code{\link{kldstudent}}
#' }
#' \item Tools for MTD:
#' \itemize{
#' \item Probability density: \code{\link{dmtd}}
#' \item Simulation from a MTD: \code{\link{rmtd}}
#' \item Plot of the density of a MTD with 2 variables: \code{\link{plotmtd}}, \code{\link{contourmtd}}
#' }
#' }
#'
#' @name mstudentd-package
#' @aliases mstudentd-package mstudentd
#' @docType package
#' @author Pierre Santagostini <pierre.santagostini@agrocampus-ouest.fr>,
#' Nizar Bouhlel <nizar.bouhlel@agrocampus-ouest.fr>
#' @references
#' S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
#'
#' N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
#' \doi{10.1109/LSP.2023.3324594}
#' #' @keywords internal
"_PACKAGE"

NULL
