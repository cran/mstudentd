contourmtd <- function(nu, mu, Sigma,
                        xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                        ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]),
                        zlim = NULL, npt = 30, nx = npt, ny = npt,
                        main = "Multivariate t density",
                        sub = NULL, nlevels = 10,
                        levels = pretty(zlim, nlevels),
                        tol = 1e-6, ...) {
  #' Contour Plot of the Bivariate \eqn{t} Density
  #'
  #' Draws the contour plot of the probability density of the multivariate \eqn{t} distribution with 2 variables
  #' with location parameter \code{mu} and scatter matrix \code{Sigma}.
  #'
  #' @aliases contourmtd
  #'
  #' @usage contourmtd(nu, mu, Sigma,
  #'                    xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
  #'                    ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]),
  #'                    zlim = NULL, npt = 30, nx = npt, ny = npt,
  #'                    main = "Multivariate t density",
  #'                    sub = NULL, nlevels = 10,
  #'                    levels = pretty(zlim, nlevels), tol = 1e-6, ...)
  #' @param nu numeric. The degrees of freedom.
  #' @param mu length 2 numeric vector.
  #' @param Sigma symmetric, positive-definite square matrix of order 2. The scatter matrix.
  #' @param main,sub main and sub title, as for \code{\link{title}}.
  #' @param xlim,ylim x-and y- limits.
  #' @param zlim z- limits. If NULL, it is the range of the values of the density on the x and y values within `xlim` and `ylim`.
  #' @param npt number of points for the discretisation.
  #' @param nx,ny number of points for the discretisation among the x- and y- axes.
  #' @param nlevels,levels arguments to be passed to the \code{\link{contour}} function.
  #' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma, for the estimation of the density. see \code{\link{dmtd}}.
  #' @param ... additional arguments to \code{\link{plot.window}}, \code{\link{title}}, \code{\link{Axis}} and \code{\link{box}}, typically \link{graphical parameters} such as \code{cex.axis}.
  #' @return Returns invisibly the probability density function.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #'
  #' @references S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
  #'
  #' @seealso \code{\link{dmtd}}: probability density of a multivariate \eqn{t} density
  #'
  #' \code{\link{plotmtd}}: 3D plot of a bivariate \eqn{t} density.
  #'
  #' @examples
  #' nu <- 1
  #' mu <- c(1, 4)
  #' Sigma <- matrix(c(0.8, 0.2, 0.2, 0.2), nrow = 2)
  #' contourmtd(nu, mu, Sigma)
  #'
  #' @importFrom graphics contour
  #' @importFrom graphics par
  #' @export

  if (length(mu)!=2 | nrow(Sigma)!=2 | ncol(Sigma)!=2)
    stop(paste("contourmtd only allows plotting a t density with 2 variables.",
               "mu must be a length 2 numeric vector and Sigma must be a 2*2 square matrix.", sep = "\n"))

  # Estimation of the density
  f <- function(x) dmtd(x, nu = nu, mu = mu, Sigma = Sigma, tol = tol)
  ff <- function(x, y) sapply(1:length(x), function(i) as.numeric(f(c(x[i], y[i]))))

  x <- seq(xlim[1], xlim[2], length = nx)
  y <- seq(ylim[1], ylim[2], length = ny)
  z <- outer(x, y, ff)
  if (is.null(zlim)) zlim <- range(z)

  # Plot
  contour(x, y, z, nlevels = nlevels, levels = levels, labels = NULL,
          xlim = xlim, ylim = ylim, zlim = zlim,  labcex = 0.6,
          drawlabels = TRUE, method = "flattest", vfont = NULL, axes = TRUE,
          frame.plot = TRUE, col = par("fg"), lty = par("lty"),
          lwd = par("lwd"), add = FALSE, main = main, ...)

  return(invisible(f))
}
