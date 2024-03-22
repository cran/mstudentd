kldstudent <- function(nu1, Sigma1, nu2, Sigma2, eps = 1e-06) {
  #' Kullback-Leibler Divergence between Centered Multivariate \eqn{t} Distributions
  #'
  #' Computes the Kullback-Leibler divergence between two random vectors distributed
  #' according to multivariate \eqn{t} distributions (MTD) with zero location vector.
  #'
  #' @aliases kldstudent
  #'
  #' @usage kldstudent(nu1, Sigma1, nu2, Sigma2, eps = 1e-06)
  #' @param nu1 numeric. The degrees of freedom of the first distribution.
  #' @param Sigma1 symmetric, positive-definite matrix. The scatter matrix of the first distribution.
  #' @param nu2 numeric. The degrees of freedom of the second distribution.
  #' @param Sigma2 symmetric, positive-definite matrix. The scatter matrix of the second distribution.
  #' @param eps numeric. Precision for the computation of the partial derivative of the Lauricella \eqn{D}-hypergeometric function (see Details). Default: 1e-06.
  #' @return A numeric value: the Kullback-Leibler divergence between the two distributions,
  #' with two attributes \code{attr(, "epsilon")} (precision of the partial derivative of the Lauricella \eqn{D}-hypergeometric function,see Details)
  #' and \code{attr(, "k")} (number of iterations).
  #'
  #' @details Given \eqn{X_1}, a random vector of \eqn{R^p} distributed according to the centered MTD
  #' with parameters \eqn{(\nu_1, 0, \Sigma_1)}
  #' and \eqn{X_2}, a random vector of \eqn{R^p} distributed according to the MCD
  #' with parameters \eqn{(\nu_2, 0, \Sigma_2)}.
  #'
  #' Let \eqn{\lambda_1, \dots, \lambda_p} the eigenvalues of the square matrix \eqn{\Sigma_1 \Sigma_2^{-1}}
  #' sorted in increasing order: \deqn{\lambda_1 < \dots < \lambda_{p-1} < \lambda_p}
  #' The Kullback-Leibler divergence of \eqn{X_1} from \eqn{X_2} is given by:
  #' \deqn{
  #' \displaystyle{ D_{KL}(\mathbf{X}_1\|\mathbf{X}_2) = \ln\left(\frac{\Gamma\left(\frac{\nu_1+p}{2}\right) \Gamma\left(\frac{\nu_2}{2}\right) \nu_2^{\frac{p}{2}}}{\Gamma\left(\frac{\nu_2+p}{2}\right) \Gamma\left(\frac{\nu_1}{2}\right) \nu_1^{\frac{p}{2}}} \right) + \frac{\nu_2-\nu_1}{2} \left[\psi\left(\frac{\nu_1+p}{2} \right) - \psi\left(\frac{\nu_1}{2}\right)\right] - \frac{1}{2} \sum_{i=1}^p{\ln\lambda_i} - \frac{\nu_2+p}{2} \times D }
  #' }
  #' where \eqn{\psi} is the digamma function (see \link{Special})
  #' and \eqn{D} is given by:
  #' \itemize{
  #' \item If \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_1 > 1}}:
  #'
  #' \eqn{
  #' \displaystyle{ D = \prod_{i=1}^p{\left(\frac{\nu_2}{\nu_1} \frac{1}{\lambda_i}\right)^\frac{1}{2}} \frac{\partial}{\partial a}\left.\left\{ F_D^{(p)}\left(\frac{\nu_1+p}{2}, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p; a + \frac{\nu_1+p}{2}; 1-\frac{\nu_2}{\nu_1}\frac{1}{\lambda_1}, \dots, 1-\frac{\nu_2}{\nu_1}\frac{1}{\lambda_p}\right) \right\}\right|_{a=0} }
  #' }
  #' \item If \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_p < 1}}:
  #'
  #' \eqn{
  #' \displaystyle{ D = \frac{\partial}{\partial a}\left.\left\{ F_D^{(p)}\left(a, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p; a + \frac{\nu_1+p}{2}; 1-\frac{\nu_1}{\nu_2}\lambda_1, \dots, 1-\frac{\nu_1}{\nu_2}\lambda_p\right) \right\}\right|_{a=0} }
  #' }
  #' \item If \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_1 < 1}} and \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_p > 1}}:
  #'
  #' \eqn{
  #' \displaystyle{ D = -\ln\left(\frac{\nu_1}{\nu_2}\lambda_p\right) + \frac{\partial}{\partial a}\left.\left\{F_D^{(p)}\left(a, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}, a+\frac{\nu_1}{2}}_p; a+\frac{\nu_1+p}{2}; 1-\frac{\lambda_1}{\lambda_p}, \dots, 1-\frac{\lambda_{p-1}}{\lambda_p}, 1-\frac{\nu_2}{\nu_1}\frac{1}{\lambda_p}\right)\right\}\right|_{a=0} }
  #' }
  #' }
  #'
  #' \eqn{F_D^{(p)}} is the Lauricella \eqn{D}-hypergeometric function defined for \eqn{p} variables:
  #' \deqn{ \displaystyle{ F_D^{(p)}\left(a; b_1, ..., b_p; g; x_1, ..., x_p\right) = \sum\limits_{m_1 \geq 0} ... \sum\limits_{m_p \geq 0}{ \frac{ (a)_{m_1+...+m_p}(b_1)_{m_1} ... (b_p)_{m_p} }{ (g)_{m_1+...+m_p} } \frac{x_1^{m_1}}{m_1!} ... \frac{x_p^{m_p}}{m_p!} } } }
  #'
  #' The computation of the partial derivative uses the \code{\link{pochhammer}} function.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
  #' \doi{10.1109/LSP.2023.3324594}
  #'
  #' @examples
  #' nu1 <- 2
  #' Sigma1 <- matrix(c(2, 1.2, 0.4, 1.2, 2, 0.6, 0.4, 0.6, 2), nrow = 3)
  #' nu2 <- 4
  #' Sigma2 <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.4, 0.1, 0.4, 1), nrow = 3)
  #'
  #' kldstudent(nu1, Sigma1, nu2, Sigma2)
  #' kldstudent(nu2, Sigma2, nu1, Sigma1)
  #'
  #' @importFrom utils combn
  #' @export

  # Sigma1 and Sigma2 must be matrices
  if (is.numeric(Sigma1) & !is.matrix(Sigma1))
    Sigma1 <- matrix(Sigma1)
  if (is.numeric(Sigma2) & !is.matrix(Sigma2))
    Sigma2 <- matrix(Sigma2)

  # Number of variables
  p <- nrow(Sigma1)

  # Sigma1 and Sigma2 must be square matrices with the same size
  if (ncol(Sigma1) != p | nrow(Sigma2) != p | ncol(Sigma2) != p)
    stop("Sigma1 et Sigma2 must be square matrices with rank p.")

  # IS Sigma1 symmetric, positive-definite?
  if (!isSymmetric(Sigma1))
    stop("Sigma1 must be a symmetric, positive-definite matrix.")
  lambda1 <- eigen(Sigma1, only.values = TRUE)$values
  if (any(lambda1 < .Machine$double.eps))
    stop("Sigma1 must be a symmetric, positive-definite matrix.")

  # IS Sigma2 symmetric, positive-definite?
  if (!isSymmetric(Sigma2))
    stop("Sigma2 must be a symmetric, positive-definite matrix.")
  lambda2 <- eigen(Sigma2, only.values = TRUE)$values
  if (any(lambda2 < .Machine$double.eps))
    stop("Sigma2 must be a symmetric, positive-definite matrix.")

  # Eigenvalues of Sigma1 %*% inv(Sigma2)
  lambda <- sort(eigen(Sigma1 %*% solve(Sigma2), only.values = TRUE)$values, decreasing = FALSE)
  lambdanu <- lambda*nu1/nu2
  prodlambdanu <- prod(lambdanu)

  k <- 5

  # M: data.frame of the indices for the nested sums
  # (i.e. all arrangements of n elements from {0:k})
  M <- expand.grid(rep(list(0:k), p))
  M <- M[-1, , drop = FALSE]
  Msum <- apply(M, 1, sum)
  kstep <- 5

  if (lambdanu[p] < 1) {  # lambda[1] < ... < lambda[p] < 1

    # The first 5 elements of the sum
    d <- 0
    for (i in 1:length(Msum)) {
      commun <- prod(
        sapply(1:p, function(j) {
          pochhammer(0.5, M[i, j])*(1 - lambdanu[j])^M[i, j]/factorial(M[i, j])
        })
      )
      d <- d + commun * pochhammer(1, Msum[i]) / ( pochhammer((nu1 + p)/2, Msum[i]) * Msum[i] )
    }

    # Next elements of the sum, until the expected precision
    k1 <- 1:k
    derive <- 0
    while (abs(d) > eps/10 & !is.nan(d)) {
      epsret <- signif(abs(d), 1)*10
      k <- k1[length(k1)]
      k1 <- k + (1:kstep)
      derive <- derive + d

      # M: data.frame of the indices for the nested sums
      M <- expand.grid(rep(list(k1), p))
      if (p > 1) {
        for (i in 1:(p-1)) {
          indsupp <- combn(p, i)
          for (j in 1:ncol(indsupp)) {
            jsupp <- indsupp[, j]
            Mlist <- vector("list", p)
            for (l in jsupp) Mlist[[l]] <- k1
            for (l in (1:p)[-jsupp]) Mlist[[l]] <- 0:k
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }

      Msum <- apply(M, 1, sum)

      d <- 0
      for (i in 1:length(Msum)) {
        commun <- prod(
          sapply(1:p, function(j) {
            pochhammer(0.5, M[i, j])*(1 - lambdanu[j])^M[i, j]/factorial(M[i, j])
          })
        )
        d <- d + commun * pochhammer(1, Msum[i]) / ( pochhammer((nu1 + p)/2, Msum[i]) * Msum[i] )
      }

    }

    # Next elements of the sum, with step=1, while not NaN
    if (is.nan(d)) {
      k1 <- k
      d <- 0
      while (!is.nan(d)) {
        if (d > 0)
          epsret <- signif(abs(d), 1)*10
        k <- k1
        k1 <- k + 1
        derive <- derive + d

        # M: data.frame of the indices for the nested sums
        M <- expand.grid(rep(list(k1), p))
        if (p > 1) {
          for (i in 1:(p-1)) {
            indsupp <- combn(p, i)
            for (j in 1:ncol(indsupp)) {
              jsupp <- indsupp[, j]
              Mlist <- vector("list", p)
              for (l in jsupp) Mlist[[l]] <- k1
              for (l in (1:p)[-jsupp]) Mlist[[l]] <- 0:k
              M <- rbind(M, expand.grid(Mlist))
            }
          }
        }

        Msum <- apply(M, 1, sum)

        d <- 0
        for (i in 1:length(Msum)) {
          commun <- prod(
            sapply(1:p, function(j) {
              pochhammer(0.5, M[i, j])*(1 - lambdanu[j])^M[i, j]/factorial(M[i, j])
            })
          )
          d <- d + commun * pochhammer(1, Msum[i]) / ( pochhammer((nu1 + p)/2, Msum[i]) * Msum[i] )
        }

      }
    }

  } else if (lambdanu[1] > 1) { # 1 < lambda[1] < ... < lambda[p]

    # The first 5 elements of the sum
    d <- 0
    for (i in 1:length(Msum)) {
      commun <- prod(
        sapply(1:p, function(j) {
          pochhammer(0.5, M[i, j])*(1 - 1/lambdanu[j])^M[i, j]/factorial(M[i, j])
        })
      )
      A <- sum(1/(0:(Msum[i]-1) + (nu1+p)/2))
      d <- d - commun * A # / pochhammer((1 + p)/2, Msum[i])
    }

    # Next elements of the sum, until the expected precision
    k1 <- 1:k
    derive <- 0
    # vd <- vderive <- numeric()
    while (abs(d) > eps/10 & !is.nan(d)) {
      epsret <- signif(abs(d), 1)*10
      k <- k1[length(k1)]
      k1 <- k + (1:kstep)
      derive <- derive + d
      # vd <- c(vd, d); vderive <- c(vderive, derive)

      # M: data.frame of the indices for the nested sums
      M <- expand.grid(rep(list(k1), p))
      if (p > 1) {
        for (i in 1:(p-1)) {
          indsupp <- combn(p, i)
          for (j in 1:ncol(indsupp)) {
            jsupp <- indsupp[, j]
            Mlist <- vector("list", p)
            for (l in jsupp) Mlist[[l]] <- k1
            for (l in (1:p)[-jsupp]) Mlist[[l]] <- 0:k
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }

      Msum <- apply(M, 1, sum)

      d <- 0
      for (i in 1:length(Msum)) {
        commun <- prod(
          sapply(1:p, function(j) {
            pochhammer(0.5, M[i, j])*(1 - 1/lambdanu[j])^M[i, j]/factorial(M[i, j])
          })
        )
        A <- sum(1/(0:(Msum[i]-1) + (nu1+p)/2))
        d <- d - commun * A # / pochhammer((1 + p)/2, Msum[i])
      }

    }

    # Next elements of the sum, with step=1, while not NaN
    if (is.nan(d)) {
      k1 <- k
      d <- 0
      while (!is.nan(d)) {
        if (d > 0)
          epsret <- signif(abs(d), 1)*10
        k <- k1
        k1 <- k + 1
        derive <- derive + d

        # M: data.frame of the indices for the nested sums
        M <- expand.grid(rep(list(k1), p))
        if (p > 1) {
          for (i in 1:(p-1)) {
            indsupp <- combn(p, i)
            for (j in 1:ncol(indsupp)) {
              jsupp <- indsupp[, j]
              Mlist <- vector("list", p)
              for (l in jsupp) Mlist[[l]] <- k1
              for (l in (1:p)[-jsupp]) Mlist[[l]] <- 0:k
              M <- rbind(M, expand.grid(Mlist))
            }
          }
        }

        Msum <- apply(M, 1, sum)

        d <- 0
        for (i in 1:length(Msum)) {
          commun <- prod(
            sapply(1:p, function(j) {
              pochhammer(0.5, M[i, j])*(1 - 1/lambdanu[j])^M[i, j]/factorial(M[i, j])
            })
          )
          A <- sum(1/(0:(Msum[i]-1) + (nu1+p)/2))
          d <- d - commun * A # / pochhammer((1 + p)/2, Msum[i])
        }

      }
    }
    derive <- prod(1/sqrt(lambdanu)) * derive

  } else { # lambda[1] < ... < 1 < ... < lambda[p]

    # The first 5 elements of the sum
    d <- 0
    for (i in 1:length(Msum)) {
      commun <- prod(
        sapply(1:(p-1), function(j) {
          pochhammer(0.5, M[i, j])*(1 - lambda[j]/lambda[p])^M[i, j]/factorial(M[i, j])
        })
      )
      commun <- commun*(1 - 1/lambdanu[p])^M[i, p]/factorial(M[i, p])
      d <- d + commun * pochhammer(0.5, M[i, p])*pochhammer(1, Msum[i]) / ( pochhammer((nu1 + p)/2, Msum[i]) * Msum[i] )
    }

    # Next elements of the sum, until the expected precision
    k1 <- 1:k
    derive <- 0
    while (abs(d) > eps/10 & !is.nan(d)) {
      epsret <- signif(abs(d), 1)*10
      k <- k1[length(k1)]
      k1 <- k + (1:kstep)
      derive <- derive + d

      # M: data.frame of the indices for the nested sums
      M <- expand.grid(rep(list(k1), p))
      if (p > 1) {
        for (i in 1:(p-1)) {
          indsupp <- combn(p, i)
          for (j in 1:ncol(indsupp)) {
            jsupp <- indsupp[, j]
            Mlist <- vector("list", p)
            for (l in jsupp) Mlist[[l]] <- k1
            for (l in (1:p)[-jsupp]) Mlist[[l]] <- 0:k
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }

      Msum <- apply(M, 1, sum)

      d <- 0
      for (i in 1:length(Msum)) {
        commun <- prod(
          sapply(1:(p-1), function(j) {
            pochhammer(0.5, M[i, j])*(1 - lambda[j]/lambda[p])^M[i, j]/factorial(M[i, j])
          })
        )
        commun <- commun*(1 - 1/lambdanu[p])^M[i, p]/factorial(M[i, p])
        d <- d + commun * pochhammer(0.5, M[i, p])*pochhammer(1, Msum[i]) / ( pochhammer((nu1 + p)/2, Msum[i]) * Msum[i] )
      }

    }

    # Next elements of the sum, with step=1, while not NaN
    if (is.nan(d)) {
      k1 <- k
      d <- 0
      while (!is.nan(d)) {
        if (d > 0)
          epsret <- signif(abs(d), 1)*10
        k <- k1
        k1 <- k + 1
        derive <- derive + d

        # # M: data.frame of the indices for the nested sums
        # M <- as.data.frame(matrix(nrow = 0, ncol = p))
        # if (p > 1) {
        #   for (i in 1:(p-1)) {
        #     Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
        #     M <- rbind( M, expand.grid(Mlist) )
        #     for (j in 1:(p-1)) {
        #       Mlist <- Mlist[c(p, 1:(p-1))]
        #       M <- rbind(M, expand.grid(Mlist))
        #     }
        #   }
        # }
        # M <- rbind( M, rep(k1, p) )

        # M: data.frame of the indices for the nested sums
        M <- expand.grid(rep(list(k1), p))
        if (p > 1) {
          for (i in 1:(p-1)) {
            indsupp <- combn(p, i)
            for (j in 1:ncol(indsupp)) {
              jsupp <- indsupp[, j]
              Mlist <- vector("list", p)
              for (l in jsupp) Mlist[[l]] <- k1
              for (l in (1:p)[-jsupp]) Mlist[[l]] <- 0:k
              M <- rbind(M, expand.grid(Mlist))
            }
          }
        }

        Msum <- apply(M, 1, sum)

        d <- 0
        for (i in 1:length(Msum)) {
          commun <- prod(
            sapply(1:(p-1), function(j) {
              pochhammer(0.5, M[i, j])*(1 - lambda[j]/lambda[p])^M[i, j]/factorial(M[i, j])
            })
          )
          commun <- commun*(1 - 1/lambdanu[p])^M[i, p]/factorial(M[i, p])
          d <- d + commun * pochhammer(0.5, M[i, p])*pochhammer(1, Msum[i]) / ( pochhammer((nu1 + p)/2, Msum[i]) * Msum[i] )
        }

      }
    }
    derive <- -log(lambdanu[p]) + derive

  }

  result <- log(gamma((nu1+p)/2)*gamma(nu2/2)*nu2^(p/2) / (gamma((nu2+p)/2)*gamma(nu1/2)*nu1^(p/2)))
  result <- result + (nu2-nu1)/2 * (digamma((nu1+p)/2) - digamma(nu1/2)) - 0.5 * sum(log(lambda))
  result <- result - (nu2 + p)/2 * derive
  result <- as.numeric(result)

  if (is.nan(d)) {
    epsret <- signif(epsret, 1)
    warning("Cannot reach the precision ", eps, " due to NaN\n",
            "Number of iteration: ", k, "\n",
            "Precision reached:", epsret)
    attr(result, "epsilon") <- epsret
  } else {
    attr(result, "epsilon") <- eps
  }
  attr(result, "k") <- k

  return(result)
}
