#' Poisson-Lognormal Distribution
#'
#' @description
#' These functions provide density, distribution function, quantile function, and
#' random number generation for the Poisson-Lognormal (PLogN) distribution.
#'
#' @details
#' The density function is calculated using the C++ function `dpln_rcpp`, which
#' uses 100-point Gauss-Hermite quadrature for accurate numerical integration.
#'
#' The user-facing `mu` parameter represents the mean of the Poisson-Lognormal
#' distribution, \eqn{E[Y] = \mu}. This is internally converted to the mean of
#' the underlying Lognormal distribution on the log scale (\eqn{\mu_{LN}}) via
#' \eqn{\mu_{LN} = \log(\mu) - \sigma^2/2}.
#'
#' @param x A vector of quantiles (non-negative integers).
#' @param q A vector of quantiles.
#' @param p A vector of probabilities.
#' @param n The number of random numbers to generate.
#' @param mu A vector of positive means of the distribution.
#' @param sigma A vector of positive dispersion parameters.
#' @param log,log.p Logical; if `TRUE`, probabilities are given/returned as log(p).
#' @param lower.tail Logical; if `TRUE` (the default), probabilities are
#'   \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#'
#' @return
#' `dpLnorm` gives the density (PMF), `ppLnorm` gives the distribution
#' function (CDF), `qpLnorm` gives the quantile function, and `rpLnorm`
#' generates random deviates.
#'
#' @name PoissonLognormal
#' @aliases dpLnorm ppLnorm qpLnorm rpLnorm
#'
#' @examples
#' dpLnorm(c(0, 1, 2), mu = 1.5, sigma = 0.8)
#' ppLnorm(5, mu = 1.5, sigma = 0.8)
#' qpLnorm(0.95, mu = 1.5, sigma = 0.8)
#' rpLnorm(10, mu = 1.5, sigma = 0.8)
NULL

#' @rdname PoissonLognormal
#' @export
dpLnorm <- function(x, mu = 1, sigma = 1, log = FALSE) {
  if (any(mu <= 0) || any(sigma <= 0)) {
    stop("Parameters `mu` and `sigma` must be positive.")
  }
  # Convert the mean of the PLN, mu, to the log-mean of the underlying Lognormal
  mu_ln <- log(mu) - sigma^2 / 2
  p <- dpln_rcpp(x, mu = mu_ln, sigma = sigma, log_p = log)
  return(p)
}

#' @rdname PoissonLognormal
#' @export
ppLnorm <- Vectorize(function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  if (q < 0) return(0)
  y <- seq(0, floor(q), 1)
  probs <- dpLnorm(y, mu = mu, sigma = sigma)
  p <- sum(probs)
  
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    return(log(p))
  } else {
    return(p)
  }
})

#' @rdname PoissonLognormal
#' @export
qpLnorm <- Vectorize(function(p, mu = 1, sigma = 1) {
  if (p < 0 || p > 1) stop("`p` must be between 0 and 1.")
  if (p == 0) return(0)
  if (p == 1) return(Inf)
  
  y <- 0
  p_value <- ppLnorm(y, mu = mu, sigma = sigma)
  while (p_value < p) {
    y <- y + 1
    # Safety break for extreme quantiles
    if (y > 50000) {
      warning("Quantile is very large, returning Inf.")
      return(Inf)
    }
    p_value <- ppLnorm(y, mu = mu, sigma = sigma)
  }
  return(y)
})

#' @rdname PoissonLognormal
#' @export
rpLnorm <- function(n, mu = 1, sigma = 1) {
  if (any(mu <= 0) || any(sigma <= 0)) {
    stop("Parameters `mu` and `sigma` must be positive.")
  }
  u <- runif(n)
  y <- qpLnorm(u, mu = mu, sigma = sigma)
  return(y)
}
