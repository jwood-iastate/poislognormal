#' Family Function for Poisson-Lognormal Distribution
#'
#' @description
#' Creates a `family` object for the Poisson-Lognormal distribution to be used in
#' model fitting functions like `glm()`. The family uses a log-link and requires a
#' fixed dispersion parameter `sigma`.
#'
#' @details
#' The Poisson-Lognormal distribution is a compound distribution that can model
#' over-dispersed count data. The mean-variance relationship is:
#' \deqn{Var(Y) = E[Y] + (e^{\sigma^2}-1) E[Y]^2}
#' where \eqn{\sigma} is the dispersion parameter of the underlying Lognormal
#' distribution. As \eqn{\sigma \to 0}, the distribution converges to a
#' standard Poisson, and the variance converges to the mean.
#'
#' The family uses a log-link function, where the linear predictor \eqn{\eta}
#' relates to the mean \eqn{\mu + exp(\sigma^2/2) = E[Y]} as \eqn{\eta = \log(\mu)}.
#'
#' The underlying \eqn{\mu_{LN}} parameter of the Lognormal component required by
#' the internal density function (`dpln_rcpp`) is related to the mean of the
#' Poisson-Lognormal distribution \eqn{\mu} by:
#' \deqn{\mu_{LN} = \log(\mu) - \sigma^2/2}
#'
#' Since `glm` can only estimate one dispersion parameter (typically \eqn{\phi} in
#' \eqn{Var(Y) = \phi V(\mu)}), the Lognormal dispersion parameter `sigma` must be
#' provided as a fixed value to this family function. To estimate `sigma`, one
#' might fit the model over a range of `sigma` values and select the one that
#' maximizes the log-likelihood or minimizes AIC.
#'
#' @param sigma A numeric scalar specifying the fixed dispersion parameter
#'   (\eqn{\sigma > 0}) of the underlying Lognormal distribution.
#'
#' @return An object of class `family` for use with functions like `glm`.
#'
#' @export
#' @importFrom stats dpois make.link
poisLogn <- function(sigma = 1) {
  
  if (!is.numeric(sigma) || sigma <= 0) {
    stop("Dispersion parameter 'sigma' must be a positive number.")
  }
  
  # Store sigma^2 for repeated use
  sigma2 <- sigma^2
  
  # Environment to store sigma, so it's accessible to the functions below
  env <- new.env(parent = .GlobalEnv)
  env$sigma <- sigma
  env$sigma2 <- sigma2
  
  # Variance function
  # Var(Y) = E[Y] + (exp(sigma^2) - 1) * E[Y]^2
  variance <- function(mu) {
    mu * exp(sigma2/2) + (exp(sigma2) - 1) * mu^2 * exp(sigma2/2)^2
  }
  
  # Link functions (standard log-link)
  link_info <- stats::make.link("log")
  linkfun <- link_info$linkfun
  linkinv <- link_info$linkinv
  mu.eta <- link_info$mu.eta
  
  # Deviance residuals
  # Deviance = 2 * (logLik_saturated - logLik_model)
  # Saturated log-likelihood for a count model is based on mu=y
  # This need updated to use the correct functions
  dev.resids <- function(y, mu, wt) {
    # The mu parameter for dpln_rcpp is the log-mean of the Lognormal,
    # which is log(E[Y]) - sigma^2 / 2.
    mu_ln <- log(mu) - sigma2 / 2
    
    # Calculate the log-likelihood for the fitted model
    ll_model <- dpln_rcpp(y, mu = mu_ln, sigma = sigma, log_p = TRUE)
    
    # Calculate the log-likelihood for the saturated model (where mu = y)
    # This uses the Poisson(y) as the saturated model, a standard choice.
    ll_saturated <- stats::dpois(y, lambda = y, log = TRUE)
    
    # Deviance is 2 * (ll_saturated - ll_model)
    # The sign of the result is taken from (y - mu)
    dev <- 2 * wt * (ll_saturated - ll_model)
    ifelse(y > mu, sqrt(dev), -sqrt(dev))
  }
  
  # AIC function
  # AIC = -2 * logLik + 2 * k, where k is number of parameters.
  # The family$aic function in R should return AIC for a model with k parameters
  # (including one for the dispersion if it is estimated).
  aic <- function(y, n, mu, wt, dev) {
    # The log-likelihood can be calculated from the deviance.
    # dev = 2 * (ll_saturated - ll_model) => ll_model = ll_saturated - dev/2
    mu_ln <- log(mu) - sigma2 / 2
    ll_model <- sum(dpln_rcpp(y, mu = mu_ln, sigma = sigma, log_p = TRUE) * wt)
    
    # Add 1 to the rank (p) for the dispersion parameter sigma,
    # which is standard practice for AIC comparison with other over-dispersed
    # models like quasipoisson or negative.binomial.
    p <- attr(dev, "rank") + 1
    -2 * ll_model + 2 * p
  }
  
  # Initialize expression
  # Use the same initialization as quasipoisson
  initialize <- expression({
    if (any(y < 0))
      warning("negative values not allowed for the Poisson-Lognormal family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  
  # Return the family object
  structure(
    list(
      family = "poisLogn",
      link = "log",
      linkfun = linkfun,
      linkinv = linkinv,
      variance = variance,
      dev.resids = dev.resids,
      aic = aic,
      mu.eta = mu.eta,
      initialize = initialize,
      validmu = function(mu) all(is.finite(mu)) && all(mu > 0),
      valideta = function(eta) all(is.finite(eta)),
      # Pass the environment so functions can access sigma
      environment = env
    ),
    class = "family"
  )
}
