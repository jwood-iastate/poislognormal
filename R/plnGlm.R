#' Fit a Poisson-Lognormal Regression Model
#'
#' @description
#' Fits a Poisson-Lognormal regression model using direct maximum likelihood
#' estimation via the `maxLik` package.
#'
#' @details
#' This function directly optimizes the Poisson-Lognormal log-likelihood function
#' to estimate both the regression coefficients (betas) and the dispersion
#' parameter (`sigma`) simultaneously. This provides a direct estimate and
#' standard error for `sigma`.
#'
#' The `sigma` parameter is estimated on the log scale (`log(sigma)`) to ensure
#' it remains positive during optimization. The returned estimate and standard
#' error are transformed back to the original scale.
#'
#' @param formula A formula specifying the model.
#' @param data A data frame containing the variables in the model.
#' @param offset An optional character string specifying an offset variable in `data`.
#' @param weights An optional character string specifying a weights variable in `data`.
#' @param method A character string specifying the optimization method for `maxLik`.
#'   See `?maxLik::maxLik` for options.
#' @param verbose A logical value. If `TRUE`, detailed optimization output is printed.
#' @param ... Additional arguments passed to `maxLik::maxLik` (e.g., `iterlim`).
#'
#' @return An object of class `poislogn_fit`.
#'
#' @export
#' @importFrom maxLik maxLik
#' @importFrom stats .getXlevels model.frame model.matrix model.offset model.response model.weights poisson glm.fit pnorm printCoefmat delete.response terms
#'
#' @examples
#' ## --- Example with washington_roads dataset ---
#' data(washington_roads)
#'
#' # Create a binary indicator for high AADT
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' # Fit a model for animal-related crashes, using road length as an exposure offset
#' # The offset term must be on the log scale.
#' pln_model_roads <- pln.glm(
#'   Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k + offset(lnlength),
#'   data = washington_roads
#' )
#'
#' # View the model summary
#' summary(pln_model_roads)
#' 
#' predict(pln_model_roads, type = "response", newdata = washington_roads)
#' 
pln.glm <- function(formula, data, offset = NULL, weights = NULL,
                    method = "BFGS", verbose = FALSE, ...) {
  
  # --- 1. Data Preparation ---
  cl <- match.call()
  mf <- stats::model.frame(formula, data)
  X <- stats::model.matrix(formula, mf)
  y <- stats::model.response(mf)
  
  # Handle offset
  off <- stats::model.offset(mf)
  if (is.null(off)) {
    off <- rep(0, length(y))
  }
  
  # Handle weights
  wts <- stats::model.weights(mf)
  if (is.null(wts)) {
    wts <- rep(1, length(y))
  }
  
  # --- 2. Define the Log-Likelihood Function ---
  log_likelihood <- function(params) {
    # Separate parameters
    n_betas <- ncol(X)
    betas <- params[1:n_betas]
    log_sigma <- params[n_betas + 1]
    sigma <- exp(log_sigma)
    
    # Calculate mu_ln for the underlying Lognormal distribution
    eta <- as.vector(X %*% betas + off)
    mu_y <- exp(eta)
    mu_ln <- log(mu_y) - sigma^2 / 2
    
    # Calculate the log-likelihood for each observation using the C++ function
    log_lik_vec <- dpln_rcpp(y, mu = mu_ln, sigma = sigma, log_p = TRUE)
    
    # Return the weighted sum
    sum(wts * log_lik_vec)
  }
  
  # --- 3. Get Starting Values ---
  message("Fitting initial Poisson model for starting values...")
  initial_fit <- stats::glm.fit(
    x = X, y = y, family = stats::poisson(),
    weights = wts, offset = off
  )
  start_betas <- initial_fit$coefficients
  start_log_sigma <- 0 # Corresponds to sigma = 1
  
  start_params <- c(start_betas, start_log_sigma)
  param_names <- c(colnames(X), "log_sigma")
  names(start_params) <- param_names
  
  # --- 4. Perform Maximum Likelihood Estimation ---
  message("Optimizing Poisson-Lognormal log-likelihood...")
  
  # Set up control list for maxLik
  control_list <- list(...)
  if (verbose) {
    control_list$printLevel <- 2
  }
  
  fit <- maxLik::maxLik(
    logLik = log_likelihood,
    start = start_params,
    method = method,
    control = control_list
  )
  
  message("Estimation complete.")
  
  # --- 5. Structure and Return the Output ---
  final_params <- fit$estimate
  n_betas <- ncol(X)
  final_betas <- final_params[1:n_betas]
  final_log_sigma <- final_params[n_betas + 1]
  
  final_sigma <- exp(final_log_sigma)
  names(final_sigma) <- "sigma"
  
  # Calculate final fitted values
  final_eta <- as.vector(X %*% final_betas + off)
  final_mu <- exp(final_eta)
  
  # Create a structured result object
  result <- list(
    fit = fit,
    coefficients = final_betas,
    sigma = final_sigma,
    vcov = tryCatch(solve(-fit$hessian), error = function(e) NULL),
    logLik = fit$maximum,
    fitted.values = final_mu,
    residuals = y - final_mu,
    call = cl,
    formula = formula,
    terms = stats::terms(mf),
    xlevels = stats::.getXlevels(stats::terms(mf), mf),
    contrasts = attr(X, "contrasts"),
    nobs = length(y)
  )
  
  class(result) <- "poislogn_fit"
  return(result)
}

#' @export
print.poislogn_fit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nDispersion parameter (sigma):\n")
  print(x$sigma)
}

#' @export
summary.poislogn_fit <- function(object, ...) {
  # Get parameter estimates and vcov from maxLik object
  est <- object$fit$estimate
  vcov_log <- tryCatch(solve(-object$fit$hessian), error = function(e) NULL)
  
  if (is.null(vcov_log)) {
    warning("Could not invert Hessian matrix to calculate standard errors.")
  }
  
  # Use delta method to get SE for sigma (from log_sigma)
  # Var(f(x)) approx (f'(x))^2 * Var(x)
  # f(log_sigma) = exp(log_sigma) = sigma. f' = exp(log_sigma) = sigma
  se_log_sigma <- sqrt(diag(vcov_log)[length(est)])
  sigma_val <- exp(est[length(est)])
  se_sigma <- sigma_val * se_log_sigma
  
  # Combine SEs
  se <- sqrt(diag(vcov_log))
  se[length(se)] <- se_sigma
  
  # Combine estimates
  estimates <- c(object$coefficients, object$sigma)
  names(se) <- names(estimates)
  
  # Calculate z-values and p-values
  z_val <- estimates / se
  p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)
  
  coef_table <- cbind(
    Estimate = estimates,
    `Std. Error` = se,
    `z value` = z_val,
    `Pr(>|z|)` = p_val
  )
  
  result <- list(
    call = object$call,
    coefficients = coef_table,
    logLik = object$logLik,
    nobs = object$nobs
  )
  
  class(result) <- "summary.poislogn_fit"
  return(result)
}

#' @export
print.summary.poislogn_fit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nLog-Likelihood:", x$logLik, "(", nrow(x$coefficients), "df )")
  cat("\nNumber of observations:", x$nobs, "\n")
}

#' @export
predict.poislogn_fit <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  
  # If no new data, predict on the original data
  if (is.null(newdata)) {
    mu <- object$fitted.values
    if (type == "link") {
      # The link is log(mu), which is the linear predictor eta
      # We need to add back the offset to get the full eta
      eta <- log(mu)
      return(eta)
    } else {
      return(mu)
    }
  }
  
  # Prepare new data using the standard R approach
  tt <- object$terms
  Terms <- stats::delete.response(tt)
  
  # Create the model frame from the new data
  mf <- stats::model.frame(Terms, newdata, xlev = object$xlevels)
  
  # Create the model matrix
  X <- stats::model.matrix(Terms, mf, contrasts.arg = object$contrasts)
  
  # Extract the offset from the new model frame
  off <- stats::model.offset(mf)
  if (is.null(off)) {
    off <- rep(0, nrow(X))
  }
  
  # Calculate linear predictor
  eta <- as.vector(X %*% object$coefficients) + off
  
  if (type == "link") {
    return(eta)
  } else {
    return(exp(eta))
  }
}