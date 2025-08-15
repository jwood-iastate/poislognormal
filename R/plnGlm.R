#' Fit a Poisson-Lognormal Regression Model
#'
#' @description
#' Fits a Poisson-Lognormal regression model using direct maximum likelihood
#' estimation via the `maxLik` package.
#'
#' @details
#' This function directly optimizes the Poisson-Lognormal log-likelihood 
#' function to estimate both the regression coefficients (betas) and the 
#' dispersion parameter (`sigma`) simultaneously. This provides a direct 
#' estimate and standard error for `sigma`.
#'
#' The `sigma` parameter is estimated on the log scale (`log(sigma)`) to ensure
#' it remains positive during optimization. The returned estimate and standard
#' error are transformed back to the original scale.
#'
#' @param formula A formula specifying the model.
#' @param data A data frame containing the variables in the model.
#' @param offset An optional character string specifying an offset variable in 
#'                `data`.
#' @param weights An optional character string specifying a weights variable in 
#'                `data`.
#' @param method A character string specifying the optimization method for 
#'              `maxLik`. See `?maxLik::maxLik` for options.
#' @param verbose A logical value. If `TRUE`, detailed optimization output is 
#'                printed.
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
#' # Fit a model for total crashes
#' # The offset term must be on the log scale.
#' pln_model_roads <- pln.glm(
#'   Total_crashes ~ lnaadt  + lnlength,
#'   data = washington_roads
#' )
#'
#' # View the model summary
#' summary(pln_model_roads)
#' 
#' predict(pln_model_roads, type = "response", newdata = washington_roads)
#' 
pln.glm <- function(formula, data, offset = NULL, weights = NULL,
                    method = "NM", verbose = FALSE, ...) {
  
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
    
    # Calculate mu for the underlying Lognormal distribution
    eta <- as.vector(X %*% betas + off)
    eta <- ifelse(is.null(eta),1e-10, eta)
    mu<- exp(eta)

    
    # Calculate the log-likelihood for each observation
    log_lik_vec <- dpLnorm(y, mu, sigma, log = TRUE)
    
    # Return the weighted sum
    return(wts * log_lik_vec)
  }
  
  # --- 3. Get Starting Values ---
  initial_fit <- stats::glm.fit(
    x = X, y = y, family = stats::poisson(),
    weights = wts, offset = off
  )
  start_betas <- initial_fit$coefficients
  start_log_sigma <- 0 # Corresponds to sigma = 1
  
  start_params <- c(start_betas, start_log_sigma)
  param_names <- c(colnames(X), "ln(sigma)")
  names(start_params) <- param_names
  
  # --- 4. Perform Maximum Likelihood Estimation ---
  
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
  vcov_log <- -object$fit$hessian
  
  if (is.null(vcov_log)) {
    warning("Could not invert Hessian matrix to calculate standard errors.")
  }
  
  # Combine SEs
  se <- sqrt(1/diag(vcov_log))
  
  
  # Calculate z-values and p-values
  z_val <- est / se
  p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)
  
  coef_table <- cbind(
    Estimate = est,
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