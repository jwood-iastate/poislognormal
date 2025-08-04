#' Fit a Poisson-Lognormal GEE Model (Fixed Sigma)
#'
#' @description
#' Fits a Poisson-Lognormal model for correlated count data using Generalized
#' Estimating Equations (GEE). This version requires a fixed, known sigma.
#'
#' @details
#' This function acts as a wrapper around `geeM::geem` to fit a GEE model
#' with a custom variance function for the Poisson-Lognormal distribution.
#'
#' The variance function is defined as:
#' \deqn{Var(Y) = \mu + (e^{\sigma^2}-1) \mu^2}
#'
#' The dispersion parameter `sigma` is **not** estimated by this function.
#' It must be provided as a fixed, known value. For iterative estimation of
#' sigma, see `pln.gee.est`.
#'
#' @param formula A formula specifying the model.
#' @param id A vector that identifies the clusters.
#' @param data A data frame containing the variables.
#' @param sigma A single, positive numeric value for the fixed dispersion parameter.
#' @param corstr A character string specifying the correlation structure. See
#'   `?geeM::geem` for options (e.g., "independence", "exchangeable", "ar1").
#' @param ... Additional arguments passed to `geeM::geem`.
#'
#' @return An object of class `pln_gee`, which inherits from `geem`.
#'
#' @export
#' @importFrom geeM geem
#' @importFrom stats terms model.frame make.link
pln.gee <- function(formula, id, data, sigma, corstr = "independence", ...) {
  
  # --- 1. Input Validation ---
  if (missing(sigma) || !is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
    stop("`sigma` must be provided as a single, positive number.", call. = FALSE)
  }
  cl <- match.call()
  
  # --- 2. Define the Custom GEE Family ---
  # Create a complete family object with all necessary components for geeM
  link_obj <- stats::make.link("log")
  pln_family <- list(
    link = "log",
    linkfun = link_obj$linkfun,
    linkinv = link_obj$linkinv,
    mu.eta = link_obj$mu.eta, # This component was missing
    variance = function(mu) {
      mu + (exp(sigma^2) - 1) * mu^2
    }
  )
  class(pln_family) <- "family"
  
  # --- 3. Fit the GEE model using geeM ---
  fit <- geeM::geem(
    formula = formula,
    id = id,
    data = data,
    family = pln_family,
    corstr = corstr,
    ...
  )
  
  # --- 4. Structure and Return the Output ---
  fit$call <- cl
  fit$pln.sigma <- sigma
  class(fit) <- c("pln_gee", class(fit))
  return(fit)
}


#' Fit a Poisson-Lognormal GEE Model with Iterative Sigma Estimation
#'
#' @description
#' Fits a Poisson-Lognormal GEE model and iteratively estimates the dispersion
#' parameter `sigma` using the method of moments.
#'
#' @details
#' This function provides a convenient wrapper that automates the estimation of
#' the dispersion parameter `sigma`. It uses an iterative algorithm that
#' alternates between:
#' 1. Fitting the GEE model with a fixed `sigma` to estimate coefficients.
#' 2. Using the Pearson residuals from the fit to find an updated `sigma` that
#'    satisfies the moment condition (i.e., sum of squared Pearson residuals
#'    equals the residual degrees of freedom).
#'
#' This process repeats until `sigma` converges.
#'
#' @param formula A formula specifying the model.
#' @param id A character string specifying the name of the cluster/subject ID column.
#' @param data A data frame containing the variables.
#' @param corstr A character string specifying the correlation structure.
#' @param tol Tolerance for convergence of `sigma`.
#' @param max_iter Maximum number of iterations.
#' @param ... Additional arguments passed to `geeM::geem`.
#'
#' @return An object of class `pln_gee`, which inherits from `geem`, with the
#'   estimated `sigma` included.
#'
#' @export
#' @importFrom stats poisson glm optimize
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
#' gee_model_est <- pln.gee.est(
#'   formula = Total_crashes ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k + offset(lnlength),
#'   id = "ID",
#'   data = washington_roads,
#'   corstr = "exchangeable"
#' )
#'
#' # View the summary
#' summary(gee_model_est)
pln.gee.est <- function(formula, id, data, corstr = "independence", tol = 1e-4, max_iter = 50, ...) {
  
  cl <- match.call()
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  nobs <- length(y)
  
  # Robustly extract the ID vector from the data frame
  if (!is.character(id) || length(id) != 1) {
    stop("`id` must be a character string specifying the ID column name.", call. = FALSE)
  }
  if (!(id %in% names(data))) {
    stop(paste0("ID column '", id, "' not found in data."), call. = FALSE)
  }
  id_vec <- data[[id]]
  
  objective_fun <- function(sigma, y, mu, p) {
    V <- mu + (exp(sigma^2) - 1) * mu^2
    pearson_chi2 <- sum(((y - mu)^2) / V, na.rm = TRUE)
    df_resid <- nobs - p
    (pearson_chi2 - df_resid)^2
  }
  
  message("Fitting initial Poisson model for starting values...")
  init_fit <- stats::glm(formula, data = data, family = stats::poisson())
  mu <- init_fit$fitted.values
  p <- length(init_fit$coefficients)
  init_disp <- sum(((y - mu)^2) / mu, na.rm = TRUE) / (nobs - p)
  sigma_current <- if (init_disp > 1) sqrt(log(init_disp)) else 0.5
  
  message(paste("Initial sigma estimate:", round(sigma_current, 4)))
  
  iter <- 0
  converged <- FALSE
  while (iter < max_iter) {
    iter <- iter + 1
    sigma_old <- sigma_current
    gee_fit <- try(
      pln.gee(formula, id = id_vec, data = data, sigma = sigma_current, corstr = corstr, ...),
      silent = TRUE
    )
    if (inherits(gee_fit, "try-error")) {
      warning("GEE fit failed at iteration ", iter, ". Stopping.")
      break
    }
    mu <- gee_fit$fitted.values
    p <- length(gee_fit$beta)
    opt <- optimize(
      f = objective_fun, interval = c(1e-4, 20),
      y = y, mu = mu, p = p
    )
    sigma_current <- opt$minimum
    message(paste("Iteration:", iter, " Sigma:", round(sigma_current, 4)))
    if (abs(sigma_current - sigma_old) < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning("Sigma estimation did not converge within ", max_iter, " iterations.")
  }
  
  final_fit <- pln.gee(formula, id = id_vec, data = data, sigma = sigma_current, corstr = corstr, ...)
  final_fit$call <- cl
  final_fit$pln.sigma.estimated <- TRUE
  final_fit$pln.iterations <- iter
  return(final_fit)
}

#' @export
print.pln_gee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients (Naive):\n")
  print(x$beta)
  sigma_val <- round(x$pln.sigma, 4)
  if (!is.null(x$pln.sigma.estimated) && x$pln.sigma.estimated) {
    cat("\nDispersion parameter (sigma) was estimated to be:", sigma_val, "\n")
  } else {
    cat("\nDispersion parameter (sigma) was fixed at:", sigma_val, "\n")
  }
}

#' @export
summary.pln_gee <- function(object, ...) {
  res <- NextMethod()
  res$pln.sigma <- object$pln.sigma
  res$pln.sigma.estimated <- object$pln.sigma.estimated
  class(res) <- c("summary.pln_gee", class(res))
  return(res)
}

#' @export
print.summary.pln_gee <- function(x, ...) {
  NextMethod()
  sigma_val <- round(x$pln.sigma, 4)
  if (!is.null(x$pln.sigma.estimated) && x$pln.sigma.estimated) {
    cat("\nNote: Poisson-Lognormal dispersion 'sigma' was estimated to be", sigma_val, "\n")
  } else {
    cat("\nNote: Poisson-Lognormal dispersion 'sigma' was fixed at", sigma_val, "\n")
  }
}
