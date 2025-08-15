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
#' @examples
#' ## --- Example with washington_roads dataset ---
#' data(washington_roads)
#' # Create a binary indicator for high AADT
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' 
#' # Fit a model for Total crashes
#' # Note: We now pass the ID column directly, which is best practice.
#' total_crash_model <- pln.gee(Total_crashes ~ lnaadt + speed50 +
#'  ShouldWidth04 + AADTover10k + lnlength,
#'   sigma = 0.5,
#'   id = ID, 
#'   data = washington_roads
#'   )
#'
#' summary(total_crash_model)
#' 
#' @note Using `pln.gee` requires specifying the value of `sigma`. 
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
  
  # --- 2. Evaluate the id argument to get the vector ---
  # This takes the `id` argument (e.g., ID) and evaluates it in the `data` frame
  # to get the actual vector of cluster IDs.
  id_vec <- eval(substitute(id), data, parent.frame())
  
  # --- 3. Fit the GEE model using geeM ---
  fit <- geeM::geem(
    formula = formula,
    id = id_vec, # Pass the resolved vector here
    data = data,
    family = poisLogn(sigma),
    corstr = corstr,
    ...
  )
  
  # --- 4. Structure and Return the Output ---
  fit$call <- cl
  fit$pln.sigma <- sigma
  class(fit) <- c("pln_gee", class(fit))
  return(fit)
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