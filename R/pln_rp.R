#' Fit a Random Parameters Poisson-Lognormal Model
#'
#' @description
#' Estimates a random parameters Poisson-Lognormal (RP-PLN) model using maximum
#' simulated likelihood. Supports correlated or uncorrelated random parameters,
#' various parameter distributions, and panel data structures.
#'
#' @param formula A formula for the model's fixed and random parameters.
#' @param rpar_formula A one-sided formula specifying the random parameters
#'   (e.g., `~ var1 + var2`).
#' @param data A data frame containing the variables for the model.
#' @param panel A character string specifying the column in `data` that
#'   identifies panel groups. If `NULL`, each observation is treated as a
#'   separate panel (cross-sectional data).
#' @param rpardists An optional named character vector specifying the distribution
#'   for each random parameter. Options: `"n"` (normal), `"ln"` (log-normal),
#'   `"t"` (triangular), `"u"` (uniform), `"g"` (gamma). Defaults to normal.
#' @param ndraws The number of Halton draws for simulation.
#' @param correlated Logical. If `TRUE`, estimates correlated random parameters
#'   (forces all to be normally distributed).
#' @param method The optimization algorithm for `maxLik::maxLik`.
#' @param ... Additional arguments passed to `maxLik::maxLik`.
#'
#' @return A list object of class `pln.rp` containing the model results.
#'
#' @export
#' @importFrom stats update qnorm qlnorm qunif qgamma coef glm pnorm vcov cov2cor
#' @importFrom hardhat mold forge
#'
#' @examples
#' \dontrun{
#' data("washington_roads")
#'
#' # Create a binary indicator for high AADT
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' # Fit a model with uncorrelated random parameters
#' rp_model <- pln.rp(
#'   formula = Animal ~ lnaadt + speed50,
#'   rpar_formula = ~ ShouldWidth04 + AADTover10k,
#'   data = washington_roads,
#'   panel = "ID",
#'   rpardists = c(ShouldWidth04 = "n", AADTover10k = "g"),
#'   ndraws = 200,
#'   method = "BFGS"
#' )
#'
#' summary(rp_model)
#' }
pln.rp <- function(formula, rpar_formula, data, panel = NULL,
                   rpardists = NULL, ndraws = 500, correlated = FALSE,
                   method = "BFGS", ...) {
  
  cl <- match.call()
  
  # --- 1. Data Preparation with hardhat ---
  # Combine formulas for a single mold
  full_formula <- update(formula, paste("~ . +", paste(all.vars(rpar_formula), collapse = " + ")))
  processed <- hardhat::mold(full_formula, data)
  
  y <- processed$outcomes[[1]]
  X_all <- processed$predictors
  
  rand_vars <- all.vars(rpar_formula)
  fixed_vars <- setdiff(colnames(X_all), rand_vars)
  
  X_fixed <- X_all[, fixed_vars, drop = FALSE]
  X_rand <- X_all[, rand_vars, drop = FALSE]
  
  fixed_names <- colnames(X_fixed)
  rpar_names <- colnames(X_rand)
  
  if (is.null(panel)) {
    panel_index <- 1:nrow(data)
  } else {
    panel_index <- as.numeric(as.factor(data[[panel]]))
  }
  n_panels <- length(unique(panel_index))
  
  # --- 2. Validate Parameters ---
  if (correlated) {
    message("Correlated parameters requested. All random parameters will be normal.")
    rpardists <- rep("n", length(rpar_names))
    names(rpardists) <- rpar_names
  } else if (is.null(rpardists)) {
    rpardists <- rep("n", length(rpar_names))
    names(rpardists) <- rpar_names
  }
  
  # --- 3. Get Starting Values ---
  message("Fitting base model for starting values...")
  start_fit <- pln.glm(formula, data)
  start_coefs <- coef(start_fit)
  
  start_means_fixed <- start_coefs[fixed_names]
  start_means_rand <- start_coefs[rpar_names]
  
  start_params <- c(start_means_fixed, start_means_rand)
  
  if (correlated) {
    chol_mat <- diag(0.1, nrow = length(rpar_names))
    chol_vals <- chol_mat[lower.tri(chol_mat, diag = TRUE)]
    start_params <- c(start_params, chol_vals)
    names(start_params) <- c(paste0("mean_", fixed_names),
                             paste0("mean_", rpar_names),
                             paste0("chol_", 1:length(chol_vals)))
  } else {
    start_sds <- rep(0.1, length(rpar_names))
    start_params <- c(start_params, start_sds)
    names(start_params) <- c(paste0("mean_", fixed_names),
                             paste0("mean_", rpar_names),
                             paste0("sd_", rpar_names))
  }
  
  start_params <- c(start_params, "log_sigma" = log(start_fit$sigma))
  
  # --- 4. Define the Log-Likelihood Function ---
  log_likelihood_rp <- function(params) {
    n_fixed <- ncol(X_fixed)
    n_rand <- ncol(X_rand)
    
    betas_fixed <- params[paste0("mean_", fixed_names)]
    betas_rand_mean <- params[paste0("mean_", rpar_names)]
    
    # Generate Halton draws for simulation
    halton_draws <- make_draws(ndraws, n_rand)
    
    # Generate parameter draws for each random parameter
    draws <- matrix(NA, nrow = ndraws, ncol = n_rand)
    if (correlated) {
      chol_vals <- params[grep("^chol_", names(params))]
      chol_mat <- matrix(0, n_rand, n_rand)
      chol_mat[lower.tri(chol_mat, diag = TRUE)] <- chol_vals
      norm_draws <- qnorm(halton_draws)
      draws <- (norm_draws %*% t(chol_mat)) + 
        matrix(betas_rand_mean, nrow = ndraws, ncol = n_rand, byrow = TRUE)
    } else {
      rand_sds <- abs(params[grep("^sd_", names(params))])
      for (i in 1:n_rand) {
        draws[, i] <- switch(
          rpardists[i],
          "n"  = qnorm(halton_draws[, i], mean = betas_rand_mean[i], sd = rand_sds[i]),
          "ln" = qlnorm(halton_draws[, i], meanlog = betas_rand_mean[i], sdlog = rand_sds[i]),
          "t"  = qtri(halton_draws[, i], lower = betas_rand_mean[i] - rand_sds[i], upper = betas_rand_mean[i] + rand_sds[i], mode = betas_rand_mean[i]),
          "u"  = qunif(halton_draws[, i], min = betas_rand_mean[i] - rand_sds[i], max = betas_rand_mean[i] + rand_sds[i]),
          "g"  = qgamma(halton_draws[, i], shape = betas_rand_mean[i]^2 / rand_sds[i]^2, rate = betas_rand_mean[i] / rand_sds[i]^2)
        )
      }
    }
    
    # Calculate simulated probabilities
    eta_fixed <- as.vector(X_fixed %*% betas_fixed)
    eta_rand_sim <- X_rand %*% t(draws)
    
    eta_sim <- eta_fixed + eta_rand_sim
    mu_y_sim <- exp(eta_sim)
    
    sigma <- exp(params["log_sigma"])
    mu_ln_sim <- log(mu_y_sim) - sigma^2 / 2
    
    # Probability for each observation and each draw
    prob_mat <- dpln_rcpp(y, mu = mu_ln_sim, sigma = sigma, log_p = FALSE)
    prob_mat[prob_mat <= 0] <- 1e-300
    
    # Calculate panel-level likelihoods
    panel_log_probs <- rowsum(log(prob_mat), group = panel_index, reorder = FALSE)
    max_log_prob <- apply(panel_log_probs, 1, max)
    panel_probs <- exp(max_log_prob) * rowMeans(exp(panel_log_probs - max_log_prob))
    panel_probs[panel_probs <= 0] <- 1e-300
    
    log(panel_probs)
  }
  
  # --- 5. Model Fitting ---
  message("Optimizing random parameters log-likelihood...")
  fit <- maxLik::maxLik(
    logLik = log_likelihood_rp,
    start = start_params,
    method = method,
    ...
  )
  message("Estimation complete.")
  
  # --- 6. Post-Estimation Processing ---
  # Re-calculate Cholesky for correlated model to store VCV matrix
  if (correlated) {
    chol_vals <- fit$estimate[grep("^chol\\.", names(fit$estimate))]
    chol_mat <- matrix(0, length(rpar_names), length(rpar_names))
    chol_mat[lower.tri(chol_mat, diag = TRUE)] <- chol_vals
    rpar_vcov <- t(chol_mat) %*% chol_mat
    dimnames(rpar_vcov) <- list(rpar_names, rpar_names)
  } else {
    rpar_vcov <- NULL
  }
  
  result <- list(
    fit = fit,
    coefficients = fit$estimate,
    vcov = tryCatch(solve(-fit$hessian), error = function(e) NULL),
    rpar_vcov = rpar_vcov,
    logLik = sum(fit$maximum),
    nobs = length(y),
    n_panels = n_panels,
    n_params = length(start_params),
    call = cl,
    formula = formula,
    blueprint = processed$blueprint, # Store the hardhat blueprint
    rpar_formula = rpar_formula,
    correlated = correlated,
    rpardists = rpardists,
    fixed_names = fixed_names,
    rpar_names = rpar_names,
    panel_id_name = panel # Store panel column name for prediction
  )
  class(result) <- "pln.rp"
  return(result)
}


#' @export
print.pln.rp <- function(x, ...) {
  cat("Random Parameters Poisson-Lognormal Model\n")
  cat("-----------------------------------------\n")
  cat("Log-Likelihood:", round(x$logLik, 3), "\n")
  cat("Observations:", x$nobs, "\n")
  if (x$n_panels < x$nobs) {
    cat("Panels:", x$n_panels, "\n")
  }
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#' @export
summary.pln.rp <- function(object, ...) {
  est <- object$coefficients
  se <- tryCatch(sqrt(diag(object$vcov)), error = function(e) rep(NA, length(est)))
  z_val <- est / se
  p_val <- 2 * pnorm(-abs(z_val))
  
  results_df <- data.frame(
    Estimate = est,
    `Std. Error` = se,
    `z-value` = z_val,
    `Pr(>|z|)` = p_val,
    check.names = FALSE
  )
  
  cat("Random Parameters Poisson-Lognormal Model\n")
  cat("-----------------------------------------\n")
  cat("Log-Likelihood:", round(object$logLik, 3), "\n")
  
  n_obs_for_bic <- if (object$n_panels < object$nobs) object$n_panels else object$nobs
  cat("AIC:", round(-2 * object$logLik + 2 * object$n_params, 3), "\n")
  cat("BIC:", round(-2 * object$logLik + log(n_obs_for_bic) * object$n_params, 3), "\n\n")
  
  cat("Parameter Estimates:\n")
  print(results_df)
  
  if (object$correlated) {
    cat("\nRandom Parameter Correlation Matrix:\n")
    print(cov2cor(object$rpar_vcov))
  }
  
  invisible(object)
}


# Internal helper function to generate draws based on model estimates
generate_rpar_draws_pln <- function(object, ndraws) {
  params <- object$coefficients
  n_rand <- length(object$rpar_names)
  
  mean_param_names <- paste0("mean_", object$rpar_names)
  betas_rand_mean <- params[mean_param_names]
  
  halton_draws <- make_draws(ndraws, n_rand, type = "scrambled-halton-rand-perm")
  draws <- matrix(NA, nrow = ndraws, ncol = n_rand)
  
  if (object$correlated) {
    chol_mat <- t(chol(object$rpar_vcov))
    norm_draws <- qnorm(halton_draws)
    draws <- (norm_draws %*% t(chol_mat)) +
      matrix(betas_rand_mean, nrow = ndraws, ncol = n_rand, byrow = TRUE)
  } else {
    sd_param_names <- paste0("sd_", object$rpar_names)
    rand_sds <- abs(params[sd_param_names])
    
    for (i in 1:n_rand) {
      dist <- object$rpardists[i]
      draws[, i] <- switch(
        dist,
        "n"  = qnorm(halton_draws[, i], mean = betas_rand_mean[i], sd = rand_sds[i]),
        "ln" = qlnorm(halton_draws[, i], meanlog = betas_rand_mean[i], sdlog = rand_sds[i]),
        "t"  = qtri(halton_draws[, i], lower = betas_rand_mean[i] - rand_sds[i], upper = betas_rand_mean[i] + rand_sds[i], mode = betas_rand_mean[i]),
        "u"  = qunif(halton_draws[, i], min = betas_rand_mean[i] - rand_sds[i], max = betas_rand_mean[i] + rand_sds[i]),
        "g"  = qgamma(halton_draws[, i], shape = betas_rand_mean[i]^2 / rand_sds[i]^2, rate = betas_rand_mean[i] / rand_sds[i]^2)
      )
    }
  }
  return(draws)
}


#' Predict from a Random Parameters Poisson-Lognormal Model
#'
#' @name predict.pln.rp
#' @param object A model object of class `pln.rp` estimated by `pln.rp()`.
#' @param newdata A data frame containing variables for prediction.
#' @param method The prediction method: `"Simulated"` (default), `"Exact"`, or `"Individual"`.
#' @param ndraws The number of draws to use for simulation-based methods.
#' @param ... Additional arguments (currently ignored).
#'
#' @note
#' The `"Individual"` method performs a Bayesian update to find the conditional
#' parameters for each panel and thus requires the outcome variable to be present
#' in `newdata`. The `"Exact"` method is only available for normal, uniform, and
#' gamma random parameter distributions.
#'
#' @return A numeric vector of predicted values.
#' @export
#' @importFrom stats model.matrix
predict.pln.rp <- function(object, newdata, method = "Simulated", ndraws = 2000, ...) {
  
  method <- match.arg(method, c("Simulated", "Exact", "Individual"))
  
  # --- 1. Data Preparation using hardhat blueprint ---
  processed_new <- hardhat::forge(newdata, object$blueprint)
  X_all <- processed_new$predictors
  
  X_fixed <- X_all[, object$fixed_names, drop = FALSE]
  X_rand <- X_all[, object$rpar_names, drop = FALSE]
  
  params <- object$coefficients
  mean_fixed_names <- paste0("mean_", object$fixed_names)
  betas_fixed <- params[mean_fixed_names]
  eta_fixed <- as.vector(X_fixed %*% betas_fixed)
  
  # --- 2. Simulated Prediction (Population Average) ---
  if (method == "Simulated") {
    draws <- generate_rpar_draws_pln(object, ndraws)
    eta_rand_sim <- X_rand %*% t(draws)
    eta_sim <- eta_fixed + eta_rand_sim
    mu_y_sim <- exp(eta_sim)
    predictions <- rowMeans(mu_y_sim)
    return(predictions)
  }
  
  # --- 3. "Exact" Prediction (using MGF) ---
  if (method == "Exact") {
    mean_rand_names <- paste0("mean_", object$rpar_names)
    betas_rand_mean <- params[mean_rand_names]
    
    if (object$correlated) {
      # MGF for correlated normal: exp(t'μ + 0.5 * t'Σt)
      quad_term <- rowSums((X_rand %*% object$rpar_vcov) * X_rand)
      rand_factor <- exp(X_rand %*% betas_rand_mean + 0.5 * quad_term)
    } else {
      # Check for unsupported distributions
      if (any(!object$rpardists %in% c("n", "u", "g"))) {
        warning("`method = 'Exact'` is only supported for normal ('n'), uniform ('u'), and gamma ('g') distributions. Falling back to `method = 'Simulated'`.")
        return(predict(object, newdata, method = "Simulated", ndraws = ndraws, ...))
      }
      
      sd_rand_names <- paste0("sd_", object$rpar_names)
      rand_sds <- abs(params[sd_rand_names])
      rand_factors <- matrix(1, nrow = nrow(newdata), ncol = length(object$rpar_names))
      
      for(i in seq_along(object$rpar_names)) {
        mu <- betas_rand_mean[i]
        sigma <- rand_sds[i]
        x_val <- X_rand[, i]
        
        rand_factors[, i] <- switch(
          object$rpardists[i],
          "n" = exp(x_val * mu + 0.5 * x_val^2 * sigma^2),
          "u" = ifelse(x_val == 0, 1, exp(x_val * mu) * sinh(x_val * sigma) / (x_val * sigma)),
          "g" = (1 - x_val * sigma^2 / mu)^(-mu^2 / sigma^2)
        )
      }
      rand_factor <- apply(rand_factors, 1, prod)
    }
    predictions <- exp(eta_fixed) * rand_factor
    return(as.vector(predictions))
  }
  
  # --- 4. Individual Prediction (Panel or Cross-section) ---
  if (method == "Individual") {
    y_name <- all.vars(object$formula)[1]
    if (!y_name %in% names(newdata)) {
      stop("Method 'Individual' requires the outcome variable '", y_name, "' in newdata.")
    }
    y <- newdata[[y_name]]
    
    panel_id_col <- object$panel_id_name
    panel_index <- if (is.null(panel_id_col)) 1:nrow(newdata) else as.numeric(as.factor(newdata[[panel_id_col]]))
    
    draws <- generate_rpar_draws_pln(object, ndraws)
    eta_rand_sim <- X_rand %*% t(draws)
    eta_sim <- eta_fixed + eta_rand_sim
    mu_y_sim <- exp(eta_sim)
    
    sigma <- exp(params["log_sigma"])
    mu_ln_sim <- log(mu_y_sim) - sigma^2 / 2
    
    prob_mat <- dpln_rcpp(y, mu = mu_ln_sim, sigma = sigma, log_p = FALSE)
    prob_mat[prob_mat <= 0] <- 1e-300
    
    panel_log_probs <- rowsum(log(prob_mat), group = panel_index, reorder = FALSE)
    max_log_prob <- apply(panel_log_probs, 1, max)
    panel_likelihoods <- exp(panel_log_probs - max_log_prob)
    
    numerator <- panel_likelihoods %*% draws
    denominator <- rowSums(panel_likelihoods)
    
    cond_betas_panel <- numerator / denominator
    cond_betas_obs <- cond_betas_panel[panel_index, ]
    
    predictions <- exp(eta_fixed + rowSums(X_rand * cond_betas_obs))
    return(predictions)
  }
}