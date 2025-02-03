#' Flash Both Types: Fit Independent and Smooth Factors Using Flashier
#'
#' @description This function fits both independent and smooth factors to a dataset
#' using the `flashier` package. The order in which factors are added can be specified
#' via the \code{which_first} argument. Optionally, backfitting can be performed
#' to refine the final solution.
#'
#' @param data A matrix of observed data.
#' @param S (Optional) A matrix of standard errors corresponding to the data in \code{data}. Defaults to \code{NULL}.
#' @param k_smooth An integer specifying the number of smooth factors to fit. If \code{0}, this step is skipped. Defaults to \code{5}.
#' @param k_ind An integer specifying the number of independent factors to fit. If \code{0}, this step is skipped. Defaults to \code{5}.
#' @param ebnm_smooth A function for empirical Bayes smoothing of factors.
#' @param ebnm_indep (Optional) A function for empirical Bayes estimation of independent factors. Defaults to \code{ebnm::ebnm_point_normal}.
#' @param var_type (Optional) An integer specifying the variance type. Defaults to \code{NULL}.
#' @param which_first A string specifying whether to fit independent factors (\code{"ind"}) or smooth factors (\code{"smooth"}) first. Defaults to \code{"ind"}.
#' @param verbose An integer controlling verbosity. Set to \code{0} for no messages, or higher values for more detailed messages. Defaults to \code{0}.
#' @param backfitting_steps An integer specifying the number of backfitting iterations. If \code{0}, backfitting is skipped. Defaults to \code{0}.
#' @param convergence_fn A function defining the convergence criterion. Defaults to \code{flashier::flash_conv_crit_max_chg}.
#' @param convergence_tol A numeric value specifying the convergence tolerance for smooth factor fitting. Defaults to \code{1e-4}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: A `flashier` fit object containing the results of the fitted independent and smooth factors.
#'   \item \code{metadata}: A list summarizing the fitting process, including:
#'     \itemize{
#'       \item Number of independent and smooth factors fitted.
#'       \item Total runtime for all steps.
#'       \item Whether backfitting was performed.
#'     }
#' }
#'
#' @details
#' This function uses the following steps:
#' 1. Fits either independent or smooth factors first, as determined by \code{which_first}.
#' 2. Adds the other type of factor to the model using \code{flashier::flash_greedy()}.
#' 3. Optionally performs backfitting to refine the solution if \code{backfitting_steps > 0}.
#'
#' Verbose messages provide details about the runtime and number of factors fitted at each stage.
#'
#'
#' @export
flash_both_types <- function(data, S = NULL, k_smooth = 5, k_ind = 5, ebnm_smooth, ebnm_indep = NULL, var_type = NULL, which_first = "ind", verbose = 0, backfitting_steps = 0, convergence_fn = flashier::flash_conv_crit_max_chg, convergence_tol = 1e-4) {

  # Input validation
  if (!is.matrix(data)) {
    stop("`data` must be a matrix.")
  }
  if (!is.numeric(k_ind) || k_ind < 0) {
    stop("`k_ind` must be a non-negative integer.")
  }
  if (!is.numeric(k_smooth) || k_smooth < 0) {
    stop("`k_smooth` must be a non-negative integer.")
  }
  if (!is.numeric(backfitting_steps) || backfitting_steps < 0) {
    stop("`backfitting_steps` must be a non-negative integer.")
  }
  if (!which_first %in% c("ind", "smooth")) {
    stop("`which_first` must be either 'ind' or 'smooth'.")
  }

  # Set default for ebnm_indep if not provided
  if (is.null(ebnm_indep)) {
    ebnm_indep <- ebnm::ebnm_point_normal
  }

  # Initialize fit object and start tracking time
  overall_start_time <- Sys.time()
  fit <- flashier::flash_init(data = data, S = S, var_type = var_type)

  # Helper function for verbose messaging
  verbose_message <- function(message) {
    if (verbose > 0) cat(message, "\n")
  }

  # Step 1: Fit factors based on `which_first`
  if (which_first == "ind" && k_ind > 0) {
    verbose_message("Starting fit for independent factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose) |>
      flashier::flash_greedy(Kmax = k_ind, ebnm_fn = ebnm_indep)
    verbose_message(sprintf("Independent factors fit in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    k_ind_effective <- fit$n_factors
    verbose_message(sprintf("Number of independent factors: %d", k_ind_effective))
  } else if (which_first == "ind") {
    verbose_message("Skipping independent factor fitting as k_ind = 0.")
  }

  if (which_first == "smooth" && k_smooth > 0) {
    verbose_message("Starting fit for smooth factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose) |>
      flashier::flash_set_conv_crit(fn = convergence_fn, tol = convergence_tol) |>
      flashier::flash_greedy(Kmax = k_smooth, ebnm_fn = ebnm_smooth)
    verbose_message(sprintf("Smooth factors fit in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    k_smooth_effective <- fit$n_factors
    verbose_message(sprintf("Number of smooth factors: %d", k_smooth_effective))
  } else if (which_first == "smooth") {
    verbose_message("Skipping smooth factor fitting as k_smooth = 0.")
  }

  # Step 2: Add the other type of factors
  if (which_first == "ind" && k_smooth > 0) {
    verbose_message("Adding smooth factors after independent factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose) |>
      flashier::flash_set_conv_crit(fn = convergence_fn, tol = convergence_tol) |>
      flashier::flash_greedy(Kmax = k_smooth, ebnm_fn = ebnm_smooth)
    k_smooth_effective <- fit$n_factors - k_ind_effective
    verbose_message(sprintf("Smooth factors added in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    verbose_message(sprintf("Total number of smooth factors: %d", k_smooth_effective))
  } else if (which_first == "smooth" && k_ind > 0) {
    verbose_message("Adding independent factors after smooth factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose) |>
      flashier::flash_greedy(Kmax = k_ind, ebnm_fn = ebnm_indep)
    k_ind_effective <- fit$n_factors - k_smooth_effective
    verbose_message(sprintf("Independent factors added in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    verbose_message(sprintf("Total number of independent factors: %d", k_ind_effective))
  }

  # Step 3: Perform backfitting if backfitting_steps > 0
  if (backfitting_steps > 0) {
    verbose_message("Starting backfitting...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_backfit(maxiter = backfitting_steps) |>
      flashier::flash_nullcheck()
    verbose_message(sprintf("Backfitting completed in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
  } else {
    verbose_message("Skipping backfitting as backfitting_steps = 0.")
  }

  # Calculate total runtime
  total_runtime <- Sys.time() - overall_start_time

  # Prepare metadata
  metadata <- list(
    k_ind_effective = k_ind_effective,
    k_smooth_effective = k_smooth_effective,
    which_first = which_first,
    backfitting_steps = backfitting_steps,
    total_runtime = total_runtime
  )

  # Return fit object and metadata
  return(list(fit = fit, metadata = metadata))
}




#' Initialize Rank with Smoothness Estimation
#'
#' @description This function initializes the rank of factors or loadings in a dataset
#' using the `flashier` package and evaluates their smoothness. The function first fits
#' a specified number of independent factors and then applies a smoothing function to
#' estimate the smoothness of either the factors or loadings.
#'
#' @param data A matrix of observed data.
#' @param S (Optional) A matrix of standard errors corresponding to the data in \code{data}. Defaults to \code{NULL}.
#' @param k_total An integer specifying the total number of factors to fit. Defaults to \code{5}.
#' @param ebnm_indep (Optional) A function for empirical Bayes estimation of independent factors. Defaults to \code{ebnm::ebnm_point_normal}.
#' @param var_type (Optional) An integer specifying the variance type. Defaults to \code{NULL}.
#' @param verbose An integer controlling verbosity. Set to \code{0} for no messages, or higher values for more detailed messages. Defaults to \code{0}.
#' @param backfitting_steps An integer specifying the number of backfitting iterations. If \code{0}, backfitting is skipped. Defaults to \code{0}.
#' @param convergence_fn A function defining the convergence criterion. Defaults to \code{flashier::flash_conv_crit_max_chg}.
#' @param convergence_tol A numeric value specifying the convergence tolerance for smooth factor fitting. Defaults to \code{1e-4}.
#' @param smooth_term A string indicating whether to estimate smoothness for factors (\code{"F"}) or loadings (\code{"L"}). Defaults to \code{"F"}.
#' @param smooth_fun A smoothing function to estimate the smoothness of the selected terms.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: A `flashier` fit object containing the fitted factors and loadings.
#'   \item \code{smoothness}: A numeric vector of estimated smoothness values for the selected terms (factors or loadings).
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Fits a specified number of independent factors using \code{flashier::flash_greedy()}.
#' 2. Optionally performs backfitting to refine the solution.
#' 3. Applies a smoothing function to estimate the smoothness of the factors or loadings, as specified by \code{smooth_term}.
#'
#' @export
init_rank <- function(data, S = NULL, k_total = 5, ebnm_indep = NULL, var_type = NULL, verbose = 0, backfitting_steps = 0, convergence_fn = flashier::flash_conv_crit_max_chg, convergence_tol = 1e-4, smooth_term = "F", smooth_fun) {

  verbose_message <- function(message) {
    if (verbose > 0) cat(message, "\n")
  }

  # Input validation
  if (!is.matrix(data)) {
    stop("`data` must be a matrix.")
  }

  # Set default for ebnm_indep if not provided
  if (is.null(ebnm_indep)) {
    ebnm_indep <- ebnm::ebnm_point_normal
  }

  # Initialize fit object and start tracking time
  overall_start_time <- Sys.time()
  fit <- flashier::flash_init(data = data, S = S, var_type = var_type)

  verbose_message("Starting fit for independent factors...")
  start_time <- Sys.time()
  fit <- fit |>
    flashier::flash_set_verbose(verbose = verbose) |>
    flashier::flash_greedy(Kmax = k_total, ebnm_fn = ebnm_indep)
  verbose_message(sprintf(
    "Independent factors fit in %.2f seconds",
    as.numeric(Sys.time() - start_time, units = "secs")
  ))
  verbose_message(sprintf("Number of independent factors: %d", fit$n_factors))

  if (backfitting_steps > 0) {
    verbose_message("Starting backfitting...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_backfit(maxiter = backfitting_steps) |>
      flashier::flash_nullcheck()
    verbose_message(sprintf("Backfitting completed in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
  } else {
    verbose_message("Skipping backfitting as backfitting_steps = 0.")
  }

  # Apply smoothing to each factor or loading
  verbose_message("Estimating smoothness of factors or loadings...")
  start_time <- Sys.time()
  k_effective <- fit$n_factors
  smoothness <- rep(0, k_effective)
  for (k in 1:k_effective) {
    if (smooth_term == "F") {
      smoothness_result <- smooth_fun(x = fit$F_pm[, k], s = fit$F_psd[, k])
      if(class(smoothness_result$g_init) == "Matern"){
        smoothness[k] <- as.numeric(smoothness_result$fitted_g$theta[1])
      } else {
        smoothness[k] <- 1/as.numeric(smoothness_result$fitted_g$theta)
      }

    } else if (smooth_term == "L") {
      smoothness_result <- smooth_fun(x = fit$L_pm[, k], s = fit$L_psd[, k])
      if(class(smoothness_result$g_init) == "Matern"){
        smoothness[k] <- as.numeric(smoothness_result$fitted_g$theta[1])
      } else {
        smoothness[k] <- 1/as.numeric(smoothness_result$fitted_g$theta)
      }
    } else {
      stop("Invalid value for `smooth_term`.")
    }
  }

  verbose_message(sprintf(
    "Smoothness estimated in %.2f seconds",
    as.numeric(Sys.time() - start_time, units = "secs")
  ))

  # Return fit object and metadata
  return(list(fit = fit, smoothness = smoothness))
}



#' Fit Independent and Smooth Factors with Automatic Initialization
#'
#' @description This function fits both independent and smooth factors to a dataset
#' using the `flashier` package. The initialization of factors is based on their
#' correlation, with smoothness estimated using a smoothing function.
#'
#' @param data A matrix of observed data.
#' @param S (Optional) A matrix of standard errors corresponding to the data in \code{data}. Defaults to \code{NULL}.
#' @param k_smooth An integer specifying the number of smooth factors to fit. Defaults to \code{5}.
#' @param k_ind An integer specifying the number of independent factors to fit. Defaults to \code{5}.
#' @param ebnm_smooth A list containing functions for empirical Bayes smoothing of factors (\code{[2]}) and loadings (\code{[1]}).
#' @param ebnm_indep (Optional) A function for empirical Bayes estimation of independent factors. Defaults to \code{ebnm::ebnm_point_normal}.
#' @param var_type (Optional) An integer specifying the variance type. Defaults to \code{NULL}.
#' @param which_first A string specifying whether to fit independent factors (\code{"ind"}) or smooth factors (\code{"smooth"}) first. Defaults to \code{"ind"}.
#' @param verbose An integer controlling verbosity. Set to \code{0} for no messages, or higher values for more detailed messages. Defaults to \code{0}.
#' @param backfitting_steps An integer specifying the number of backfitting iterations. If \code{0}, backfitting is skipped. Defaults to \code{0}.
#' @param convergence_fn A function defining the convergence criterion. Defaults to \code{flashier::flash_conv_crit_max_chg}.
#' @param convergence_tol A numeric value specifying the convergence tolerance for smooth factor fitting. Defaults to \code{1e-4}.
#' @param smooth_term A string indicating whether to estimate smoothness for factors (\code{"F"}) or loadings (\code{"L"}). Defaults to \code{"F"}.
#' @param smooth_fun (Optional) A user-defined smoothing function. If \code{NULL}, it is inferred from \code{ebnm_smooth} based on \code{smooth_term}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: A `flashier` fit object containing the results of the fitted independent and smooth factors.
#'   \item \code{metadata}: A list summarizing the fitting process, including:
#'     \itemize{
#'       \item \code{k_ind_effective}: Number of independent factors fitted.
#'       \item \code{k_smooth_effective}: Number of smooth factors fitted.
#'       \item \code{which_first}: The order of factor fitting (\code{"ind"} or \code{"smooth"}).
#'       \item \code{backfitting_steps}: Number of backfitting steps performed.
#'       \item \code{initialization_result}: The result of the initialization step, including smoothness estimates.
#'     }
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Initializes rank by fitting all factors and estimating their smoothness.
#' 2. Orders factors by smoothness (increasing or decreasing) based on \code{which_first}.
#' 3. Fits factors iteratively, starting with the type specified in \code{which_first}.
#' 4. Optionally performs backfitting to refine the final solution.
#'
#'
#' @export
flash_both_types_init <- function(data, S = NULL, k_smooth = 5, k_ind = 5, ebnm_smooth, ebnm_indep = NULL, var_type = NULL, which_first = "ind", verbose = 0, backfitting_steps = 0, convergence_fn = flashier::flash_conv_crit_max_chg, convergence_tol = 1e-4, smooth_term = "F", smooth_fun = NULL) {

  # Set default for ebnm_indep if not provided
  if (is.null(ebnm_indep)) {
    ebnm_indep <- ebnm::ebnm_point_normal
  }

  # Initialize fit object and start tracking time
  fit <- flashier::flash_init(data = data, S = S, var_type = var_type)


  verbose_message <- function(message) {
    if (verbose > 0) cat(message, "\n")
  }

  if(is.null(smooth_fun)){
    if(smooth_term == "F"){
      smooth_fun = ebnm_smooth[[2]]
    } else if(smooth_term == "L"){
      smooth_fun = ebnm_smooth[[1]]
    }
  }

  # initialize the rank based on the correlation
  verbose_message("Initializing rank based on correlation...")
  start_time <- Sys.time()
  initialization_result <- init_rank(data = data, S = S, k_total = (k_smooth + k_ind), ebnm_indep = ebnm_indep, var_type = var_type, verbose = 0, backfitting_steps = 0, convergence_fn = convergence_fn, convergence_tol = convergence_tol, smooth_term = smooth_term, smooth_fun = smooth_fun)
  verbose_message(sprintf(
    "Rank initialized in %.2f seconds",
    as.numeric(Sys.time() - start_time, units = "secs")
  ))
  k_total_effective <- initialization_result$fit$n_factors

  # if which_first == "smooth", order the initialization from the most smooth to the least smooth
  # if which_first == "ind", order the initialization from the least smooth to the most smooth
  smoothness_rank <- if(which_first == "smooth") order(initialization_result$smoothness, decreasing = TRUE) else order(initialization_result$smoothness, decreasing = FALSE)

  customized_init_fn_generator <- function(which_term){
    func <- function(flash){
      # stop if which_term is not of the right length
      if(length(which_term) != 1){
        stop(paste("which_term is not of the right length: ", length(which_term)))
      }

      # stop if which_term is not in the right range
      if(which_term > k_total_effective){
        # say it will be set to k_total_effective
        warning(paste("which_term is not in the right range: ", which_term, " > ", k_total_effective, ". It will be set to ", k_total_effective))
        which_term = k_total_effective
      }

      # extract the num_term th factor/loading
      return(list(loading = as.vector(initialization_result$fit$L_pm[, smoothness_rank[which_term]]),
                  factor = as.vector(initialization_result$fit$F_pm[, smoothness_rank[which_term]]))
      )
    }
    return(func)
  }

  # Step 1: Fit factors based on `which_first`
  if (which_first == "ind" && k_ind > 0) {
    verbose_message("Starting fit for independent factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose)
    for (i in 1:k_ind) {
      fit <- fit |>
        flashier::flash_greedy(Kmax = 1, ebnm_fn = ebnm_indep, init_fn = customized_init_fn_generator(which_term = i))
    }
    verbose_message(sprintf("Independent factors fit in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    k_ind_effective <- fit$n_factors
    verbose_message(sprintf("Number of independent factors: %d", k_ind_effective))
  } else if (which_first == "ind") {
    verbose_message("Skipping independent factor fitting as k_ind = 0.")
  }
  if (which_first == "smooth" && k_smooth > 0) {
    verbose_message("Starting fit for smooth factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose)
    for (i in 1:k_smooth) {
      fit <- fit |>
        flashier::flash_greedy(Kmax = 1, ebnm_fn = ebnm_smooth, init_fn = customized_init_fn_generator(which_term = i))
    }
    verbose_message(sprintf("Smooth factors fit in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    k_smooth_effective <- fit$n_factors
    verbose_message(sprintf("Number of smooth factors: %d", k_smooth_effective))
  } else if (which_first == "smooth") {
    verbose_message("Skipping smooth factor fitting as k_smooth = 0.")
  }

  # Step 2: Add the other type of factors
  if (which_first == "ind" && k_smooth > 0) {
    verbose_message("Adding smooth factors after independent factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose)

    for (i in 1:k_smooth) {
      fit <- fit |>
        flashier::flash_greedy(Kmax = 1, ebnm_fn = ebnm_smooth, init_fn = customized_init_fn_generator(which_term = (k_ind_effective + i)))
    }

    k_smooth_effective <- fit$n_factors - k_ind_effective
    verbose_message(sprintf("Smooth factors added in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    verbose_message(sprintf("Total number of smooth factors: %d", k_smooth_effective))
  } else if (which_first == "smooth" && k_ind > 0) {
    verbose_message("Adding independent factors after smooth factors...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_set_verbose(verbose = verbose)
    for (i in 1:k_ind) {
      fit <- fit |>
        flashier::flash_greedy(Kmax = 1, ebnm_fn = ebnm_indep, init_fn = customized_init_fn_generator(which_term = (k_smooth_effective + i)))
    }
    k_ind_effective <- fit$n_factors - k_smooth_effective
    verbose_message(sprintf("Independent factors added in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
    verbose_message(sprintf("Total number of independent factors: %d", k_ind_effective))
  }

  # Step 3: Perform backfitting if backfitting_steps > 0
  if (backfitting_steps > 0) {
    verbose_message("Starting backfitting...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_backfit(maxiter = backfitting_steps) |>
      flashier::flash_nullcheck()
    verbose_message(sprintf("Backfitting completed in %.2f seconds", as.numeric(Sys.time() - start_time, units = "secs")))
  } else {
    verbose_message("Skipping backfitting as backfitting_steps = 0.")
  }

  # Prepare metadata
  metadata <- list(
    k_ind_effective = k_ind_effective,
    k_smooth_effective = k_smooth_effective,
    which_first = which_first,
    backfitting_steps = backfitting_steps,
    initialization_result = initialization_result
  )

  return(list(fit = fit, metadata = metadata))

}



#' Flash Both Types Sequential: Fit Smooth and Independent Factors Alternately
#'
#' @description This function fits smooth and independent factors sequentially to a dataset
#' using the `flashier` package. It alternates between fitting smooth and independent factors,
#' selecting the one with the better ELBO (evidence lower bound optimization) at each step.
#' Backfitting can optionally be performed to refine the final solution.
#'
#' @param data A matrix of observed data.
#' @param S (Optional) A matrix of standard errors corresponding to the data in \code{data}. Defaults to \code{NULL}.
#' @param k_total An integer specifying the total number of factors to fit. Defaults to \code{10}.
#' @param ebnm_smooth A function for empirical Bayes smoothing of factors.
#' @param ebnm_indep (Optional) A function for empirical Bayes estimation of independent factors. Defaults to \code{ebnm::ebnm_point_normal}.
#' @param var_type (Optional) An integer specifying the variance type. Defaults to \code{NULL}.
#' @param verbose An integer controlling verbosity. Set to \code{0} for no messages, or higher values for more detailed messages. Defaults to \code{0}.
#' @param backfitting_steps An integer specifying the number of backfitting iterations. If \code{0}, backfitting is skipped. Defaults to \code{0}.
#' @param convergence_fn A function defining the convergence criterion. Defaults to \code{flashier::flash_conv_crit_max_chg}.
#' @param convergence_tol A numeric value specifying the convergence tolerance for smooth factor fitting. Defaults to \code{1e-4}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: A `flashier` fit object containing the fitted factors and loadings.
#'   \item \code{type_factor}: A character vector indicating the type of each fitted factor ("smooth" or "independent").
#' }
#'
#' @details
#' This function performs the following steps:
#' 1. Alternates between fitting smooth and independent factors using greedy optimization.
#'    - At each step, the factor type with the better ELBO is retained.
#'    - Verbose messages indicate the progress and results of each step.
#' 2. Optionally performs backfitting to refine the solution if \code{backfitting_steps > 0}.
#'
#' @examples
#' library(flashier)
#' library(ebnm)
#' data <- matrix(rnorm(100), nrow = 10)
#' result <- flash_both_types_sequential(
#'   data = data,
#'   k_total = 5,
#'   ebnm_smooth = ebnm::ebnm_point_normal,
#'   verbose = 1
#' )
#' print(result$type_factor)
#'
#' @export
flash_both_types_sequential <- function(data, S = NULL, k_total = 10, ebnm_smooth, ebnm_indep = NULL, var_type = NULL, verbose = 0, backfitting_steps = 0, convergence_fn = flashier::flash_conv_crit_max_chg, convergence_tol = 1e-5) {

  # Set default for ebnm_indep if not provided
  if (is.null(ebnm_indep)) {
    ebnm_indep <- ebnm::ebnm_point_normal
  }

  # Initialize fit object and start tracking time
  fit <- flashier::flash_init(data = data, S = S, var_type = var_type) |>
    flashier::flash_set_verbose(verbose = verbose) |>
    flashier::flash_set_conv_crit(fn = convergence_fn, tol = convergence_tol)


  verbose_message <- function(message) {
    if (verbose > 0) cat(message, "\n")
  }

  type_factor <- character()

  # Step 1: Fit factors sequentially
  verbose_message(sprintf("Starting sequential greedy-fitting for %d factors...", k_total))
  for (i in seq_len(k_total)) {
    # Fitting the smooth factor
    verbose_message(sprintf("Fitting smooth factor %d...", i))
    start_time <- Sys.time()
    fit_smooth <- fit |>
      flashier::flash_greedy(Kmax = 1, ebnm_fn = ebnm_smooth)
    smooth_time <- Sys.time() - start_time
    verbose_message(sprintf("Smooth factor %d fit in %.2f seconds", i, as.numeric(smooth_time, units = "secs")))

    # Fitting the independent factor
    verbose_message(sprintf("Fitting independent factor %d...", i))
    start_time <- Sys.time()
    fit_indep <- fit |>
      flashier::flash_greedy(Kmax = 1, ebnm_fn = ebnm_indep)
    indep_time <- Sys.time() - start_time
    verbose_message(sprintf("Independent factor %d fit in %.2f seconds", i, as.numeric(indep_time, units = "secs")))

    # Compare ELBO values and update the fit
    if (fit_smooth$elbo > fit_indep$elbo) {
      verbose_message(sprintf("Smooth factor has better ELBO: %.2f compared to independent factor: %.2f", fit_smooth$elbo, fit_indep$elbo))
      if (fit_smooth$n_factors > fit$n_factors) {
        fit <- fit_smooth
        type_factor <- c(type_factor, "smooth")
        verbose_message(sprintf("Smooth factor added to the model. Total factors: %d", fit$n_factors))
      } else {
        verbose_message("Smooth factor not added. Total factors remain unchanged.")
      }
    } else {
      verbose_message(sprintf("Independent factor has better ELBO: %.2f compared to smooth factor: %.2f", fit_indep$elbo, fit_smooth$elbo))
      if (fit_indep$n_factors > fit$n_factors) {
        fit <- fit_indep
        type_factor <- c(type_factor, "independent")
        verbose_message(sprintf("Independent factor added to the model. Total factors: %d", fit$n_factors))
      } else {
        verbose_message("Independent factor not added. Total factors remain unchanged.")
      }
    }
  }
  verbose_message(sprintf("Sequential greedy-fitting for %d factors completed.", k_total))
  verbose_message(sprintf("Total independent factors: %d", sum(type_factor == "independent")))
  verbose_message(sprintf("Total smooth factors: %d", sum(type_factor == "smooth")))

  # Step 2: Perform backfitting if backfitting_steps > 0
  if (backfitting_steps > 0) {
    verbose_message("Starting backfitting...")
    start_time <- Sys.time()
    fit <- fit |>
      flashier::flash_backfit(maxiter = backfitting_steps) |>
      flashier::flash_nullcheck()
    backfit_time <- Sys.time() - start_time
    verbose_message(sprintf("Backfitting completed in %.2f seconds", as.numeric(backfit_time, units = "secs")))
  } else {
    verbose_message("Skipping backfitting as backfitting_steps = 0.")
  }

  return(list(fit = fit, type_factor = type_factor))
}



#' Plot Factors and Loadings from a Flashier Model fitted with LGP prior
#'
#' @description This function plots the normalized factors and loadings from a fitted
#' flashier model alongside the true normalized factors and loadings.
#'
#' @param f_mod A flashier model object containing \code{L_pm} (loadings) and \code{F_pm} (factors).
#' @param t_vec A numeric vector representing the time points for the smooth loadings.
#' @param true_L A matrix of true loadings (optional). Defaults to \code{NULL}.
#' @param true_F A matrix of true factors (optional). Defaults to \code{NULL}.
#'
#' @return Generates plots showing the comparison between estimated and true loadings/factors for all components.
#'
#' @examples
#' set.seed(123)
#' # Simulated data
#' sim <- simulate_matrix(n = 200, p = 100, k_smooth = 2, k_ind = 0, sigma_E = 1)
#' f.mod <- flashier::flash_init(data = sim$Y) |>
#'   flashier::flash_greedy(Kmax = 2, ebnm_fn = ebnm::ebnm_point_normal)
#' plot_factors_and_loadings_1d(f_mod, sim$t_vec, sim$L_mat[, 1:2], sim$F_mat[, 1:2])
#'
#' @export
plot_factors_and_loadings_1d <- function(f_mod, t_vec, true_L = NULL, true_F = NULL) {
  # Normalize estimated loadings and factors
  L_norm <- apply(f_mod$L_pm, 2, function(x) x / max(x))
  F_norm <- sweep(f_mod$F_pm, 2, apply(f_mod$L_pm, 2, max), `*`)

  # Normalize true loadings and factors (if provided) using the same normalization
  if (!is.null(true_L)) {
    true_L_orig <- true_L
    true_L <- sweep(true_L, 2, apply(true_L, 2, max), `/`)
  }
  if (!is.null(true_F)) {
    if(!is.null(true_L)){
      true_F <- sweep(true_F, 2, apply(true_L_orig, 2, max), `*`)
    }
  }

  # Determine number of components and set plot layout
  n_components <- ncol(L_norm)
  n_rows <- n_components
  n_cols <- 2
  par(mfrow = c(n_rows, n_cols))

  # Plot loadings and factors for all components
  for (i in seq_len(n_components)) {
    # Plot loadings
    plot(t_vec, L_norm[, i], type = "l", col = "red", lty = "dashed",
         ylim = range(c(L_norm[, i], if (!is.null(true_L)) true_L[, i] else NULL), na.rm = TRUE),
         xlab = "t", ylab = "Loadings", main = paste("Loading", i))
    if (!is.null(true_L)) {
      lines(t_vec, true_L[, i], col = "black", lty = "solid")
    }

    # Plot factors
    plot(1:nrow(F_norm), F_norm[, i], type = "l", col = "red", lty = "dashed",
         ylim = range(c(F_norm[, i], if (!is.null(true_F)) true_F[, i] else NULL), na.rm = TRUE),
         xlab = "Features", ylab = "Factors", main = paste("Factor", i))
    if (!is.null(true_F)) {
      lines(1:nrow(F_norm), true_F[, i], col = "black", lty = "solid", type = "o")
    }
  }

  # Reset plotting layout
  par(mfrow = c(1, 1))
}
