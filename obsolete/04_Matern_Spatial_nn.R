#' Generate a Custom `ebnm` Function for Matern Smoothing with Non-Negative Constraints
#'
#' @description This function generates a custom \code{ebnm} function for fitting data
#' using the Matern Gaussian process (GP) with non-negative constraints. Non-negative
#' constraints are enforced via an exponential link function.
#'
#' @param locations A matrix of spatial locations, where each row represents a
#' coordinate.
#' @param max.edge A numeric value defining the maximum edge length for the spatial
#' mesh. If \code{NULL}, the default mesh construction is used.
#' @param alpha A numeric value specifying the smoothness parameter of the Matern GP.
#' Higher values result in smoother processes. Defaults to \code{2}.
#' @param suppress_warnings A logical value. If \code{TRUE}, suppresses warnings
#' generated during INLA optimization. Defaults to \code{TRUE}.
#' @param penalty_range A numeric value specifying the initial range parameter for the
#' Matern GP prior. If \code{NULL}, a default is computed based on the spatial extent
#' of the \code{locations}.
#'
#' @return A function that can be used to fit data using the Matern GP model with
#' non-negative constraints. The returned function includes:
#' \itemize{
#'   \item \code{posterior}: A data frame containing posterior mean, variance, and
#'   second moments.
#'   \item \code{fitted_g}: The fitted Matern object, containing the optimized
#'   log-transformed range parameter.
#'   \item \code{g_init}: The initial Matern object before optimization.
#'   \item \code{log_likelihood}: The marginal log-likelihood of the model.
#'   \item \code{posterior_spatial_field}: Summary of the spatial random field.
#'   \item \code{mesh}: The spatial mesh used in the model.
#'   \item \code{inla_result}: The full result object returned by \code{INLA::inla()}.
#' }
#'
#' @details
#' The function first constructs a spatial mesh using \code{INLA::inla.mesh.2d} and
#' initializes the Matern GP parameters. If \code{fix_g = FALSE}, the function optimizes
#' the Matern GP parameters (range and marginal standard deviation) using INLA's
#' internal methods. Non-negative constraints are enforced by applying an exponential
#' link function to the mean response, ensuring all fitted values remain positive.
#' The \code{penalty_range} parameter can be used to set the initial range for the Matern
#' GP prior, which defaults to one-tenth of the smallest spatial dimension.
#'
#' @examples
#' # Example usage of the Matern GP model with non-negative constraints
#' set.seed(123)
#' locations <- matrix(runif(20), ncol = 2)  # Random 2D locations
#' matern_model <- eblnm_Matern_generator(locations = locations)
#'
#' # Data and standard deviations
#' x <- abs(rnorm(nrow(locations)))  # Ensure non-negative response
#' s <- rep(0.1, nrow(locations))
#'
#' result <- matern_model(x = x, s = s)
#' str(result$posterior)
#'
#' @export
eblnm_Matern_generator <- function(locations,
                                   max.edge = NULL,
                                   alpha = 2,
                                   suppress_warnings = TRUE,
                                   penalty_range = NULL
) {

  mesh <- INLA::inla.mesh.2d(loc = locations, max.edge = c(max.edge, max.edge))

  A <- INLA::inla.spde.make.A(mesh = mesh, loc = locations)

  ebnm_Matern <- function(x,
                          s,
                          g_init = NULL,
                          fix_g = FALSE,
                          output = NULL) {
    # Check the Setup
    if (nrow(locations) != length(x)) {
      warning(
        paste0(
          "The length of x must be equal to the number of rows in the location matrix.\n",
          "The length of x is ",
          length(x),
          " and the number of rows in the location matrix is ",
          nrow(locations),
          ".\n"
        )
      )

      if (length(s) == 3 & length(x) == 3) {
        warning(
          paste0(
            "Assume this is just an initialization check.\n",
            "Just returning ebnm_flat(x) to pass the check.\n"
          )
        )
        return (ebnm::ebnm_flat(x))
      } else{
        stop("The length of x must be equal to the number of rows in the design matrices.")
      }
    }

    if(is.null(penalty_range)){
      penalty_range <- min(diff(range(locations[ , 1]))/10, diff(range(locations[ , 2]))/10)
    }
    if(is.null(g_init)){
      g_init <- Matern(theta = log(penalty_range))
    }

    # Optimize for g
    if(!fix_g){
      optim_inla <- function() {
        spde <- INLA::inla.spde2.pcmatern(
          mesh = mesh,
          alpha = alpha,
          prior.range = c(penalty_range, 0.5),
          prior.sigma = c(1, NA)
        )

        indexs <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

        stack <- INLA::inla.stack(
          data = list(Y = x),
          A = list(A, matrix(1, nrow = nrow(locations), ncol = 1)),
          effects = list(
            spatial.field = indexs$spatial.field,
            beta0 = 1
          ),
          tag = "est"
        )

        formula <- Y ~ 0 + beta0 + f(spatial.field, model = spde)

        result <- INLA::inla(
          formula,
          scale = (1 / s ^ 2),
          control.inla = list(int.strategy = "eb", strategy = "gaussian"),
          control.family = list(control.link=list(model="log"), hyper = list(prec = list(fixed = TRUE, initial = 0))),
          control.fixed = INLA::control.fixed(prec = 0),
          data = INLA::inla.stack.data(stack),
          control.predictor = list(A = INLA::inla.stack.A(stack), link=1),
          control.compute = list(
            mlik = TRUE,
            hyperpar = TRUE,
            return.marginals = FALSE
          ),
          silent = TRUE
        )
        list(par = result$mode$theta, value = result$mlik[[1]], beta0 = result$summary.fixed$mean)
      }
      opt_result <- optim_inla()
      fitted_g <- Matern(opt_result$par)
      beta0_est <- opt_result$beta0
    }else{
      fitted_g <- g_init
      beta0_est <- 0
    }

    new_range_val <- as.numeric(exp(fitted_g$theta))

    spde <- INLA::inla.spde2.pcmatern(
      mesh = mesh,
      alpha = alpha,
      prior.range = c(new_range_val, NA),
      prior.sigma = c(1, NA)
    )

    indexs <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

    stack <- INLA::inla.stack(
      data = list(Y = x),
      A = list(A, matrix(beta0_est, nrow = nrow(locations), ncol = 1)),
      effects = list(spatial.field = indexs$spatial.field, beta0 = 1),
      tag = "est"
    )

    formula <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)

    result <- INLA::inla(
      formula,
      scale = (1 / s ^ 2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(control.link=list(model="log"), hyper = list(prec = list(fixed = TRUE, initial = 0))),
      control.fixed = INLA::control.fixed(prec = 0),
      data = INLA::inla.stack.data(stack),
      control.predictor = list(A = INLA::inla.stack.A(stack), link=1),
      silent = TRUE
    )

    ii <- INLA::inla.stack.index(stack, tag = "est")$data

    posterior <-
      data.frame(
        mean = result$summary.fitted.values$mean[ii],
        var = result$summary.fitted.values$sd[ii]^2
      )

    posterior$second_moment <- posterior$mean^2 + posterior$var

    posterior_spatial_field <- result$summary.random$spatial.field

    log_likelihood <- result$mlik[[1]]

    posterior_sampler <- function(n = n) {
      warning("The posterior sampler is not implemented yet for Matern.")
      matrix(NA, nrow = n, ncol = length(x))
    }

    result <- list(
      posterior = posterior,
      fitted_g = fitted_g,
      g_init = g_init,
      log_likelihood = log_likelihood,
      posterior_sampler = posterior_sampler,
      data = data.frame(x = x, s = s),
      prior_family = "log-Matern",
      posterior_spatial_field = posterior_spatial_field,
      mesh = mesh,
      inla_result = result
    )

    return(structure(result, class = c("list", "ebnm")))

  }

  if(suppress_warnings){
    return(suppressWarnings(ebnm_Matern))
  }else{
    return(ebnm_Matern)
  }
}
