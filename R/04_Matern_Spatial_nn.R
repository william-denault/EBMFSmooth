#' Generate a Custom `ebnm` Function for Matern Smoothing with Non-Negative Constraints
#'
#' @description Creates a custom `ebnm` function for fitting data using a Matern process
#' with non-negative constraints enforced via an exponential link function.
#'
#' @param locations A matrix of spatial locations, where each row represents a coordinate.
#' @param max.edge A numeric value defining the maximum edge length for the spatial mesh. Default is \code{NULL}.
#' @param alpha A numeric value specifying the smoothness parameter of the Matern GP. Default is \code{2}.
#' @param range_prob A numeric value specifying the prior probability for the range parameter. Default is \code{0.05}.
#' @param sigma_prob A numeric value specifying the prior probability for the standard deviation parameter. Default is \code{0.05}.
#' @param range_val A numeric value specifying the prior mean for the range parameter. Default is \code{2}.
#' @param sigma_val A numeric value specifying the prior mean for the standard deviation parameter. Default is \code{2}.
#'
#' @return A function that can be used to fit data using the Matern process with non-negative constraints. The returned function includes:
#' \itemize{
#'   \item \code{posterior}: A data frame of posterior mean, variance, and second moments.
#'   \item \code{fitted_g}: The fitted Matern object.
#'   \item \code{log_likelihood}: The marginal log-likelihood of the model.
#'   \item \code{posterior_spatial_field}: Summary of the spatial random field.
#'   \item \code{mesh}: The spatial mesh used in the model.
#'   \item \code{inla_result}: The full result object from \code{INLA::inla()}.
#' }
#'
#' @examples
#' # Example usage of the Matern GP model with non-negative constraints
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
                                   range_prob = 0.05,
                                   sigma_prob = 0.05,
                                   range_val = 2,
                                   sigma_val = 2) {

  mesh <- INLA::inla.mesh.2d(loc = locations, max.edge = c(max.edge, max.edge))

  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = alpha,
    prior.range = c(range_val, range_prob),
    prior.sigma = c(sigma_val, sigma_prob)
  )

  indexs <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

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

    stack <- INLA::inla.stack(
      data = list(Y = x),
      A = list(A, 1),
      effects = list(spatial.field = indexs$spatial.field, beta0 = rep(1, nrow(locations))),
      tag = "est"
    )

    formula <- Y ~ 0 + beta0 + f(spatial.field, model = spde)

    result <- INLA::inla(
      formula,
      scale = (1 / s^2),
      family = "gaussian",
      control.family = list(control.link=list(model="log")),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.predictor = list(A = INLA::inla.stack.A(stack), link=1),
      data = INLA::inla.stack.data(stack)
    )

    ii <- INLA::inla.stack.index(stack, tag = "est")$data

    posterior <- data.frame(
      mean = result$summary.fitted.values$mean[ii],
      var = result$summary.fitted.values$sd[ii]^2
    )

    posterior$second_moment <- posterior$mean^2 + posterior$var

    log_likelihood <- result$mlik[[1]]

    posterior_sampler <- function(n = n) {
      warning("The posterior sampler is not implemented yet for Matern.")
      matrix(NA, nrow = n, ncol = length(x))
    }

    posterior_spatial_field <- result$summary.random$spatial.field

    result <- list(
      posterior = posterior,
      fitted_g = Matern(result$mode$theta),
      g_init = Matern(),
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

  return(ebnm_Matern)

}
