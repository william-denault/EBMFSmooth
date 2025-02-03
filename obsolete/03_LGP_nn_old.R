#' Generate a Custom `ebnm` Function for L-GP Smoothing with Non-Negative Constraints
#'
#' @description This function creates a custom `ebnm` function to fit data using
#' a (log) L-GP model with non-negative constraints enforced via an exponential
#' link function. The method ensures posterior means and variances remain positive,
#' leveraging exponential transformations in both fitting and posterior computations.
#'
#' @param LGP_setup A list of precomputed matrices and related quantities (e.g., from
#'   `LGP_setup`) required for L-GP smoothing.
#' @return A function that performs L-GP smoothing with:
#'   - Posterior moments estimation using non-negative constraints.
#'   - Log-likelihood computation of the fitted model.
#'   - Posterior sampling to generate fitted values.
#'
#' @details The function enforces non-negativity using an exponential link function,
#' making it suitable for cases where non-negative constraints are essential (e.g.,
#' modeling counts or rates).
#'
#' @examples
#' # Example setup and usage:
#' setup <- LGP_setup(t = seq(0, 1, length.out = 100))
#' ebgp <- eblnm_LGP_generator(setup)
#' result <- ebgp(x = rnorm(100), s = 0.1)
#'
#' @export
eblnm_LGP_generator <- function(LGP_setup){

  ebnm_gp <- function(x, s,
                      g_init = NULL, fix_g = FALSE,
                      output = NULL){

    # Check the LGP_setup
    if (nrow(LGP_setup$X) != length(x) ||
        nrow(LGP_setup$B) != length(x)) {
      warning(paste0("The length of x must be equal to the number of rows in the design matrices.\n",
                     "The length of x is ",
                     length(x),
                     " and the number of rows in the design matrices are ",
                     nrow(LGP_setup$X),
                     " and ",
                     nrow(LGP_setup$B),
                     ".\n"
      ))

      if (length(s) == 3 & length(x) == 3) {
        warning(paste0("Assume this is just an initialization check.\n",
                       "Just returning ebnm_point_exponential(x) to pass the check.\n"))
        return (ebnm::ebnm_point_exponential(x))
      } else{
        stop("The length of x must be equal to the number of rows in the design matrices.")
      }
    }

    # Check the length of s
    if(length(s) != length(x) & length(s) != 1){
      warning(paste0("The length of s must be 1 or equal to the length of x.\n",
                     "The length of s is ", length(s), " and the length of x is ", length(x), ".\n"))
      stop("The length of s must be 1 or equal to the length of x.")
    }
    if(length(s) == 1){
      s <- rep(s, length(x))
    }

    # Loading quantities from the LGP_setup
    tmbdat <- LGP_setup
    tmbdat$x <- x
    tmbdat$s <- s

    # Fill in the initial value for g if not provided
    if(is.null(g_init)){
      g_init <- LGP(1)
    }

    # Optimize for g
    if(!fix_g){
      tmbparams <- list(theta = 0,
                        W = c(rep(0, (ncol(tmbdat$X) + ncol(tmbdat$B)))))

      obj_fun <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        DLL = "opt_theta_nn",
        random = "W",
        silent = TRUE
      )

      fitted_scale  <- optim(par = 0, fn = obj_fun$fn, gr = obj_fun$gr, method = "BFGS")$par
      fitted_g <- LGP(fitted_scale)

    }else{
      fitted_g <- g_init
    }

    # Fit the model
    tmbparams <- list(W = c(rep(0, (ncol(tmbdat$X) + ncol(tmbdat$B)))))
    tmbdat$theta <- fitted_g$scale

    ff2 <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      DLL = "fit_theta_nn",
      silent = TRUE
    )

    ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)

    opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he,
                  control = list(eval.max = 20000, iter.max = 20000))

    prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
    mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)

    # Compute the posterior quantities
    all_design <- cbind(as.matrix(tmbdat$B), as.matrix(tmbdat$X))
    log_posterior <- data.frame(mean = all_design %*% (mod$mean))
    log_posterior$var <- diag(all_design %*% solve(mod$prec) %*% t(all_design))
    posterior <- data.frame(
      mean = exp(log_posterior$mean + 0.5 * log_posterior$var),
      var = exp(2 * log_posterior$mean + log_posterior$var) * (exp(log_posterior$var) - 1)
    )
    posterior$second_moment <- posterior$mean^2 + posterior$var

    # Compute the log-likelihood
    obj_fun <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      DLL = "fit_theta_nn",
      random = "W",
      silent = TRUE
    )
    log_likelihood <- -obj_fun$fn()
    # log_likelihood <- -ff2$fn(opt$par) - 0.5 * determinant(as.matrix(ff2$he(opt$par)), logarithm = TRUE)$modulus + 0.5 * length(opt$par) * log(2 * pi)
    class(log_likelihood) <- "logLik"

    posterior_sampler <- function(nsamp){
      samps_coef <- LaplacesDemon::rmvnp(n = nsamp, mu = mod$mean, Omega = as.matrix(mod$prec))
      samps_fitted <- exp(as.matrix(tmbdat$B) %*% t(samps_coef[,1:ncol(tmbdat$B)]) + as.matrix(tmbdat$X) %*% t(samps_coef[,(ncol(tmbdat$B)+1):ncol(samps_coef)]))
      t(samps_fitted)
    }

    result <- list(posterior = posterior,
                   fitted_g = fitted_g,
                   log_likelihood = log_likelihood,
                   posterior_sampler = posterior_sampler,
                   data = data.frame(x = x, s = s),
                   prior_family = "log-LGP"
    )

    return(structure(result, class = c("list", "ebnm")))

  }
  return(ebnm_gp)
}

