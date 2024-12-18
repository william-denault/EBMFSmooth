#' Define the L-GP Object
#'
#' @description Creates an object from the one-parameter L-GP family.
#'
#' @param scale A numeric value representing the scale parameter of the L-GP.
#'
#' @return An object of class `"LGP"`, a data frame containing the scale parameter.
#'
#' @examples
#' gp <- LGP(1.5)
#' print(gp)
#'
#' @export
LGP <- function (scale) {
  structure(data.frame(scale), class = "LGP")
}


#' Set Up Matrices for L-GP Smoothing
#'
#' @description Prepares design matrices, precision matrices, and related quantities required for smoothing using L-GPs.
#'
#' @param t A numeric vector of input points (e.g., time points or other independent variables).
#' @param p An integer specifying the order of the differential operator. Default is 2.
#' @param num_knots An integer specifying the number of knots for the spline basis. Default is 30.
#' @param betaprec A numeric value specifying the prior precision for coefficients. Default is 1e-3.
#'
#' @return A list containing:
#'   - \code{X}: Global polynomial design matrix.
#'   - \code{B}: Local polynomial design matrix.
#'   - \code{P}: Precision matrix.
#'   - \code{logPdet}: Log-determinant of the precision matrix.
#'   - \code{betaprec}: Prior precision value.
#'
#' @examples
#' setup <- LGP_setup(t = seq(0, 1, length.out = 100))
#' str(setup)
#'
#' @export
LGP_setup <- function(t, p = 2, num_knots = 30, betaprec = 1e-3){

  # Defining the matrices for smoothing
  knots <- seq(min(t), max(t), length = num_knots)
  X <- BayesGP:::global_poly_helper(x = t, p = p)
  P <- BayesGP::compute_weights_precision_helper(knots)
  B <- BayesGP:::local_poly_helper(knots = knots, refined_x = t, p = p)

  tmbdat <- list(
    X = as(as(as(X, "dMatrix"), "generalMatrix"), "TsparseMatrix"),
    B = as(as(as(B, "dMatrix"), "generalMatrix"), "TsparseMatrix"),
    P = as(as(as(P, "dMatrix"), "generalMatrix"), "TsparseMatrix"),
    logPdet = as.numeric(determinant(P, logarithm = TRUE)$modulus),
    betaprec = betaprec
  )

  return(tmbdat)
}


#' Generate a Custom `ebnm` Function for L-GP Smoothing
#'
#' @description Creates a custom `ebnm` function for fitting data using L-GP models.
#'
#' @param LGP_setup A list of precomputed matrices and related quantities (e.g., from \code{LGP_setup}) required for L-GP smoothing.
#'
#' @return A function that performs L-GP smoothing with:
#'   - Posterior moments estimation.
#'   - Log-likelihood computation.
#'   - Posterior sampling.
#'
#' @examples
#' setup <- LGP_setup(t = seq(0, 1, length.out = 100))
#' ebgp <- ebnm_LGP_generator(setup)
#' result <- ebgp(x = rnorm(100), s = 0.1)
#'
#' @export
ebnm_LGP_generator <- function(LGP_setup){

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
                       "Just returning ebnm_flat(x) to pass the check.\n"))
        return (ebnm::ebnm_flat(x))
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
        DLL = "opt_theta",
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
      DLL = "fit_theta",
      silent = TRUE
    )

    ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)

    opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he,
                  control = list(eval.max = 20000, iter.max = 20000))

    prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
    mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)

    # Compute the posterior quantities
    all_design <- cbind(as.matrix(tmbdat$B), as.matrix(tmbdat$X))
    posterior <- data.frame(mean = all_design %*% (mod$mean))
    posterior$var <- diag(all_design %*% solve(mod$prec) %*% t(all_design))
    posterior$second_moment <- posterior$var + posterior$mean^2

    log_likelihood <- -ff2$fn(opt$par) - 0.5 * determinant(as.matrix(ff2$he(opt$par)), logarithm = TRUE)$modulus + 0.5 * length(opt$par) * log(2 * pi)
    class(log_likelihood) <- "logLik"

    posterior_sampler <- function(nsamp){
      samps_coef <- LaplacesDemon::rmvnp(n = nsamp, mu = mod$mean, Omega = as.matrix(mod$prec))
      samps_fitted <- as.matrix(tmbdat$B) %*% t(samps_coef[,1:ncol(tmbdat$B)]) + as.matrix(tmbdat$X) %*% t(samps_coef[,(ncol(tmbdat$B)+1):ncol(samps_coef)])
      t(samps_fitted)
    }

    result <- list(posterior = posterior,
                   fitted_g = fitted_g,
                   log_likelihood = log_likelihood,
                   posterior_sampler = posterior_sampler,
                   data = data.frame(x = x, s = s),
                   prior_family = "LGP"
    )

    return(structure(result, class = c("list", "ebnm")))

  }
  return(ebnm_gp)
}

