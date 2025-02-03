#' Simulate a Random Matrix with Smooth and Independent Factors
#'
#' @description This function generates a random matrix \eqn{Y} of size \eqn{n \times p} as a sum of
#' \eqn{k_smooth} smooth factors and \eqn{k_ind} independent factors, with sparsity enforced and optional Gaussian noise.
#'
#' @param n Number of rows in the generated matrix.
#' @param p Number of columns in the generated matrix.
#' @param k_smooth Number of smooth factors to include.
#' @param k_ind Number of independent factors to include.
#' @param sigma_E Standard deviation of the Gaussian noise. Defaults to \code{1}.
#' @param seed (Optional) Seed for random number generation. Defaults to \code{NULL}.
#' @param positive Logical; if \code{TRUE}, factors and loadings will be non-negative. Defaults to \code{FALSE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Y}: The simulated matrix of size \eqn{n \times p}.
#'   \item \code{L_mat}: The loading matrix of size \eqn{n \times (k_smooth + k_ind)}.
#'   \item \code{F_mat}: The factor matrix of size \eqn{p \times (k_smooth + k_ind)}.
#'   \item \code{E_mat}: The Gaussian noise matrix of size \eqn{n \times p}.
#'   \item \code{t_vec}: The time vector used for generating smooth loadings.
#' }
#'
#' @examples
#' set.seed(123)
#' sim <- simulate_matrix(n = 200, p = 100, k_smooth = 2, k_ind = 3, sigma_E = 1)
#' str(sim)
#'
#' @export
simulate_matrix <- function(n, p, k_smooth, k_ind, sigma_E = 1, seed = NULL, positive = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  # Generate smooth factors
  t_vec <- seq(0, 5, length.out = n)
  L_smooth <- matrix(0, nrow = n, ncol = k_smooth)
  F_smooth <- matrix(0, nrow = p, ncol = k_smooth)

  for (i in seq_len(k_smooth)) {
    # Generate smooth loadings
    smooth_loading <- runif(1, -pi, pi) * sin(t_vec) + runif(1, -pi, pi) * cos(t_vec)
    if (positive) smooth_loading <- exp(smooth_loading)
    L_smooth[, i] <- smooth_loading

    # Generate sparse smooth factors
    sparsity <- runif(1, 0.3, 0.7)  # Random sparsity level
    smooth_factor <- rexp(p) * rbinom(p, 1, sparsity)
    if (positive) smooth_factor <- exp(smooth_factor)  # Ensure non-negative values
    F_smooth[, i] <- smooth_factor
  }

  # Generate independent factors
  L_ind <- matrix(0, nrow = n, ncol = k_ind)
  for (i in seq_len(k_ind)) {
    sparsity_L <- runif(1, 0.3, 0.7)  # Random sparsity level for loadings
    if (positive) {
      L_ind[, i] <- rexp(n) * rbinom(n, 1, sparsity_L)  # Sparse, positive loadings
    } else {
      L_ind[, i] <- rnorm(n) * rbinom(n, 1, sparsity_L)  # Sparse, real-valued loadings
    }
  }

  F_ind <- matrix(0, nrow = p, ncol = k_ind)
  for (i in seq_len(k_ind)) {
    sparsity <- runif(1, 0.3, 0.7)  # Random sparsity level
    independent_factor <- rexp(p) * rbinom(p, 1, sparsity)
    F_ind[, i] <- independent_factor
  }

  # Combine smooth and independent factors
  L_mat <- cbind(L_smooth, L_ind)
  F_mat <- cbind(F_smooth, F_ind)

  # Generate Gaussian noise
  E_mat <- matrix(rnorm(n * p, sd = sigma_E), nrow = n, ncol = p)

  # Compute the resulting matrix
  Y <- L_mat %*% t(F_mat) + E_mat

  return(list(Y = Y, L_mat = L_mat, F_mat = F_mat, E_mat = E_mat, t_vec = t_vec))
}
