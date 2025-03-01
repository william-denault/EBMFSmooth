% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_simulation.R
\name{simulate_matrix}
\alias{simulate_matrix}
\title{Simulate a Random Matrix with Smooth and Independent Factors}
\usage{
simulate_matrix(
  n,
  p,
  k_smooth,
  k_ind,
  sigma_E = 1,
  seed = NULL,
  positive = FALSE
)
}
\arguments{
\item{n}{Number of rows in the generated matrix.}

\item{p}{Number of columns in the generated matrix.}

\item{k_smooth}{Number of smooth factors to include.}

\item{k_ind}{Number of independent factors to include.}

\item{sigma_E}{Standard deviation of the Gaussian noise. Defaults to \code{1}.}

\item{seed}{(Optional) Seed for random number generation. Defaults to \code{NULL}.}

\item{positive}{Logical; if \code{TRUE}, factors and loadings will be non-negative. Defaults to \code{FALSE}.}
}
\value{
A list containing:
\itemize{
  \item \code{Y}: The simulated matrix of size \eqn{n \times p}.
  \item \code{L_mat}: The loading matrix of size \eqn{n \times (k_smooth + k_ind)}.
  \item \code{F_mat}: The factor matrix of size \eqn{p \times (k_smooth + k_ind)}.
  \item \code{E_mat}: The Gaussian noise matrix of size \eqn{n \times p}.
  \item \code{t_vec}: The time vector used for generating smooth loadings.
}
}
\description{
This function generates a random matrix \eqn{Y} of size \eqn{n \times p} as a sum of
\eqn{k_smooth} smooth factors and \eqn{k_ind} independent factors, with sparsity enforced and optional Gaussian noise.
}
\examples{
set.seed(123)
sim <- simulate_matrix(n = 200, p = 100, k_smooth = 2, k_ind = 3, sigma_E = 1)
str(sim)

}
