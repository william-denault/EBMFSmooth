% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_utility.R
\name{flash_both_types_sequential}
\alias{flash_both_types_sequential}
\title{Flash Both Types Sequential: Fit Smooth and Independent Factors Alternately}
\usage{
flash_both_types_sequential(
  data,
  S = NULL,
  k_total = 10,
  ebnm_smooth,
  ebnm_indep = NULL,
  var_type = NULL,
  verbose = 0,
  backfitting_steps = 0,
  convergence_fn = flashier::flash_conv_crit_max_chg,
  convergence_tol = 1e-05
)
}
\arguments{
\item{data}{A matrix of observed data.}

\item{S}{(Optional) A matrix of standard errors corresponding to the data in \code{data}. Defaults to \code{NULL}.}

\item{k_total}{An integer specifying the total number of factors to fit. Defaults to \code{10}.}

\item{ebnm_smooth}{A function for empirical Bayes smoothing of factors.}

\item{ebnm_indep}{(Optional) A function for empirical Bayes estimation of independent factors. Defaults to \code{ebnm::ebnm_point_normal}.}

\item{var_type}{(Optional) An integer specifying the variance type. Defaults to \code{NULL}.}

\item{verbose}{An integer controlling verbosity. Set to \code{0} for no messages, or higher values for more detailed messages. Defaults to \code{0}.}

\item{backfitting_steps}{An integer specifying the number of backfitting iterations. If \code{0}, backfitting is skipped. Defaults to \code{0}.}

\item{convergence_fn}{A function defining the convergence criterion. Defaults to \code{flashier::flash_conv_crit_max_chg}.}

\item{convergence_tol}{A numeric value specifying the convergence tolerance for smooth factor fitting. Defaults to \code{1e-4}.}
}
\value{
A list containing:
\itemize{
  \item \code{fit}: A `flashier` fit object containing the fitted factors and loadings.
  \item \code{type_factor}: A character vector indicating the type of each fitted factor ("smooth" or "independent").
}
}
\description{
This function fits smooth and independent factors sequentially to a dataset
using the `flashier` package. It alternates between fitting smooth and independent factors,
selecting the one with the better ELBO (evidence lower bound optimization) at each step.
Backfitting can optionally be performed to refine the final solution.
}
\details{
This function performs the following steps:
1. Alternates between fitting smooth and independent factors using greedy optimization.
   - At each step, the factor type with the better ELBO is retained.
   - Verbose messages indicate the progress and results of each step.
2. Optionally performs backfitting to refine the solution if \code{backfitting_steps > 0}.
}
\examples{
library(flashier)
library(ebnm)
data <- matrix(rnorm(100), nrow = 10)
result <- flash_both_types_sequential(
  data = data,
  k_total = 5,
  ebnm_smooth = ebnm::ebnm_point_normal,
  verbose = 1
)
print(result$type_factor)

}
