#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x); //response variable
  DATA_VECTOR(s); //standard error
  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(B); // Design matrix (for random effects)
  DATA_SPARSE_MATRIX(P); // Penalty matrix

  int d = P.cols(); // Number of Spline coefficients

  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  PARAMETER_VECTOR(beta); // beta estimates
  PARAMETER(theta);

  // Parameter
  PARAMETER_VECTOR(U); // eta = B * U + X * beta
  Type sigmaIWP = exp(-0.5*theta);

  // Transformations
  vector<Type> eta = X * beta + B * U;

  // Log likelihood
  Type ll = 0;
  ll = sum(dnorm(x, eta, s, TRUE));
  REPORT(ll);

  // Log prior on U
  Type lpW = 0;

  // Cross product
  if(sigmaIWP != 0){
    vector<Type> PU = P*U;
    Type UPU = (U * PU).sum();
    lpW += -0.5 * exp(-2*log(sigmaIWP)) * UPU; // U part
  }

  // Log determinant
  if(sigmaIWP != 0){
    Type logdet1 = d * (-2*log(sigmaIWP)) + logPdet;
    lpW += 0.5 * logdet1; // P part
    // wrt the dimension
    lpW += -0.5 * d * log(2*M_PI);
  }

  REPORT(lpW);

  // Final result!
  Type logpost = -1 * (ll + lpW);

  return logpost;
}

