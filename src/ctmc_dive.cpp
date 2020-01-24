// Author: Richard Glennie
// Date: Dec 2019
// Continuous-time Markov chain with SPDE 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_MATRIX(Xdive); // dive design matrix
  DATA_MATRIX(Xsurf); // surface design matrix
  DATA_SPARSE_MATRIX(S_dive); // smoothing matrix
  DATA_SPARSE_MATRIX(S_surface); // smoothing matrix
  DATA_MATRIX(A_dive); // projector for dive times
  DATA_MATRIX(A_surf); // projector for surface times
  DATA_MATRIX(A_grid); // projector for integration
  DATA_SPARSE_MATRIX(indD); // integration points within surfacings
  DATA_SPARSE_MATRIX(indS); // integration point within dives
  DATA_SCALAR(dt); // time step in integration

  PARAMETER_VECTOR(par_dive); // dive parameters
  PARAMETER_VECTOR(par_surf); // surface parameters
  PARAMETER(log_lambda_dive); // dive log smoothing parameter
  PARAMETER(log_lambda_surf); // surface log smoothing parameter
  PARAMETER_VECTOR(s_dive); // dive random effects
  PARAMETER_VECTOR(s_surf); // surface random effects

  Type lambda_dive = exp(log_lambda_dive); // dive smoothing parameter
  Type lambda_surf = exp(log_lambda_surf); // surface smoothing parameter

  // Negative log-likelihood is nll
  Type nll = 0;

  // add smoothing penalties
  nll -= Type(0.5) * S_dive.cols() * log_lambda_dive - 0.5 * lambda_dive * GMRF(S_dive).Quadform(s_dive);
  nll -= Type(0.5) * S_surface.cols() * log_lambda_surf - 0.5 * lambda_surf * GMRF(S_surface).Quadform(s_surf);

  // Linear predictors with and without random effect
  vector<Type> leta_dive = Xdive * par_dive;
  vector<Type> le_dive = leta_dive + A_dive * s_dive;
  vector<Type> leta_surf = Xsurf * par_surf;
  vector<Type> le_surf  = leta_surf + A_surf * s_surf;
  vector<Type> e_dive = exp(le_dive);
  vector<Type> e_surf = exp(le_surf);
  vector<Type> eta_dive = exp(leta_dive);
  vector<Type> eta_surf = exp(leta_surf);

  // Integral of dive intensity
  vector<Type> int_dive = A_grid * s_dive;
  int_dive = exp(int_dive);
  int_dive = indD * int_dive;
  vector<Type> subint_dive = int_dive.head(int_dive.size() - 1);
  vector<Type> subeta_dive = eta_dive.head(eta_dive.size() - 1);
  subint_dive = (subint_dive.array() * subeta_dive.array()).matrix();
  subint_dive *= dt;

  // Integral of surface intensity
  vector<Type> int_surf = A_grid * s_surf;
  int_surf = exp(int_surf);
  int_surf = indS * int_surf;
  int_surf = (int_surf.array() * eta_surf.array()).matrix();
  int_surf *= dt; 

  // Likelihood contributions
  for(int i = 0; i < e_dive.size(); i++) {
    nll -= le_dive(i) + le_surf(i) - int_surf(i);
    if (i < e_dive.size() - 1) nll -= -subint_dive(i);
  }

  return nll;
}
