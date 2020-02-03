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
  DATA_MATRIX(Xs_grid_surface); // projector for integration (fixed effects)
  DATA_MATRIX(Xs_grid_dive); // projector for integration (fixed effects)
  DATA_MATRIX(A_grid_surface); // projector for integration (random effects)
  DATA_MATRIX(A_grid_dive); // projector for integration (random effects)
  DATA_SPARSE_MATRIX(indD); // integration points within surfacings
  DATA_SPARSE_MATRIX(indS); // integration point within dives
  DATA_SCALAR(dt); // time step in integration
  DATA_IVECTOR(S_dive_n); // sizes of sparse matrices
  DATA_IVECTOR(S_surface_n); // sizes of sparse matrices

  PARAMETER_VECTOR(par_dive); // dive parameters
  PARAMETER_VECTOR(par_surf); // surface parameters
  PARAMETER_VECTOR(log_lambda_dive); // dive log smoothing parameter
  PARAMETER_VECTOR(log_lambda_surf); // surface log smoothing parameter
  PARAMETER_VECTOR(s_dive); // dive random effects
  PARAMETER_VECTOR(s_surf); // surface random effects

  vector<Type> lambda_dive = exp(log_lambda_dive); // dive smoothing parameter
  vector<Type> lambda_surf = exp(log_lambda_surf); // surface smoothing parameter

  // Negative log-likelihood is nll
  Type nll = 0;

  // add smoothing penalties, need to do this wonky bit to unblock S in each case
  // data setup
  int Sn = 0;
  int S_start = 0;
  vector<Type> this_beta = s_dive;
  vector<Type> this_s;
  SparseMatrix<Type> this_S;

  // dive bit
  for(int i = 0; i < S_dive_n.size(); i++) {
    Sn = S_dive_n(i);
    this_S = S_dive.block(S_start, S_start, Sn, Sn);
    this_s = s_dive.segment(S_start, Sn);
    nll -= Type(0.5) * Type(Sn) * log_lambda_dive(i) - 0.5 * lambda_dive(i) * GMRF(this_S).Quadform(this_s);
    S_start += Sn;
  }

  // surface bit
  Sn = 0;
  S_start = 0;
  for(int i = 0; i < S_surface_n.size(); i++) {
    Sn = S_surface_n(i);
    this_S = S_surface.block(S_start, S_start, Sn, Sn);
    this_s = s_surf.segment(S_start, Sn);
    nll -= Type(0.5) * Type(Sn) * log_lambda_surf(i) - 0.5 * lambda_surf(i) * GMRF(this_S).Quadform(this_s);
    S_start += Sn;
  }

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
  vector<Type> int_dive = A_grid_dive * s_dive + Xs_grid_dive * par_dive;
  int_dive = exp(int_dive);
  int_dive = indD * int_dive;
  vector<Type> subint_dive = int_dive.head(int_dive.size() - 1);
  subint_dive *= dt;

  // Integral of surface intensity
  vector<Type> int_surf = A_grid_surface * s_surf + Xs_grid_surface * par_surf;
  int_surf = exp(int_surf);
  int_surf = indS * int_surf;
  int_surf *= dt;

  // Likelihood contributions
  for(int i = 0; i < e_dive.size(); i++) {
    nll -= le_dive(i) + le_surf(i) - int_surf(i);
    if (i < e_dive.size() - 1) nll -= -subint_dive(i);
  }

  return nll;
}
