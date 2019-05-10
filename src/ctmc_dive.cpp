// Author: Richard Glennie
// Date: Feb 2019 
// Continuous-time Markov chain with SPDE 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace R_inla;
  using namespace density; 
  using namespace Eigen;  
  
  DATA_VECTOR(ldive); // log-dive durations 
  DATA_VECTOR(lsurf); // log-surface durations 
  DATA_MATRIX(Xdive); // dive design matrix 
  DATA_MATRIX(Xsurf); // surface design matrix 
  DATA_SPARSE_MATRIX(A_dive); // projector for dives 
  DATA_SPARSE_MATRIX(A_surf); // projector for surface
  DATA_SPARSE_MATRIX(C); // FEM <psi_i, psi_j> matrix
  DATA_SPARSE_MATRIX(G1); // FEM <dpsi_i, psi_j> matrix 
  DATA_SPARSE_MATRIX(G2); // FEM G2 = G1C^{-1}G1 
  
  PARAMETER_VECTOR(par_dive); // dive parameters 
  PARAMETER_VECTOR(par_surf); // surface parameters 
  PARAMETER(log_sigma_surf); // surface variability 
  PARAMETER(log_sigma_dive); // dive variability 
  PARAMETER(cor_surfdive); // covariance of surface and dive 
  PARAMETER(log_kappa_dive); // dive correlation parameter 
  PARAMETER(log_kappa_surf); // dive correlation parameter 
  PARAMETER(log_tau_dive); // dive rf variance 
  PARAMETER(log_tau_surf); // surface rf variance 
  PARAMETER_VECTOR(s_dive); // dive random field 
  PARAMETER_VECTOR(s_surf); // surface random field 
  
  Type sigma_dive = exp(log_sigma_dive); 
  Type sigma_surf = exp(log_sigma_surf); 
  Type kappa_dive = exp(log_kappa_dive); 
  Type kappa_surf = exp(log_kappa_surf); 
  Type tau_dive = exp(log_tau_dive); 
  Type tau_surf = exp(log_tau_surf);
  SparseMatrix<Type> Q_dive = pow(kappa_dive, 4) * C + 2.0 * pow(kappa_dive, 2) * G1 + G2; 
  SparseMatrix<Type> Q_surf = pow(kappa_surf, 4) * C + 2.0 * pow(kappa_surf, 2) * G1 + G2; 
  
  vector<Type> e_dive = A_dive * s_dive; 
  vector<Type> e_surf = A_surf * s_surf;  
  
  Type nll = 0 ; 
  
  nll += SCALE(GMRF(Q_dive), 1 / tau_dive)(s_dive); 
  nll += SCALE(GMRF(Q_surf), 1 / tau_surf)(s_surf); 
   
  vector<Type> eta_dive = Xdive * par_dive + e_dive;
  vector<Type> eta_surf = Xsurf * par_surf + e_surf;
  matrix<Type> Sigma(2, 2); 
  Sigma(0, 0) = sigma_dive; 
  Sigma(1, 1) = sigma_surf; 
  Sigma(0, 1) = cor_surfdive; 
  Sigma(1, 0) = cor_surfdive; 

  for(int i = 0; i < ldive.size(); i++) {
    vector<Type> dat(2); 
    dat(0) = ldive(i) - eta_dive(i); 
    dat(1) = lsurf(i) - eta_surf(i); 
    nll += MVNORM(Sigma)(dat); 
    //nll -= dnorm(lsurf(i), eta_surf(i), sigma_surf, true); 
    //nll -= dnorm(ldive(i), eta_dive(i), sigma_dive, true); 
  }

  double nu = 2 - 1.0 / 2; // nu = alpha - d / 2
  Type mu_dive = e_dive.mean(); 
  Type mu_surf = e_surf.mean(); 
  ADREPORT(mu_dive);
  ADREPORT(mu_surf); 
  return nll;
}
