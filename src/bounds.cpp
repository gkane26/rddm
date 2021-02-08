#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//' Collapsing boundary functions
//' 
//' @description evaluate diffusion model boundary using the hyperbolic ratio or weibull functions
//'
//' @param t vector; time points to evaluate boundary
//' @param a numeric or vector; initial boundary
//' @param aprime numeric or vector; degree of collapse (weibull only)
//' @param kappa numeric or vector; slope of collapse
//' @param tc numeric or vector; time constant of collapse
//' 
//' @return 
//' column vector with boundary at time t, 
//' or a time point x parameter vector matrix of boundaries 
//' (each column represents a time varying boundary for a given parameter set)
//' 
//' @name bounds

//' @rdname bounds
//' @export
// [[Rcpp::export]]
arma::vec hyperbolic_ratio_bound(arma::vec t, double a, double kappa=0, double tc=.25){
  return a / 2 * (1 - kappa * (t / (t + tc)));
}

//' @rdname bounds
//' @export
// [[Rcpp::export]]
arma::mat hyperbolic_ratio_bound_vec(arma::vec t, arma::vec a, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true){

  if (check_pars) {
    arma::uvec lens = {a.n_elem, kappa.n_elem, tc.n_elem};
    if (a.n_elem < lens.max()) {
      a = arma::zeros(lens.max()) + a(0);  
    }
    if (kappa.n_elem < lens.max()) {
      if (kappa(0) == 0) {
        kappa(0) = 1;
      }
      kappa = arma::zeros(lens.max()) + kappa(0);
    }
    if (tc.n_elem < lens.max()) {
      if (tc(0) == 0) {
        tc(0) = 0.25;
      }
      tc = arma::zeros(lens.max()) + tc(0);
    }
  }
  
  arma::mat bmat(t.n_elem, a.n_elem);
  for (unsigned int i=0; i < a.n_elem; i++) {
    bmat.col(i) = a(i) / 2 * (1 - kappa(i) * (t / (t + tc(i))));
  }

  return bmat;

}

//' @rdname bounds
//' @export
// [[Rcpp::export]]
arma::vec weibull_bound(arma::vec t, double a, double aprime=0, double kappa=1, double tc=.25){
  return a/2 - (1 - exp(-1 * arma::pow(t/tc,kappa))) * (a/2 * aprime);
}

//' @rdname bounds
//' @export
// [[Rcpp::export]]
arma::mat weibull_bound_vec(arma::vec t, arma::vec a, arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true){

  if (check_pars) {
    arma::uvec lens = {a.n_elem, aprime.n_elem, kappa.n_elem, tc.n_elem};
    if (a.n_elem < lens.max()) {
      a = arma::zeros(lens.max()) + a(0);  
    }
    if (aprime.n_elem < lens.max()) {
      aprime = arma::zeros(lens.max()) + aprime(0);
    }
    if (kappa.n_elem < lens.max()) {
      if (kappa(0) == 0) {
        kappa(0) = 1;
      }
      kappa = arma::zeros(lens.max()) + kappa(0);
    }
    if (tc.n_elem < lens.max()) {
      if (tc(0) == 0) {
        tc(0) = 0.25;
      }
      tc = arma::zeros(lens.max()) + tc(0);
    }
  }
  
  arma::mat bmat(t.n_elem, a.n_elem);
  for (unsigned int i=0; i < a.n_elem; i++) {
    bmat.col(i) = a(i) / 2 - (1 - exp(-1 * arma::pow(t / tc(i), kappa(i)))) * (a(i) / 2 * aprime(i));
  }

  return bmat;
}
