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
//' @param a numeric; initial boundary
//' @param a_prime numeric; degree of collapse (weibull only)
//' @param kappa numeric; slope of collapse
//' @param tc numeric; time constant of collapse
//' 
//' @return boundary at time t
//' 
//' @name bounds

//' @rdname bounds
//' @export
// [[Rcpp::export]]
arma::vec hyperbolic_ratio_bound(arma::vec t, double a, double kappa=0, double tc=.25){
  // hyperbolic ratio function for time-varying boundary
  // t = time
  // a = initial boundary
  // kappa = degree of collapse
  // tc = time constant of collapse (time at which boundary has collapsed halfway)
  
  return a/2 * (1 - kappa * (t/(t+tc)));
}

//' @rdname bounds
//' @export
// [[Rcpp::export]]
arma::vec weibull_bound(arma::vec t, double a, double a_prime=0, double kappa=1, double tc=.25){
  // weibull time-varying boundary
  // t = time
  // a = initial boundary
  // a_prime = degree of collapse (0 = no collapse; 1 = full collapse)
  // tc = time constant of collapse
  // kappa = slope of collapse
  
  return a/2 - (1 - exp(-1 * arma::pow(t/tc,kappa))) * (a/2 * a_prime);
}
