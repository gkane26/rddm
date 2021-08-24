#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//' Urgency functions
//' 
//' @description evaluate diffusion model urgency signals (linear and logistic)
//'
//' @param t vector; time points to evaluate boundary
//' @param uslope; slope of urgency signal
//' @param udelay; delay to start of rising urgency
//' @param umag; magnitude of urgency signal (only for logistic)
//' 
//' @return column vector with urgency at time t, 
//' 
//' @name urgency

//' @rdname urgency
//' @export
// [[Rcpp::export]]
arma::vec linear_urgency(arma::vec t, double uslope=0, double udelay=0, double umag=0){
  arma::vec gamma = uslope * (t - udelay) + 1;
  gamma.clamp(1, arma::datum::inf);
  return gamma;
}

//' @rdname urgency
//' @export
// [[Rcpp::export]]
arma::vec logistic_urgency(arma::vec t, double uslope=0, double udelay=0, double umag=0){
  arma::vec s1 = arma::exp(uslope * (t - udelay));
  double s2 = exp(-uslope * udelay);
  arma::vec gamma = (umag*s1 / (1 + s1)) + (1 + (1-umag)*s2) / (1 + s2);
  return gamma;
}
