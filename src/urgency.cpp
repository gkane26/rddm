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
arma::mat linear_urgency_vec(arma::vec t, arma::vec uslope=0, arma::vec udelay=0, arma::vec umag=0, bool check_pars=true){
  
  if(check_pars) {
    arma::uvec lens = {uslope.n_elem, udelay.n_elem, umag.n_elem};
    if (uslope.n_elem < lens.max()) {
      uslope = arma::zeros(lens.max()) + uslope(0);  
    }
    if (udelay.n_elem < lens.max()) {
      udelay = arma::zeros(lens.max()) + udelay(0);
    }
    if (umag.n_elem < lens.max()) {
      umag = arma::zeros(lens.max()) + umag(0);
    }
  }
  
  arma::mat gamma(t.n_elem, uslope.n_elem);
  for (unsigned int i=0; i<uslope.n_elem; i++) {
    gamma.col(i) = uslope(i) * (t - udelay(i)) + 1;
  }
  
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

//' @rdname urgency
//' @export
// [[Rcpp::export]]
arma::mat logistic_urgency_vec(arma::vec t, arma::vec uslope=0, arma::vec udelay=0, arma::vec umag=0, bool check_pars=true){
  
  if(check_pars) {
    arma::uvec lens = {uslope.n_elem, udelay.n_elem, umag.n_elem};
    if (uslope.n_elem < lens.max()) {
      uslope = arma::zeros(lens.max()) + uslope(0);  
    }
    if (udelay.n_elem < lens.max()) {
      udelay = arma::zeros(lens.max()) + udelay(0);
    }
    if (umag.n_elem < lens.max()) {
      umag = arma::zeros(lens.max()) + umag(0);
    }
  }
  
  arma::mat gamma(t.n_elem, uslope.n_elem);
  for (unsigned int i=0; i<uslope.n_elem; i++) {
    arma::vec s1 = arma::exp(uslope(i) * (t - udelay(i)));
    double s2 = exp(-uslope(i) * udelay(i));
    gamma.col(i) = (umag(i)*s1 / (1 + s1)) + (1 + (1-umag(i))*s2) / (1 + s2);
  }
  
  return gamma;
  
}
