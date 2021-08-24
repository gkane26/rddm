#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec linear_urgency(arma::vec t, double uslope=0, double udelay=0, double umag=0);
arma::vec logistic_urgency(arma::vec t, double uslope=0, double udelay=0, double umag=0);