#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec hyperbolic_ratio_bound(arma::vec t, double a, double kappa=0, double tc=.25);
arma::vec weibull_bound(arma::vec t, double a, double a_prime=0, double kappa=1, double tc=.25);