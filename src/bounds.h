#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec hyperbolic_ratio_bound(arma::vec t, double a, double kappa=0, double tc=.25);
arma::mat hyperbolic_ratio_bound_vec(arma::vec t, arma::vec a, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true);
arma::vec weibull_bound(arma::vec t, double a, double aprime=0, double kappa=1, double tc=.25);
arma::mat weibull_bound_vec(arma::vec t, arma::vec a, arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true);