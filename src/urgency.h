#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec linear_urgency(arma::vec t, double uslope=0, double udelay=0, double umag=0);
arma::mat linear_urgency_vec(arma::vec t, arma::vec uslope=0, arma::vec udelay=0, arma::vec umag=0, bool check_pars=true);
arma::vec logistic_urgency(arma::vec t, double uslope=0, double udelay=0, double umag=0);
arma::mat logistic_urgency_vec(arma::vec t, arma::vec uslope=0, arma::vec udelay=0, arma::vec umag=0, bool check_pars=true);
  