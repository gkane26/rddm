#include "RcppArmadillo.h"
#include "RcppZiggurat.h"
#include "Ziggurat.h"
// [[Rcpp::depends(RcppArmadillo, RcppZiggurat)]]

using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

//' random normal
//'
//' @description random normal
//'
//' @param n int; number of random draws
//'
//' @return random numbers
//'
//' @export
// [[Rcpp::export]]
arma::vec zrandn(int n) {
  arma::vec x(n);
  for(int i=0; i<n; i++) {
    x(i) = zigg.norm();
  }
  return x;
}
