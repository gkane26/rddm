#include <RcppArmadillo.h>
#include <RcppZiggurat.h>
#include <Ziggurat.h>

// [[Rcpp::depends(RcppArmadillo, RcppZiggurat)]]

arma::vec zrandn(int n);
void zrandseed(unsigned long int s);