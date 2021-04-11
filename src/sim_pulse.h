#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

List sim_pulse(int n, arma::mat stimulus, double v, double a, double t0,
               double z=.5, double dc=0, double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
               double aprime=0, double kappa=0, double tc=.25, 
               double v_scale=1, double dt=.001, int bounds=0, int n_threads=1);