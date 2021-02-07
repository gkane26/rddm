#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

DataFrame sim_pulse(int n, arma::ivec up_stimulus, arma::ivec down_stimulus, double v, double a, double t0,
                    double z=.5, double d=0, double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                    double aprime=0, double kappa=0, double tc=.25, 
                    double dt=.001, double v_scale=100, bool use_weibull_bound=false, int n_threads=1);