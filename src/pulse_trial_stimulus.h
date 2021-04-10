#include "RcppArmadillo.h"
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

arma::vec pulse_trial_stimulus(arma::vec stim_seq, double dur=0.01, double isi=0.1, double pre_stim=0, double dt=0.001);