#include "RcppArmadillo.h"
#include "bounds.h"
#include "sim_pulse.h"
#include "pulse_nll.h"
#include "omp.h"
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp)]]

////////////////////////////////////////////////
/* *** pulse model predictions *** */
////////////////////////////////////////////////

//' Get predicted behavior from pulse model
//'
//' @param n int; number of predicted samples to take per stimulus
//' @param stimuli list; list of stimulus matrices
//' @param v numeric; drift rate, either single value or vector for each trial
//' @param a numeric; initial boundary, either single value or vector for each trial
//' @param t0 numeric; non-decision time, either single value or vector for each trial
//' @param z numeric; starting point, , either single value or vector for each trial, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
//' @param sv numeric; standard deviation of variability in drift rate, either single value or vector for each trial, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time, either single value or vector for each trial. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point, either single value or vector for each trial. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric; standard deviation in wiener diffusion noise, either single value or vector for each trial, default = 1
//' @param lambda numeric; O-U process slope, either single value or vector for each trial
//' @param aprime numeric; degree of collapse, either single value or vector for each trial, default = 0
//' @param kappa numeric; slope of collapse, either single value or vector for each trial, default = 0
//' @param tc numeric; time constant of collapse, either single value or vector for each trial, default = .25
//' @param check_pars logical; if True, check that parameters are vectors of the same length as choices and rts. Must be true if providing scalar parameters. default = true
//' @param dt numeric; time step of simulation, default = .001
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param n_threads int; number of threads (trials) to run in parallel
//'
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), response time
//'
//' @export
// [[Rcpp::export]]
DataFrame pulse_predict(int n, List stimuli,
                        arma::vec v, arma::vec a, arma::vec t0,
                        arma::vec z=0, arma::vec dc=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0,
                        arma::vec lambda=0, arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true,
                        double dt=.001, double dx=.05, int bounds=0, int n_threads=1){
  
  omp_set_num_threads(n_threads);
  
  // check parameter vectors
  unsigned int stim_length = stimuli.length();
  if(check_pars){
    z(arma::find(z==0)).fill(0.5);
    s(arma::find(s==0)).fill(1);
    tc(arma::find(tc==0)).fill(0.25);
    if(v.n_elem != stim_length)
      v = arma::zeros(stim_length)+v(0);
    if(a.n_elem != stim_length)
      a = arma::zeros(stim_length)+a(0);
    if(t0.n_elem != stim_length)
      t0 = arma::zeros(stim_length)+t0(0);
    if(z.n_elem != stim_length)
      z = arma::zeros(stim_length)+z(0);
    if(dc.n_elem != stim_length)
      dc = arma::zeros(stim_length) + dc(0);
    if(sv.n_elem != stim_length)
      sv = arma::zeros(stim_length)+sv(0);
    if(sz.n_elem != stim_length)
      sz = arma::zeros(stim_length)+sz(0);
    if(st0.n_elem != stim_length)
      st0 = arma::zeros(stim_length)+st0(0);
    if(s.n_elem != stim_length)
      s = arma::zeros(stim_length)+s(0);
    if(lambda.n_elem != stim_length)
      lambda = arma::zeros(stim_length)+lambda(0);
    if(aprime.n_elem != stim_length)
      aprime = arma::zeros(stim_length)+aprime(0);
    if(kappa.n_elem != stim_length)
      kappa = arma::zeros(stim_length)+kappa(0);
    if(tc.n_elem != stim_length)
      tc = arma::zeros(stim_length)+tc(0);
  }
  
  arma::ivec choices(n*stim_length), trials(n*stim_length);
  arma::vec rt(n*stim_length);
  
#pragma omp parallel for
  for(unsigned int i=0; i<stim_length; i++){
    
    Rcout << i << "\n";
    
    DataFrame this_sim = sim_pulse(n, stimuli[i], v(i), a(i), t0(i), z(i), dc(i), sv(i), st0(i), sz(i), s(i),
                                   lambda(i), aprime(i), kappa(i), tc(i),
                                   dt, bounds, n_threads);
    
    choices.subvec(i*n,i*n+n-1) = as<arma::ivec>(this_sim[0]);
    rt.subvec(i*n, i*n+n-1) = as<arma::vec>(this_sim[1]);
    trials.subvec(i*n, i*n+n-1) = arma::zeros<arma::ivec>(n) + i;
  }
  
  return DataFrame::create(Named("trial")=trials, Named("response")=choices, Named("rt")=rt);
  
}