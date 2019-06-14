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
//' @param stim_seq vector of strings; vector of stimulus sequence for all trials, each element should be a string of 0s and 1s (0 for evidence to lower, 1 for evidence to upper)
//' @param v numeric; drift rate, either single value or vector for each trial
//' @param a numeric; initial boundary, either single value or vector for each trial
//' @param t0 numeric; non-decision time, either single value or vector for each trial
//' @param z numeric; starting point, , either single value or vector for each trial, 0 < z < 1, default = .5
//' @param sv numeric; standard deviation of variability in drift rate, either single value or vector for each trial, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time, either single value or vector for each trial. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point, either single value or vector for each trial. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric; standard deviation in wiener diffusion noise, either single value or vector for each trial, default = 1
//' @param lambda numeric; O-U process slope, either single value or vector for each trial
//' @param a_prime numeric; degree of collapse, either single value or vector for each trial, default = 0
//' @param kappa numeric; slope of collapse, either single value or vector for each trial, default = 0
//' @param tc numeric; time constant of collapse, either single value or vector for each trial, default = .25
//' @param check_pars logical; if True, check that parameters are vectors of the same length as choices and rts. Must be true if providing scalar parameters. default = true
//' @param dt numeric; time step of simulation, default = .001
//' @param v_scale numeric; scale drift rate to be similar to boundary separation a, default = 100
//' @param use_weibull_bound logical; if True, use weibull function for collapsing bounds, if False, use hyperbolic ratio function
//' @param dur numeric; duration of stimulus
//' @param isi numeric; interstimulus interval
//' @param n_threads int; number of threads (trials) to run in parallel
//'
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), response time
//'
//' @export
// [[Rcpp::export]]
DataFrame pulse_predict(int n, std::vector<std::string> stim_seq,
                        arma::vec v, arma::vec a, arma::vec t0,
                        arma::vec z=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0,
                        arma::vec lambda=0, arma::vec a_prime=0, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true,
                        double dt=.002, double dx=.05, double v_scale=100, bool use_weibull_bound=false, double dur=.01, double isi=.1, int n_threads=1){
  
  // check parameter vectors
  int stim_length = stim_seq.size();
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
      z = arma::zeros(stim_length)+.5;
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
    if(a_prime.n_elem != stim_length)
      a_prime = arma::zeros(stim_length)+a_prime(0);
    if(kappa.n_elem != stim_length)
      kappa = arma::zeros(stim_length)+kappa(0);
    if(tc.n_elem != stim_length)
      tc = arma::zeros(stim_length)+tc(0);
  }
  
  arma::ivec choices(n*stim_length), trials(n*stim_length);
  arma::vec rt(n*stim_length);
  
  for(unsigned int i=0; i<stim_length; i++){
    arma::ivec this_stim = get_stimulus(stim_seq[i], dur, isi);
    DataFrame this_sim = sim_pulse(n, this_stim, v(i), a(i), t0(i), z(i), sv(i), st0(i), sz(i), s(i),
                                   lambda(i), a_prime(i), kappa(i), tc(i),
                                   dt, v_scale, use_weibull_bound, n_threads);
    
    choices.subvec(i*n,i*n+n-1) = as<arma::ivec>(this_sim[0]);
    rt.subvec(i*n, i*n+n-1) = as<arma::vec>(this_sim[1]);
    trials.subvec(i*n, i*n+n-1) = arma::zeros<arma::ivec>(n) + i;
  }
  
  return DataFrame::create(Named("trial")=trials, Named("response")=choices, Named("rt")=rt);
}