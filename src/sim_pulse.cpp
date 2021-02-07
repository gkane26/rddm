#include "RcppArmadillo.h"
#include "omp.h"
#include "bounds.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


//' Simulate drift diffusion model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param up_stimulus float vector; stimulus to the positive boundary 
//' @param down_stimulus float vector; stimulus to the negative boundary
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param d numeric; drift criterion, upper drift = v, lower drift = v-d
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param lambda numeric; O-U process slope
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 1
//' @param tc numeric; time constant of collapse, default = .25
//' @param dt numeric; time step of simulation, default = .001
//' @param v_scale numeric; scale drift rate to be similar to boundary separation a, default = 100
//' @param use_weibull_bound logical; if True, use weibull function for collapsing bounds, if False, use hyperbolic ratio function
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' 
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//' 
//' @export
// [[Rcpp::export]]
DataFrame sim_pulse(int n, arma::ivec up_stimulus, arma::ivec down_stimulus, double v, double a, double t0,
                    double z=.5, double d=0, double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                    double aprime=0, double kappa=0, double tc=.25, 
                    double dt=.001, double v_scale=100, bool use_weibull_bound=false, int n_threads=1){
  
  omp_set_num_threads(n_threads);
  int n_on_thread = n / n_threads;
  
  arma::vec rt_full = arma::zeros(n),
    response_full = arma::zeros(n),
    x_full = arma::zeros(n);
  
  // get bounds
  arma::vec bound;
  if(use_weibull_bound)
    bound = weibull_bound(arma::regspace(dt, dt, up_stimulus.n_elem*dt), a, aprime, kappa, tc);
  else
    bound = hyperbolic_ratio_bound(arma::regspace(dt, dt, up_stimulus.n_elem*dt), a, kappa, tc);
  
  if (down_stimulus.n_elem != up_stimulus.n_elem) {
    if (down_stimulus.n_elem == 1) {
      down_stimulus = arma::zeros<arma::ivec>(up_stimulus.n_elem) + down_stimulus[0];
    }
  }
  
#pragma omp parallel for
  for (int i=0; i<n_threads; i++){
    
    arma::vec t0_var = t0 + st0 * (arma::randu(n_on_thread)-.5),
      x = (2*z-1)*a/2 + sqrt(sz) * arma::randn(n_on_thread),
      v_var = arma::zeros(n_on_thread),
      d_var = arma::zeros(n_on_thread),
      rt = arma::zeros(n_on_thread);
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
    unsigned int step = 0, n_drift = n_on_thread;
    
    while((n_drift > 0) & (step < up_stimulus.n_elem)){
      v_var = v + sqrt(sv)*arma::randn(n_drift);
      d_var = d + sqrt(sv)*arma::randn(n_drift);
      x(still_drift) += sqrt(s*dt)*arma::randn(n_drift) + up_stimulus[step]*v_var*dt - down_stimulus[step]*(v_var-d_var)*dt + lambda*x(still_drift)*dt;
      rt(still_drift) += dt;
      still_drift = arma::find(abs(x) < bound(step));
      n_drift = still_drift.n_elem;
      step++;
    }
    
    double final_bound = bound(bound.n_elem-1);
    arma::vec response = -arma::ones(n_on_thread);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    
#pragma omp critical
{
  rt_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = rt;
  response_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = response;
  x_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = x;
}

  }
  
  DataFrame sim = DataFrame::create(Named("response") = response_full, Named("rt") = rt_full+t0+st0*(arma::randu(n)-.5), Named("evidence") = x_full);
  return sim;
}