#include "RcppArmadillo.h"
#include "omp.h"
#include "bounds.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


//' Simulate drift diffusion model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param a_prime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 0
//' @param tc numeric; time constant of collapse, default = .25
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param dt numeric; time step of simulation, default = .001
//' @param max_time numeric; max time of simulation, default = 10
//' @param use_weibull_bound logical; if True, use weibull function for collapsing bounds, if False, use hyperbolic ratio function
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' 
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
//' 
//' @export
// [[Rcpp::export]]
DataFrame sim_ddm(int n, double v, double a, double t0, double z=.5,
                  double sv=0, double st0=0, double sz=0,
                  double a_prime=0, double kappa=0, double tc=.25,
                  double s=1, double dt=.001, double max_time=10, bool use_weibull_bound=false, int n_threads=1){
  
  omp_set_num_threads(n_threads);
  int n_on_thread = n / n_threads;
  
  arma::vec rt_full = arma::zeros(n),
    response_full = arma::zeros(n);
  
  arma::vec tvec = arma::regspace(dt, dt, max_time),
    bound;
  if(use_weibull_bound){
    bound = weibull_bound(tvec, a, a_prime, kappa, tc);
  }else{
    bound = hyperbolic_ratio_bound(tvec, a, kappa, tc);
  }
  
  double dW = s*sqrt(dt);
  max_time = max_time - t0;
  
#pragma omp parallel for
  for (int i=0; i<n_threads; i++){
    
    double step = 0,
      t = dt;
    
    arma::vec v_var = v + sv * arma::randn(n_on_thread),
      z_var = z + sz*(arma::randu(n_on_thread)-.5),
      x = a/2 * (2*z_var-1),
      rt = arma::zeros(n_on_thread);
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
    
    while((still_drift.n_elem > 0) & (t < max_time)){
      x(still_drift) +=  v_var(still_drift)*dt + dW*arma::randn(still_drift.size());
      rt(still_drift) += dt;
      still_drift = arma::find(arma::abs(x) < bound(step));
      t+=dt;
      step++;
    }
    
    double final_bound = bound(bound.n_elem-1);
    arma::vec response = -arma::ones(n_on_thread);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    
    rt_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = rt;
    response_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = response;

  }
  
  DataFrame sim = DataFrame::create(Named("response")=response_full, Named("rt")=rt_full+t0+st0*(arma::randu(n)-.5));
  return sim;
}

