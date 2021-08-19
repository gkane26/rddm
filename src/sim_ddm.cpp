#include "RcppArmadillo.h"
#include "omp.h"
#include "bounds.h"
#include "urgency.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


//' Simulate drift diffusion model with fixed or collapsing boundary, with or without urgency
//'
//' @param n integer; number of decisions to simulate
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 0
//' @param tc numeric; time constant of collapse, default = .25
//' @param uslope numeric; urgency scaling factor, default = 0;
//' @param umag numeric; urgency magnitude, default = 0;
//' @param udelay numeric; urgency delay, default = 0;
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param dt numeric; time step of simulation, default = .001
//' @param max_time numeric; max time of simulation, default = 10
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
//' @param urgency int: 0 for none, 1 for linear, 2 for logistic
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' @param return_accu bool; if True, return full trajectory of accumulators
//'
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
//' 
//' @export
// [[Rcpp::export]]
List sim_ddm(int n, double v, double a, double t0, double z=.5, double dc=0,
             double sv=0, double st0=0, double sz=0,
             double aprime=0, double kappa=0, double tc=.25,
             double uslope=0, double umag=0, double udelay=0,
             double s=1, double dt=.001, double max_time=10,
             int bounds=0, int urgency=0, int n_threads=1, bool return_accu=false){
  
  omp_set_num_threads(n_threads);
  int n_on_thread = n / n_threads;
  
  arma::vec rt_full = arma::zeros(n),
    response_full = arma::zeros(n);
  
  arma::vec tvec = arma::regspace(dt, dt, max_time),
    bound;
  if (bounds == 1) {
    bound = hyperbolic_ratio_bound(tvec, a, kappa, tc);
  } else if (bounds == 2) {
    bound = weibull_bound(tvec, a, aprime, kappa, tc);
  } else if (bounds == 3) {
    bound = linear_bound(tvec, a, kappa, tc);
  } else {
    bound = rep(a/2, tvec.n_elem);
  }
  
  // get urgency signal
  arma::vec gamma(tvec.n_elem);
  if (urgency == 1) {
    gamma = linear_urgency(tvec, uslope, udelay, umag);
  } else if (urgency == 2) {
    gamma = logistic_urgency(tvec, uslope, udelay, umag);
  } else {
    gamma.fill(1);
  }
  
  arma::mat  accumulators_full(n, tvec.n_elem+1);
  if (return_accu) {
    accumulators_full.fill(arma::datum::nan);
  }
  
  double dW = s*sqrt(dt);
  
#pragma omp parallel for
  for (int i=0; i<n_threads; i++){
    
    double step = 0,
      t = dt;
    
    arma::vec v_var = v + sv * arma::randn(n_on_thread) + dc,
      z_var = z + sz*(arma::randu(n_on_thread)-.5),
      x = a/2 * (z_var - 0.5),
      rt = arma::zeros(n_on_thread);
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
    
    arma::mat accumulators(n_on_thread, tvec.n_elem+1);
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0) = x;
    }
    
    while((still_drift.n_elem > 0) & (t < max_time)){
      x(still_drift) +=  gamma(step) * (v_var(still_drift)*dt + dW*arma::randn(still_drift.size()));
      rt(still_drift) += dt;
      still_drift = arma::find(arma::abs(x) < bound(step));
      t+=dt;
      step++;
      
      if (return_accu) {
        accumulators.col(step) = x;
      }
    }
    
    double final_bound = bound(bound.n_elem-1);
    arma::vec response = -arma::ones(n_on_thread);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    
    rt_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = rt;
    response_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = response;
    
    if (return_accu) {
      accumulators_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1), arma::span::all) = accumulators;
    }
  }
  
  DataFrame sim = DataFrame::create(Named("response")=response_full, Named("rt")=rt_full+t0+st0*(arma::randu(n)-.5));
  if (return_accu) {
    return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  } else {
    return List::create(Named("behavior") = sim);
  }
  
}

//' Simulate drift diffusion model with fixed or collapsing boundary
//'
//' @param v vector; drift rate
//' @param a vector; initial boundary
//' @param t0 vector; non-decision time
//' @param z vector; starting point, 0 < z < 1, default = .5
//' @param dc vector; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
//' @param sv vector; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 vector; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz vector; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param aprime vector; degree of collapse, default = 0
//' @param kappa vector; slope of collapse, default = 0
//' @param tc vector; time constant of collapse, default = .25
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param dt numeric; time step of simulation, default = .001
//' @param max_time numeric; max time of simulation, default = 10
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param check_pars bool; if True (default) check parameter vector lengths and default values
//' @param n_threads integer; number of threads to run in parallel, default = 1
//'
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
//'
//' @export
// [[Rcpp::export]]
DataFrame sim_ddm_vec(arma::vec v, arma::vec a, arma::vec t0, arma::vec z=0, arma::vec dc=0,
                      arma::vec sv=0, arma::vec st0=0, arma::vec sz=0,
                      arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0,
                      double s=1, double dt=.001, double max_time=10,
                      int bounds=0, bool check_pars=true, int n_threads=1){
  
  omp_set_num_threads(n_threads);
  
  if (check_pars) {
    arma::uvec lens = {v.n_elem, a.n_elem, t0.n_elem, z.n_elem, dc.n_elem, sv.n_elem, st0.n_elem, sz.n_elem, aprime.n_elem, kappa.n_elem, tc.n_elem};
    if (v.n_elem < lens.max()) {
      v = arma::zeros(lens.max()) + v(0);
    }
    if (a.n_elem < lens.max()) {
      a = arma::zeros(lens.max()) + a(0);
    }
    if (t0.n_elem < lens.max()) {
      t0 = arma::zeros(lens.max()) + t0(0);
    }
    if (z.n_elem < lens.max()) {
      if (z(0) == 0)
        z(0) = 0.5;
      z = arma::zeros(lens.max()) + z(0);
    }
    if (dc.n_elem < lens.max()) {
      dc = arma::zeros(lens.max()) + dc(0);
    }
    if (sv.n_elem < lens.max()) {
      sv = arma::zeros(lens.max()) + sv(0);
    }
    if (st0.n_elem < lens.max()) {
      st0 = arma::zeros(lens.max()) + st0(0);
    }
    if (sz.n_elem < lens.max()) {
      sz = arma::zeros(lens.max()) + sz(0);
    }
    if (aprime.n_elem < lens.max()) {
      aprime = arma::zeros(lens.max()) + aprime(0);
    }
    if (kappa.n_elem < lens.max()) {
      kappa = arma::zeros(lens.max()) + kappa(0);
    }
    if (tc.n_elem < lens.max()) {
      if (tc(0) == 0)
        tc(0) = 0.25;
      tc = arma::zeros(lens.max()) + tc(0);
    }
  }
  
  int n_on_thread = v.n_elem / n_threads;
  
  arma::vec rt_full = arma::zeros(v.n_elem),
    response_full = arma::zeros(v.n_elem);
  
  arma::vec tvec = arma::regspace(dt, dt, max_time);
  arma::mat bound(tvec.n_elem, a.n_elem);
  if(bounds == 1){
    bound = hyperbolic_ratio_bound_vec(tvec, a, kappa, tc);
  }else if(bounds == 2){
    bound = weibull_bound_vec(tvec, a, aprime, kappa, tc);
  } else if(bounds == 3) {
    bound = linear_bound_vec(tvec, a, kappa, tc);
  }else {
    bound.each_col() = arma::zeros(tvec.n_elem) + a/2;
  }
  
  double dW = s*sqrt(dt);
  
#pragma omp parallel for
  for (int i=0; i<n_threads; i++) {
    
    double step = 0,
      t = dt;
    
    arma::span thread_span = arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1);
    
    arma::vec v_var = v(thread_span) + sv(thread_span) % arma::randn(n_on_thread) + dc(thread_span),
      z_var = z(thread_span) + sz(thread_span) % (arma::randu(n_on_thread) - 0.5),
      x = a(thread_span) / 2 % (z_var(thread_span) - 0.5),
      rt = arma::zeros(n_on_thread);
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
    arma::mat bnd = bound.cols(i*n_on_thread, i*n_on_thread+n_on_thread-1);
    arma::mat bnd_t = bnd.t();
    
    while((still_drift.n_elem > 0) & (t < max_time)){
      x(still_drift) +=  v_var(still_drift)*dt + dW*arma::randn(still_drift.size());
      rt(still_drift) += dt;
      still_drift = arma::find(arma::abs(x) < bnd_t.col(step));
      t+=dt;
      step++;
    }
    
    double final_bound = bound(bound.n_rows-1);
    arma::vec response = -arma::ones(n_on_thread);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    
    rt_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = rt;
    response_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1)) = response;
    
  }
  
  DataFrame sim = DataFrame::create(Named("response")=response_full,
                                    Named("rt")=rt_full+t0+st0 % (arma::randu(v.n_elem)-.5));
  
  return sim;
  
}


