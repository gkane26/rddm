#include <RcppArmadillo.h>
#include <omp.h>
#include "bounds.h"
#include "urgency.h"
#include "fast_rand.h"
#include <chrono>
#include <thread>
using namespace Rcpp;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds
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
//' @param sz numeric; variability in starting point (percentage of boundary a). Uniform from [z-sz/2, z+sz/2], 0 < sz < 1, default = 0
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
//' @param seed integer; set random seed
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
             int bounds=0, int urgency=0, int n_threads=1, bool return_accu=false, int seed=-1) {
  
  if (seed >= 0) {
    zrandseed(seed);
  }
  
  int thread_leftover = n % n_threads;
  arma::uvec n_on_thread(n_threads);
  n_on_thread.fill(n / n_threads);
  if (thread_leftover > 0) {
    n_on_thread.rows(0, thread_leftover-1) += 1;
  }
  arma::uvec thread_starts = arma::cumsum(n_on_thread);
  
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
  
#pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<n_threads; i++){
    
    int thread_start = thread_starts(i) - n_on_thread(i);
    int thread_end = thread_starts(i) - 1;
    
    double step = 0,
      t = dt;
    
    arma::vec v_var = v + sv * zrandn(n_on_thread(i)) + dc,
      x = ((2 * z - 1) * a / 2) + (sz * a * (arma::randu(n_on_thread(i)) - 0.5)),
      rt = arma::zeros(n_on_thread(i));
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread(i)-1, n_on_thread(i));
    
    arma::mat accumulators(n_on_thread(i), tvec.n_elem+1);
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0) = x;
    }
    
    while((still_drift.n_elem > 0) & (t < max_time)){
      x(still_drift) +=  gamma(step) * (v_var(still_drift) * dt + dW * zrandn(still_drift.n_elem));
      rt(still_drift) += dt;
      still_drift = arma::find(arma::abs(x) < bound(step));
      t+=dt;
      step++;
      
      if (return_accu) {
        accumulators.col(step) = x;
      }
    }
    
    double final_bound = bound(bound.n_elem-1);
    arma::vec response(n_on_thread(i));
    response.fill(arma::datum::nan);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    
    arma::span thread_span = arma::span(thread_start, thread_end);
    rt_full(thread_span) = rt;
    response_full(thread_span) = response;
    rt(arma::find_nonfinite(response)).fill(arma::datum::nan);
    
    if (return_accu) {
      accumulators_full(thread_span, arma::span::all) = accumulators;
    }
  }
  
  arma::vec t0_vec = t0 + (t0 * st0 * (arma::randu(n) - 0.5));
  t0_vec.clamp(0, arma::datum::inf);
  
  DataFrame sim = DataFrame::create(Named("response")=response_full,
                                    Named("rt")=rt_full + t0_vec);
  
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
//' @param sz vector; variability in starting point (percentage of boundary a). Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param aprime vector; degree of collapse, default = 0
//' @param kappa vector; slope of collapse, default = 0
//' @param tc vector; time constant of collapse, default = .25
//' @param uslope numeric; urgency scaling factor, default = 0;
//' @param umag numeric; urgency magnitude, default = 0;
//' @param udelay numeric; urgency delay, default = 0;
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param dt numeric; time step of simulation, default = .001
//' @param max_time numeric; max time of simulation, default = 10
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param urgency int: 0 for none, 1 for linear, 2 for logistic
//' @param check_pars bool; if True (default) check parameter vector lengths and default values
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' @param return_accu bool; if True, return full trajectory of accumulators
//' @param seed integer; set random seed
//'
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
//'
//' @export
// [[Rcpp::export]]
List sim_ddm_vec(arma::vec v, arma::vec a, arma::vec t0,
                 arma::vec z=0, arma::vec dc=0,
                 arma::vec sv=0, arma::vec st0=0, arma::vec sz=0,
                 arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0,
                 arma::vec uslope=0, arma::vec udelay=0, arma::vec umag=0,
                 arma::vec s=0, double dt=.001, double max_time=10,
                 int bounds=0, int urgency=0, bool check_pars=true, int n_threads=1, bool return_accu=false, int seed=-1) {
  
  if (seed >= 0) {
    zrandseed(seed);
  }
  
  if (check_pars) {
    arma::uvec lens = {v.n_elem, a.n_elem, t0.n_elem,
                       z.n_elem, dc.n_elem,
                       sv.n_elem, st0.n_elem, sz.n_elem,
                       aprime.n_elem, kappa.n_elem, tc.n_elem,
                       uslope.n_elem, udelay.n_elem, umag.n_elem,
                       s.n_elem};
    
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
    if (uslope.n_elem < lens.max()) {
      uslope = arma::zeros(lens.max()) + uslope(0);
    }
    if (udelay.n_elem < lens.max()) {
      udelay = arma::zeros(lens.max()) + udelay(0);
    }
    if (umag.n_elem < lens.max()) {
      umag = arma::zeros(lens.max()) + umag(0);
    }
    if (s.n_elem < lens.max()) {
      if (s(0) == 0)
        s(0) = 1;
      s = arma::zeros(lens.max()) + s(0);
    }
  }
  
  int thread_leftover = v.n_elem % n_threads;
  arma::uvec n_on_thread(n_threads);
  n_on_thread.fill(v.n_elem / n_threads);
  if (thread_leftover > 0) {
    n_on_thread.rows(0, thread_leftover-1) += 1;
  }
  arma::uvec thread_cum = arma::cumsum(n_on_thread);
  
  arma::vec rt_full = arma::zeros(v.n_elem),
    response_full = arma::zeros(v.n_elem);
  
  // get boundary
  arma::vec tvec = arma::regspace(dt, dt, max_time);
  arma::mat bound(tvec.n_elem, v.n_elem);
  if (bounds == 1){
    bound = hyperbolic_ratio_bound_vec(tvec, a, kappa, tc);
  }else if (bounds == 2){
    bound = weibull_bound_vec(tvec, a, aprime, kappa, tc);
  } else if (bounds == 3) {
    bound = linear_bound_vec(tvec, a, kappa, tc);
  }else {
    for (unsigned int i = 0; i < v.n_elem; i++) {
      bound.col(i).fill(a(i)/2);
    }
  }
  
  // get urgency signal
  arma::mat gamma(tvec.n_elem, v.n_elem);
  if (urgency == 1) {
    gamma = linear_urgency_vec(tvec, uslope, udelay, umag);
  } else if (urgency == 2) {
    gamma = logistic_urgency_vec(tvec, uslope, udelay, umag);
  } else {
    gamma.fill(1);
  }
  
  arma::mat  accumulators_full(v.n_elem, tvec.n_elem+1);
  if (return_accu) {
    accumulators_full.fill(arma::datum::nan);
  }
  
#pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<n_threads; i++) {
    
    int thread_start = thread_cum(i) - n_on_thread(i);
    int thread_end = thread_cum(i) - 1;
    arma::span thread_span = arma::span(thread_start, thread_end);
    
    double step = 0,
      t = dt;
    
    arma::vec v_var = v(thread_span) + sv(thread_span) % zrandn(n_on_thread(i)) + dc(thread_span),
      x = ((2 * z(thread_span) - 1) * a(thread_span) / 2) + (sz(thread_span) * a(thread_span) * (arma::randu(n_on_thread(i)) - 0.5)),
      dW = s(thread_span) * sqrt(dt),
      rt = arma::zeros(n_on_thread(i));
    
    arma::uvec still_drift = arma::regspace<arma::uvec>(0, n_on_thread(i)-1);
    arma::mat gamma_thread = gamma.cols(thread_start, thread_end);
    arma::mat gamma_thread_t = gamma_thread.t();
    arma::mat bnd = bound.cols(thread_start, thread_end);
    arma::mat bnd_t = bnd.t();
    
    arma::mat accumulators(n_on_thread(i), tvec.n_elem+1);
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0) = x;
    }
    
    while((still_drift.n_elem > 0) & (t < max_time)){
      
      arma::vec gamma_t = gamma_thread_t.col(step);
      x(still_drift) += gamma_t(still_drift) % (v_var(still_drift) * dt + dW(still_drift) % zrandn(still_drift.n_elem));
      rt(still_drift) += dt;
      still_drift = arma::find(bnd_t.col(step) - arma::abs(x) > 0);
      t+=dt;
      step++;
      
      if (return_accu) {
        accumulators.col(step) = x;
      }
      
    }
    
    double final_bound = bound(bound.n_rows-1);
    arma::vec response(n_on_thread(i));
    response.fill(arma::datum::nan);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    rt(arma::find_nonfinite(response)).fill(arma::datum::nan);
    
    rt_full(thread_span) = rt;
    response_full(thread_span) = response;
    
    if (return_accu) {
      accumulators_full(thread_span, arma::span::all) = accumulators;
    }
    
  }
  
  arma::vec t0_vec = t0 + (t0 * st0 % (arma::randu(v.n_elem) - 0.5));
  t0_vec.clamp(0, arma::datum::inf);
  
  DataFrame sim = DataFrame::create(Named("response")=response_full,
                                    Named("rt")=rt_full + t0_vec);
  
  if (return_accu) {
    return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  } else {
    return List::create(Named("behavior") = sim);
  }
  
}


