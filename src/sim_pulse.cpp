#include "RcppArmadillo.h"
#include "omp.h"
#include "bounds.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


//' Simulate drift diffusion model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param stimulus matrix; stimulus to simulate, 1 or 2 rows X timepoints. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, upper drift = v, lower drift = v-d
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param lambda numeric; O-U process slope
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 1
//' @param tc numeric; time constant of collapse, default = .25
//' @param dt numeric; time step of simulation, default = .001
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' 
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//' 
//' @export
// [[Rcpp::export]]
List sim_pulse(int n, arma::mat stimulus, double v, double a, double t0,
                    double z=.5, double dc=0, double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                    double aprime=0, double kappa=0, double tc=.25, 
                    double dt=.001, int bounds=0, int n_threads=1){
  
  omp_set_num_threads(n_threads);
  int n_on_thread = n / n_threads;
  
  arma::vec rt_full = arma::zeros(n),
    response_full = arma::zeros(n),
    x_full = arma::zeros(n);
  arma::mat accumulators_full(n, stimulus.n_cols+1);
  accumulators_full.fill(arma::datum::nan);
  
  // get boundary vector
  arma::vec tvec = arma::regspace(dt, dt, stimulus.n_cols*dt+dt);
  arma::vec bound;
  if(bounds == 2) {
    bound = weibull_bound(tvec, a, aprime, kappa, tc);
  } else if(bounds == 1) {
    bound = hyperbolic_ratio_bound(tvec, a, kappa, tc);
  } else {
    bound = rep(a/2, stimulus.n_cols); 
  }
  
  if (stimulus.n_rows == 1) {
    stimulus = arma::join_vert(stimulus, arma::zeros(stimulus.n_cols));
  }
  
  double dW = s * sqrt(dt);
  
#pragma omp parallel for
  for (int i=0; i<n_threads; i++){
    
    arma::vec t0_var = t0 + st0 * (arma::randu(n_on_thread)-.5),
      x = a/2*(z-0.5) + sqrt(sz) * arma::randn(n_on_thread),
      v_var = arma::zeros(n_on_thread),
      rt = arma::zeros(n_on_thread);
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
    arma::mat accumulators(n_on_thread, stimulus.n_cols+1);
    accumulators.fill(arma::datum::nan); //= arma::zeros(n_on_thread, stimulus.n_cols+1);
    accumulators.col(0) = x;
    unsigned int step = 0, n_drift = n_on_thread;
    
    while((n_drift > 0) & (step < stimulus.n_cols)){
      v_var = v + sqrt(sv)*arma::randn(n_drift);
      x(still_drift) += dW*arma::randn(n_drift) + ((v_var+dc)*stimulus(0, step) + (-v_var+dc)*stimulus(1, step) + lambda*x(still_drift)) * dt;
      rt(still_drift) += dt;
      still_drift = arma::find(abs(x) < bound(step));
      n_drift = still_drift.n_elem;
      step++;
      accumulators.col(step) = x;
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
  accumulators_full(arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1), arma::span::all) = accumulators;
}

  }
  
  DataFrame sim = DataFrame::create(Named("response") = response_full, Named("rt") = rt_full+t0+st0*(arma::randu(n)-.5), Named("evidence") = x_full);
  return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  
}


//' Simulate drift diffusion model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param stimulus matrix; stimulus to simulate, 1 or 2 rows X timepoints. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
//' @param v numeric vector; drift rate
//' @param a numeric vector; initial boundary
//' @param t0 numeric vector; non-decision time
//' @param z numeric vector; starting point, 0 < z < 1, default = .5
//' @param dc numeric vector; drift criterion, upper drift = v, lower drift = v-d
//' @param sv numeric vector; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric vector; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric vector; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric vector; standard deviation in wiener diffusion noise, default = 1
//' @param lambda numeric vector; O-U process slope
//' @param aprime numeric vector; degree of collapse, default = 0
//' @param kappa numeric vector; slope of collapse, default = 1
//' @param tc numeric vector; time constant of collapse, default = .25
//' @param dt numeric; time step of simulation, default = .001
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param check_pars bool; if True (default) check parameter vector lengths and default values
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' 
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//' 
//' @export
// [[Rcpp::export]]
DataFrame sim_pulse_vec(int n, List stimuli, arma::vec v, arma::vec a, arma::vec t0,
                        arma::vec z=0, arma::vec dc=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0, arma::vec lambda=0,
                        arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0, 
                        double dt=.001, int bounds=0, bool check_pars=true, int n_threads=1){
  
  omp_set_num_threads(n_threads);

  if (check_pars) {
    arma::uvec lens = {stimuli.length(), v.n_elem, a.n_elem, t0.n_elem, z.n_elem, dc.n_elem, sv.n_elem, st0.n_elem, sz.n_elem, s.n_elem, lambda.n_elem, aprime.n_elem, kappa.n_elem, tc.n_elem};
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
    if (s.n_elem < lens.max()) {
      if(s(0) == 0)
        s(0) = 1;
      s = arma::zeros(lens.max()) + s(0);
    }
    if (lambda.n_elem < lens.max()) {
      lambda = arma::zeros(lens.max()) + lambda(0);
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
  
  int total_n = n * stimuli.length();
  arma::vec rt = arma::zeros(total_n),
    response = arma::zeros(total_n),
    evidence = arma::zeros(total_n);
  
#pragma omp parallel for
  for (int i=0; i<stimuli.length(); i++) {
    
    DataFrame sim_results = sim_pulse(n, stimuli[i], v(i), a(i), t0(i),
                                      z(i), dc(i), sv(i), st0(i), sz(i), s(i), lambda(i),
                                      aprime(i), kappa(i), tc(i), 
                                      dt, bounds);
    
    arma::span this_span = arma::span(i*n, i*n+n-1);
    
    response(this_span) = as<arma::vec>(sim_results["response"]);
    rt(this_span) = as<arma::vec>(sim_results["rt"]);
    evidence(this_span) = as<arma::vec>(sim_results["evidence"]);
    
  }
  
  return DataFrame::create(Named("response") = response,
                           Named("rt") = rt,
                           Named("evidence") = evidence);
  
}


/*** NOT IMPLEMENTED ********************
 
 // convert list of stimulus matrices to one large stimulus matrix + a vector of the length of each stimulus matrix
 List stimulus_list_to_matrix() {
 
 }
 
 //' Simulate drift diffusion model with fixed or collapsing boundary
 //'
 //' @param n integer; number of decisions to simulate
 //' @param stimulus matrix; stimulus to simulate, 1 or 2 rows X timepoints. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
 //' @param v numeric vector; drift rate
 //' @param a numeric vector; initial boundary
 //' @param t0 numeric vector; non-decision time
 //' @param z numeric vector; starting point, 0 < z < 1, default = .5
 //' @param dc numeric vector; drift criterion, upper drift = v, lower drift = v-d
 //' @param sv numeric vector; standard deviation of variability in drift rate, sv >= 0, default = 0
 //' @param st0 numeric vector; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
 //' @param sz numeric vector; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
 //' @param s numeric vector; standard deviation in wiener diffusion noise, default = 1
 //' @param lambda numeric vector; O-U process slope
 //' @param aprime numeric vector; degree of collapse, default = 0
 //' @param kappa numeric vector; slope of collapse, default = 1
 //' @param tc numeric vector; time constant of collapse, default = .25
 //' @param dt numeric; time step of simulation, default = .001
 //' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
 //' @param check_pars bool; if True (default) check parameter vector lengths and default values
 //' @param n_threads integer; number of threads to run in parallel, default = 1
 //' 
 //' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
 //' 
 //' @export
 // [[Rcpp::export]]
 DataFrame sim_pulse_vec(int n, List stimuli, arma::vec v, arma::vec a, arma::vec t0,
 arma::vec z=0, arma::vec dc=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0, arma::vec lambda=0,
 arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0, 
 double dt=.001, int bounds=0, bool check_pars=true, int n_threads=1){
 
 omp_set_num_threads(n_threads);
 int n_on_thread = n / n_threads;
 
 if (check_pars) {
 arma::uvec lens = {stimuli.length(), v.n_elem, a.n_elem, t0.n_elem, z.n_elem, dc.n_elem, sv.n_elem, st0.n_elem, sz.n_elem, s.n_elem, lambda.n_elem, aprime.n_elem, kappa.n_elem, tc.n_elem};
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
 if (s.n_elem < lens.max()) {
 if(s(0) == 0)
 s(0) = 1;
 s = arma::zeros(lens.max()) + s(0);
 }
 if (lambda.n_elem < lens.max()) {
 lambda = arma::zeros(lens.max()) + lambda(0);
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
 
 arma::vec rt_full = arma::zeros(n),
 response_full = arma::zeros(n),
 x_full = arma::zeros(n);
 
 // get boundary vector
 
 arma::vec tvec = arma::regspace(dt, dt, stimulus.n_cols*dt+dt);
 arma::mat bound(tvec.n_elem, a.n_elem);
 if(bounds == 2) {
 bound = weibull_bound_vec(tvec, a, aprime, kappa, tc);
 } else if(bounds == 1) {
 bound = hyperbolic_ratio_bound_vec(tvec, a, kappa, tc);
 } else {
 for (unsigned int i=0; i < a.n_elem; i++) {
 bound.col(i) = arma::zeros(tvec.n_elem) + a(i)/2;
 }
 }
 
 // if (stimulus.n_rows == 1) {
 //   stimulus = arma::join_vert(stimulus, arma::zeros(stimulus.n_cols));
 // }
 
 arma::vec dW = s * sqrt(dt);
 
#pragma omp parallel for
 for (int i=0; i<n_threads; i++){
 
 for (unsigned int s=i*n_on_thread; s<i*n_on_thread+n_on_thread; s++) {
 if (stimuli[s].nrow() == 1) {
 stimuli[s] = rbind(stimuli[s], rep(0, stimuli[s].ncol()));
 }
 }
 
 arma::span thread_span = arma::span(i*n_on_thread, i*n_on_thread+n_on_thread-1);
 
 Rcout << "vector lengths\n";
 Rcout << t0.n_elem << " "<< st0.n_elem << " " << n_on_thread << "\n";
 
 arma::vec t0_var = t0(thread_span) + st0(thread_span) % (arma::randu(n_on_thread)-.5);
 
 Rcout << "added t0\n";
 
 arma::vec x = a(thread_span) / 2 % (z(thread_span) - 0.5) + sqrt(sz(thread_span)) % arma::randn(n_on_thread),
 v_sub = v(thread_span),
 sv_sub = sv(thread_span),
 v_var = arma::zeros(n_on_thread),
 dc_sub = dc(thread_span),
 rt = arma::zeros(n_on_thread);
 
 Rcout << "sub bounds!\n";
 arma::mat bnd = bound.cols(i*n_on_thread, i*n_on_thread+n_on_thread-1);
 arma::mat bnd_t = bnd.t();
 Rcout << "sub bounds complete!\n";
 
 arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
 unsigned int step = 0, n_drift = n_on_thread;
 
 while((n_drift > 0) & (step < stimulus.n_cols)){
 
 lapply(stimuli[still_drift], )
 
 v_var = v_sub(still_drift) + sqrt(sv_sub(still_drift))*arma::randn(n_drift);
 x(still_drift) += dW%arma::randn(n_drift) + ((v_var+dc_sub(still_drift))*stimulus(0, step) + (-v_var+dc_sub(still_drift))*stimulus(1, step) + lambda*x(still_drift)) * dt;
 rt(still_drift) += dt;
 still_drift = arma::find(abs(x) < bnd_t.col(step));
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
 
 ***********************   **/