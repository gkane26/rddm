#include "RcppArmadillo.h"
#include "dqrng.h"
#include "RcppZiggurat.h"
#include "omp.h"
#include "bounds.h"
#include "urgency.h"
#include "fast_rand.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo, dqrng, RcppZiggurat)]]
// [[Rcpp::plugins(openmp)]]


//' Simulate pulse diffusion model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param stimuli cube; stimuli to simulate, 2 rows X timepoints x trials.
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param s numeric; standard deviation in wiener diffusion noise
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, upper drift = v, lower drift = v-d
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param lambda numeric; O-U process slope
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 1
//' @param tc numeric; time constant of collapse, default = .25
//' @param uslope numeric; urgency scaling factor, default = 0;
//' @param umag numeric; urgency magnitude, default = 0;
//' @param udelay numeric; urgency delay, default = 0;
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
//' @param dt numeric; time step of simulation, default = .001
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
//' @param urgency int: 0 for none, 1 for linear, 2 for logistic
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' @param return_accu bool; if True, return full trajectory of accumulators
//'
//' @return List containing 1) data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence and 2) matrix with full accumulator trajectories
//' 
//' @export
// [[Rcpp::export]]
List sim_pulse(int n, arma::cube stimuli, double a, double t0, double s,
               double z=.5, double dc=0,
               double sv=0, double st0=0, double sz=0,
               double lambda=0, double aprime=0, double kappa=0, double tc=.25, 
               double uslope=0, double umag=0, double udelay=0,
               double v_scale=1, double dt=.001, int bounds=0, int urgency=0,
               int n_threads=1, bool return_accu=false){
  
  int n_total = n * stimuli.n_slices;
  n_threads = std::min(n_threads, n_total);
  int n_on_thread = n_total / n_threads;
  
  sv *= v_scale;
  dc *= v_scale;
  
  arma::uvec trial_vector = arma::regspace<arma::uvec>(0, stimuli.n_slices-1);
  trial_vector = arma::repmat(trial_vector, n, 1);
  
  arma::vec rt_full = arma::zeros(n_total),
    response_full = arma::zeros(n_total),
    x_full = arma::zeros(n_total);
  arma::mat accumulators_full(n_total, stimuli.n_cols+1);
  if (return_accu) {
    accumulators_full.fill(arma::datum::nan);
  }
  
  // get boundary vector
  arma::vec tvec = arma::regspace(dt, dt, stimuli.n_cols*dt+dt),
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
  
  double dW = s * sqrt(dt);
  
#pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<n_threads; i++){
    
    arma::uvec thread_trials = trial_vector.rows(i*n_on_thread, i*n_on_thread+n_on_thread-1);

    arma::vec t0_var = t0 + st0 * (arma::randu(n_on_thread)-.5),
      x = a/2*(z-0.5) + sqrt(sz) * zrandn(n_on_thread),
      v_var = arma::zeros(n_on_thread),
      rt = arma::zeros(n_on_thread);
    
    arma::uvec still_drift = arma::linspace<arma::uvec>(0, n_on_thread-1, n_on_thread);
    
    arma::mat accumulators(n_on_thread, stimuli.n_cols+1);
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0) = x;
    }
    
    unsigned int step = 0;
    
    while((still_drift.n_elem > 0) & (step < stimuli.n_cols)){
      
      arma::uvec trials_drifting = thread_trials(still_drift);
      arma::vec stim0 = stimuli.tube(0, step),
        stim1 = stimuli.tube(1, step);
      stim0 = stim0(trials_drifting);
      stim1 = stim1(trials_drifting);
      
      // v_var(still_drift) = 1 + sqrt(sv)*arma::randn(still_drift.n_elem);
      v_var(still_drift) = v_scale + sqrt(sv) * zrandn(still_drift.n_elem);
      
      // x(still_drift) += (gamma(step) * dW * arma::randn(still_drift.n_elem)) +
      //   (gamma(step) * ((v_var(still_drift) + dc) % stim0 * dt)) +
      //   (gamma(step) * ((-v_var(still_drift) + dc) % stim1 * dt)) +
      //   (lambda * x(still_drift) * dt);
      x(still_drift) += (gamma(step) * dW * zrandn(still_drift.n_elem)) +
        (gamma(step) * ((v_var(still_drift) + dc) % stim0 * dt)) +
        (gamma(step) * ((-v_var(still_drift) + dc) % stim1 * dt)) +
        (lambda * x(still_drift) * dt);

      rt(still_drift) += dt;
      still_drift = arma::find(abs(x) < bound(step));
      step++;
      if (return_accu) {
        accumulators.col(step) = x;
      }
      
    }
    
    double final_bound = bound(bound.n_elem-1);
    arma::vec response(n_on_thread);
    response.fill(arma::datum::nan);
    response(arma::find(x >= final_bound)).fill(1);
    response(arma::find(x <= -final_bound)).fill(0);
    rt(arma::find_nonfinite(response)).fill(arma::datum::nan);
    
    arma::span thread_span = arma::span(i * n_on_thread, i * n_on_thread + n_on_thread - 1);
    rt_full(thread_span) = rt;
    response_full(thread_span) = response;
    x_full(thread_span) = x;
    
    if (return_accu) {
      accumulators_full(thread_span, arma::span::all) = accumulators;
    }
    
  }
  
  DataFrame sim = DataFrame::create(Named("trial") = trial_vector,
                                    Named("response") = response_full,
                                    Named("rt") = rt_full+t0+st0*(arma::randu(n_total)-.5),
                                    Named("evidence") = x_full);
  if (return_accu) {
    return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  } else {
    return List::create(Named("behavior") = sim);
  }
  
}


/* ***
 //' Simulate drift diffusion model with fixed or collapsing boundary
 //'
 //' @param n integer; number of decisions to simulate
 //' @param stimuli list of matrices; each element is a stimulus to simulate, 1 or 2 rows X timepoints. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
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
 //' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
 //' @param dt numeric; time step of simulation, default = .001
 //' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
 //' @param check_pars bool; if True (default) check parameter vector lengths and default values
 //' @param n_threads integer; number of threads to run in parallel, default = 1
 //' @param return_accu bool; if True, return full trajectory of accumulators
 //'
 //' @return List containing 1) data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence and 2) matrix with full accumulator trajectories
 //'
 //' @export
 // [[Rcpp::export]]
 List sim_pulse_vec(int n, List stimuli, arma::vec v, arma::vec a, arma::vec t0,
 arma::vec z=0, arma::vec dc=0,
 arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0,
 arma::vec lambda=0, arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0,
 double v_scale=1, double dt=.001, int bounds=0,
 bool check_pars=true, int n_threads=1, bool return_accu=false){
 
 v *= v_scale;
 sv *= v_scale;
 dc *= v_scale;
 
 if (check_pars) {
 arma::uvec lens = {stimuli.length(), v.n_elem, a.n_elem, t0.n_elem, z.n_elem,
 dc.n_elem, sv.n_elem, st0.n_elem, sz.n_elem, s.n_elem,
 lambda.n_elem, aprime.n_elem, kappa.n_elem, tc.n_elem};
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
 arma::vec trial_full = arma::zeros(total_n),
 rt_full = arma::zeros(total_n),
 response_full = arma::zeros(total_n),
 x_full = arma::zeros(total_n);
 
 arma::field<arma::mat> stimuli_field = as<arma::field<arma::mat>>(stimuli);
 int max_stimulus_len = 0;
 for (int i=0; i<stimuli_field.n_elem; i++) {
 int this_cols = stimuli_field(i).n_cols;
 max_stimulus_len = std::max(max_stimulus_len, this_cols);
 }
 
 arma::mat accumulators_full(total_n, max_stimulus_len+1);
 if (return_accu) {
 accumulators_full.fill(arma::datum::nan);
 }
 
#pragma omp parallel for num_threads(n_threads)
 for (int i=0; i<stimuli_field.n_elem; i++) {
 
 arma::mat stimulus = stimuli_field(i);
 
 // get boundary vector
 arma::vec tvec = arma::regspace(dt, dt, stimulus.n_cols*dt+dt);
 arma::vec bound;
 if(bounds == 2) {
 bound = weibull_bound(tvec, a(i), aprime(i), kappa(i), tc(i));
 } else if(bounds == 1) {
 bound = hyperbolic_ratio_bound(tvec, a(i), kappa(i), tc(i));
 } else {
 bound = rep(a(i)/2, stimulus.n_cols);
 }
 
 double dW = s(i) * sqrt(dt);
 
 arma::vec t0_var = t0(i) + st0(i) * (arma::randu(n)-.5),
 x = a(i)/2*(z(i)-0.5) + sqrt(sz(i)) * arma::randn(n),
 v_var = arma::zeros(n),
 rt = arma::zeros(n);
 arma::uvec still_drift = arma::linspace<arma::uvec>(0, n-1, n);
 unsigned int step = 0, n_drift = n;
 arma::mat accumulators(n, stimulus.n_cols+1);
 
 if (return_accu) {
 accumulators.fill(arma::datum::nan);
 accumulators.col(0) = x;
 }
 
 while((n_drift > 0) & (step < stimulus.n_cols)){
 v_var(still_drift) = v(i) + sqrt(sv(i)) * arma::randn(n_drift);
 x(still_drift) += dW*arma::randn(n_drift) + ((v_var(still_drift)+dc(i))*stimulus(0, step) + (-v_var(still_drift)+dc(i))*stimulus(1, step) + lambda(i)*x(still_drift)) * dt;
 rt(still_drift) += dt;
 still_drift = arma::find(abs(x) < bound(step));
 n_drift = still_drift.n_elem;
 step++;
 if (return_accu) {
 accumulators.col(step) = x;
 }
 }
 
 double final_bound = bound(bound.n_elem-1);
 arma::vec response = -arma::ones(n);
 response(arma::find(x >= final_bound)).fill(1);
 response(arma::find(x <= -final_bound)).fill(0);
 
 arma::span this_span = arma::span(i*n, i*n+n-1);
 trial_full(this_span).fill(i+1);
 rt_full(this_span) = rt + t0(i) + st0(i) * (arma::randu(n)-0.5);
 response_full(this_span) = response;
 x_full(this_span) = x;
 if (return_accu) {
 accumulators_full(this_span, arma::span(0, stimulus.n_cols)) = accumulators;
 }
 
 }
 
 DataFrame sim = DataFrame::create(Named("trial") = trial_full,
 Named("response") = response_full,
 Named("rt") = rt_full,
 Named("evidence") = x_full);
 
 if (return_accu) {
 return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
 } else {
 return List::create(Named("behavior") = sim);
 }
 
 }
 
 *** */