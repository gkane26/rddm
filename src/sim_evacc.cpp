#include "RcppArmadillo.h"
#include "RcppZiggurat.h"
#include "omp.h"
#include "bounds.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppZiggurat)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Simulate EvAcc (2-accumulator) model with fixed or collapsing boundary
//'
//' @param n integer; number of simulations for each stimulus
//' @param stimulus array; stimulus to simulate, 1 or 2 rows X timepoints x number of stimuli. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; boundary separation bias, left boundary = a - z, right = a + z
//' @param dc numeric; drift criterion, left = v + dc, right = v - dc
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting bias (drawn from gaussian), sz > 0, default = 0
//' @param s numeric; diffusion noise
//' @param lambda numeric; accumulator recurrent weight
//' @param mi numeric; mutual inhibition
//' @param sv2 numeric; sampling noise for second layer
//' @param s2 numeric; diffusion noise for layer 2
//' @param lambda2 numeric; second layer accumulator recurrent weight
//' @param mi2 numeric; mutual inhibition in second layer. Set mi=0 and mi2 > 0 to implement static sample model from Scott et al., 2015
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 1
//' @param tc numeric; time constant of collapse, default = .25
//' @param v numeric; drift rate. Default = 1, recommended to leave this parameter fixed
//' @param two_layer bool; if true, use two layer accumulator model
//' @param accumulator_gain; if true, gain of first accumulator layer scales with value
//' @param scalar_stimulus_noise bool; if true, drift rate variability scales with accumulator value (see Koay et al., 2020)
//' @param scalar_diffusion_noise bool; if true, diffusion noise scales with accumulator value
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv and criterion dc are multiplied by this number
//' @param dt numeric; time step of simulation, default = .001
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' @param return_accu bool; if True, return full trajectory of accumulators
//'
//' @return List containing 1) data frame with three columns: response (1 for right, 0 for left), response time, and evidence and 2) matrix with full accumulator trajectories
//' 
//' @export
// [[Rcpp::export]]
List sim_evacc(int n, arma::cube stimuli,
                double a, double t0,
                double z=0, double dc=0,
                double sv=0, double st0=0, double sz=0, double s=1,
                double lambda=0, double mi=0,
                double sv2=0, double s2=1,
                double lambda2=1, double mi2=0,
                double aprime=0, double kappa=0, double tc=.25,
                double v=1,
                bool accumulator_gain=false, bool two_layer=false,
                bool scalar_stimulus_noise=false, bool scalar_diffusion_noise=false,
                double v_scale=1, double dt=.001, int bounds=0,
                int n_threads=1, bool return_accu=false){

  int n_total = n * stimuli.n_slices,
    n_on_thread = n_total / n_threads;
  
  arma::uvec trial_vector = arma::regspace<arma::uvec>(0, stimuli.n_slices-1);
  trial_vector = arma::repmat(trial_vector, n, 1);
  
  v *= v_scale;
  sv *= v_scale;
  dc *= v_scale;
  
  arma::vec
    rt_full(n_total),
    response_full(n_total);
  rt_full.fill(arma::datum::nan);
  response_full.fill(arma::datum::nan);
  arma::cube accumulators_full(2, stimuli.n_cols+1, n_total);
  if (return_accu) {
    accumulators_full.fill(arma::datum::nan);
  }
  
  // get time-varying boundary vector
  arma::vec tvec = arma::regspace(dt, dt, stimuli.n_cols*dt+dt);
  arma::vec bound;
  if(bounds == 2) {
    bound = weibull_bound(tvec, 2*a, aprime, kappa, tc);
  } else if(bounds == 1) {
    bound = hyperbolic_ratio_bound(tvec, 2*a, kappa, tc);
  } else {
    bound = rep(a, stimuli.n_cols); 
  }
  
  double dt_sr = sqrt(dt);

#pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<n_threads; i++) {
    
    arma::uvec thread_trials = trial_vector.rows(i*n_on_thread, i*n_on_thread+n_on_thread-1);
    
    arma::mat
      x(n_on_thread, 2, arma::fill::zeros),
      y(n_on_thread, 2, arma::fill::zeros),
      gain(n_on_thread, 2, arma::fill::ones),
      sv_mat(n_on_thread, 2),
      dw_mat(n_on_thread, 2),
      dw2_mat(n_on_thread, 2),
      sv2_mat(n_on_thread, 2);
    sv_mat.fill(sv);
    dw_mat.fill(s * dt_sr);
    dw2_mat.fill(s2 * dt_sr);
    sv2_mat.fill(sv2);
    
    arma::vec
      z_var = z + sqrt(sz) * RcppZiggurat::zrnorm(n_on_thread),
        response,
        rt;
    
    arma::cube accumulators(2, stimuli.n_cols+1, n_on_thread),
    accumulators_done;
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0).fill(0);
    }
    
    double step = 0;
    
    while ((x.n_rows > 0) & (step < stimuli.n_cols)) {
      
      if (scalar_stimulus_noise) {
        sv_mat = sv * x;
        sv2_mat = sv2 * y;
      }

      if (scalar_diffusion_noise) {
        dw_mat = s * x * dt_sr;
        dw2_mat = s2 * y * dt_sr;
      }
      
      if (accumulator_gain) {
        gain = x;
      }
      
      arma::mat x_temp = x;
      arma::vec stim0 = stimuli.tube(0, step);
      stim0 = stim0(thread_trials);
      arma::vec stim1;
      if (stimuli.n_rows == 1) {
        stim1 = arma::zeros(stim0.n_elem);
      } else{
        stim1 = stimuli.tube(1, step);
        stim1 = stim1(thread_trials);
      }
      
      x.col(0) += gain.col(0) % 
        (((v + dc + sv_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))) % stim0 * dt) -
        (lambda * x_temp.col(0) * dt) -
        (mi * x_temp.col(1) * dt) +
        (dw_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))));
      x.col(1) += gain.col(1) % 
        (((v - dc + sv_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))) % stim1 * dt) -
        (lambda * x_temp.col(1) * dt) -
        (mi * x_temp.col(0) * dt) +
        (dw_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))));
      x(find(x < 0)).fill(0);
      
      if (two_layer) {
        arma::mat y_temp = y;
        y.col(0) += (x_temp.col(0) + sv2_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)) * dt) - 
          (lambda2 * y_temp.col(0) * dt) - 
          (mi2 * y_temp.col(1) * dt) +
          (dw2_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(y_temp.n_rows)));
        y.col(1) += (x_temp.col(1) + sv2_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(y_temp.n_rows)) * dt) - 
          (lambda2 * y_temp.col(0) * dt) - 
          (mi2 * y_temp.col(0) * dt) +
          (dw2_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(y_temp.n_rows)));
        y(find(y < 0)).fill(0);
      } else {
        y = x;  
      }
      
      if (return_accu) {
        accumulators.col(step+1) = y.t();
      }
      
      arma::uvec lefts = arma::find(y.col(0) >= bound(step) - z_var);
      response = arma::join_cols(response, arma::zeros(lefts.n_elem));
      arma::vec temp_rt_left = rep((step + 1) * dt, lefts.n_elem);
      rt = arma::join_cols(rt, temp_rt_left);
      x.shed_rows(lefts);
      y.shed_rows(lefts);
      gain.shed_rows(lefts);
      sv_mat.shed_rows(lefts);
      dw_mat.shed_rows(lefts);
      sv2_mat.shed_rows(lefts);
      z_var.shed_rows(lefts);
      thread_trials.shed_rows(lefts);
      if (return_accu) {
        accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(lefts));
        accumulators.shed_slices(lefts);
      }
      
      arma::uvec rights = arma::find(y.col(1) >= bound(step) + z_var);
      response = arma::join_cols(response, arma::ones(rights.n_elem));
      arma::vec temp_rt_right = rep((step + 1) * dt, rights.n_elem);
      rt = arma::join_cols(rt, temp_rt_right);
      x.shed_rows(rights);
      y.shed_rows(rights);
      gain.shed_rows(rights);
      sv_mat.shed_rows(rights);
      dw_mat.shed_rows(rights);
      sv2_mat.shed_rows(rights);
      z_var.shed_rows(rights);
      thread_trials.shed_rows(rights);
      if (return_accu) {
        accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(rights));
        accumulators.shed_slices(rights);
      }
      
      step++;
      
    }
    
    if (rt.n_elem > 0) {
      rt_full(arma::span(i*n_on_thread, i*n_on_thread+rt.n_elem-1)) = rt + t0 + st0 * (arma::randu(rt.n_elem) - 0.5);
    }
    if (response.n_elem > 0) {
      response_full(arma::span(i*n_on_thread, i*n_on_thread+response.n_elem-1)) = response;
    }
    
    if (return_accu) {
      accumulators_done = arma::join_slices(accumulators_done, accumulators);
      accumulators_full(arma::span(), arma::span(), arma::span(i*n_on_thread, i*n_on_thread+accumulators_done.n_slices-1)) = accumulators_done;
    }
    
  }
  
  DataFrame sim = DataFrame::create(
    Named("trial") = trial_vector,
    Named("response") = response_full,
    Named("rt") = rt_full);
  
  if (return_accu) {
    return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  } else {
    return List::create(Named("behavior") = sim);
  }
  
}

/*** EVERYTHING BELOW COMMENTED OUT ***
 
//' Simulate EvAcc (2-accumulator) model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param stimulus matrix; stimulus to simulate, 1 or 2 rows X timepoints. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; boundary separation bias, left boundary = a - z, right = a + z
//' @param dc numeric; drift criterion, left = v + dc, right = v - dc
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting bias (drawn from gaussian), sz > 0, default = 0
//' @param s numeric; diffusion noise
//' @param lambda numeric; accumulator recurrent weight
//' @param mi numeric; mutual inhibition
//' @param mi2 numeric; mutual inhibition in second layer. Set mi=0 and mi2 > 0 to implement static sample model from Scott et al., 2015
//' @param sv2 numeric; sampling noise for second layer
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 1
//' @param tc numeric; time constant of collapse, default = .25
//' @param scalar_diffusion_noise bool; if true, diffusion noise scales with accumulator value
//' @param scalar_drift_noise bool; if true, drift rate variability scales with accumulator value (see Koay et al., 2020)
//' @param two_layer bool; if true, implement two layer accumulator (i.e. static sample model from Scott et al., 2015)
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv and criterion dc are multiplied by this number
//' @param dt numeric; time step of simulation, default = .001
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' @param return_accu bool; if True, return full trajectory of accumulators
//'
//' @return List containing 1) data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence and 2) matrix with full accumulator trajectories
//' 
//' @export
// [[Rcpp::export]]
List sim_evacc1(int n, arma::mat stimulus,
               double v, double a, double t0,
               double z=0, double dc=0,
               double sv=0, double st0=0, double sz=0, double s=1,
               double lambda=0, double mi=0, double mi2=0, double sv2=0,
               double aprime=0, double kappa=0, double tc=.25,
               bool scalar_diffusion_noise=false, bool scalar_drift_noise=false, 
               double v_scale=1, double dt=.001, int bounds=0,
               int n_threads=1, bool return_accu=false){
  
  int n_on_thread = n / n_threads;
  
  v *= v_scale;
  sv *= v_scale;
  dc *= v_scale;
  
  arma::vec
    rt_full(n),
    response_full(n);
  rt_full.fill(arma::datum::nan);
  response_full.fill(arma::datum::nan);
  arma::cube accumulators_full(2, stimulus.n_cols+1, n);
  if (return_accu) {
    accumulators_full.fill(arma::datum::nan);
  }
  
  // get time-varying boundary vector
  arma::vec tvec = arma::regspace(dt, dt, stimulus.n_cols*dt+dt);
  arma::vec bound;
  if(bounds == 2) {
    bound = weibull_bound(tvec, 2*a, aprime, kappa, tc);
  } else if(bounds == 1) {
    bound = hyperbolic_ratio_bound(tvec, 2*a, kappa, tc);
  } else {
    bound = rep(a, stimulus.n_cols); 
  }
  
  if (stimulus.n_rows == 1) {
    stimulus = arma::join_vert(stimulus, arma::zeros(stimulus.n_cols));
  }
  
  double dt_sr = sqrt(dt);
  
  // #pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<n_threads; i++) {
    
    arma::mat
    x(n_on_thread, 2, arma::fill::zeros),
    y(n_on_thread, 2, arma::fill::zeros),
    sv_mat(n_on_thread, 2),
    dw_mat(n_on_thread, 2),
    sv2_mat(n_on_thread, 2);
    sv_mat.fill(sv);
    dw_mat.fill(s * dt_sr);
    sv2_mat.fill(sv2);
    
    arma::vec
      z_var = z + sqrt(sz) * RcppZiggurat::zrnorm(n_on_thread),
      response,
      rt;
    
    arma::cube accumulators(2, stimulus.n_cols+1, n_on_thread),
    accumulators_done;
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0).fill(0);
    }
    
    double step = 0;
    
    while ((x.n_rows > 0) & (step < stimulus.n_cols)) {
      
      if (scalar_drift_noise) {
        sv_mat = sv * x;
        sv2_mat = sv2 * y;
      }
      
      if (scalar_diffusion_noise) {
        dw_mat = s * x * dt_sr;
      }
      
      arma::mat x_temp = x;
      
      x.col(0) += ((v + dc + sv_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))) * stimulus(0, step) * dt) -
        (lambda * x_temp.col(0) * dt) -
        (mi * x_temp.col(1) * dt) +
        (dw_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
      x.col(1) += ((v - dc + sv_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))) * stimulus(1, step) * dt) -
        (lambda * x_temp.col(1) * dt) -
        (mi * x_temp.col(0) * dt) +
        (dw_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
      x(find(x < 0)).fill(0);
      
      if (sv2 == 0) {
        y = x_temp;
      } else {
        y.col(0) = (x_temp.col(0) + sv2_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
        y.col(1) = (x_temp.col(1) + sv2_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
      }
      
      if(mi2 > 0) {
        arma::mat y_temp = y;
        y.col(0) -= mi2 * y_temp.col(1) * dt;
        y.col(1) -= mi2 * y_temp.col(0) * dt;
      }
      
      if (return_accu) {
        accumulators.col(step+1) = y.t();
      }
      
      arma::uvec lefts = arma::find(y.col(0) >= bound(step) - z_var);
      response = arma::join_cols(response, arma::zeros(lefts.n_elem));
      arma::vec temp_rt_left = rep((step + 1) * dt, lefts.n_elem);
      rt = arma::join_cols(rt, temp_rt_left);
      x.shed_rows(lefts);
      y.shed_rows(lefts);
      sv_mat.shed_rows(lefts);
      dw_mat.shed_rows(lefts);
      sv2_mat.shed_rows(lefts);
      z_var.shed_rows(lefts);
      if (return_accu) {
        accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(lefts));
        accumulators.shed_slices(lefts);
      }
      
      arma::uvec rights = arma::find(y.col(1) >= bound(step) + z_var);
      response = arma::join_cols(response, arma::ones(rights.n_elem));
      arma::vec temp_rt_right = rep((step + 1) * dt, rights.n_elem);
      rt = arma::join_cols(rt, temp_rt_right);
      x.shed_rows(rights);
      y.shed_rows(rights);
      sv_mat.shed_rows(rights);
      dw_mat.shed_rows(rights);
      sv2_mat.shed_rows(rights);
      z_var.shed_rows(rights);
      if (return_accu) {
        accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(rights));
        accumulators.shed_slices(rights);
      }
      
      step++;
      
    }
    
    if (rt.n_elem > 0) {
      rt_full(arma::span(i*n_on_thread, i*n_on_thread+rt.n_elem-1)) = rt + t0 + st0 * (arma::randu(rt.n_elem) - 0.5);
    }
    if (response.n_elem > 0) {
      response_full(arma::span(i*n_on_thread, i*n_on_thread+response.n_elem-1)) = response;
    }
    
    if (return_accu) {
      accumulators_done = arma::join_slices(accumulators_done, accumulators);
      accumulators_full(arma::span(), arma::span(), arma::span(i*n_on_thread, i*n_on_thread+accumulators_done.n_slices-1)) = accumulators_done;
    }
    
  }
  
  DataFrame sim = DataFrame::create(
    Named("response") = response_full,
    Named("rt") = rt_full);
  
  if (return_accu) {
    return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  } else {
    return List::create(Named("behavior") = sim);
  }
  
}

//' Simulate EvAcc (2-accumulator) model with fixed or collapsing boundary
//'
//' @param n integer; number of decisions to simulate
//' @param stimuli list of matrices; each element is a stimulus to simulate, 1 or 2 rows X timepoints. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
//' @param v numeric vector; drift rate
//' @param a numeric vector; initial boundary
//' @param t0 numeric vector; non-decision time
//' @param z numeric vector; boundary separation bias, left boundary = a - z, right = a + z
//' @param dc numeric vector; drift criterion, left = v + dc, right = v - dc
//' @param sv numeric vector; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric vector; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric vector; variability in starting bias (drawn from gaussian), sz > 0, default = 0
//' @param s numeric vector; diffusion noise
//' @param lambda numeric vector; accumulator recurrent weight
//' @param mi numeric vector; mutual inhibition
//' @param mi2 numeric vector; mutual inhibition in second layer. Set mi=0 and mi2 > 0 to implement static sample model from Scott et al., 2015
//' @param sv2 numeric vector; sampling noise for second layer
//' @param aprime numeric vector; degree of collapse, default = 0
//' @param kappa numeric vector; slope of collapse, default = 1
//' @param tc numeric vector; time constant of collapse, default = .25
//' @param scalar_diffusion_noise bool; if true, diffusion noise scales with accumulator value
//' @param scalar_drift_noise bool; if true, drift rate variability scales with accumulator value (see Koay et al., 2020)
//' @param two_layer bool; if true, implement two layer accumulator (i.e. static sample model from Scott et al., 2015)
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv and criterion dc are multiplied by this number
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
List sim_evacc_vec(int n, arma::cube stimuli,
                   arma::vec v, arma::vec a, arma::vec t0,
                   arma::vec z=0, arma::vec dc=0,
                   arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0,
                   arma::vec lambda=0, arma::vec mi=0, arma::vec mi2=0, arma::vec sv2=0,
                   arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0,
                   bool scalar_diffusion_noise=false, bool scalar_drift_noise=false, 
                   double v_scale=1, double dt=.001, int bounds=0,
                   bool check_pars=true, int n_threads=1, bool return_accu=false) {
  
  int n_total = n * stimuli.n_slices,
    n_on_thread = (n * stimuli.n_slices) / n_threads;

  arma::uvec trial_vector = arma::regspace<arma::uvec>(0, stimuli.n_slices-1);
  trial_vector = arma::repmat(trial_vector, n, 1);
  
  v *= v_scale;
  sv *= v_scale;
  dc *= v_scale;
  
  if (check_pars) {
    arma::uvec lens = {stimuli.n_slices, v.n_elem, a.n_elem, t0.n_elem, z.n_elem,
                       dc.n_elem, sv.n_elem, st0.n_elem, sz.n_elem, s.n_elem,
                       lambda.n_elem, mi.n_elem, mi2.n_elem, sv2.n_elem,
                       aprime.n_elem, kappa.n_elem, tc.n_elem};
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
      s = arma::zeros(lens.max()) + s(0);
    }
    if (lambda.n_elem < lens.max()) {
      lambda = arma::zeros(lens.max()) + lambda(0);
    }
    if (mi.n_elem < lens.max()) {
      mi = arma::zeros(lens.max()) + mi(0);
    }
    if (mi2.n_elem < lens.max()) {
      mi2 = arma::zeros(lens.max()) + mi2(0);  
    }
    if (sv2.n_elem < lens.max()) {
      sv2 = arma::zeros(lens.max()) + sv2(0);
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
  
  arma::vec
    rt_full = arma::zeros(n_total),
    response_full = arma::zeros(n_total);
  
  arma::cube accumulators_full(2, stimuli.n_slices, n_total);
  if (return_accu) {
    accumulators_full.fill(arma::datum::nan);
  }
  
  double dt_sr = sqrt(dt);
  
#pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<n_threads; i++) {
    
    arma::uvec thread_trials = trial_vector.rows(i*n_on_thread, i*n_on_thread+n_on_thread-1);

    // get boundary vector
    arma::vec tvec = arma::regspace(dt, dt, stimuli.n_cols*dt+dt);
    arma::mat bound;
    if(bounds == 2) {
      bound = weibull_bound_vec(tvec, 2*a(thread_trials), aprime(thread_trials), kappa(thread_trials), tc(thread_trials));
    } else if(bounds == 1) {
      bound = hyperbolic_ratio_bound_vec(tvec, 2*a(thread_trials), kappa(thread_trials), tc(thread_trials));
    } else {
      bound = repmat(arma::trans(a(thread_trials)), stimuli.n_cols, 1);
    }
    bound = bound.t();
    
    arma::mat
      x(n_on_thread, 2, arma::fill::zeros),
      y(n_on_thread, 2, arma::fill::zeros),
      sv_mat = repmat(sv(thread_trials), 1, 2),
      dw_mat = repmat(s(thread_trials) * dt_sr, 1, 2),
      sv2_mat = repmat(sv2(thread_trials), 1, 2);
    
    arma::vec z_var = z(thread_trials) + sqrt(sz(thread_trials)) % as<arma::vec>(RcppZiggurat::zrnorm(n_on_thread)),
      response,
      rt;
    
    arma::cube accumulators(2, stimuli.n_cols+1, n_on_thread),
    accumulators_done;
    if (return_accu) {
      accumulators.fill(arma::datum::nan);
      accumulators.col(0).fill(0);
    }
    
    double step = 0;
    
    while ((x.n_rows > 0) & (step < stimuli.n_cols)) {
      
      if (scalar_drift_noise) {
        sv_mat = x.each_col() % sv(thread_trials);
        sv2_mat = y.each_col() % sv2(thread_trials);
      }
      
      if (scalar_diffusion_noise) {
        dw_mat = x.each_col() % s(thread_trials) * dt_sr;
      }
      
      arma::mat x_temp = x;
      arma::vec stim0 = stimuli.tube(0, step),
        stim1 = stimuli.tube(1, step);
      stim0 = stim0(thread_trials);
      stim1 = stim1(thread_trials);
      
      x.col(0) += ((v(thread_trials) + dc(thread_trials) + sv_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))) % stim0 * dt) -
        (lambda(thread_trials) % x_temp.col(0) * dt) -
        (mi(thread_trials) % x_temp.col(1) * dt) +
        (dw_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
      x.col(1) += ((v(thread_trials) - dc(thread_trials) + sv_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows))) % stim1 * dt) -
        (lambda(thread_trials) % x_temp.col(1) * dt) -
        (mi(thread_trials) % x_temp.col(0) * dt) +
        (dw_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
      x(find(x < 0)).fill(0);
      
      if (sv2(i) == 0) {
        y = x_temp;
      } else {
        y.col(0) = (x_temp.col(0) + sv2_mat.col(0) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
        y.col(1) = (x_temp.col(1) + sv2_mat.col(1) % as<arma::vec>(RcppZiggurat::zrnorm(x_temp.n_rows)));
      }
      
      if(mi2(i) > 0) {
        arma::mat y_temp = y;
        y.col(0) -= mi2(thread_trials) % y_temp.col(1) * dt;
        y.col(1) -= mi2(thread_trials) % y_temp.col(0) * dt;
      }
      
      if (return_accu) {
        accumulators.col(step+1) = y.t();
      }
      
      arma::uvec lefts = arma::find(y.col(0) >= bound.col(step) - z_var);
      response = arma::join_cols(response, arma::ones(lefts.n_elem));
      arma::vec temp_rt_left(lefts.n_elem);
      temp_rt_left.fill((step + 1) * dt);
      rt = arma::join_cols(rt, temp_rt_left);
      x.shed_rows(lefts);
      y.shed_rows(lefts);
      sv_mat.shed_rows(lefts);
      dw_mat.shed_rows(lefts);
      sv2_mat.shed_rows(lefts);
      z_var.shed_rows(lefts);
      thread_trials.shed_rows(lefts);
      bound.shed_rows(lefts);
      if (return_accu) {
        accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(lefts));
        accumulators.shed_slices(lefts);
      }
      
      arma::uvec rights = arma::find(y.col(1) >= bound.col(step) + z_var);
      response = arma::join_cols(response, arma::zeros(rights.n_elem));
      arma::vec temp_rt_right(rights.n_elem);
      temp_rt_right.fill((step + 1) * dt);
      rt = arma::join_cols(rt, temp_rt_right);
      x.shed_rows(rights);
      y.shed_rows(rights);
      sv_mat.shed_rows(rights);
      dw_mat.shed_rows(rights);
      sv2_mat.shed_rows(rights);
      z_var.shed_rows(rights);
      thread_trials.shed_rows(rights);
      bound.shed_rows(rights);
      if (return_accu) {
        accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(rights));
        accumulators.shed_slices(rights);
      }
      
      // arma::uvec stim_complete = arma::find_nonfinite(stim0);
      // x.shed_rows(stim_complete);
      // y.shed_rows(stim_complete);
      // sv_mat.shed_rows(stim_complete);
      // dw_mat.shed_rows(stim_complete);
      // sv2_mat.shed_rows(stim_complete);
      // z_var.shed_rows(stim_complete);
      // thread_trials.shed_rows(stim_complete);
      // bound.shed_cols(stim_complete);
      // if (return_accu) {
      //   accumulators_done = arma::join_slices(accumulators_done, accumulators.slices(stim_complete));
      //   accumulators.shed_slices(stim_complete);
      // }
      
      step++;
      
    }
    
    if (rt.n_elem > 0) {
      rt_full(arma::span(i*n_on_thread, i*n_on_thread+rt.n_elem-1)) = rt + t0(i) + st0(i) * (arma::randu(rt.n_elem) - 0.5);
    }
    if (response.n_elem > 0) {
      response_full(arma::span(i*n, i*n+response.n_elem-1)) = response;
    }
    
    if (return_accu) {
      accumulators_done = arma::join_slices(accumulators_done, accumulators);
      accumulators_full.slices(i*n, i*n+accumulators_done.n_slices-1) = accumulators_done;
    }
    
  }
  
  DataFrame sim = DataFrame::create(Named("trial") = trial_vector,
                                    Named("response") = response_full,
                                    Named("rt") = rt_full);
  
  if (return_accu) {
    return List::create(Named("behavior") = sim, Named("accumulators") = accumulators_full);
  } else {
    return List::create(Named("behavior") = sim);
  }
  
}

*** END COMMENT ***/
