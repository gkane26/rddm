#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <unistd.h>
#include <omp.h>
#include "bounds.h"
#include "urgency.h"
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp)]]


//' Get first passage time distribution of pulse diffusion model by simulating probability mass
//'
//' @param stimulus matrix; stimulus to simulate (row 1 is evidence to upper boundary, row 2 to lower boundary)
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param s numeric; standard deviation in wiener diffusion noise
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param lambda numeric; O-U process slope
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 0
//' @param tc numeric; time constant of collapse, default = .25
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
//' @param uslope numeric; urgency scaling factor, default = 0;
//' @param umag numeric; urgency magnitude, default = 0;
//' @param udelay numeric; urgency delay, default = 0;
//' @param dt numeric; time step of simulation, default = .001
//' @param xbins numeric; number of evidence bins, default = 100
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
//' @param urgency int; 0 for none, 1 for linear, 2 for logistic
//'
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//'
//' @export
// [[Rcpp::export]]
arma::mat pulse_fp_fpt(arma::mat stimulus, double v, double a, double t0,
                       double s=1, double z=0.5, double dc=0,
                       double sv=0, double st0=0, double sz=0, double lambda=0,
                       double aprime=0, double kappa=0, double tc=.25,
                       double uslope=0, double umag=0, double udelay=0,
                       double v_scale=1, double dt=.001, int xbins=200, int bounds=0, int urgency=0){
  
  // scale drift rate
  v *= v_scale;
  sv *= v_scale;
  dc *= v_scale;
  
  // make bins
  int n_x_breaks = xbins;
  if(n_x_breaks % 2 == 1)
    n_x_breaks++;
  
  arma::vec bin_breaks = arma::linspace(-a/2, a/2, n_x_breaks);
  arma::vec bin_centers = bin_breaks.head(bin_breaks.n_elem-1) + diff(bin_breaks)/2;
  
  // initialize prob mass matrix and transition matrix
  arma::mat p_x = arma::zeros(bin_breaks.n_elem+1, stimulus.n_cols+1),
    p_rt = arma::zeros(stimulus.n_cols, 2),
    t_mat = arma::zeros(bin_breaks.n_elem+1, bin_breaks.n_elem+1);
  t_mat(0,0)=1;
  t_mat(t_mat.n_rows-1,t_mat.n_cols-1)=1;
  
  // populate first bin of prob mass matrix: gaussian with mean=z, sd=sqrt(sz)
  if (sz < 1e-10)
    sz = 1e-10;
  arma::vec p_breaks = arma::normcdf(bin_breaks, a/2*(z-0.5), sqrt(sz));
  p_x(0, 0) = p_breaks(0);
  p_x(arma::span(1,p_x.n_rows-2), 0) = arma::diff(p_breaks);
  p_x(p_x.n_rows-1, 0) = 1-p_breaks(p_breaks.n_elem-1);
  p_rt(0, 0) = p_x(p_x.n_rows-1, 1);
  p_rt(0, 1) = p_x(0, 1);
  
  if (stimulus.n_rows == 1) {
    stimulus = arma::join_vert(stimulus, arma::zeros(stimulus.n_cols));
  }
  
  // get boundary vector
  arma::vec tvec = arma::regspace(dt, dt, stimulus.n_cols*dt+dt),
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
  
  unsigned int low_bound_index = 0,
    up_bound_index = t_mat.n_rows-1;
  
  //iterate through time and over bins
  for(unsigned int t=0; t<stimulus.n_cols; t++){
    
    // adjust transition matrix for dynamic bounds
    low_bound_index = arma::max(arma::find(bin_breaks <= -bound(t)));
    up_bound_index = arma::min(arma::find(bin_breaks >= bound(t)))+1;
    
    // get transition matrix
    double sigma2 = gamma(t)*gamma(t)*s*s*dt + gamma(t)*sv*abs(stimulus(0, t)*(v+dc))*dt*dt + gamma(t)*sv*abs(stimulus(1, t)*(-v+dc))*dt*dt;
    for(unsigned int j=1; j<t_mat.n_rows-1; j++){
      double mu = exp(-lambda*dt)*bin_centers(j-1) + gamma(t)*(((v+dc)*stimulus(0, t) + (-v+dc)*stimulus(1, t)) * dt);
      p_breaks = arma::normcdf(bin_breaks, mu, sqrt(sigma2));
      t_mat(0, j) = p_breaks(0);
      t_mat(arma::span(1, t_mat.n_rows-2), j) = arma::diff(p_breaks);
      t_mat(t_mat.n_rows-1, j) = 1-p_breaks(p_breaks.n_elem-1);
    }
    
    t_mat.cols(arma::span(0, low_bound_index)).fill(0);
    t_mat.cols(arma::span(up_bound_index, t_mat.n_cols-1)).fill(0);
    t_mat(0, arma::span(0, low_bound_index)).fill(1);
    t_mat(t_mat.n_rows-1, arma::span(up_bound_index, t_mat.n_cols-1)).fill(1);
    
    p_x(arma::span::all, t+1) = t_mat * p_x(arma::span::all, t);
    p_rt(t, 0) += p_x(p_x.n_rows-1, t+1) - p_x(p_x.n_rows-1, t);
    p_rt(t, 1) += p_x(0, t+1) - p_x(0, t);
    
  }
  
  // add t0 to response time distribution
  
  int t0_bin_min = round((t0-st0/2)/dt);
  if (t0_bin_min < 0)
    t0_bin_min = 0;
  int t0_bin_max = round((t0+st0/2)/dt);
  if (t0_bin_max < 0)
    t0_bin_max = 0;
  arma::mat t0_rt = arma::zeros(t0_bin_max+p_rt.n_rows,2);
  for(int i=t0_bin_min; i<=t0_bin_max; i++){
    t0_rt(arma::span(i,i+p_rt.n_rows-1), 0) += p_rt.col(0);
    t0_rt(arma::span(i,i+p_rt.n_rows-1), 1) += p_rt.col(1);
  }
  
  return t0_rt / (t0_bin_max-t0_bin_min+1);
  
}


////////////////////////////////////////////////
/* *** likelihood function via simulating probability mass *** */
////////////////////////////////////////////////


//' Get pulse model likelihood for a given trial
//'
//' @param choice int; decision on trial, 0 for lower boundary, 1 for upper
//' @param rt numeric; response time on trial
//' @param stimulus matrix; stimulus to simulate (row 1 is evidence to upper boundary, row 2 to lower boundary)
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param s numeric; standard deviation in wiener diffusion noise
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param lambda numeric; O-U process slope
//' @param aprime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 0
//' @param tc numeric; time constant of collapse, default = .25
//' @param uslope numeric; urgency scaling factor, default = 0;
//' @param umag numeric; urgency magnitude, default = 0;
//' @param udelay numeric; urgency delay, default = 0;
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
//' @param dt numeric; time step of simulation, default = .001
//' @param xbins numeric; number of evidence bins, default = 100
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
//' @param urgency int; 0 for none, 1 for linear, 2 for logistic
//'
//' @return probability of choice and rt for trial given pulse model parameters
//'
//' @export
// [[Rcpp::export]]
double pulse_trial_lik(int choice, double rt, arma::mat stimulus,
                       double v, double a, double t0,
                       double s=1, double z=0.5, double dc=0,
                       double sv=0, double st0=0, double sz=0, double lambda=0,
                       double aprime=0, double kappa=0, double tc=.25,
                       double uslope=0, double umag=0, double udelay=0,
                       double v_scale=1, double dt=.001, int xbins=200, int bounds=0, int urgency=0){
  
  arma::mat fpt_density = pulse_fp_fpt(stimulus,
                                       v, a, t0,
                                       s, z, dc,
                                       sv, st0, sz,
                                       lambda, aprime, kappa, tc,
                                       uslope, umag, udelay,
                                       v_scale, dt, xbins, bounds, urgency);
  
  return fpt_density(rt/dt-1, abs(1-choice));
  
}


//' Get pulse model negative log likelihood for a set of trials
//'
//' @param choice integer; vector of decisions, 0 for lower boundary, 1 for upper
//' @param rt numeric; vector of response times
//' @param stimuli array; 2 x timepoints x trials array of pulse stimuli
//' @param up_sequence vector of strings; string of 0s, and 1s of stimulus values (0 no evidence, 1 to upper). If down_sequence not specified, (0 to lower, 1 to upper).
//' @param down_sequence vector of strings; string of 0s, and 1s of stimulus values (0 is no evidence, 1 to lower). If not specified, up_sequence is (0 to lower, 1 to upper)
//' @param v numeric; drift rate, either single value or vector for each trial
//' @param a numeric; initial boundary, either single value or vector for each trial
//' @param s numeric; standard deviation in wiener diffusion noise, either single value or vector for each trial
//' @param t0 numeric; non-decision time, either single value or vector for each trial
//' @param z numeric; starting point, , either single value or vector for each trial, 0 < z < 1, default = .5
//' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
//' @param sv numeric; standard deviation of variability in drift rate, either single value or vector for each trial, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time, either single value or vector for each trial. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point, either single value or vector for each trial. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param lambda numeric; O-U process slope, either single value or vector for each trial
//' @param aprime numeric; degree of collapse, either single value or vector for each trial, default = 0
//' @param kappa numeric; slope of collapse, either single value or vector for each trial, default = 0
//' @param tc numeric; time constant of collapse, either single value or vector for each trial, default = .25
//' @param uslope numeric; urgency scaling factor, default = 0;
//' @param umag numeric; urgency magnitude, default = 0;
//' @param udelay numeric; urgency delay, default = 0;
//' @param check_pars logical; if True, check that parameters are vectors of the same length as choices and rts. Must be true if providing scalar parameters. default = true
//' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
//' @param dt numeric; time step of simulation, default = .002
//' @param xbins integer; number of evidence bins, default = 100
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
//' @param urgency int: 0 for none, 1 for linear, 2 for logistic
//' @param n_threads int; number of threads (trials) to run in parallel
//'
//' @return negative log likelihood of all choices and rts given pulse model parameters
//'
//' @export
// [[Rcpp::export]]
double pulse_nll(arma::vec choices, arma::vec rt, arma::cube stimuli,
                 arma::vec v, arma::vec a, arma::vec t0,
                 arma::vec s=0, arma::vec z=0, arma::vec dc=0,
                 arma::vec sv=0, arma::vec st0=0, arma::vec sz=0,
                 arma::vec lambda=0, arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0,
                 arma::vec uslope=0, arma::vec umag=0, arma::vec udelay=0, bool check_pars=true,
                 double v_scale=1, double dt=.001, int xbins=200, int bounds=0, int urgency=0, int n_threads=1){
  
  omp_set_num_threads(n_threads);
  
  // check parameter vectors
  if(check_pars){
    s(arma::find(s==0)).fill(1);
    z(arma::find(z==0)).fill(0.5);
    tc(arma::find(tc==0)).fill(0.25);
    if(v.n_elem != choices.n_elem)
      v = arma::zeros(choices.n_elem)+v(0);
    if(a.n_elem != choices.n_elem)
      a = arma::zeros(choices.n_elem)+a(0);
    if(t0.n_elem != choices.n_elem)
      t0 = arma::zeros(choices.n_elem)+t0(0);
    if(s.n_elem != choices.n_elem)
      s = arma::zeros(choices.n_elem)+s(0);
    if(z.n_elem != choices.n_elem)
      z = arma::zeros(choices.n_elem)+z(0);
    if(dc.n_elem != choices.n_elem)
      dc = arma::zeros(choices.n_elem)+ dc(0);
    if(sv.n_elem != choices.n_elem)
      sv = arma::zeros(choices.n_elem)+sv(0);
    if(sz.n_elem != choices.n_elem)
      sz = arma::zeros(choices.n_elem)+sz(0);
    if(st0.n_elem != choices.n_elem)
      st0 = arma::zeros(choices.n_elem)+st0(0);
    if(lambda.n_elem != choices.n_elem)
      lambda = arma::zeros(choices.n_elem)+lambda(0);
    if(aprime.n_elem != choices.n_elem)
      aprime = arma::zeros(choices.n_elem)+aprime(0);
    if(kappa.n_elem != choices.n_elem)
      kappa = arma::zeros(choices.n_elem)+kappa(0);
    if(tc.n_elem != choices.n_elem)
      tc = arma::zeros(choices.n_elem)+tc(0);
    if(uslope.n_elem != choices.n_elem)
      uslope = arma::zeros(choices.n_elem)+uslope(0);
    if(umag.n_elem != choices.n_elem)
      umag = arma::zeros(choices.n_elem)+umag(0);
    if(udelay.n_elem != choices.n_elem)
      udelay = arma::zeros(choices.n_elem)+udelay(0);
  }
  
  double p_rt=0;
  
  // omp parallel over decisions
#pragma omp parallel for num_threads(n_threads)
  for(unsigned int i=0; i<choices.n_elem; i++){
    
    double p_local = pulse_trial_lik(choices(i), rt(i), stimuli.slice(i),
                                     v(i), a(i), t0(i),
                                     s(i), z(i), dc(i),
                                     sv(i), st0(i), sz(i),
                                     lambda(i), aprime(i), kappa(i), tc(i),
                                     uslope(i), umag(i), udelay(i),
                                     v_scale, dt, xbins, bounds, urgency);
    
    if(p_local<1e-10)
      p_local=1e-10;
    
#pragma omp critical
{
  p_rt += log(p_local);
}   
  }
  
  return -p_rt;
  
}

