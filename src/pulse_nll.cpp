#include "RcppArmadillo.h"
#include "bounds.h"
#include "omp.h"
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp)]]

////////////////////////////////////////////////
/* *** helper functions *** */
////////////////////////////////////////////////

//' Get pulse stimulus
//'
//' @description get full pulse stimulus
//'
//' @param sequence string; string of 0s, and 1s of stimulus values (0 is evidence to lower, 1 is evidence to upper)
//' @param intensity numeric; intensity of each pulse, default = 1
//' @param dur numeric; duration of pulse, default = .01
//' @param isi numeric; inter-pulse interval, default = .1
//' @param dt numeric; timestep for simulation, default = .001
//'
//' @return stimulus train
//'
//' @export
// [[Rcpp::export]]
arma::ivec get_stimulus(std::string sequence, double intensity=1, double dur=0.01, double isi=0.1, double dt=0.001){

  int stim_bins = round(dur/dt),
    blink_bins = round(isi/dt);

  arma::ivec stimulus = arma::zeros<arma::ivec>(blink_bins*(sequence.size()+1));
  for(unsigned int i=0; i<sequence.size(); i++){
    int this_stim = std::stoi(sequence.substr(i,1));
    stimulus(arma::span((i+1)*blink_bins, (i+1)*blink_bins+stim_bins-1)).fill(2*this_stim-1);
  }

  return stimulus;
}

//' Get first passage time distribution of pulse diffusion model by simulating probability mass
//'
//' @param stimulus integer vector; stimulus to simulate (1 for evidence to upper boundary, -1 to lower boundary, 0 for no evidence)
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param lambda numeric; O-U process slope
//' @param a_prime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 0
//' @param tc numeric; time constant of collapse, default = .25
//' @param dt numeric; time step of simulation, default = .002
//' @param dx numeric; size of evidence bins, default = .05
//' @param v_scale numeric; scale drift rate to be similar to boundary separation a, default = 10
//' @param use_weibull_bound logical; if True, use weibull function for collapsing bounds, if False, use hyperbolic ratio function
//'
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//'
//' @export
// [[Rcpp::export]]
arma::mat pulse_pmass_fpt(arma::ivec stimulus, double v, double a, double t0, double z=0,
                      double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                      double a_prime=0, double kappa=0, double tc=.25,
                      double dt=.002, double dx=.05, double v_scale=10, bool use_weibull_bound=false){

  // make bins
  int n_x_breaks = round(a / dx);
  if(n_x_breaks % 2 == 1)
    n_x_breaks++;
  arma::vec bin_breaks = arma::linspace(-a/2, a/2, n_x_breaks);
  arma::vec bin_centers = bin_breaks.head(bin_breaks.n_elem-1) + diff(bin_breaks)/2;

  // initialize prob mass matrix and transition matrix
  arma::mat p_x = arma::zeros(bin_breaks.n_elem+1, stimulus.n_elem+1),
    p_rt = arma::zeros(stimulus.n_elem, 2),
    t_mat = arma::zeros(bin_breaks.n_elem+1, bin_breaks.n_elem+1);
  t_mat(0,0)=1;
  t_mat(t_mat.n_rows-1,t_mat.n_cols-1)=1;

  // populate first bin of prob mass matrix: gaussian with mean=z, sd=sqrt(sz)
  if (sz < 1e-10)
    sz = 1e-10;
  arma::vec p_breaks = arma::normcdf(bin_breaks, (2*z-1)*a/2, sqrt(sz));
  p_x(0, 0) = p_breaks(0);
  p_x(arma::span(1,p_x.n_rows-2), 0) = arma::diff(p_breaks);
  p_x(p_x.n_rows-1, 0) = 1-p_breaks(p_breaks.n_elem-1);

  // get boundary vector
  arma::vec bound;
  if(use_weibull_bound)
    bound = weibull_bound(arma::regspace(dt, dt, stimulus.n_elem*dt+dt), a, a_prime, kappa, tc);
  else
    bound = hyperbolic_ratio_bound(arma::regspace(dt, dt, stimulus.n_elem*dt+dt), a, kappa, tc);

  unsigned int low_bound_index = 0,
    up_bound_index = t_mat.n_rows-1;
  
  //iterate through time and over bins
  for(unsigned int t=1; t<stimulus.n_elem; t++){
    
    // adjust transition matrix for dynamic bounds
    low_bound_index = arma::max(arma::find(bin_breaks <= -bound(t)));
    up_bound_index = arma::min(arma::find(bin_breaks >= bound(t)))+1;
    
    // get transition matrix
    double sigma2 = s*dt + sv*abs(stimulus(t))*v*dt*dt;
    for(unsigned int j=1; j<t_mat.n_rows-1; j++){
      double mu = exp(lambda*dt)*bin_centers(j-1) + v*stimulus(t)*dt;
      p_breaks = arma::normcdf(bin_breaks, mu, sqrt(sigma2));
      t_mat(0, j) = p_breaks(0);
      t_mat(arma::span(1, t_mat.n_rows-2), j) = arma::diff(p_breaks);
      t_mat(t_mat.n_rows-1, j) = 1-p_breaks(p_breaks.n_elem-1);
    }

    t_mat.cols(arma::span(0, low_bound_index)).fill(0);
    t_mat.cols(arma::span(up_bound_index, t_mat.n_cols-1)).fill(0);
    t_mat(0, arma::span(0, low_bound_index)).fill(1);
    t_mat(t_mat.n_rows-1, arma::span(up_bound_index, t_mat.n_cols-1)).fill(1);

    p_x(arma::span::all, t) = t_mat * p_x(arma::span::all, t-1);
    p_rt(t-1, 0) = p_x(p_x.n_rows-1, t) - p_x(p_x.n_rows-1, t-1);
    p_rt(t-1, 1) = p_x(0, t) - p_x(0, t-1);

  }

  p_rt(0, 0) = p_x(p_x.n_rows-1, 1);
  p_rt(0, 1) = p_x(0, 1);

  // add t0 to response time distribution

  int t0_bin_min = round((t0-st0/2)/dt), t0_bin_max = round((t0+st0/2)/dt);
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
//' @param stimulus string; stimulus sequence for trial, a string of 0s and 1s (0 for evidence to lower, 1 for evidence to upper)
//' @param v numeric; drift rate
//' @param a numeric; initial boundary
//' @param t0 numeric; non-decision time
//' @param z numeric; starting point, 0 < z < 1, default = .5
//' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
//' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
//' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param lambda numeric; O-U process slope
//' @param a_prime numeric; degree of collapse, default = 0
//' @param kappa numeric; slope of collapse, default = 0
//' @param tc numeric; time constant of collapse, default = .25
//' @param dt numeric; time step of simulation, default = .002
//' @param dx numeric; size of evidence bins, default = .05
//' @param v_scale numeric; scale drift rate to be similar to boundary separation a, default = 10
//' @param use_weibull_bound logical; if True, use weibull function for collapsing bounds, if False, use hyperbolic ratio function
//' @param dur numeric; duration of stimulus
//' @param isi numeric; interstimulus interval
//'
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//'
//' @export
// [[Rcpp::export]]
double pulse_trial_lik(int choice, double rt, std::string blink_seq,
                       double v, double a, double t0, double z=0,
                       double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                       double a_prime=0, double kappa=0, double tc=.25,
                       double dt=.002, double dx=.05, double v_scale=10, bool use_weibull_bound=false, double dur=.01, double isi=.1){
  
  arma::ivec stimulus = get_stimulus(blink_seq, 1, dur, isi);
  arma::mat fpt_density = pulse_pmass_fpt(stimulus, v, a, t0, z, sv, st0, sz, s, lambda, a_prime, kappa, tc, dt, dx, v_scale, use_weibull_bound);
  return fpt_density(rt/dt-1, abs(1-choice));
}


//' Get pulse model negative log likelihood for a set of trials
//'
//' @param choice integer; vector of decisions, 0 for lower boundary, 1 for upper
//' @param rt numeric; vector of response times
//' @param stimulus vector of strings; vector of stimulus sequence for all trials, each element should be a string of 0s and 1s (0 for evidence to lower, 1 for evidence to upper)
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
//' @param dt numeric; time step of simulation, default = .002
//' @param dx numeric; size of evidence bins, default = .05
//' @param v_scale numeric; scale drift rate to be similar to boundary separation a, default = 10
//' @param use_weibull_bound logical; if True, use weibull function for collapsing bounds, if False, use hyperbolic ratio function
//' @param dur numeric; duration of stimulus
//' @param isi numeric; interstimulus interval
//' @param n_threads int; number of threads (trials) to run in parallel
//'
//' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
//'
//' @export
// [[Rcpp::export]]
double pulse_nll(arma::vec choices, arma::vec rt, std::vector<std::string> blink_seq,
                 arma::vec v, arma::vec a, arma::vec t0,
                 arma::vec z=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0,
                 arma::vec lambda=0, arma::vec a_prime=0, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true,
                 double dt=.002, double dx=.05, double v_scale=10, bool use_weibull_bound=false, double dur=.01, double isi=.1, int n_threads=1){

  omp_set_num_threads(n_threads);
  
  // check parameter vectors
  if(check_pars){
    z(arma::find(z==0)).fill(0.5);
    s(arma::find(s==0)).fill(1);
    tc(arma::find(tc==0)).fill(0.25);
    if(v.n_elem != choices.n_elem)
      v = arma::zeros(choices.n_elem)+v(0);
    if(a.n_elem != choices.n_elem)
      a = arma::zeros(choices.n_elem)+a(0);
    if(t0.n_elem != choices.n_elem)
      t0 = arma::zeros(choices.n_elem)+t0(0);
    if(z.n_elem != choices.n_elem)
      z = arma::zeros(choices.n_elem)+.5;
    if(sv.n_elem != choices.n_elem)
      sv = arma::zeros(choices.n_elem)+sv(0);
    if(sz.n_elem != choices.n_elem)
      sz = arma::zeros(choices.n_elem)+sz(0);
    if(st0.n_elem != choices.n_elem)
      st0 = arma::zeros(choices.n_elem)+st0(0);
    if(s.n_elem != choices.n_elem)
      s = arma::zeros(choices.n_elem)+s(0);
    if(lambda.n_elem != choices.n_elem)
      lambda = arma::zeros(choices.n_elem)+lambda(0);
    if(a_prime.n_elem != choices.n_elem)
      a_prime = arma::zeros(choices.n_elem)+a_prime(0);
    if(kappa.n_elem != choices.n_elem)
      kappa = arma::zeros(choices.n_elem)+kappa(0);
    if(tc.n_elem != choices.n_elem)
      tc = arma::zeros(choices.n_elem)+tc(0);
  }
  
  double p_rt=0;
  
  
  // omp parallel over decisions
#pragma omp parallel for
  for(unsigned int i=0; i<choices.n_elem; i++){
    
    double p_local = pulse_trial_lik(choices(i), rt(i), blink_seq[i],
                                     v(i), a(i), t0(i), z(i), sv(i), st0(i), sz(i), s(i),
                                     lambda(i), a_prime(i), kappa(i), tc(i),
                                     dt, dx, v_scale, use_weibull_bound, dur, isi);

    if(p_local<1e-10)
      p_local=1e-10;
    
#pragma omp critical
{
    p_rt += log(p_local);
}   
  }
  
    return -p_rt;
    
}
