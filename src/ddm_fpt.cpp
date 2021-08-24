/*
 
 functions to get analytical first passage time distribution for
diffusion model with time varying boundary
from Voskuilen, Ratcliff & Smith (2016). J Math Psychol

*/

#include <RcppArmadillo.h>
#include <omp.h>
#include "bounds.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/* helper functions */

double free_transition_density(double x, double t, double y, double tau, double v=0, double s=1){
  // probability that a particle starting at y at time tau will be at x at time t
  //  v = drift rate
  
  return (1 / sqrt(2*PI*pow(s,2)*(t-tau))) * exp(-pow(x - y - v*(t-tau), 2) / (2 * pow(s,2) * (t-tau)));
}

arma::vec free_transition_density_vec(double x, double t, arma::vec y, double tau, double v=0, double s=1){
  // probability that a particle starting at y at time tau will be at x at time t
  //  v = drift rate
  
  return (1 / sqrt(2*PI*pow(s,2)*(t-tau))) * exp(-pow(x - y - v*(t-tau), 2) / (2 * pow(s,2) * (t-tau)));
}

double kernel_function(double x, double t, double y, double tau, double dx=0, double v=0, double s=1){
  // probability that a particle starting at 'y' at time 'tau' will be at 'x' at time 't'
  
  return (free_transition_density(x, t, y, tau, v, s) / 2) * (dx - (x-y)/(t-tau));
}

arma::vec kernel_function_vec(double x, double t, arma::vec y, double tau, double dx=0, double v=0, double s=1){
  // probability that a particle starting at 'y' at time 'tau' will be at 'x' at time 't'
  
  return (free_transition_density_vec(x, t, y, tau, v, s) / 2) % (dx - (x-y)/(t-tau));
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//' Simulate drift diffusion model with fixed or collapsing boundary
//'
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
//' @param s numeric; standard deviation in wiener diffusion noise, default = 1
//' @param sv_points integer; number of points to approximate integral over drift rate variability, default = 19
//' @param sz_points integer; number of points to approximate integral over starting point variability, default = 11
//' @param st0_points integer; number of points to approximate integral over non-decision time variability, default = 11
//' @param dt numeric; time step of simulation, default = .01
//' @param max_time numeric; max time of simulation, default = 10
//' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
//' @param n_threads integer; number of threads to run in parallel, default = 1
//' 
//' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
//' 
//' @export
// [[Rcpp::export]]
arma::mat ddm_integral_fpt(double v, double a, double t0, double z=.5, double dc=0,
                       double sv=0, double sz=0, double st0=0,
                       double aprime=0, double kappa=0, double tc=.25, double s=1,
                       double sv_points=19, double sz_points=11, double st0_points=11,
                       double dt=.01, double max_time=10, int bounds=0, int n_threads=1){
  
  t0 = round(t0/dt)*dt;
  z = (2*z-1)*a/2;
  int nBins = round(max_time / dt)+1;
  
  arma::vec v_vec, p_v, z_vec;
  v += dc;
  if(sv==0){
    v_vec = v;
    p_v = 1;
  }else{
    v_vec = arma::linspace(v-3*sv, v+3*sv, sv_points);
    p_v = normpdf(v_vec, v, sv);
    p_v = p_v / sum(p_v);
  }
  
  if(sz==0)
    z_vec = z;
  else
    z_vec = arma::linspace(z-sz/2, z+sz/2, sz_points);
  
  arma::vec tvec = arma::regspace(dt, dt, max_time+dt);
  arma::vec bound;
  if(bounds == 2){
    bound = weibull_bound(tvec, a, aprime, kappa, tc);
  }else if(bounds == 1){
    bound = hyperbolic_ratio_bound(tvec, a, kappa, tc);
  } else {
    bound = rep(a/2, tvec.n_elem); 
  }
  
  arma::mat density_shared = arma::zeros(nBins, 2);
  
  if(n_threads==0)
    n_threads = v_vec.n_elem * z_vec.n_elem;
  
  //parallelize over variable drift rates (v_vec) and variable starting points
#pragma omp parallel for num_threads(n_threads) collapse(2)
  
  for(unsigned int i=0; i<v_vec.n_elem; i++){
    
    // loop through variable starting points (z_vec)
    for(unsigned int j=0; j<z_vec.n_elem; j++){
      
      // initialize density on local thread
      arma::mat density_thread = arma::zeros(nBins, 2);
      density_thread(0,0) = -2 * kernel_function(bound[0], dt, z_vec(j), 0, 0, v_vec(i), s);
      density_thread(0,1) = 2 * kernel_function(-bound[0], dt, z_vec(j), 0, 0, v_vec(i), s);
      
      double dk=0, up_a, up_b, down_a, down_b;
      
      for(int k=1; k<nBins; k++){
        
        dk = (bound[k]-bound[k-1])/dt;
        density_thread(k,0) = -2 * kernel_function(bound[k], (k+1)*dt, z_vec(j), 0, dk, v_vec(i), s);
        density_thread(k,1) = 2 * kernel_function(-bound[k], (k+1)*dt, z_vec(j), 0, -dk, v_vec(i), s);
        
        up_a=0, up_b=0, down_a=0, down_b=0;
        for(int l=0; l<k; l++){
          up_a += density_thread(l,0) * kernel_function(bound[k], (k+1)*dt, bound[l], (l+1)*dt, dk, v_vec(i), s);
          down_a += density_thread(l,1) * kernel_function(bound[k], (k+1)*dt, -bound[l], (l+1)*dt, dk, v_vec(i), s);
          up_b += density_thread(l,0) * kernel_function(-bound[k], (k+1)*dt, bound[l], (l+1)*dt, -dk, v_vec(i), s);
          down_b += density_thread(l,1) * kernel_function(-bound[k], (k+1)*dt, -bound[l], (l+1)*dt, -dk, v_vec(i), s);
        }
        
        density_thread(k,0) += 2*dt*(up_a + down_a);
        if(density_thread(k,0) < 0)
          density_thread(k,0) = 0;
        density_thread(k,1) += -2*dt*(up_b + down_b);
        if(density_thread(k,1) < 0)
          density_thread(k,1) = 0;
      }
      
#pragma omp critical
      density_shared += density_thread * p_v(i) * (1.0/z_vec.n_elem);
    }
  }
  
  // add t0 to get rt_density
  double min_t0_bins = round((t0-st0/2)/dt),
    max_t0_bins = round((t0+st0/2)/dt),
    p_t0 = 1 / (max_t0_bins - min_t0_bins + 1);
  arma::mat rt_density = arma::zeros(density_shared.n_rows + max_t0_bins + 1, 2);
  
  for(int i=min_t0_bins; i<=max_t0_bins; i++){
    rt_density(arma::span(i, i+density_shared.n_rows-1), 0) += density_shared.col(0) * p_t0;
    rt_density(arma::span(i, i+density_shared.n_rows-1), 1) += density_shared.col(1) * p_t0;
  }
  
  return rt_density * dt;
  
}