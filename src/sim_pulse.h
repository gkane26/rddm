#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

List sim_pulse(int n, arma::cube stimuli, double a, double t0, double s,
               double z=.5, double dc=0,
               double sv=0, double st0=0, double sz=0,
               double lambda=0, double aprime=0, double kappa=0, double tc=.25, 
               double uslope=0, double umag=0, double udelay=0,
               double v_scale=1, double dt=.001, int bounds=0, int urgency=0,
               int n_threads=1, bool return_accu=false, int seed=-1);

// List sim_pulse_vec(int n, List stimuli, arma::vec v, arma::vec a, arma::vec t0,
//                    arma::vec z=0, arma::vec dc=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0, arma::vec lambda=0,
//                    arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0, 
//                    double v_scale=1, double dt=.001, int bounds=0, bool check_pars=true, int n_threads=1, bool return_accu=false);
