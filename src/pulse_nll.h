#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

arma::mat pulse_fp_fpt(arma::mat stimulus, double v, double a, double t0,
                       double s=1, double z=0, double dc=0,
                       double sv=0, double st0=0, double sz=0, double lambda=0,
                       double aprime=0, double kappa=0, double tc=.25,
                       double uslope=0, double umag=0, double udelay=0,
                       double v_scale=1, double dt=.001, int xbins=200, int bounds=0, int urgency=0);
double pulse_trial_lik(int choice, double rt, arma::mat stimulus,
                       double v, double a, double t0,
                       double s=1, double z=0, double dc=0,
                       double sv=0, double st0=0, double sz=0, double lambda=0,
                       double aprime=0, double kappa=0, double tc=.25,
                       double uslope=0, double umag=0, double udelay=0,
                       double v_scale=1, double dt=.001, int xbins=200, int bounds=0, double urgency=0);
double pulse_nll(arma::vec choices, arma::vec rt, List stimuli,
                 arma::vec v, arma::vec a, arma::vec t0,
                 arma::vec s=0, arma::vec z=0, arma::vec dc=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0,
                 arma::vec lambda=0, arma::vec aprime=0, arma::vec kappa=0, arma::vec tc=0,
                 arma::vec uslope=0, arma::vec umag=0, arma::vec udelay=0, bool check_pars=true,
                 double v_scale=1, double dt=.001, int xbins=200, int bounds=0, int urgency=0, int n_threads=1);
