#include "RcppArmadillo.h"
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

arma::ivec get_stimulus(std::string sequence, double intensity=1, double dur=0.01, double isi=0.1, double dt=0.001);
arma::mat pulse_pmass_fpt(arma::ivec stimulus, double v, double a, double t0, double z=0,
                          double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                          double a_prime=0, double kappa=0, double tc=.25,
                          double dt=.002, double dx=.05, double v_scale=100, bool use_weibull_bound=false);
double pulse_trial_lik(int choice, double rt, std::string stim_seq,
                       double v, double a, double t0, double z=0,
                       double sv=0, double st0=0, double sz=0, double s=1, double lambda=0,
                       double a_prime=0, double kappa=0, double tc=.25,
                       double dt=.002, double dx=.05, double v_scale=100, bool use_weibull_bound=false, double dur=.01, double isi=.1);
double pulse_nll(arma::vec choices, arma::vec rt, std::vector<std::string> stim_seq,
                 arma::vec v, arma::vec a, arma::vec t0,
                 arma::vec z=0, arma::vec sv=0, arma::vec st0=0, arma::vec sz=0, arma::vec s=0,
                 arma::vec lambda=0, arma::vec a_prime=0, arma::vec kappa=0, arma::vec tc=0, bool check_pars=true,
                 double dt=.002, double dx=.05, double v_scale=100, bool use_weibull_bound=false, double dur=.01, double isi=.1, int n_threads=1);
