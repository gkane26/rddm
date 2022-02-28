#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//' Get pulse stimulus
//'
//' @description get full pulse stimulus
//'
//' @param stim_seq numeric vector; vecotr of timepoint by timepoint evidence
//' @param dur numeric; duration of pulse, default = .01
//' @param isi numeric; inter-pulse interval, default = .1
//' @param pre_stim numeric; time before pulse within stimulus bin
//' @param dt numeric; timestep for simulation, default = .001
//'
//' @return stimulus train vector
//'
//' @export
// [[Rcpp::export]]
arma::rowvec pulse_trial_stimulus(arma::vec stim_seq, double dur=0.01, double isi=0.1, double pre_stim=0, double dt=0.001){
  
  int pre_bins = round(pre_stim / dt),
    pulse_bins = round(dur / dt),
    stim_bins = round(isi / dt);
  
  arma::rowvec stimulus(stim_bins * stim_seq.n_elem + pre_bins, arma::fill::zeros);
  
  for(unsigned int i=0; i<stim_seq.n_elem; i++){
    
    arma::span this_pulse = arma::span(i * stim_bins + pre_bins,
                                       i * stim_bins + pre_bins + pulse_bins - 1);
    stimulus(this_pulse).fill(stim_seq(i));
    
  }
  
  return stimulus;
  
}
