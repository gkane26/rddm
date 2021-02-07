// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// hyperbolic_ratio_bound
arma::vec hyperbolic_ratio_bound(arma::vec t, double a, double kappa, double tc);
RcppExport SEXP _rddm_hyperbolic_ratio_bound(SEXP tSEXP, SEXP aSEXP, SEXP kappaSEXP, SEXP tcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    rcpp_result_gen = Rcpp::wrap(hyperbolic_ratio_bound(t, a, kappa, tc));
    return rcpp_result_gen;
END_RCPP
}
// hyperbolic_ratio_bound_vec
arma::mat hyperbolic_ratio_bound_vec(arma::vec t, arma::vec a, arma::vec kappa, arma::vec tc, bool check_pars);
RcppExport SEXP _rddm_hyperbolic_ratio_bound_vec(SEXP tSEXP, SEXP aSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    rcpp_result_gen = Rcpp::wrap(hyperbolic_ratio_bound_vec(t, a, kappa, tc, check_pars));
    return rcpp_result_gen;
END_RCPP
}
// weibull_bound
arma::vec weibull_bound(arma::vec t, double a, double aprime, double kappa, double tc);
RcppExport SEXP _rddm_weibull_bound(SEXP tSEXP, SEXP aSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    rcpp_result_gen = Rcpp::wrap(weibull_bound(t, a, aprime, kappa, tc));
    return rcpp_result_gen;
END_RCPP
}
// weibull_bound_vec
arma::mat weibull_bound_vec(arma::vec t, arma::vec a, arma::vec aprime, arma::vec kappa, arma::vec tc, bool check_pars);
RcppExport SEXP _rddm_weibull_bound_vec(SEXP tSEXP, SEXP aSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    rcpp_result_gen = Rcpp::wrap(weibull_bound_vec(t, a, aprime, kappa, tc, check_pars));
    return rcpp_result_gen;
END_RCPP
}
// ddm_integral_fpt
arma::mat ddm_integral_fpt(double v, double a, double t0, double z, double dc, double sv, double sz, double st0, double aprime, double kappa, double tc, double s, double sv_points, double sz_points, double st0_points, double dt, double max_time, bool use_weibull_bound, int n_threads);
RcppExport SEXP _rddm_ddm_integral_fpt(SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP szSEXP, SEXP st0SEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP sSEXP, SEXP sv_pointsSEXP, SEXP sz_pointsSEXP, SEXP st0_pointsSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP use_weibull_boundSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type sv_points(sv_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type sz_points(sz_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type st0_points(st0_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ddm_integral_fpt(v, a, t0, z, dc, sv, sz, st0, aprime, kappa, tc, s, sv_points, sz_points, st0_points, dt, max_time, use_weibull_bound, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// get_stimulus
arma::ivec get_stimulus(std::string sequence, double intensity, double dur, double isi, double dt);
RcppExport SEXP _rddm_get_stimulus(SEXP sequenceSEXP, SEXP intensitySEXP, SEXP durSEXP, SEXP isiSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< double >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< double >::type dur(durSEXP);
    Rcpp::traits::input_parameter< double >::type isi(isiSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(get_stimulus(sequence, intensity, dur, isi, dt));
    return rcpp_result_gen;
END_RCPP
}
// pulse_pmass_fpt
arma::mat pulse_pmass_fpt(arma::ivec stimulus, double v, double a, double t0, double z, double sv, double st0, double sz, double s, double lambda, double aprime, double kappa, double tc, double dt, double dx, double v_scale, bool use_weibull_bound);
RcppExport SEXP _rddm_pulse_pmass_fpt(SEXP stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP v_scaleSEXP, SEXP use_weibull_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec >::type stimulus(stimulusSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_pmass_fpt(stimulus, v, a, t0, z, sv, st0, sz, s, lambda, aprime, kappa, tc, dt, dx, v_scale, use_weibull_bound));
    return rcpp_result_gen;
END_RCPP
}
// pulse_trial_lik
double pulse_trial_lik(int choice, double rt, std::string stim_seq, double v, double a, double t0, double z, double sv, double st0, double sz, double s, double lambda, double aprime, double kappa, double tc, double dt, double dx, double v_scale, bool use_weibull_bound, double dur, double isi);
RcppExport SEXP _rddm_pulse_trial_lik(SEXP choiceSEXP, SEXP rtSEXP, SEXP stim_seqSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP v_scaleSEXP, SEXP use_weibull_boundSEXP, SEXP durSEXP, SEXP isiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type choice(choiceSEXP);
    Rcpp::traits::input_parameter< double >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< std::string >::type stim_seq(stim_seqSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< double >::type dur(durSEXP);
    Rcpp::traits::input_parameter< double >::type isi(isiSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_trial_lik(choice, rt, stim_seq, v, a, t0, z, sv, st0, sz, s, lambda, aprime, kappa, tc, dt, dx, v_scale, use_weibull_bound, dur, isi));
    return rcpp_result_gen;
END_RCPP
}
// pulse_nll
double pulse_nll(arma::vec choices, arma::vec rt, std::vector<std::string> stim_seq, arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec s, arma::vec lambda, arma::vec aprime, arma::vec kappa, arma::vec tc, bool check_pars, double dt, double dx, double v_scale, bool use_weibull_bound, double dur, double isi, int n_threads);
RcppExport SEXP _rddm_pulse_nll(SEXP choicesSEXP, SEXP rtSEXP, SEXP stim_seqSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP v_scaleSEXP, SEXP use_weibull_boundSEXP, SEXP durSEXP, SEXP isiSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type stim_seq(stim_seqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz(szSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< double >::type dur(durSEXP);
    Rcpp::traits::input_parameter< double >::type isi(isiSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_nll(choices, rt, stim_seq, v, a, t0, z, sv, st0, sz, s, lambda, aprime, kappa, tc, check_pars, dt, dx, v_scale, use_weibull_bound, dur, isi, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// pulse_predict
DataFrame pulse_predict(int n, std::vector<std::string> stim_seq, arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec d, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec s, arma::vec lambda, arma::vec aprime, arma::vec kappa, arma::vec tc, bool check_pars, double dt, double dx, double v_scale, bool use_weibull_bound, double dur, double isi, int n_threads);
RcppExport SEXP _rddm_pulse_predict(SEXP nSEXP, SEXP stim_seqSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP v_scaleSEXP, SEXP use_weibull_boundSEXP, SEXP durSEXP, SEXP isiSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type stim_seq(stim_seqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz(szSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< double >::type dur(durSEXP);
    Rcpp::traits::input_parameter< double >::type isi(isiSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_predict(n, stim_seq, v, a, t0, z, d, sv, st0, sz, s, lambda, aprime, kappa, tc, check_pars, dt, dx, v_scale, use_weibull_bound, dur, isi, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_ddm
DataFrame sim_ddm(int n, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double aprime, double kappa, double tc, double s, double dt, double max_time, bool use_weibull_bound, int n_threads);
RcppExport SEXP _rddm_sim_ddm(SEXP nSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP sSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP use_weibull_boundSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ddm(n, v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, s, dt, max_time, use_weibull_bound, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_ddm_vec
DataFrame sim_ddm_vec(arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec aprime, arma::vec kappa, arma::vec tc, double s, double dt, double max_time, bool use_weibull_bound, bool check_pars, int n_threads);
RcppExport SEXP _rddm_sim_ddm_vec(SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP sSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP use_weibull_boundSEXP, SEXP check_parsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz(szSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ddm_vec(v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, s, dt, max_time, use_weibull_bound, check_pars, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_pulse
DataFrame sim_pulse(int n, arma::ivec up_stimulus, arma::ivec down_stimulus, double v, double a, double t0, double z, double d, double sv, double st0, double sz, double s, double lambda, double aprime, double kappa, double tc, double dt, double v_scale, bool use_weibull_bound, int n_threads);
RcppExport SEXP _rddm_sim_pulse(SEXP nSEXP, SEXP up_stimulusSEXP, SEXP down_stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP dtSEXP, SEXP v_scaleSEXP, SEXP use_weibull_boundSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type up_stimulus(up_stimulusSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type down_stimulus(down_stimulusSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type use_weibull_bound(use_weibull_boundSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pulse(n, up_stimulus, down_stimulus, v, a, t0, z, d, sv, st0, sz, s, lambda, aprime, kappa, tc, dt, v_scale, use_weibull_bound, n_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rddm_hyperbolic_ratio_bound", (DL_FUNC) &_rddm_hyperbolic_ratio_bound, 4},
    {"_rddm_hyperbolic_ratio_bound_vec", (DL_FUNC) &_rddm_hyperbolic_ratio_bound_vec, 5},
    {"_rddm_weibull_bound", (DL_FUNC) &_rddm_weibull_bound, 5},
    {"_rddm_weibull_bound_vec", (DL_FUNC) &_rddm_weibull_bound_vec, 6},
    {"_rddm_ddm_integral_fpt", (DL_FUNC) &_rddm_ddm_integral_fpt, 19},
    {"_rddm_get_stimulus", (DL_FUNC) &_rddm_get_stimulus, 5},
    {"_rddm_pulse_pmass_fpt", (DL_FUNC) &_rddm_pulse_pmass_fpt, 17},
    {"_rddm_pulse_trial_lik", (DL_FUNC) &_rddm_pulse_trial_lik, 21},
    {"_rddm_pulse_nll", (DL_FUNC) &_rddm_pulse_nll, 23},
    {"_rddm_pulse_predict", (DL_FUNC) &_rddm_pulse_predict, 23},
    {"_rddm_sim_ddm", (DL_FUNC) &_rddm_sim_ddm, 17},
    {"_rddm_sim_ddm_vec", (DL_FUNC) &_rddm_sim_ddm_vec, 17},
    {"_rddm_sim_pulse", (DL_FUNC) &_rddm_sim_pulse, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_rddm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
