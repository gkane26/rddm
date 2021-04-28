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
arma::mat ddm_integral_fpt(double v, double a, double t0, double z, double dc, double sv, double sz, double st0, double aprime, double kappa, double tc, double s, double sv_points, double sz_points, double st0_points, double dt, double max_time, int bounds, int n_threads);
RcppExport SEXP _rddm_ddm_integral_fpt(SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP szSEXP, SEXP st0SEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP sSEXP, SEXP sv_pointsSEXP, SEXP sz_pointsSEXP, SEXP st0_pointsSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP boundsSEXP, SEXP n_threadsSEXP) {
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
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ddm_integral_fpt(v, a, t0, z, dc, sv, sz, st0, aprime, kappa, tc, s, sv_points, sz_points, st0_points, dt, max_time, bounds, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// pulse_fp_fpt
arma::mat pulse_fp_fpt(arma::mat stimulus, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double s, double lambda, double aprime, double kappa, double tc, double v_scale, double dt, double dx, int bounds);
RcppExport SEXP _rddm_pulse_fp_fpt(SEXP stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type stimulus(stimulusSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_fp_fpt(stimulus, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, v_scale, dt, dx, bounds));
    return rcpp_result_gen;
END_RCPP
}
// pulse_trial_lik
double pulse_trial_lik(int choice, double rt, arma::mat stimulus, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double s, double lambda, double aprime, double kappa, double tc, double v_scale, double dt, double dx, int bounds);
RcppExport SEXP _rddm_pulse_trial_lik(SEXP choiceSEXP, SEXP rtSEXP, SEXP stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type choice(choiceSEXP);
    Rcpp::traits::input_parameter< double >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type stimulus(stimulusSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_trial_lik(choice, rt, stimulus, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, v_scale, dt, dx, bounds));
    return rcpp_result_gen;
END_RCPP
}
// pulse_nll
double pulse_nll(arma::vec choices, arma::vec rt, List stimuli, arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec s, arma::vec lambda, arma::vec aprime, arma::vec kappa, arma::vec tc, bool check_pars, double v_scale, double dt, double dx, int bounds, int n_threads);
RcppExport SEXP _rddm_pulse_nll(SEXP choicesSEXP, SEXP rtSEXP, SEXP stimuliSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP boundsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< List >::type stimuli(stimuliSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz(szSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_nll(choices, rt, stimuli, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, check_pars, v_scale, dt, dx, bounds, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// pulse_predict
DataFrame pulse_predict(int n, List stimuli, arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec s, arma::vec lambda, arma::vec aprime, arma::vec kappa, arma::vec tc, bool check_pars, double dt, double dx, int bounds, int n_threads);
RcppExport SEXP _rddm_pulse_predict(SEXP nSEXP, SEXP stimuliSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP, SEXP dtSEXP, SEXP dxSEXP, SEXP boundsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type stimuli(stimuliSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dc(dcSEXP);
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
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_predict(n, stimuli, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, check_pars, dt, dx, bounds, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// pulse_trial_stimulus
arma::rowvec pulse_trial_stimulus(arma::vec stim_seq, double dur, double isi, double pre_stim, double dt);
RcppExport SEXP _rddm_pulse_trial_stimulus(SEXP stim_seqSEXP, SEXP durSEXP, SEXP isiSEXP, SEXP pre_stimSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type stim_seq(stim_seqSEXP);
    Rcpp::traits::input_parameter< double >::type dur(durSEXP);
    Rcpp::traits::input_parameter< double >::type isi(isiSEXP);
    Rcpp::traits::input_parameter< double >::type pre_stim(pre_stimSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_trial_stimulus(stim_seq, dur, isi, pre_stim, dt));
    return rcpp_result_gen;
END_RCPP
}
// sim_ddm
List sim_ddm(int n, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double aprime, double kappa, double tc, double s, double dt, double max_time, int bounds, int n_threads);
RcppExport SEXP _rddm_sim_ddm(SEXP nSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP sSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP boundsSEXP, SEXP n_threadsSEXP) {
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
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ddm(n, v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, s, dt, max_time, bounds, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_ddm_vec
DataFrame sim_ddm_vec(arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec aprime, arma::vec kappa, arma::vec tc, double s, double dt, double max_time, int bounds, bool check_pars, int n_threads);
RcppExport SEXP _rddm_sim_ddm_vec(SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP sSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP boundsSEXP, SEXP check_parsSEXP, SEXP n_threadsSEXP) {
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
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ddm_vec(v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, s, dt, max_time, bounds, check_pars, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_pulse
List sim_pulse(int n, arma::mat stimulus, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double s, double lambda, double aprime, double kappa, double tc, double v_scale, double dt, int bounds, int n_threads);
RcppExport SEXP _rddm_sim_pulse(SEXP nSEXP, SEXP stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP boundsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type stimulus(stimulusSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pulse(n, stimulus, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, v_scale, dt, bounds, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_pulse_vec
List sim_pulse_vec(int n, List stimuli, arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec s, arma::vec lambda, arma::vec aprime, arma::vec kappa, arma::vec tc, double v_scale, double dt, int bounds, bool check_pars, int n_threads);
RcppExport SEXP _rddm_sim_pulse_vec(SEXP nSEXP, SEXP stimuliSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP boundsSEXP, SEXP check_parsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type stimuli(stimuliSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz(szSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pulse_vec(n, stimuli, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, v_scale, dt, bounds, check_pars, n_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rddm_hyperbolic_ratio_bound", (DL_FUNC) &_rddm_hyperbolic_ratio_bound, 4},
    {"_rddm_hyperbolic_ratio_bound_vec", (DL_FUNC) &_rddm_hyperbolic_ratio_bound_vec, 5},
    {"_rddm_weibull_bound", (DL_FUNC) &_rddm_weibull_bound, 5},
    {"_rddm_weibull_bound_vec", (DL_FUNC) &_rddm_weibull_bound_vec, 6},
    {"_rddm_ddm_integral_fpt", (DL_FUNC) &_rddm_ddm_integral_fpt, 19},
    {"_rddm_pulse_fp_fpt", (DL_FUNC) &_rddm_pulse_fp_fpt, 18},
    {"_rddm_pulse_trial_lik", (DL_FUNC) &_rddm_pulse_trial_lik, 20},
    {"_rddm_pulse_nll", (DL_FUNC) &_rddm_pulse_nll, 22},
    {"_rddm_pulse_predict", (DL_FUNC) &_rddm_pulse_predict, 20},
    {"_rddm_pulse_trial_stimulus", (DL_FUNC) &_rddm_pulse_trial_stimulus, 5},
    {"_rddm_sim_ddm", (DL_FUNC) &_rddm_sim_ddm, 17},
    {"_rddm_sim_ddm_vec", (DL_FUNC) &_rddm_sim_ddm_vec, 17},
    {"_rddm_sim_pulse", (DL_FUNC) &_rddm_sim_pulse, 19},
    {"_rddm_sim_pulse_vec", (DL_FUNC) &_rddm_sim_pulse_vec, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_rddm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
