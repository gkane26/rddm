// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// linear_bound
arma::vec linear_bound(arma::vec t, double a, double kappa, double tc);
RcppExport SEXP _rddm_linear_bound(SEXP tSEXP, SEXP aSEXP, SEXP kappaSEXP, SEXP tcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_bound(t, a, kappa, tc));
    return rcpp_result_gen;
END_RCPP
}
// linear_bound_vec
arma::mat linear_bound_vec(arma::vec t, arma::vec a, arma::vec kappa, arma::vec tc, bool check_pars);
RcppExport SEXP _rddm_linear_bound_vec(SEXP tSEXP, SEXP aSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP check_parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_bound_vec(t, a, kappa, tc, check_pars));
    return rcpp_result_gen;
END_RCPP
}
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
// zrandn
arma::vec zrandn(int n);
RcppExport SEXP _rddm_zrandn(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(zrandn(n));
    return rcpp_result_gen;
END_RCPP
}
// zrandseed
void zrandseed(unsigned long int s);
RcppExport SEXP _rddm_zrandseed(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned long int >::type s(sSEXP);
    zrandseed(s);
    return R_NilValue;
END_RCPP
}
// pulse_fp_fpt
arma::mat pulse_fp_fpt(arma::mat stimulus, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double lambda, double aprime, double kappa, double tc, double uslope, double umag, double udelay, double s, double v_scale, double dt, int xbins, int bounds, int urgency);
RcppExport SEXP _rddm_pulse_fp_fpt(SEXP stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP uslopeSEXP, SEXP umagSEXP, SEXP udelaySEXP, SEXP sSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP xbinsSEXP, SEXP boundsSEXP, SEXP urgencySEXP) {
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
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< double >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< double >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type xbins(xbinsSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type urgency(urgencySEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_fp_fpt(stimulus, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, v_scale, dt, xbins, bounds, urgency));
    return rcpp_result_gen;
END_RCPP
}
// pulse_trial_lik
double pulse_trial_lik(int choice, double rt, arma::mat stimulus, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double lambda, double aprime, double kappa, double tc, double uslope, double umag, double udelay, double s, double v_scale, double dt, int xbins, int bounds, int urgency);
RcppExport SEXP _rddm_pulse_trial_lik(SEXP choiceSEXP, SEXP rtSEXP, SEXP stimulusSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP uslopeSEXP, SEXP umagSEXP, SEXP udelaySEXP, SEXP sSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP xbinsSEXP, SEXP boundsSEXP, SEXP urgencySEXP) {
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
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< double >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< double >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type xbins(xbinsSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type urgency(urgencySEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_trial_lik(choice, rt, stimulus, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, v_scale, dt, xbins, bounds, urgency));
    return rcpp_result_gen;
END_RCPP
}
// pulse_nll
double pulse_nll(arma::vec choices, arma::vec rt, arma::cube stimuli, arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec lambda, arma::vec aprime, arma::vec kappa, arma::vec tc, arma::vec uslope, arma::vec umag, arma::vec udelay, arma::vec s, bool check_pars, double v_scale, double dt, int xbins, int bounds, int urgency, int n_threads);
RcppExport SEXP _rddm_pulse_nll(SEXP choicesSEXP, SEXP rtSEXP, SEXP stimuliSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP uslopeSEXP, SEXP umagSEXP, SEXP udelaySEXP, SEXP sSEXP, SEXP check_parsSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP xbinsSEXP, SEXP boundsSEXP, SEXP urgencySEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type stimuli(stimuliSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz(szSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type xbins(xbinsSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type urgency(urgencySEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(pulse_nll(choices, rt, stimuli, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, check_pars, v_scale, dt, xbins, bounds, urgency, n_threads));
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
List sim_ddm(int n, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double aprime, double kappa, double tc, double uslope, double umag, double udelay, double s, double dt, double max_time, int bounds, int urgency, int n_threads, bool return_accu, int seed);
RcppExport SEXP _rddm_sim_ddm(SEXP nSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP uslopeSEXP, SEXP umagSEXP, SEXP udelaySEXP, SEXP sSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP boundsSEXP, SEXP urgencySEXP, SEXP n_threadsSEXP, SEXP return_accuSEXP, SEXP seedSEXP) {
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
    Rcpp::traits::input_parameter< double >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< double >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< double >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type urgency(urgencySEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_accu(return_accuSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ddm(n, v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, uslope, umag, udelay, s, dt, max_time, bounds, urgency, n_threads, return_accu, seed));
    return rcpp_result_gen;
END_RCPP
}
// sim_ddm_vec
List sim_ddm_vec(arma::vec v, arma::vec a, arma::vec t0, arma::vec z, arma::vec dc, arma::vec sv, arma::vec st0, arma::vec sz, arma::vec aprime, arma::vec kappa, arma::vec tc, arma::vec uslope, arma::vec udelay, arma::vec umag, arma::vec s, double dt, double max_time, int bounds, int urgency, bool check_pars, int n_threads, bool return_accu, int seed);
RcppExport SEXP _rddm_sim_ddm_vec(SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP uslopeSEXP, SEXP udelaySEXP, SEXP umagSEXP, SEXP sSEXP, SEXP dtSEXP, SEXP max_timeSEXP, SEXP boundsSEXP, SEXP urgencySEXP, SEXP check_parsSEXP, SEXP n_threadsSEXP, SEXP return_accuSEXP, SEXP seedSEXP) {
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
    Rcpp::traits::input_parameter< arma::vec >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type urgency(urgencySEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_accu(return_accuSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ddm_vec(v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, uslope, udelay, umag, s, dt, max_time, bounds, urgency, check_pars, n_threads, return_accu, seed));
    return rcpp_result_gen;
END_RCPP
}
// sim_evacc
List sim_evacc(int n, arma::cube stimuli, double a, double t0, double z, double dc, double sv, double st0, double sz, double s, double lambda, double mi, double sv2, double s2, double lambda2, double mi2, double aprime, double kappa, double tc, double v, bool accumulator_gain, bool two_layer, bool scalar_stimulus_noise, bool scalar_diffusion_noise, double v_scale, double dt, int bounds, int n_threads, bool return_accu);
RcppExport SEXP _rddm_sim_evacc(SEXP nSEXP, SEXP stimuliSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP sSEXP, SEXP lambdaSEXP, SEXP miSEXP, SEXP sv2SEXP, SEXP s2SEXP, SEXP lambda2SEXP, SEXP mi2SEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP vSEXP, SEXP accumulator_gainSEXP, SEXP two_layerSEXP, SEXP scalar_stimulus_noiseSEXP, SEXP scalar_diffusion_noiseSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP boundsSEXP, SEXP n_threadsSEXP, SEXP return_accuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type stimuli(stimuliSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type mi(miSEXP);
    Rcpp::traits::input_parameter< double >::type sv2(sv2SEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< double >::type mi2(mi2SEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type accumulator_gain(accumulator_gainSEXP);
    Rcpp::traits::input_parameter< bool >::type two_layer(two_layerSEXP);
    Rcpp::traits::input_parameter< bool >::type scalar_stimulus_noise(scalar_stimulus_noiseSEXP);
    Rcpp::traits::input_parameter< bool >::type scalar_diffusion_noise(scalar_diffusion_noiseSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_accu(return_accuSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_evacc(n, stimuli, a, t0, z, dc, sv, st0, sz, s, lambda, mi, sv2, s2, lambda2, mi2, aprime, kappa, tc, v, accumulator_gain, two_layer, scalar_stimulus_noise, scalar_diffusion_noise, v_scale, dt, bounds, n_threads, return_accu));
    return rcpp_result_gen;
END_RCPP
}
// sim_pulse
List sim_pulse(int n, arma::cube stimuli, double v, double a, double t0, double z, double dc, double sv, double st0, double sz, double lambda, double aprime, double kappa, double tc, double uslope, double umag, double udelay, double s, double v_scale, double dt, int bounds, int urgency, int n_threads, bool return_accu, int seed);
RcppExport SEXP _rddm_sim_pulse(SEXP nSEXP, SEXP stimuliSEXP, SEXP vSEXP, SEXP aSEXP, SEXP t0SEXP, SEXP zSEXP, SEXP dcSEXP, SEXP svSEXP, SEXP st0SEXP, SEXP szSEXP, SEXP lambdaSEXP, SEXP aprimeSEXP, SEXP kappaSEXP, SEXP tcSEXP, SEXP uslopeSEXP, SEXP umagSEXP, SEXP udelaySEXP, SEXP sSEXP, SEXP v_scaleSEXP, SEXP dtSEXP, SEXP boundsSEXP, SEXP urgencySEXP, SEXP n_threadsSEXP, SEXP return_accuSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type stimuli(stimuliSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type st0(st0SEXP);
    Rcpp::traits::input_parameter< double >::type sz(szSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type aprime(aprimeSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type tc(tcSEXP);
    Rcpp::traits::input_parameter< double >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< double >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< double >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type v_scale(v_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int >::type urgency(urgencySEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_accu(return_accuSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pulse(n, stimuli, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, v_scale, dt, bounds, urgency, n_threads, return_accu, seed));
    return rcpp_result_gen;
END_RCPP
}
// linear_urgency
arma::vec linear_urgency(arma::vec t, double uslope, double udelay, double umag);
RcppExport SEXP _rddm_linear_urgency(SEXP tSEXP, SEXP uslopeSEXP, SEXP udelaySEXP, SEXP umagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< double >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< double >::type umag(umagSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_urgency(t, uslope, udelay, umag));
    return rcpp_result_gen;
END_RCPP
}
// linear_urgency_vec
arma::mat linear_urgency_vec(arma::vec t, arma::vec uslope, arma::vec udelay, arma::vec umag, bool check_pars);
RcppExport SEXP _rddm_linear_urgency_vec(SEXP tSEXP, SEXP uslopeSEXP, SEXP udelaySEXP, SEXP umagSEXP, SEXP check_parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_urgency_vec(t, uslope, udelay, umag, check_pars));
    return rcpp_result_gen;
END_RCPP
}
// logistic_urgency
arma::vec logistic_urgency(arma::vec t, double uslope, double udelay, double umag);
RcppExport SEXP _rddm_logistic_urgency(SEXP tSEXP, SEXP uslopeSEXP, SEXP udelaySEXP, SEXP umagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< double >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< double >::type umag(umagSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_urgency(t, uslope, udelay, umag));
    return rcpp_result_gen;
END_RCPP
}
// logistic_urgency_vec
arma::mat logistic_urgency_vec(arma::vec t, arma::vec uslope, arma::vec udelay, arma::vec umag, bool check_pars);
RcppExport SEXP _rddm_logistic_urgency_vec(SEXP tSEXP, SEXP uslopeSEXP, SEXP udelaySEXP, SEXP umagSEXP, SEXP check_parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type uslope(uslopeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type udelay(udelaySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type umag(umagSEXP);
    Rcpp::traits::input_parameter< bool >::type check_pars(check_parsSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_urgency_vec(t, uslope, udelay, umag, check_pars));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rddm_linear_bound", (DL_FUNC) &_rddm_linear_bound, 4},
    {"_rddm_linear_bound_vec", (DL_FUNC) &_rddm_linear_bound_vec, 5},
    {"_rddm_hyperbolic_ratio_bound", (DL_FUNC) &_rddm_hyperbolic_ratio_bound, 4},
    {"_rddm_hyperbolic_ratio_bound_vec", (DL_FUNC) &_rddm_hyperbolic_ratio_bound_vec, 5},
    {"_rddm_weibull_bound", (DL_FUNC) &_rddm_weibull_bound, 5},
    {"_rddm_weibull_bound_vec", (DL_FUNC) &_rddm_weibull_bound_vec, 6},
    {"_rddm_ddm_integral_fpt", (DL_FUNC) &_rddm_ddm_integral_fpt, 19},
    {"_rddm_zrandn", (DL_FUNC) &_rddm_zrandn, 1},
    {"_rddm_zrandseed", (DL_FUNC) &_rddm_zrandseed, 1},
    {"_rddm_pulse_fp_fpt", (DL_FUNC) &_rddm_pulse_fp_fpt, 22},
    {"_rddm_pulse_trial_lik", (DL_FUNC) &_rddm_pulse_trial_lik, 24},
    {"_rddm_pulse_nll", (DL_FUNC) &_rddm_pulse_nll, 26},
    {"_rddm_pulse_predict", (DL_FUNC) &_rddm_pulse_predict, 20},
    {"_rddm_pulse_trial_stimulus", (DL_FUNC) &_rddm_pulse_trial_stimulus, 5},
    {"_rddm_sim_ddm", (DL_FUNC) &_rddm_sim_ddm, 23},
    {"_rddm_sim_ddm_vec", (DL_FUNC) &_rddm_sim_ddm_vec, 23},
    {"_rddm_sim_evacc", (DL_FUNC) &_rddm_sim_evacc, 29},
    {"_rddm_sim_pulse", (DL_FUNC) &_rddm_sim_pulse, 25},
    {"_rddm_linear_urgency", (DL_FUNC) &_rddm_linear_urgency, 4},
    {"_rddm_linear_urgency_vec", (DL_FUNC) &_rddm_linear_urgency_vec, 5},
    {"_rddm_logistic_urgency", (DL_FUNC) &_rddm_logistic_urgency, 4},
    {"_rddm_logistic_urgency_vec", (DL_FUNC) &_rddm_logistic_urgency_vec, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_rddm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
