# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Collapsing boundary functions
#' 
#' @description evaluate diffusion model boundary using the hyperbolic ratio or weibull functions
#'
#' @param t vector; time points to evaluate boundary
#' @param a numeric or vector; initial boundary
#' @param aprime numeric or vector; degree of collapse (weibull only)
#' @param kappa numeric or vector; slope of collapse
#' @param tc numeric or vector; time constant of collapse
#' 
#' @return 
#' column vector with boundary at time t, 
#' or a time point x parameter vector matrix of boundaries 
#' (each column represents a time varying boundary for a given parameter set)
#' 
#' @name bounds
NULL

#' @rdname bounds
#' @export
linear_bound <- function(t, a, kappa = 0, tc = 0) {
    .Call('_rddm_linear_bound', PACKAGE = 'rddm', t, a, kappa, tc)
}

#' @rdname bounds
#' @export
linear_bound_vec <- function(t, a, kappa = 0L, tc = 0L, check_pars = TRUE) {
    .Call('_rddm_linear_bound_vec', PACKAGE = 'rddm', t, a, kappa, tc, check_pars)
}

#' @rdname bounds
#' @export
hyperbolic_ratio_bound <- function(t, a, kappa = 0, tc = .25) {
    .Call('_rddm_hyperbolic_ratio_bound', PACKAGE = 'rddm', t, a, kappa, tc)
}

#' @rdname bounds
#' @export
hyperbolic_ratio_bound_vec <- function(t, a, kappa = 0L, tc = 0L, check_pars = TRUE) {
    .Call('_rddm_hyperbolic_ratio_bound_vec', PACKAGE = 'rddm', t, a, kappa, tc, check_pars)
}

#' @rdname bounds
#' @export
weibull_bound <- function(t, a, aprime = 0, kappa = 1, tc = .25) {
    .Call('_rddm_weibull_bound', PACKAGE = 'rddm', t, a, aprime, kappa, tc)
}

#' @rdname bounds
#' @export
weibull_bound_vec <- function(t, a, aprime = 0L, kappa = 0L, tc = 0L, check_pars = TRUE) {
    .Call('_rddm_weibull_bound_vec', PACKAGE = 'rddm', t, a, aprime, kappa, tc, check_pars)
}

#' Simulate drift diffusion model with fixed or collapsing boundary
#'
#' @param v numeric; drift rate
#' @param a numeric; initial boundary
#' @param t0 numeric; non-decision time
#' @param z numeric; starting point, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
#' @param aprime numeric; degree of collapse, default = 0
#' @param kappa numeric; slope of collapse, default = 0
#' @param tc numeric; time constant of collapse, default = .25
#' @param s numeric; standard deviation in wiener diffusion noise, default = 1
#' @param sv_points integer; number of points to approximate integral over drift rate variability, default = 19
#' @param sz_points integer; number of points to approximate integral over starting point variability, default = 11
#' @param st0_points integer; number of points to approximate integral over non-decision time variability, default = 11
#' @param dt numeric; time step of simulation, default = .01
#' @param max_time numeric; max time of simulation, default = 10
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
#' @param n_threads integer; number of threads to run in parallel, default = 1
#' 
#' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
#' 
#' @export
ddm_integral_fpt <- function(v, a, t0, z = .5, dc = 0, sv = 0, sz = 0, st0 = 0, aprime = 0, kappa = 0, tc = .25, s = 1, sv_points = 19, sz_points = 11, st0_points = 11, dt = .01, max_time = 10, bounds = 0L, n_threads = 1L) {
    .Call('_rddm_ddm_integral_fpt', PACKAGE = 'rddm', v, a, t0, z, dc, sv, sz, st0, aprime, kappa, tc, s, sv_points, sz_points, st0_points, dt, max_time, bounds, n_threads)
}

#' random normal
#'
#' @description random normal
#'
#' @param n int; number of random draws
#'
#' @return random numbers
#'
#' @name zrand
#' @rdname zrand
#' 
#' @export
zrandn <- function(n) {
    .Call('_rddm_zrandn', PACKAGE = 'rddm', n)
}

#' @rdname zrand
#' @export
zrandseed <- function(s) {
    invisible(.Call('_rddm_zrandseed', PACKAGE = 'rddm', s))
}

#' Get first passage time distribution of pulse diffusion model by simulating probability mass
#'
#' @param stimulus matrix; stimulus to simulate (row 1 is evidence to upper boundary, row 2 to lower boundary)
#' @param v numeric; drift rate
#' @param a numeric; initial boundary
#' @param t0 numeric; non-decision time
#' @param z numeric; starting point, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point (percentage of boundary a). Normal with mean=z and sd = a*sz. 0 < sz < 1, default = 0
#' @param lambda numeric; O-U process slope
#' @param aprime numeric; degree of collapse, default = 0
#' @param kappa numeric; slope of collapse, default = 0
#' @param tc numeric; time constant of collapse, default = .25
#' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
#' @param uslope numeric; urgency scaling factor, default = 0;
#' @param umag numeric; urgency magnitude, default = 0;
#' @param udelay numeric; urgency delay, default = 0;
#' @param s numeric; standard deviation in wiener diffusion noise
#' @param dt numeric; time step of simulation, default = .001
#' @param xbins numeric; number of evidence bins, default = 100
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
#' @param urgency int; 0 for none, 1 for linear, 2 for logistic
#'
#' @return data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence
#'
#' @export
pulse_fp_fpt <- function(stimulus, v, a, t0, z = 0.5, dc = 0, sv = 0, st0 = 0, sz = 0, lambda = 0, aprime = 0, kappa = 0, tc = .25, uslope = 0, umag = 0, udelay = 0, s = 1, v_scale = 1, dt = .001, xbins = 200L, bounds = 0L, urgency = 0L) {
    .Call('_rddm_pulse_fp_fpt', PACKAGE = 'rddm', stimulus, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, v_scale, dt, xbins, bounds, urgency)
}

#' Get pulse model likelihood for a given trial
#'
#' @param choice int; decision on trial, 0 for lower boundary, 1 for upper
#' @param rt numeric; response time on trial
#' @param stimulus matrix; stimulus to simulate (row 1 is evidence to upper boundary, row 2 to lower boundary)
#' @param v numeric; drift rate
#' @param a numeric; initial boundary
#' @param t0 numeric; non-decision time
#' @param z numeric; starting point, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point (percentage of boundary a). Normal with mean=z and sd = a*sz. 0 < sz < 1, default = 0
#' @param lambda numeric; O-U process slope
#' @param aprime numeric; degree of collapse, default = 0
#' @param kappa numeric; slope of collapse, default = 0
#' @param tc numeric; time constant of collapse, default = .25
#' @param uslope numeric; urgency scaling factor, default = 0;
#' @param umag numeric; urgency magnitude, default = 0;
#' @param udelay numeric; urgency delay, default = 0;
#' @param s numeric; standard deviation in wiener diffusion noise
#' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
#' @param dt numeric; time step of simulation, default = .001
#' @param xbins numeric; number of evidence bins, default = 100
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
#' @param urgency int; 0 for none, 1 for linear, 2 for logistic
#'
#' @return probability of choice and rt for trial given pulse model parameters
#'
#' @export
pulse_trial_lik <- function(choice, rt, stimulus, v, a, t0, z = 0.5, dc = 0, sv = 0, st0 = 0, sz = 0, lambda = 0, aprime = 0, kappa = 0, tc = .25, uslope = 0, umag = 0, udelay = 0, s = 1, v_scale = 1, dt = .001, xbins = 200L, bounds = 0L, urgency = 0L) {
    .Call('_rddm_pulse_trial_lik', PACKAGE = 'rddm', choice, rt, stimulus, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, v_scale, dt, xbins, bounds, urgency)
}

#' Get pulse model negative log likelihood for a set of trials
#'
#' @param choice integer; vector of decisions, 0 for lower boundary, 1 for upper
#' @param rt numeric; vector of response times
#' @param stimuli array; 2 x timepoints x trials array of pulse stimuli
#' @param up_sequence vector of strings; string of 0s, and 1s of stimulus values (0 no evidence, 1 to upper). If down_sequence not specified, (0 to lower, 1 to upper).
#' @param down_sequence vector of strings; string of 0s, and 1s of stimulus values (0 is no evidence, 1 to lower). If not specified, up_sequence is (0 to lower, 1 to upper)
#' @param v numeric; drift rate, either single value or vector for each trial
#' @param a numeric; initial boundary, either single value or vector for each trial
#' @param t0 numeric; non-decision time, either single value or vector for each trial
#' @param z numeric; starting point, , either single value or vector for each trial, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv numeric; standard deviation of variability in drift rate, either single value or vector for each trial, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time, either single value or vector for each trial. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point (percentage of boundary a), either single value or vector for each trial. Normal with mean=z and sd = a*sz. 0 < sz < 1, default = 0
#' @param lambda numeric; O-U process slope, either single value or vector for each trial
#' @param aprime numeric; degree of collapse, either single value or vector for each trial, default = 0
#' @param kappa numeric; slope of collapse, either single value or vector for each trial, default = 0
#' @param tc numeric; time constant of collapse, either single value or vector for each trial, default = .25
#' @param uslope numeric; urgency scaling factor, default = 0;
#' @param umag numeric; urgency magnitude, default = 0;
#' @param udelay numeric; urgency delay, default = 0;
#' @param s numeric; standard deviation in wiener diffusion noise, either single value or vector for each trial
#' @param check_pars logical; if True, check that parameters are vectors of the same length as choices and rts. Must be true if providing scalar parameters. default = true
#' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
#' @param dt numeric; time step of simulation, default = .002
#' @param xbins integer; number of evidence bins, default = 100
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
#' @param urgency int: 0 for none, 1 for linear, 2 for logistic
#' @param n_threads int; number of threads (trials) to run in parallel
#'
#' @return negative log likelihood of all choices and rts given pulse model parameters
#'
#' @export
pulse_nll <- function(choices, rt, stimuli, v, a, t0, z = 0L, dc = 0L, sv = 0L, st0 = 0L, sz = 0L, lambda = 0L, aprime = 0L, kappa = 0L, tc = 0L, uslope = 0L, umag = 0L, udelay = 0L, s = 0L, check_pars = TRUE, v_scale = 1, dt = .001, xbins = 200L, bounds = 0L, urgency = 0L, n_threads = 1L) {
    .Call('_rddm_pulse_nll', PACKAGE = 'rddm', choices, rt, stimuli, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, check_pars, v_scale, dt, xbins, bounds, urgency, n_threads)
}

#' Get predicted behavior from pulse model
#'
#' @param n int; number of predicted samples to take per stimulus
#' @param stimuli list; list of stimulus matrices
#' @param v numeric; drift rate, either single value or vector for each trial
#' @param a numeric; initial boundary, either single value or vector for each trial
#' @param t0 numeric; non-decision time, either single value or vector for each trial
#' @param z numeric; starting point, , either single value or vector for each trial, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv numeric; standard deviation of variability in drift rate, either single value or vector for each trial, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time, either single value or vector for each trial. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point, either single value or vector for each trial. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
#' @param s numeric; standard deviation in wiener diffusion noise, either single value or vector for each trial, default = 1
#' @param lambda numeric; O-U process slope, either single value or vector for each trial
#' @param aprime numeric; degree of collapse, either single value or vector for each trial, default = 0
#' @param kappa numeric; slope of collapse, either single value or vector for each trial, default = 0
#' @param tc numeric; time constant of collapse, either single value or vector for each trial, default = .25
#' @param check_pars logical; if True, check that parameters are vectors of the same length as choices and rts. Must be true if providing scalar parameters. default = true
#' @param dt numeric; time step of simulation, default = .001
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
#' @param n_threads int; number of threads (trials) to run in parallel
#'
#' @return data frame with two columns: response (1 for upper boundary, 0 for lower), response time
#'
#' @export
pulse_predict <- function(n, stimuli, v, a, t0, z = 0L, dc = 0L, sv = 0L, st0 = 0L, sz = 0L, s = 0L, lambda = 0L, aprime = 0L, kappa = 0L, tc = 0L, check_pars = TRUE, dt = .001, dx = .05, bounds = 0L, n_threads = 1L) {
    .Call('_rddm_pulse_predict', PACKAGE = 'rddm', n, stimuli, v, a, t0, z, dc, sv, st0, sz, s, lambda, aprime, kappa, tc, check_pars, dt, dx, bounds, n_threads)
}

#' Get pulse stimulus
#'
#' @description get full pulse stimulus
#'
#' @param stim_seq numeric vector; vecotr of timepoint by timepoint evidence
#' @param dur numeric; duration of pulse, default = .01
#' @param isi numeric; inter-pulse interval, default = .1
#' @param pre_stim numeric; time before pulse within stimulus bin
#' @param dt numeric; timestep for simulation, default = .001
#'
#' @return stimulus train vector
#'
#' @export
pulse_trial_stimulus <- function(stim_seq, dur = 0.01, isi = 0.1, pre_stim = 0, dt = 0.001) {
    .Call('_rddm_pulse_trial_stimulus', PACKAGE = 'rddm', stim_seq, dur, isi, pre_stim, dt)
}

#' Simulate drift diffusion model with fixed or collapsing boundary, with or without urgency
#'
#' @param n integer; number of decisions to simulate
#' @param v numeric; drift rate
#' @param a numeric; initial boundary
#' @param t0 numeric; non-decision time
#' @param z numeric; starting point, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point (percentage of boundary a). Uniform from [z-sz/2, z+sz/2], 0 < sz < 1, default = 0
#' @param aprime numeric; degree of collapse, default = 0
#' @param kappa numeric; slope of collapse, default = 0
#' @param tc numeric; time constant of collapse, default = .25
#' @param uslope numeric; urgency scaling factor, default = 0;
#' @param umag numeric; urgency magnitude, default = 0;
#' @param udelay numeric; urgency delay, default = 0;
#' @param s numeric; standard deviation in wiener diffusion noise, default = 1
#' @param dt numeric; time step of simulation, default = .001
#' @param max_time numeric; max time of simulation, default = 10
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
#' @param urgency int: 0 for none, 1 for linear, 2 for logistic
#' @param n_threads integer; number of threads to run in parallel, default = 1
#' @param return_accu bool; if True, return full trajectory of accumulators
#' @param seed integer; set random seed
#'
#' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
#' 
#' @export
sim_ddm <- function(n, v, a, t0, z = .5, dc = 0, sv = 0, st0 = 0, sz = 0, aprime = 0, kappa = 0, tc = .25, uslope = 0, umag = 0, udelay = 0, s = 1, dt = .001, max_time = 10, bounds = 0L, urgency = 0L, n_threads = 1L, return_accu = FALSE, seed = -1L) {
    .Call('_rddm_sim_ddm', PACKAGE = 'rddm', n, v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, uslope, umag, udelay, s, dt, max_time, bounds, urgency, n_threads, return_accu, seed)
}

#' Simulate drift diffusion model with fixed or collapsing boundary
#'
#' @param v vector; drift rate
#' @param a vector; initial boundary
#' @param t0 vector; non-decision time
#' @param z vector; starting point, 0 < z < 1, default = .5
#' @param dc vector; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0
#' @param sv vector; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 vector; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz vector; variability in starting point (percentage of boundary a). Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0
#' @param aprime vector; degree of collapse, default = 0
#' @param kappa vector; slope of collapse, default = 0
#' @param tc vector; time constant of collapse, default = .25
#' @param uslope numeric; urgency scaling factor, default = 0;
#' @param umag numeric; urgency magnitude, default = 0;
#' @param udelay numeric; urgency delay, default = 0;
#' @param s numeric; standard deviation in wiener diffusion noise, default = 1
#' @param dt numeric; time step of simulation, default = .001
#' @param max_time numeric; max time of simulation, default = 10
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
#' @param urgency int: 0 for none, 1 for linear, 2 for logistic
#' @param check_pars bool; if True (default) check parameter vector lengths and default values
#' @param n_threads integer; number of threads to run in parallel, default = 1
#' @param return_accu bool; if True, return full trajectory of accumulators
#' @param seed integer; set random seed
#'
#' @return data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
#'
#' @export
sim_ddm_vec <- function(v, a, t0, z = 0L, dc = 0L, sv = 0L, st0 = 0L, sz = 0L, aprime = 0L, kappa = 0L, tc = 0L, uslope = 0L, udelay = 0L, umag = 0L, s = 0L, dt = .001, max_time = 10, bounds = 0L, urgency = 0L, check_pars = TRUE, n_threads = 1L, return_accu = FALSE, seed = -1L) {
    .Call('_rddm_sim_ddm_vec', PACKAGE = 'rddm', v, a, t0, z, dc, sv, st0, sz, aprime, kappa, tc, uslope, udelay, umag, s, dt, max_time, bounds, urgency, check_pars, n_threads, return_accu, seed)
}

#' Simulate EvAcc (2-accumulator) model with fixed or collapsing boundary
#'
#' @param n integer; number of simulations for each stimulus
#' @param stimulus array; stimulus to simulate, 1 or 2 rows X timepoints x number of stimuli. If 2 rows, row 1 is evidence to upper boundary, row 2 to lower boundary.
#' @param a numeric; initial boundary
#' @param t0 numeric; non-decision time
#' @param z numeric; boundary separation bias, left boundary = a - z, right = a + z
#' @param dc numeric; drift criterion, left = v + dc, right = v - dc
#' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting bias (drawn from gaussian), sz > 0, default = 0
#' @param s numeric; diffusion noise
#' @param lambda numeric; accumulator recurrent weight
#' @param mi numeric; mutual inhibition
#' @param sv2 numeric; sampling noise for second layer
#' @param s2 numeric; diffusion noise for layer 2
#' @param lambda2 numeric; second layer accumulator recurrent weight
#' @param mi2 numeric; mutual inhibition in second layer. Set mi=0 and mi2 > 0 to implement static sample model from Scott et al., 2015
#' @param aprime numeric; degree of collapse, default = 0
#' @param kappa numeric; slope of collapse, default = 1
#' @param tc numeric; time constant of collapse, default = .25
#' @param v numeric; drift rate. Default = 1, recommended to leave this parameter fixed
#' @param two_layer bool; if true, use two layer accumulator model
#' @param accumulator_gain; if true, gain of first accumulator layer scales with value
#' @param scalar_stimulus_noise bool; if true, drift rate variability scales with accumulator value (see Koay et al., 2020)
#' @param scalar_diffusion_noise bool; if true, diffusion noise scales with accumulator value
#' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv and criterion dc are multiplied by this number
#' @param dt numeric; time step of simulation, default = .001
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds
#' @param n_threads integer; number of threads to run in parallel, default = 1
#' @param return_accu bool; if True, return full trajectory of accumulators
#'
#' @return List containing 1) data frame with three columns: response (1 for right, 0 for left), response time, and evidence and 2) matrix with full accumulator trajectories
#' 
#' @export
sim_evacc <- function(n, stimuli, a, t0, z = 0, dc = 0, sv = 0, st0 = 0, sz = 0, s = 1, lambda = 0, mi = 0, sv2 = 0, s2 = 1, lambda2 = 1, mi2 = 0, aprime = 0, kappa = 0, tc = .25, v = 1, accumulator_gain = FALSE, two_layer = FALSE, scalar_stimulus_noise = FALSE, scalar_diffusion_noise = FALSE, v_scale = 1, dt = .001, bounds = 0L, n_threads = 1L, return_accu = FALSE) {
    .Call('_rddm_sim_evacc', PACKAGE = 'rddm', n, stimuli, a, t0, z, dc, sv, st0, sz, s, lambda, mi, sv2, s2, lambda2, mi2, aprime, kappa, tc, v, accumulator_gain, two_layer, scalar_stimulus_noise, scalar_diffusion_noise, v_scale, dt, bounds, n_threads, return_accu)
}

#' Simulate pulse diffusion model with fixed or collapsing boundary
#'
#' @param n integer; number of decisions to simulate
#' @param stimuli cube; stimuli to simulate, 2 rows X timepoints x trials.
#' @param v numeric; drift rate
#' @param a numeric; initial boundary
#' @param t0 numeric; non-decision time
#' @param z numeric; starting point, 0 < z < 1, default = .5
#' @param dc numeric; drift criterion, upper drift = v, lower drift = v-d
#' @param sv numeric; standard deviation of variability in drift rate, sv >= 0, default = 0
#' @param st0 numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0
#' @param sz numeric; variability in starting point (percentage of boundary a). Normal with mean=z and sd = a*sz. 0 < sz < 1, default = 0
#' @param lambda numeric; O-U process slope
#' @param aprime numeric; degree of collapse, default = 0
#' @param kappa numeric; slope of collapse, default = 1
#' @param tc numeric; time constant of collapse, default = .25
#' @param uslope numeric; urgency scaling factor, default = 0;
#' @param umag numeric; urgency magnitude, default = 0;
#' @param udelay numeric; urgency delay, default = 0;
#' @param s numeric; standard deviation in wiener diffusion noise
#' @param v_scale numeric; scale for the drift rate. drift rate v and variability sv are multiplied by this number
#' @param dt numeric; time step of simulation, default = .001
#' @param bounds int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds, 3 for linear
#' @param urgency int: 0 for none, 1 for linear, 2 for logistic
#' @param n_threads integer; number of threads to run in parallel, default = 1
#' @param return_accu bool; if True, return full trajectory of accumulators
#' @param seed int; set random seed (for zrandn)
#'
#' @return List containing 1) data frame with three columns: response (1 for upper boundary, 0 for lower), response time, and evidence and 2) matrix with full accumulator trajectories
#' 
#' @export
sim_pulse <- function(n, stimuli, v, a, t0, z = .5, dc = 0, sv = 0, st0 = 0, sz = 0, lambda = 0, aprime = 0, kappa = 0, tc = .25, uslope = 0, umag = 0, udelay = 0, s = 1, v_scale = 1, dt = .001, bounds = 0L, urgency = 0L, n_threads = 1L, return_accu = FALSE, seed = -1L) {
    .Call('_rddm_sim_pulse', PACKAGE = 'rddm', n, stimuli, v, a, t0, z, dc, sv, st0, sz, lambda, aprime, kappa, tc, uslope, umag, udelay, s, v_scale, dt, bounds, urgency, n_threads, return_accu, seed)
}

#' Urgency functions
#' 
#' @description evaluate diffusion model urgency signals (linear and logistic)
#'
#' @param t vector; time points to evaluate boundary
#' @param uslope; slope of urgency signal
#' @param udelay; delay to start of rising urgency
#' @param umag; magnitude of urgency signal (only for logistic)
#' 
#' @return column vector with urgency at time t, 
#' 
#' @name urgency
NULL

#' @rdname urgency
#' @export
linear_urgency <- function(t, uslope = 0, udelay = 0, umag = 0) {
    .Call('_rddm_linear_urgency', PACKAGE = 'rddm', t, uslope, udelay, umag)
}

#' @rdname urgency
#' @export
linear_urgency_vec <- function(t, uslope = 0L, udelay = 0L, umag = 0L, check_pars = TRUE) {
    .Call('_rddm_linear_urgency_vec', PACKAGE = 'rddm', t, uslope, udelay, umag, check_pars)
}

#' @rdname urgency
#' @export
logistic_urgency <- function(t, uslope = 0, udelay = 0, umag = 0) {
    .Call('_rddm_logistic_urgency', PACKAGE = 'rddm', t, uslope, udelay, umag)
}

#' @rdname urgency
#' @export
logistic_urgency_vec <- function(t, uslope = 0L, udelay = 0L, umag = 0L, check_pars = TRUE) {
    .Call('_rddm_logistic_urgency_vec', PACKAGE = 'rddm', t, uslope, udelay, umag, check_pars)
}

