#####################################
# Helper functions for pulse model object

pulse_fp_obj <- function(pars,
                         dat=NULL,
                         transform_pars=F,
                         check_constraints=T,
                         debug=F,
                         ...){
  
  ### check constraints
  
  checks = private$objective_checks(pars,
                                    transform_pars,
                                    check_constraints,
                                    debug)
  pars = checks[[1]]
  pass = checks[[2]]
  
  if (!is.na(pass) & !pass) {
    nll = 1e10
    if(debug) cat("nll =", nll, "\n")
    return(nll)
  }
  
  if(is.null(dat)) {
    dat = self$data
  }
  
  ### loop through conditions to get likelihood
  
  nll = 0
  if (private$par_matrix[, .N] < dat[, .N]) {
    for(i in 1:private$par_matrix[, .N]){
      sub_dat = copy(dat)
      for(j in 1:length(private$sim_cond)){
        sub_dat = sub_dat[get(private$sim_cond[j]) == private$par_matrix[i, get(private$sim_cond[j])]]
      }
      sub_nll = do.call(pulse_nll, c(list(choice = sub_dat[, response],
                                          rt = sub_dat[, rt],
                                          stimuli = private$stim_list[[i]]),
                                     as.list(private$par_matrix[i, -(1:length(private$sim_cond))]),
                                     bounds=private$bounds,
                                     urgency=private$urgency,
                                     ...))
      nll = nll + sub_nll
    }
  } else {
    
    nll = do.call(pulse_nll, c(list(choice=dat[,response]),
                               list(rt=dat[, rt]),
                               stimuli=list(private$stim_list[[1]]),
                               as.list(private$par_matrix),
                               bounds=private$bounds,
                               ...))
  }
  
  if (is.nan(nll)) browser()
  
  if(debug) cat("nll =", round(nll, 3), "\n")
  
  nll
  
}


#' @noRd
#' @importFrom foreach %do% %dopar%
pulse_x2_obj = function(pars,
                        data_q=NULL,
                        n_sim=1,
                        transform_pars=F,
                        check_constraints=T,
                        debug=F, ...) {
  
  ### check constraints
  
  checks = private$objective_checks(pars,
                                    transform_pars,
                                    check_constraints,
                                    debug)
  pars = checks[[1]]
  pass = checks[[2]]
  
  if (!is.na(pass) & !pass) {
    nll = 1e10
    if(debug) cat("nll =", nll, "\n")
    return(nll)
  }
  
  if (is.null(data_q)) {
    data_q = copy(self$data_q)
  }
  
  
  ### loop through conditions to get chisquare
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, -(1:(length(private$sim_cond)))]
  rt_q_cols = (length(data_q)-length(private$p_q)+2):length(data_q)
  
  chisq = 0
  
  chisq = foreach::foreach(i=1:private$par_transform[, .N], .combine=sum) %do% {
    
    # simulate trials
    par_list = as.list(pars_only_mat[i])
    this_sim = setDT(self$simulate(n_sim,
                                   private$stim_list[[i]],
                                   as.numeric(par_list),
                                   names(par_list),
                                   bounds=private$bounds,
                                   urgency=private$urgency,
                                   ...)$behavior)
    
    # get rt quantile matrix
    sub_q = copy(data_q)
    for(j in 1:length(private$sim_cond)) {
      sub_q = sub_q[get(private$sim_cond[j]) == private$par_matrix[i, get(private$sim_cond[j])]]
    }
    rt_q_mat = as.matrix(sub_q[, rt_q_cols, .(response), with=F])
    n_rt = sub_q[, n_response, .(response)][, n_response]
    sim_rts = list(this_sim[response == 0, rt],
                   this_sim[response == 1, rt],
                   this_sim[is.na(response), rt])
    
    quantile_chisquare(sim_rts, rt_q_mat, private$p_q, n_rt)
    
  }
  
  if(is.na(chisq)) chisq = 1e10
  if(debug) cat(chisq, "\n")
  chisq
  
}


check_pdm_constraints <- function(){
  par_matrix_names = names(private$par_matrix)
  
  checks = sum(private$par_matrix[, (a <= 0) | (a > 10)]) # a
  checks = checks + sum(private$par_matrix[, (t0 < 0) | (t0 > 1)]) # t0
  checks = checks + sum(private$par_matrix[, s <= 0]) # s
  if ("z" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (z <= 0) | (z >= 1)]) # z
  else
    z = 0.5
  if ("sv" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (sv < 0) | (sv > 100)]) # sv
  if ("sz" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (sz < 0) | (sz >= z)]) # sz
  if ("st0" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (st0 < 0) | (st0 >= t0)]) # st0
  if ("lambda" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (lambda < 0) | (lambda > 100)]) # lambda
  if ("aprime" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (aprime < 0) | (aprime > 1)]) # aprime
  if ("kappa" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, kappa < 0]) # kappa
  if ("tc" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, tc <= 0]) # tc
  if ("uslope" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, uslope <= 0]) # uslope
  if ("udelay" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, udelay <= 0]) # udelay
  if ("umag" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, umag <= 0]) # umag
  
  return(checks == 0)
}


#' set diffusion model objective function (for internal use)
#' 
#' Change the objective function used to fit the diffusion model. Options include:
#' \describe{
#' \item{fp}{use the fast-dm method implemented in the rtdists package}
#' \item{chisq}{simulate the first passage time distribution using the Euler-Maruyama method, and compare simulated and observed distrubtions using the chisquare statistic}
#' }
#' 
#' @usage model$set_objective(objective)
#' 
#' @param objective string; the objective function to be used; either "fp" or "chisq"
#' 
#' @return modifies the field \code{obj}
#' 
set_pdm_objective <- function(objective="chisq") {
  
  if (objective == "fp") {
    self$obj = private$fp_obj
  } else if (objective == "chisq") {
    self$obj = private$chisq_obj
  } else {
    stop("specified objective not supported")
  }
  
  invisible(self)
  
}


#' predict pulse model  (for internal use)
#' 
#' Predict behavior with given pulse model parameters.
#' This function is only intended for use with a diffusion model object,
#' and should not be called directly outside of the diffusion model class.
#' 
#' @usage model$predict(pars=NULL, n=10000, ...)
#'
#' @param pars numeric vector; vector of parameters. If NULL, uses model$solution$pars.
#' @param n integer; number of decisions to simulate for each condition. If the number of conditions is equal to the length of the data, e.g. if using as_function with a continuous predictor, ignores \code{n} and simulates one decision per condition
#' @param method string; "euler" for euler-maruyama simulation or "fp" for Fokker-Planck method
#' @param stim_list list of arrays; 2 x timepoints x trials array of stimuli to simulate for each condition
#' @param trial_code list of integer vector; trial indexes to which the new stimuli are associated with. 
#' @param ... additional arguments passed to method (either \code{sim_pulse} or \code{pulse_fp_fpt})
#'
#' @return data.table with simulation conditions, decision (upper or lower boundary) and response time
#' 
#' @keywords internal
#' 
predict_pulse_model = function(pars=NULL, n=1, method="euler", stim_list=NULL, trial_code=NULL, ...){
  
  if (is.null(stim_list) | (length(stim_list) != private$par_transform[, .N])) {
    stop("Must provide a stimulus array for each condition, stored as a list.")
  }
  
  if(is.null(pars)) {
    if (is.null(self$solution)) {
      stop("if parameters have not been fit (i.e. no model solution), must supply pars!")
    }
    pars = self$solution$pars
  }
  private$set_params(pars)
  
  ### loop through conditions
  
  d_pred = data.table()
  
  if (private$par_matrix[, .N] < self$data[, .N]) {
    
    for (i in 1:private$par_matrix[, .N]) {
      
      this_stim_trials = dim(stim_list[[i]])[3]
      
      # check stimulus lengths and trial code
      if (is.null(trial_code)) {
        if ((this_stim_trials != dim(private$stim_list[[i]])[3])) {
          stop(paste("new stimulus array length for condition", i, "is not equal to number of trials"))
        } else {
          this_trial_code = 1:this_stim_trials
        }
      } else {
        if (this_stim_trials != length(trial_code[[i]])) {
          stop(paste("new stimulus array length for condition", i, "is not equal to length of trial code"))
        } else {
          this_trial_code = trial_code[[i]]
        }
      }
      
      # get predicted behavior
      if (method == "euler") {
        
        this_par_list = as.list(private$par_matrix[i, -(1:length(private$sim_cond))])
        this_par_values = as.numeric(this_par_list)
        this_par_names = names(this_par_list)
        
        this_sim = self$simulate(n,
                                 stim_list[[i]],
                                 this_par_values,
                                 par_names=this_par_names,
                                 bounds=private$bounds,
                                 ...)
        d_pred = rbind(d_pred, setDT(this_sim$behavior))
        
      } else {
        
        stop("method not implemented")
        
      }
      
    }
    
  }
  
  d_pred
  
}


#' simulate pulse model  (for internal use)
#' 
#' simulate pulse DDM with given model parameters..
#' This function is only intended for use with a pulse model object,
#' and should not be called directly outside of the pulse model class.
#' Please use \code{sim_pulse} as a standalone function.
#' 
#' @usage model$simulate(n, stimuli, par_values, par_names=NULL, ...)
#'
#' @param n integer; number of decisions to simulate for each condition. If the number of conditions is equal to the length of the data, e.g. if using as_function with a continuous predictor, ignores \code{n} and simulates one decision per condition
#' @param stimuli array; 2 x timepoints x trials array of stimuli (can generate using \code{pulse_stimulus} or \code{pulse_sequence})
#' @param par_values numeric vector; vector of parameters. Must be named vector or used with par_names
#' @param par_names character vector; vector of parameter names
#' @param ... additional arguments passed to \code{sim_pulse}
#'
#' @return data.table with simulation conditions, decision (upper or lower boundary) and response time
#' 
#' @keywords internal
#' 
simulate_pulse_model = function(n, stimuli, par_values, par_names=NULL, ...) {
  
  if ((missing(stimuli)) | (class(stimuli) != "array")) {
    stop("must provide stimulus as a 2 x timepoints x trials array.")
  } 
  
  if (missing(par_values)) {
    stop("No parameters! Must supply parameter vector \"par_values\" as a named vector, 
         or along with a separate vector \"par_names\" with names of paramters.")
  }
  
  if (is.null(names(par_values))) {
    if (is.null(par_names)) {
      stop("No parameter names supplied. \"par_values\" must be a named vector, 
           or the parameter \"par_names\" must be supplied.")
    } else {
      if (length(par_names) != length(par_values)) {
        stop("\"par_names\" must be the same length as \"par_values\".")
      }
      names(par_values) = par_names
    }
  }
  
  do.call(sim_pulse, c(n=n,
                       list(stimuli=stimuli),
                       as.list(par_values),
                       ...))
  
}


init_pulse_model = function(dat,
                            model_name="pdm",
                            stim_var=NULL,
                            stim_sep="",
                            stim_dur=.01,
                            stim_interval=.1,
                            stim_pre=0,
                            dt=.001,
                            include=NULL,
                            depends_on=NULL,
                            as_function=NULL,
                            start_values=NULL,
                            fixed_pars=NULL,
                            extra_condition=NULL,
                            bounds=NULL,
                            urgency=NULL,
                            objective="chisq",
                            max_time=10,
                            ...){
  
  if (is.null(stim_var)) {
    stop("must provide \"stim_var\" argument, referencing the data column that contains the stimulus train.")
  }
  
  
  # set default parameter values
  all_pars = c("a", "t0", "s",
               "z", "dc",
               "sv", "sz", "st0",
               "lambda", "aprime", "kappa", "tc",
               "uslope", "udelay", "umag")
  values = c(1, .3, 1,
             .5, 0,
             0, 0, 0,
             0, 0.5, 1, .25,
             0, 0, 0)
  lower = c(.1, 1e-10, 1e-10,
            .2, -100,
            0, 0, 0,
            0, 0, 0, 1e-10,
            0, 0, 0)
  upper = c(10, 1, 5,
            .8, 100,
            10, .5, .5,
            100, 1, 5, 2,
            10, 10, 10)
  default_pars = c("a", "t0", "s")
  start_if_include = c(sv=0.1,
                       sz=0.1,
                       st0=0.1,
                       lambda=1,
                       uslope=1,
                       umag=1,
                       udelay=1)
  
  super$initialize(dat,
                   model_name,
                   par_names=all_pars,
                   par_values=values,
                   par_lower=lower,
                   par_upper=upper,
                   default_pars=default_pars,
                   start_if_include=start_if_include,
                   include=include,
                   depends_on=depends_on,
                   as_function=as_function,
                   start_values=start_values,
                   fixed_pars=fixed_pars,
                   max_time=max_time,
                   extra_condition=extra_condition,
                   bounds=bounds,
                   urgency=urgency,
                   ...)
  
  self$set_objective(objective)
  
  private$stim_list = list()
  for(i in 1:private$par_transform[, .N]) {
    dsub = copy(self$data)
    for (c in private$sim_cond) {
      dsub = dsub[get(c) == private$par_transform[i, get(c)]]
    }
    
    up_stims = dsub[, get(stim_var[1])]
    if (length(stim_var) > 1) {
      down_stims = dsub[, get(stim_var[2])]
    } else {
      down_stims = NULL
    }
    
    private$stim_list[[i]] = pulse_stimulus(up_stims,
                                            down_stims,
                                            pattern=stim_sep,
                                            dur=stim_dur,
                                            isi=stim_interval,
                                            pre_stim=stim_pre,
                                            dt=dt,
                                            as_array=T)
    
  }
  
}


#' Pulse diffusion model R6 Class
#' 
#' @description
#' 
#' R6 Class that defines a pulse diffusion model (Brunton et al., 2012, Science) to be applied to a set of behavioral data.
#'
#' @details 
#' 
#' For details regarding other methods, see:
#' \itemize{
#' \item{pm$set_objective: \code{\link{set_pdm_objective}}}
#' \item{pm$fit: \code{\link{fit_diffusion_model}}}
#' \item{pm$predict: \code{\link{predict_pulse_model}}}
#' \item{pm$simulate: \code{\link{simulate_pulse_model}}}
#' }
#' 
#' @usage pm <- pulse_model$new(dat, stim_var, stim_sep="", stim_dur=.01, stim_interval=.1, model_name="pdm", include=NULL, depends_on=NULL, as_function=NULL, start_values=NULL, fixed_pars=NULL, extra_condition=NULL, bounds=0L, objective="fp")
#' @usage pm$set_objective(objective="fp")
#' @usage pm$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
#' @usage pm$predict(pars=NULL, n=10000, ...)
#' @usage pm$simulate(n, stimuli, par_values, par_names=NULL, ...)
#'
#' @param dat data table; contains at least 2 columns: rt - response time for trial, response - upper or lower boundary (1 or 0)
#' @param stim_var string; name of stimulus column in dat. If two names, first is upper boundary evidence and second is lower boundary evidence. If one variable name, assumes binary evidence with 1=evidence to upper boundary and  0=evidence to lower boundary
#' @param stim_sep string; pattern separating timepoint by timepoint evidence within evidence string
#' @param stim_dur numeric; stimulus pulse duration (in seconds)
#' @param stim_interval numeric; inter-stimulus interval (in seconds)
#' @param stim_pre numeric; time before stimulus pulse within stimulus bin
#' @param dt numeric; time step to simulate/solve for first passage times
#' @param model_name string; name to identify model, default = "ddm"
#' @param include character; vector of parameters to include in model. drift rate v, boundary a, and non-decision time t0 are included by default always. Can specify ddm parameters starting point z, drift rate variability sv, non decision time variability st0, starting point variability sz, and wiener diffusion noise s. Also can include collapsing bound parameters: degree of collapse aprime (weibull only), slope of collapse kappa, time constant of collapse tc.
#' @param depends_on named character; if a parameter value depends on a task condition, for example drift rate depends on task difficulty, specify here: c("v" = "task difficulty"). Can specify multiple dependent parameters as a vector
#' @param as_function list; list specifying parameters whose value is a function of other parameters
#' @param start_values named numeric; to change default starting parameter values, specify as named vector. E.g. c("v"=2, "a"=3) to set drift rate to 2 and boundary to 3
#' @param fixed_pars named numeric; to fix a parameter at a specified value, use a named vector as with start_values
#' @param extra_condition character; vector of task condition names. Will calculate first passage times for each condition. Recommended only when comparing a model without depends_on with a model that contains a depends_on parameter.
#' @param bounds string: either "fixed" for fixed bounds, or "weibull" or "hyperbolic" for collapsing bounds according to weibull or hyperbolic ratio functions
#' @param objective character: "rtdists" to use the rtdists package (pure and extended ddm only, will not work with collapsing bounds), "integral" to use the integral method from Voskuilen et al., 2016, or "chisquare" to use the difference in chisq from actual vs. simulated response times
#'
#' @return definition of pulse diffusion model object
#'
#' @export
pulse_model = R6::R6Class("pulse_model",
                          inherit=base_diffusion_model,
                          public=list(
                            initialize=init_pulse_model,
                            set_objective=set_pdm_objective,
                            predict=predict_pulse_model,
                            simulate=simulate_pulse_model
                          ), private=list(
                            dt=NULL,
                            as_function=NULL,
                            stim_list=NULL,
                            bounds=NULL,
                            check_par_constraints=check_pdm_constraints,
                            fp_obj = pulse_fp_obj,
                            chisq_obj = pulse_x2_obj
                          ))