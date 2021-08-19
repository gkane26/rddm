#####################################
# Helper functions for evacc model object

evacc_chisq_obj <- function(pars,
                            dat=NULL,
                            n_sim=1,
                            transform_pars=F,
                            check_constraints=T,
                            debug=F,
                            qs=seq(.1, .9, 0.2),
                            rt_var="rt",
                            conditions=NULL,
                            seed=NULL,
                            ...){
  if (!is.null(seed)) {
    set.seed(seed)
    RcppZiggurat::zsetseed(seed)
  }
  
  if(is.null(dat)) {
    dat_q = self$data_q
  } else {
    dat_q = get_rt_quantiles(dat, qs=qs, rt_var=rt_var, conditions=conditions)
  }
  
  #check params
  if(transform_pars) {
    pars = private$logistic_untransform(pars)
  }
  
  if(debug){
    cat("pars = ")
    cat(round(pars, 3), sep=", ")
    cat(" // ")
  }
  
  private$set_params(pars)
  if(check_constraints){
    if(!private$check_par_constraints()){
      chisq = 1e10
      if(debug) cat("chisq =", chisq, "\n")
      return(chisq)
    }
  }
  
  # simulate trials with given parameters
  names(pars) = self$par_names
  sim_data = do.call(sim_evacc, c(n=n_sim,
                                  list(stimuli=private$stim_array),
                                  as.list(pars),
                                  ...))$behavior
  n_total = n_sim * nrow(sim_data)
  
  # get cumulative response times from simulated data
  sim_fpt = sim_to_fpt(sim_data) * n_total
  cum_fpt = apply(sim_fpt, 2, cumsum)
  cum_fpt = cum_fpt[, c(2, 1)]
  total_response = cum_fpt[nrow(cum_fpt),]
  rt_q_idx = ((length(dat_q)-length(private$p_q)+2):length(dat_q))
  
  # chisq for left decisions
  exp_left = dat_q[response == 0, p_response] * private$p_q * n_total
  t_idx_left = round(as.numeric(dat_q[response==0, rt_q_idx, with=F])/.01)
  sim_hist_left = diff(c(0, cum_fpt[t_idx_left, 1], total_response[1]))
  chisq_left = sum((sim_hist_left - exp_left)^2/exp_left)
  
  # chisq for right decisions
  exp_right = dat_q[response == 1, p_response] * private$p_q * n_total
  t_idx_right = round(as.numeric(dat_q[response==1, rt_q_idx, with=F])/.01)
  sim_hist_right= diff(c(0, cum_fpt[t_idx_right, 2], total_response[2]))
  chisq_right = sum((sim_hist_right- exp_right)^2/exp_right)
  
  chisq = chisq_left + chisq_right
  
  if(is.na(chisq)) 1e10
  if(debug) cat("chisq = ", round(chisq, 3), "\n")
  chisq
  
}

#' set evacc model objective function (for internal use)
#' 
#' Change the objective function used to fit the evacc model. Options include:
#' \describe{
#' \item{chisq}{simulate the first passage time distribution using the Euler-Maruyama method, and compare simulated and observed distrubtions using the chisquare statistic}
#' }
#' 
#' @usage model$set_objective(objective)
#' 
#' @param objective string; the objective function to be used; only "chisq" implemented
#' 
#' @return modifies the field \code{obj}
#' 
set_evacc_objective <- function(objective="chisq") {
  
  if (objective != "chisq") {
    stop("specified objective not supported")
  } else {
    self$obj = private$chisq_obj
  }
  
  invisible(self)
  
}


set_evacc_parameters = function(pars){
  
  self$par_values = pars
  # private$par_matrix = copy(private$par_transform)
  # 
  # for(i in (1+length(private$sim_cond)):(length(private$par_transform))){
  #   
  #   this_par = names(private$par_matrix)[[i]]
  #   
  #   if (this_par %in% names(private$as_function)) {
  #     
  #     fun_par_index = which(self$par_corresponding == this_par)
  #     fun_pars = stringr::str_split_fixed(self$par_names[fun_par_index], pattern="_", n=2)[,2]
  #     fun_pars_list = list()
  #     for (fp in 1:length(fun_pars)) fun_pars_list[[fun_pars[fp]]] = self$par_values[fp]
  #     use_cols = which(names(private$par_matrix) %in% formalArgs(private$as_function[[this_par]]))
  #     private$par_matrix[[this_par]] = do.call(private$as_function[[this_par]], c(as.list(private$par_matrix[, ..use_cols]), fun_pars_list))
  #     
  #   } else {
  #     
  #     private$par_matrix[[this_par]] = self$par_values[private$par_matrix[[this_par]]]
  #     
  #   }
  #   
  # }
  
}


check_evacc_constraints <- function(){
  lower_check = any(self$par_values < private$lower)
  upper_check = any(self$par_values > private$upper)
  !lower_check & !upper_check
  
  # par_matrix_names = names(private$par_matrix)
  # checks = sum(private$par_matrix[, (v < 0) | (v > 100)]) # v
  # checks = checks + sum(private$par_matrix[, (a <= 0) | (a > 25)]) # a
  # checks = checks + sum(private$par_matrix[, (t0 < 0) | (t0 > 2)]) # t0
  # if ("z" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (z <= -a) | (z >= a)]) # z
  # else
  #   z = 0
  # if ("dc" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (dc < 100) | (dc > 100)])
  # if ("sv" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (sv < 0) | (sv > 100)]) # sv
  # if ("sz" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (sz < 0) | (sz >= abs(z))]) # sz
  # if ("st0" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (st0 < 0) | (st0 >= t0)]) # st0
  # if ("s" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, s <= 0]) # s
  # if ("lambda" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (lambda < 0) | (lambda > 100)]) # lambda
  # if ("mi" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (mi < 0) | (mi > 10)]) # mi
  # if ("mi2" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (mi2 < 0) | (mi2 > 10)]) # mi2
  # if ("sv2" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (sv2 < 0) | (sv2 > 100)]) # sv2
  # if ("aprime" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, (aprime < 0) | (aprime > 1)]) # aprime
  # if ("kappa" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, kappa < 0]) # kappa
  # if ("tc" %in% par_matrix_names)
  #   checks = checks + sum(private$par_matrix[, tc <= 0]) # tc
  # 
  # return(checks == 0)
}


init_evacc_model = function(dat,
                            stim_var,
                            stim_sep="",
                            stim_dur=.01,
                            stim_interval=.1,
                            stim_pre=0,
                            dt=.001,
                            model_name="evacc",
                            include=NULL,
                            start_values=NULL,
                            fixed_pars=NULL,
                            bounds=NULL,
                            ...){
  
  super$initialize(dat, model_name)
  
  # get task conditions and rt quantiles
  q_list = get_rt_quantiles(self$data, qs=seq(.1, .9, .1), conditions = NULL, ...)
  self$data_q = q_list[[1]]
  private$p_q = q_list[[2]]
  
  # set default parameter values
  all_pars = c("v", "a", "t0", "z", "dc", "sv", "sz", "st0", "s", "lambda", "mi", "mi2", "sv2", "aprime", "kappa", "tc")
  values = c(1, 1, .3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, .25, 1)
  lower = c(0, .1, 1e-10, -5, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-10, 1e-10)
  upper = c(100, 25, 2, 5, 100, 100, 1, 1, 10, 100, 10, 10, 100, 5, 2, 5)
  
  check_default_values = c("sv", "sz", "st0", "s", "sv2")
  for (c in check_default_values) {
    if (c %in% include) values[all_pars == c] = 0.1
  }
  
  # check all supplied parameters, remove if not in all_pars
  rm_fixed = fixed_pars[!(names(fixed_pars) %in% all_pars)]
  fixed_pars = fixed_pars[names(fixed_pars) %in% all_pars]
  rm_include = include[!(include %in% all_pars)]
  include = include[include %in% all_pars]
  rm = c(rm_fixed, rm_include)
  if(length(rm) >= 1)
    warning("requested variables are not supported :: ", rm, sep=c("", rep(", ", length(rm)-1)))
  
  # overwrite defaults with supplied starting values
  if(length(start_values) > 0)
    for(i in 1:length(start_values)){
      values[all_pars==names(start_values)[i]] = start_values[i]
    }
  
  # get parameters to be fit and parameters with fixed values
  default_pars = c("v", "a","t0")
  default_pars = default_pars[!(default_pars %in% names(fixed_pars))]
  include = include[!(include %in% names(fixed_pars))]
  include = c(default_pars, include)
  
  par_names = character()
  par_values = numeric()
  par_lower = numeric()
  par_upper = numeric()
  for(i in 1:length(all_pars)){
    if(all_pars[i] %in% include){
      par_names = c(par_names, all_pars[i])
      par_values = c(par_values, values[i])
      par_lower = c(par_lower, lower[i])
      par_upper = c(par_upper, upper[i])
    }else{
      if(!(all_pars[i] %in% names(fixed_pars))){
        message(paste("parameter ::", all_pars[i], "is not specified, including as fixed parameter with value =", values[i]))
        fixed_pars = c(fixed_pars, values[i])
        names(fixed_pars)[length(fixed_pars)] = all_pars[i]
      }
    }
  }
  
  self$par_names=par_names
  self$start_values=par_values
  self$set_objective("chisq")
  
  private$dt=dt
  private$lower=par_lower
  private$upper=par_upper
  private$fixed=fixed_pars
  private$bounds = 0L
  
  # get trial by trial stimuli
  left_stims = self$data[[stim_var[1]]]
  if (length(stim_var) > 1) {
    right_stims = self$data[[stim_var[2]]]
  } else {
    right_stims = NULL
  }
  
  private$stim_array = pulse_stimulus(left_stims,
                                      right_stims,
                                      pattern=stim_sep,
                                      dur=stim_dur,
                                      isi=stim_interval,
                                      pre_stim=stim_pre,
                                      dt=dt,
                                      as_array=TRUE)
  
}


#' EvAcc model R6 Class
#' 
#' @description
#' 
#' R6 Class that defines an EvAcc model to be applied to a set of behavioral data.
#'
#' @details 
#' 
#' For details regarding other methods, see:
#' \itemize{
#' \item{evacc$set_objective: \code{\link{set_evacc_objective}}}
#' \item{evacc$fit: \code{\link{fit_diffusion_model}}}
#' }
#' 
#' @usage evacc <- evacc_model$new(dat, stim_var, stim_sep="", stim_dur=.01, stim_interval=.1, model_name="evacc", include=NULL, start_values=NULL, fixed_pars=NULL, bounds=NULL)
#' @usage evacc$set_objective(objective="chisq")
#' @usage evacc$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
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
#' @param start_values named numeric; to change default starting parameter values, specify as named vector. E.g. c("v"=2, "a"=3) to set drift rate to 2 and boundary to 3
#' @param fixed_pars named numeric; to fix a parameter at a specified value, use a named vector as with start_values
#' @param bounds string: NOT IMPLEMENTED 
#'
#' @return definition of evacc model object
#'
#' @export
#' 
evacc_model = R6::R6Class("evacc_model",
                          inherit=base_diffusion_model,
                          public=list(
                            data_q=NULL,
                            initialize=init_evacc_model,
                            set_objective=set_evacc_objective
                          ), private=list(
                            p_q=NULL,
                            dt=NULL,
                            stim_array=NULL,
                            bounds=NULL,
                            set_params=set_evacc_parameters,
                            check_par_constraints=check_evacc_constraints,
                            chisq_obj = evacc_chisq_obj
                          ))