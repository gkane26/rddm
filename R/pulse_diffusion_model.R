#####################################
# Helper functions for pulse model object

pulse_fp_obj <- function(pars,
                         dat=NULL,
                         transform_pars=F,
                         check_constraints=T,
                         debug=F,
                         ...){
  
  if(is.null(dat)) dat = self$data
  
  #check params
  if(transform_pars)
    pars = private$logistic_untransform(pars)
  
  if(debug){
    cat("pars = ")
    cat(pars, sep=", ")
    cat(" // ")
  }
  
  private$set_params(pars)
  if(check_constraints){
    if(!private$check_par_constraints()){
      if(debug)
        cat("nll =", nll, "\n")
      return(1e10)
    }
  }
  
  # get likelihood
  nll = do.call(pulse_nll, c(list(choice=dat[,response]),
                             list(rt=dat[, rt]),
                             stimuli=list(private$stim_mat_list),
                             as.list(private$par_matrix),
                             bounds=private$bounds,
                             ...))
  
  if(debug) cat("nll =", round(nll, 3), "\n")
  
  nll
  
}

pulse_chisq_obj <- function(...) stop("Not Implemented!")

set_pm_parameters = function(pars){
  
  # self$par_values = pars
  # private$par_matrix = copy(private$par_transform)
  # col_names = names(private$par_matrix)
  # 
  # for(i in 1:length(private$par_matrix)){
  #   private$par_matrix[, col_names[i] := self$par_values[get(col_names[i])]]
  # }
  # 
  # for(i in 1:length(private$fixed)) private$par_matrix[, names(private$fixed[i]) := as.numeric(private$fixed[i])]
  
  self$par_values = pars
  private$par_matrix = copy(private$par_transform)
  
  for(i in (1+length(private$sim_cond)):(length(private$par_transform))){
    
    this_par = names(private$par_matrix)[[i]]
    
    if (this_par %in% names(private$as_function)) {
      
      fun_par_index = which(self$par_corresponding == this_par)
      fun_pars = stringr::str_split_fixed(self$par_names[fun_par_index], pattern="_", n=2)[,2]
      fun_pars_list = list()
      for (fp in 1:length(fun_pars)) fun_pars_list[[fun_pars[fp]]] = self$par_values[fp]
      use_cols = which(names(private$par_matrix) %in% formalArgs(private$as_function[[this_par]]))
      private$par_matrix[[this_par]] = do.call(private$as_function[[this_par]], c(as.list(private$par_matrix[, ..use_cols]), fun_pars_list))
      
    } else {
      
      private$par_matrix[[this_par]] = self$par_values[private$par_matrix[[this_par]]]
      
    }
    
  }
  
}


check_pdm_constraints <- function(){
  par_matrix_names = names(private$par_matrix)
  checks = sum(private$par_matrix[, a <= 0]) # a
  checks = checks + sum(private$par_matrix[, t0 < 0]) # t0
  if ("z" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (z <= 0) | (z >= 1)]) # z
  if ("sz" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (sz < 0) | (sz >= z)]) # sz
  if ("st0" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (st0 < 0) | (st0 >= t0)]) # st0
  if ("lambda" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (lambda < -1) | (lambda > 1)]) # lambda
  if ("aprime" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (aprime < 0) | (aprime > 1)]) # aprime
  if ("kappa" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, kappa < 0]) # kappa
  if ("tc" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, tc <= 0]) # tc
  if ("s" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, s <= 0]) # s
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
set_pdm_objective <- function(objective="fp") {
  
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
#' @param ... additional arguments passed to method (either \code{sim_pulse} or \code{pulse_fp_fpt})
#'
#' @return data.table with simulation conditions, decision (upper or lower boundary) and response time
#' 
#' @keywords internal
#' 
predict_pulse_model = function(pars=NULL, n=10000, method="euler", ...){
  
  browser()
  
  if(is.null(pars)) pars = self$solution$pars
  private$set_params(pars)
  
  if (method == "euler") {
    
    # do.call(pulse_predict, c(n=n, list(stimuli=private$stim_mat_list), as.list(private$par_matrix), bounds=private$bounds, ...))
    d_pred = data.table()
    
    for (i in 1:length(private$stim_mat_list)) {
      
      this_sim = do.call(sim_pulse, c(n=n, list(stimuli=private$stim_mat_list[[i]]), as.list(private$par_matrix[i]), bounds=private$bounds, ...))
      d_pred = rbind(d_pred, setDT(this_sim))
      
    }
    
    d_pred
    
  } else if (method == "fp") {
    
    d_pred = data.table()
    
    for (i in 1:length(private$stim_mat_list)) {
      
      this_fpt = do.call(pulse_fp_fpt, c(list(stimulus=private$stim_mat_list[[i]]), as.list(private$par_matrix[i]), bounds=private$bounds, ...))
      this_sim = setDT(fpt_to_sim(this_fpt, n=n))
      d_pred = rbind(d_pred, this_sim)
      
    }
    
    d_pred
    
    
  } else {
    
    stop("method not implemented")
  }
}


init_pulse_model = function(dat,
                            stim_var,
                            stim_sep="",
                            stim_dur=.01,
                            stim_interval=.1,
                            stim_pre=0,
                            dt=.001,
                            model_name="pdm",
                            include=NULL,
                            depends_on=NULL,
                            as_function=NULL,
                            start_values=NULL,
                            fixed_pars=NULL,
                            extra_condition=NULL,
                            bounds=NULL,
                            objective="fp"){
  
  super$initialize(dat, model_name)
  
  # get variables used for as_function parameters
  as_function_vars = NULL
  if (!is.null(as_function)) {
    for (i in 1:length(as_function)) {
      
      if (class(as_function[[i]]) == "list") {
        all_args = formals(as_function[[i]][[1]])
      } else  if (class(as_function[[i]]) == "function"){
        all_args = formals(as_function[[i]])
      } else {
        stop("as_function arguments must be a function or a list with two elements: [1] a function and [2] vector or parameter names")
      }
      
      no_default = sapply(all_args, function(x) x == "")
      no_default[names(no_default) == "..."] = FALSE
      as_function_vars = c(as_function_vars, names(all_args)[no_default])
      
    }
  }
  
  # get task conditions and rt quantiles
  sort_var = c(depends_on, extra_condition, as_function_vars, "response")
  setorderv(self$data, sort_var)
  simulate_conditions = c(unique(c(depends_on, extra_condition, as_function_vars)))
  private$as_function = lapply(as_function, function(x) ifelse(class(x) == "function", x, x[[1]]))
  par_transform = data.table(trial=1:self$data[,.N])
  if(!is.null(simulate_conditions)){
    par_transform = cbind(par_transform, self$data[, get(simulate_conditions)])
    names(par_transform)[2:length(par_transform)] = simulate_conditions
  }
  
  # set default parameter values
  all_pars = c("v", "a", "t0", "z", "dc", "sv", "sz", "st0", "lambda", "aprime", "kappa", "tc", "s")
  values = c(1, 1, .3, .5, 0, 0, 0, 0, 0, 0, 0, .25, 1)
  lower = c(-10, .1, 1e-10, .2, -10, 0, 0, 0, -1, 0, 0, 1e-10, 1e-10)
  upper = c(10, 10, 1, .8, 10, 10, .2, .2, 1, 1, 5, 2, 5)
  
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
  
  if(length(depends_on) > 0)
    self$data = self$data[order(get(depends_on))]
  par_names = character()
  par_corresponding = character()
  par_values = numeric()
  par_lower = numeric()
  par_upper = numeric()
  par_transform = cbind(par_transform, matrix(0,nrow=self$data[,.N],ncol=length(include),dimnames=list(NULL,include)))
  for(i in 1:length(all_pars)){
    if(all_pars[i] %in% include){
      if(all_pars[i] %in% names(depends_on)){
        conds = self$data[, unique(get(depends_on[all_pars[i]]))]
        for(j in conds){
          par_names = c(par_names, paste(all_pars[i], j, sep="_"))
          par_corresponding = c(par_corresponding, all_pars[i])
          par_values = c(par_values, values[i])
          par_lower = c(par_lower, lower[i])
          par_upper = c(par_upper, upper[i])
          par_transform[get(depends_on[all_pars[i]]) == j, (all_pars[i]) := length(par_values)]
        }
      }else if(all_pars[i] %in% names(as_function)){
        if(class(as_function[[all_pars[i]]]) == "function") {
          fun_args = formals(as_function[[all_pars[i]]])
          has_default = sapply(all_args, function(x) x != "")
          fun_args = fun_args[has_default]
          fun_names = names(fun_args)
          fun_vals = unname(unlist(fun_args))
          fun_lower = rep(-Inf, length(fun_vals))
          fun_upper = rep(Inf, length(fun_vals))
        } else {
          fun_args = formals(as_function[[all_pars[i]]][[1]])
          has_default = sapply(all_args, function(x) x != "")
          fun_args = fun_args[has_default]
          fun_names = names(fun_args)
          fun_vals = unname(unlist(fun_args))
          fun_lower = as_function[[all_pars[i]]][[2]][1,]
          fun_upper = as_function[[all_pars[i]]][[2]][2,]
        }
        
        par_names = c(par_names, paste(all_pars[i], fun_names, sep="_"))
        par_corresponding = c(par_corresponding, rep(all_pars[i], length(par_names)))
        par_values = c(par_values, fun_vals)
        par_lower = c(par_lower, fun_lower)
        par_upper = c(par_upper, fun_upper)
        par_transform[, all_pars[i] := NA]
        
      }else{
        par_names = c(par_names, all_pars[i])
        par_corresponding = c(par_corresponding, all_pars[i])
        par_values = c(par_values, values[i])
        par_lower = c(par_lower, lower[i])
        par_upper = c(par_upper, upper[i])
        par_transform[, (all_pars[i]) := length(par_values)]
      }
    }else{
      if(!(all_pars[i] %in% names(fixed_pars))){
        message(paste("parameter ::", all_pars[i], "is not specified, including as fixed parameter with value =", values[i]))
        fixed_pars = c(fixed_pars, values[i])
        names(fixed_pars)[length(fixed_pars)] = all_pars[i]
      }
    }
  }
  
  rm_trans_par_cols = 1 + ifelse(is.null(simulate_conditions), 0, length(simulate_conditions))
  par_transform = par_transform[, -(1:rm_trans_par_cols)]
  
  self$par_names=par_names
  self$par_corresponding=par_corresponding
  self$start_values=par_values
  self$set_objective(objective)
  
  private$dt=dt
  private$lower=par_lower
  private$upper=par_upper
  private$fixed=fixed_pars
  private$par_transform=par_transform
  private$sim_cond=simulate_conditions
  
  if (is.null(bounds)) {
    if("aprime" %in% par_corresponding) {
      private$bounds = 2L
    } else if (any(c("kappa", "tc") %in% par_corresponding)) {
      private$bounds = 1L
    } else {
      private$bounds = 0L
    }
  } else {
    if (!bounds %in% c("fixed", "hyperbolic", "weibull")) {
      warning("specified bounds not supported, using fixed bounds.")
      private$bounds = 0L
    } else {
      private$bounds = as.numeric(factor(bounds, levels=c("fixed", "hyperbolic", "weibull"))) - 1
    }
    
  }
  
  # get trial by trial stimuli
  up_stims = self$data[[stim_var[1]]]
  if (length(stim_var) > 1) {
    down_stims = self$data[[stim_var[2]]]
  } else {
    down_stims = NULL
  }
  
  private$stim_mat_list = pulse_stimulus(up_stims,
                                         down_stims,
                                         pattern=stim_sep,
                                         dur=stim_dur,
                                         isi=stim_interval,
                                         pre_stim=stim_pre,
                                         dt=dt)
  
}


#' Pulse diffusion model R6 Class
#' 
#' @description
#' 
#' R6 Class that defines a pulse diffusion model (Brunton et al., 2012, Science) to be applied to a set of behavioral data.
#'
#'#' @details 
#' 
#' For details regarding other methods, see:
#' \itemize{
#' \item{pm$set_objective: \code{\link{set_pdm_objective}}}
#' \item{pm$fit: \code{\link{fit_diffusion_model}}}
#' \item{pm$predict: \code{\link{predict_pulse_model}}}
#' }
#' 
#' @usage pdm <- pulse_model$new(dat, stim_var, stim_sep="", stim_dur=.01, stim_interval=.1, model_name="pdm", include=NULL, depends_on=NULL, as_function=NULL, start_values=NULL, fixed_pars=NULL, extra_condition=NULL, bounds=0L, objective="fp")
#' @usage pdm$set_objective(objective="fp")
#' @usage pdm$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
#' @usage dm$predict(pars=NULL, n=10000, ...)
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
                            predict=predict_pulse_model
                          ), private=list(
                            dt=NULL,
                            as_function=NULL,
                            stim_mat_list=NULL,
                            bounds=NULL,
                            set_params=set_pm_parameters,
                            check_par_constraints=check_pdm_constraints,
                            fp_obj = pulse_fp_obj,
                            chisq_obj = pulse_chisq_obj
                          ))