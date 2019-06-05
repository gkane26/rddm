#####################################
# Helper functions for pulse model object

set_pm_parameters = function(pars){
  self$par_values = pars
  private$par_matrix = copy(private$par_transform)
  col_names = names(private$par_matrix)
  for(i in 1:length(private$par_matrix)){
    private$par_matrix[, col_names[i] := self$par_values[get(col_names[i])]]
  }
  
  for(i in 1:length(private$fixed)) private$par_matrix[, names(private$fixed[i]) := as.numeric(private$fixed[i])]
}

init_pulse_model = function(dat, model_name="pdm", include=NULL, depends_on=NULL, start_values=NULL, fixed_pars=NULL, extra_condition=NULL, use_weibull_bound=F, ...){
  
  super$initialize(dat, model_name)
  
  # get task conditions and rt quantiles
  sort_var = c(depends_on, extra_condition, "response")
  setorderv(dat, sort_var)
  simulate_conditions = c(unique(c(depends_on, extra_condition)))
  par_transform = data.table(trial=1:dat[,.N])
  if(!is.null(simulate_conditions)){
    par_transform = cbind(par_transform, dat[, get(simulate_conditions)])
    names(par_transform)[2:length(par_transform)] = simulate_conditions
  }
  
  # set default parameter values
  all_pars = c("v", "a", "t0", "z", "sv", "sz", "st0", "lambda", "a_prime", "kappa", "tc", "s")
  values = c(10, 1, .3, .5, 0, 0, 0, 0, 0, 0, .25, 1)
  lower = c(-100, .01, 1e-10, .2, 0, 0, 0, -1, 0, 0, 1e-10, 1e-10)
  upper = c(100, 10, 1, .8, 10, .2, .2, 1, 1, 5, 2, 5)
  
  # check all supplied parameters, remove if not in all_pars
  rm_fixed = fixed_pars[!(names(fixed_pars) %in% all_pars)]
  fixed_pars = fixed_pars[names(fixed_pars) %in% all_pars]
  rm_include = include[!(include %in% all_pars)]
  include = include[include %in% all_pars]
  rm = c(rm_fixed, rm_include)
  if(length(rm) >= 1)
    cat("requested variables are not supported :: ", rm, sep=c("", rep(", ", length(rm)-1)))
  
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
    dat = dat[order(get(depends_on))]
  par_names = character()
  par_values = numeric()
  par_lower = numeric()
  par_upper = numeric()
  par_transform = cbind(par_transform, matrix(0,nrow=dat[,.N],ncol=length(include),dimnames=list(NULL,include)))
  for(i in 1:length(all_pars)){
    if(all_pars[i] %in% include){
      if(all_pars[i] %in% names(depends_on)){
        conds = dat[, unique(get(depends_on[all_pars[i]]))]
        for(j in conds){
          par_names = c(par_names, paste(all_pars[i], j, sep="_"))
          par_values = c(par_values, values[i])
          par_lower = c(par_lower, lower[i])
          par_upper = c(par_upper, upper[i])
          par_transform[get(depends_on[all_pars[i]]) == j, (all_pars[i]) := length(par_values)]
        }
      }else{
        par_names = c(par_names, all_pars[i])
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
  self$obj=function(pars, transform_pars=F, debug_lik=F, ...){
    
    #check params
    if(transform_pars){
      pars = logistic_untransform(pars)
    }else if(sum(pars<private$lower, pars>private$upper) > 0){
      return(1e10)
    }
    
    if(debug_lik){
      cat("pars = ")
      cat(pars, sep=", ")
      cat(" // ")
    }
    
    # get likelihood
    private$set_params(pars)
    nll = do.call(pulse_nll, c(list(choice=self$data[,response]), list(rt=self$data[, rt]), list(blink_seq=self$data[, blinkSequence]), as.list(private$par_matrix), use_weibull_bound=private$use_weibull_bound, ...))
    
    if(is.infinite(nll)) browser()
    if(is.na(nll)) browser()
    if(debug_lik) cat("nll =", nll, "\n")
    
    nll
  }
  
  self$start_values=par_values
  private$lower=par_lower
  private$upper=par_upper
  private$fixed=fixed_pars
  private$par_transform=par_transform
  private$sim_cond=simulate_conditions
  private$use_weibull_bound=use_weibull_bound
}

#' Pulse diffusion model object
#'
#' @usage pulse_model$new(dat, model_name)
#'
#' @param dat data table; contains at least 2 columns: rt - response time for trial, response - upper or lower boundary (1 or 0)
#' @param model_name string; name to identify model, default = "ddm"
#' @param include character; vector of parameters to include in model. drift rate v, boundary a, and non-decision time t0 are included by default always. Can specify ddm parameters starting point z, drift rate variability sv, non decision time variability st0, starting point variability sz, and wiener diffusion noise s. Also can include collapsing bound parameters: degree of collapse a_prime (weibull only), slope of collapse kappa, time constant of collapse tc.
#' @param depends_on named character; if a parameter value depends on a task condition, for example drift rate depends on task difficulty, specify here: c("v" = "task difficulty"). Can specify multiple dependent parameters as a vector
#' @param start_values named numeric; to change default starting parameter values, specify as named vector. E.g. c("v"=2, "a"=3) to set drift rate to 2 and boundary to 3
#' @param fixed_pars named numeric; to fix a parameter at a specified value, use a named vector as with start_values
#' @param extra_condition character; vector of task condition names. Will calculate first passage times for each condition. Recommended only when comparing a model without depends_on with a model that contains a depends_on parameter.
#' @param use_weibull_bound logical: if T, use weibull function for collapsing bounds. Default = F
#'
#' @return definition of base diffusion model object
#'#' 
#' @export
pulse_model = R6::R6Class("pulse_model",
                      inherit=base_diffusion_model,
                      public=list(
                        initialize=init_pulse_model
                      ), private=list(
                        set_params=set_pm_parameters
                      ))