#' fit diffusion model (for internal use)
#' 
#' Fit parameters of the diffusion model.
#' This function is only intended for use with a diffusion model object,
#' and should not be called directly outside of the diffusion model class.
#' 
#' @usage model$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
#'
#' @param use_bounds logical; if True, perform bounded optimization
#' @param transform_pars logical; if True, will use logistic transformation on parameters to ensure they don't violate bounds, even if using an unbounded optimization method
#' @param n_cores integer; number of cores to use for parallel processing (for sections using the foreach package)
#' @param ... additional arguments passed to \code{modelfitr::fit_model}
#'
#' @return \code{base_diffusion_model} with result saved to the field \code{solution}
#' 
#' @keywords internal
#' 
fit_diffusion_model <- function(use_bounds=FALSE, transform_pars=FALSE, n_cores=1, ...) {
  
  if (n_cores <= 0) n_cores = parallel::detectCores()
  if (n_cores == 1) {
    if (foreach::getDoParRegistered()) {
      foreach::registerDoSEQ()
    }
  } else {
    cl = parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      foreach::registerDoSEQ()
      parallel::stopCluster(cl)
    })
  }
  
  start_vals = self$start_values
  names(start_vals) = self$par_names
  
  if (use_bounds) {
    
    start_vals[start_vals < private$lower] = private$lower[start_vals < private$lower]
    start_vals[start_vals > private$upper] = private$upper[start_vals > private$upper]
    
    self$solution = modelfitr::fit_model(self$obj,
                                         start=start_vals,
                                         lower=private$lower,
                                         upper=private$upper,
                                         aic=TRUE,
                                         bic=TRUE,
                                         n_obs=self$data[, .N],
                                         ...)
    
  } else {
    
    if (transform_pars) {
      start_vals = private$inv_logistic(start_vals)
      start_vals[start_vals < -10] = -10
      start_vals[start_vals > 10] = 10
    }
    
    self$solution = modelfitr::fit_model(self$obj,
                                         start=start_vals,
                                         aic=TRUE,
                                         bic=TRUE,
                                         n_obs=self$data[, .N],
                                         transform_pars=transform_pars,
                                         ...)
    
    if (transform_pars) {
      self$solution$pars = private$logistic(self$solution$pars)
    }
    
  }
  
  private$set_params(self$solution$pars)
  
  invisible(self)
  
}


set_params = function(pars){
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
      # private$par_matrix[[names(private$par_matrix)[[i]]]] = self$par_values[private$par_matrix[[names(private$par_matrix)[i]]]]
      
    }
    
  }
  
  if("v" %in% names(private$par_matrix) & !("v" %in% names(private$as_function)))
    private$par_matrix[correctSide==0, v := -v]
  
}


objective_checks = function(pars,
                            transform_pars=F,
                            check_constraints=T,
                            debug=F) {
  
  ### check params
  
  if(transform_pars) {
    pars = private$logistic(pars)
  }
  
  if(debug){
    cat("pars = ")
    cat(round(pars, 3), sep=", ")
    cat(" // ")
  }
  
  private$set_params(pars)
  
  if(check_constraints){
    pass_checks = private$check_par_constraints()
  } else {
    pass_checks = NA
  }
  
  list(pars, pass_checks)

}


init_model = function(dat,
                      model_name,
                      par_names,
                      par_values,
                      par_lower,
                      par_upper,
                      default_pars=NULL,
                      start_if_include=NULL,
                      include=NULL,
                      depends_on=NULL,
                      as_function=NULL,
                      start_values=NULL,
                      fixed_pars=NULL,
                      max_time=10,
                      extra_condition=NULL,
                      bounds=NULL,
                      urgency=NULL,
                      ...){
  
  ### check data structure
  
  if(!("subject" %in% names(dat))){
    dat$subject = NA
  }
  
  if (!("correctSide" %in% names(dat))){
    dat$correctSide = 1
  }
  
  if ((!("response" %in% names(dat))) | (!("rt" %in% names(dat)))) {
    stop("data must contain the columns \"response\" and \"rt\".")
  }
  
  self$data=data.table::copy(data.table::setDT(dat))
  self$name=model_name
  
  ### check for supplies parameters
  
  if (missing(par_names) | missing(par_values) | missing(par_lower) | missing(par_upper)) {
    stop("must supply parameter names, values, lower and upper bounds.")
  }
  
  ### get variables used for as_function parameters
  
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
  
  ### get task conditions and rt quantiles
  
  sort_var = c(depends_on, extra_condition, as_function_vars, "correctSide", "response")
  setorderv(self$data, sort_var)
  simulate_conditions = c("correctSide", unique(c(depends_on, extra_condition, as_function_vars)))
  # self$data = self$data[rt < max_time]
  q_list = get_rt_quantiles(self$data, conditions = simulate_conditions, ...)
  self$data_q = q_list[[1]]
  private$p_q = q_list[[2]]
  private$as_function = lapply(as_function, function(x) ifelse(class(x) == "function", x, x[[1]]))
  par_transform = unique(self$data[, ..simulate_conditions])
  
  
  ### set default parameter values
  
  all_pars = par_names
  values = par_values
  lower = par_lower
  upper = par_upper
  
  check_default_values = names(start_if_include)
  for (p in names(start_if_include)) {
    if ((p %in% include) | (p %in% names(depends_on)) | (p %in% names(as_function))) values[all_pars == p] = start_if_include[p]
  }
  
  ### check all supplied parameters, remove if not in all_pars
  
  rm_fixed = fixed_pars[!(names(fixed_pars) %in% all_pars)]
  fixed_pars = fixed_pars[names(fixed_pars) %in% all_pars]
  rm_include = include[!(include %in% all_pars)]
  include = include[include %in% all_pars]
  rm = c(rm_fixed, rm_include)
  if(length(rm) >= 1)
    warning("requested variables are not supported :: ", rm, sep=c("", rep(", ", length(rm)-1)))
  
  ### overwrite defaults with supplied starting values
  
  if(length(start_values) > 0)
    for(i in 1:length(start_values)){
      values[all_pars==names(start_values)[i]] = start_values[i]
    }
  
  ### get parameters to be fit and parameters with fixed values
  
  default_pars = default_pars[!(default_pars %in% names(fixed_pars))]
  include = include[!(include %in% names(fixed_pars))]
  include = c(default_pars, include)
  
  if(length(depends_on) > 0)
    setorderv(self$data, depends_on)
  par_names = character()
  par_corresponding = character()
  par_values = numeric()
  par_lower = numeric()
  par_upper = numeric()
  par_transform = cbind(par_transform, matrix(0,nrow=par_transform[,.N],ncol=length(include), dimnames=list(NULL,include)))
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
  
  self$par_names=par_names
  self$par_corresponding=par_corresponding
  self$start_values=par_values
  
  private$max_time=max_time
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
    if (!bounds %in% c("fixed", "hyperbolic", "weibull", "linear")) {
      warning("specified bounds not supported, using fixed bounds.")
      private$bounds = 0L
    } else {
      private$bounds = as.numeric(factor(bounds, levels=c("fixed", "hyperbolic", "weibull", "linear"))) - 1
    }
    
  }
  
  if (is.null(urgency)) {
    if("umag" %in% par_corresponding) {
      private$urgency = 2L
    } else if (any(c("uslope", "udelay") %in% par_corresponding)) {
      private$urgency = 1L
    } else {
      private$urgency = 0L
    }
  } else {
    if (!(urgency %in% c("none", "linear", "logistic"))) {
      warning("specified urgency not supported, not using urgency signal.")
      private$urgency = 0L
    } else {
      private$urgency = as.numeric(factor(urgency, levels=c("none", "linear", "logistic"))) - 1
    }
  }
  
}
#' Base Diffusion Model R6 Class
#' 
#' @description 
#' 
#' This class serves as a parent object for other diffusion models.
#' 
#' @details
#' 
#' For details for the fit diffusion model method, see \code{\link{fit_diffusion_model}}
#' 
#' @usage dm <- base_diffusion_model$new(data, model_name="base_dm")
#' @usage dm$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
#' 
#' @param data data.frame or data.table with at least three columns: correctSide, which boundary is the correct answer; response, the response of the subject on that trial; and rt, the response time on the trial
#' @param model_name string; arbitrary name for the model
#' 
#' @field name arbitrary name for model
#' @field data data.frame or data.table with at least three columns: correctSide, which boundary is the correct answer; response, the response of the subject on that trial; and rt, the response time on the trial
#' @field par_names the parameter names for the model
#' @field par_corresponding which underlying diffusion model parameter does this parameter correspond to (important if as_function is used)
#' @field par_values the current parameter values
#' @field start_values starting parameter values (used for optimization)
#' @field obj the objective function used for optimization
#' @field solution information about model fit (returned from the package modelfitr). This attribute is NULL upon initiation, and is created by running the \code{fit} method
#'
#' @return \code{base_diffusion_model} object
#'
#' @import data.table
#' 
#' @keywords internal
#' 
base_diffusion_model <- R6::R6Class("base_diffusion_model",
                                    public=list(
                                      name=NULL,
                                      data=NULL,
                                      data_q=NULL,
                                      par_names=NULL,
                                      par_corresponding=NULL,
                                      par_values=NULL,
                                      start_values=NULL,
                                      obj=NULL,
                                      solution=NULL,
                                      initialize=init_model,
                                      fit=fit_diffusion_model,
                                      predict=function() stop("Not Implemented!"),
                                      simulate=function() stop("Not Implemented")
                                    ),
                                    private=list(
                                      lower=NULL,
                                      upper=NULL,
                                      fixed=NULL,
                                      par_transform=NULL,
                                      par_matrix=NULL,
                                      sim_cond=NULL,
                                      p_q=NULL,
                                      set_params=set_params,
                                      logistic=function(x) logistic(x, private$lower, private$upper),
                                      inv_logistic=function(x) inv_logistic(x, private$lower, private$upper),
                                      max_time=NULL,
                                      as_function=NULL,
                                      bounds=NULL,
                                      urgency=NULL,
                                      objective_checks=objective_checks
                                    ))
