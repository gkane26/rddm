#####################################
# Helper functions for diffusion model object

get_first_passage_density = function(pars_list, dt=.01, ...){
  density = do.call(ddm_integral_fpt, c(pars_list, dt=dt, ...))
  rownames(density) = seq(dt, dt*nrow(density), dt)
  density
}


get_rt_liks = function(dat, density_list, min_p=1e-10){
  
  min_p = as.numeric(min_p)
  
  #get time bin for each response
  tvec = as.numeric(rownames(density_list[[1]]$density))
  dt = tvec[2] - tvec[1]
  dat[, rt_bin := rt / dt]
  
  condition_idx = 2:length(density_list[[1]])
  p_response = numeric()
  new_dat = data.table()
  for(i in 1:length(density_list)){
    
    #subset only this condition
    sub_dat=dat
    for(j in 2:length(density_list[[i]])){
      sub_dat=sub_dat[get(names(density_list[[i]])[j])==density_list[[i]][j]]
    }
    
    sub_dat[response == 0, p_rt := density_list[[i]]$density[rt_bin, 2]]
    sub_dat[response == 1, p_rt := density_list[[i]]$density[rt_bin, 1]]
    new_dat = rbind(new_dat, sub_dat)
  }
  
  new_dat[p_rt < min_p, p_rt := min_p]
  return(new_dat)
}


set_dm_parameters = function(pars){
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


get_rt_quantiles = function(dat, qs=seq(.1, .9, .2), rt_var="rt", conditions = c("correctSide")){
  p_q = diff(c(0,qs,1))
  rt_qs = dat[, c("n_response" = .N, "p_response" = numeric(1), as.list(quantile(get(rt_var), qs))), c("response", conditions)]
  rt_qs[, p_response := n_response/sum(n_response), conditions]
  list(rt_qs, p_q)
}


check_ddm_constraints <- function(){
  par_matrix_names = names(private$par_matrix)
  checks = sum(private$par_matrix[, a <= 0]) # a
  checks = checks + sum(private$par_matrix[, t0 < 0]) # t0
  if ("z" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (z <= 0) & (z >= 1)]) # z
  if ("sz" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (sz < 0) & (sz >= z)]) # sz
  if ("st0" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (st0 < 0) & (st0 >= t0)]) # st0
  if ("aprime"  %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (aprime < 0) & (aprime > 1)]) # aprime
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
#' \item{rtdists}{use the fast-dm method implemented in the rtdists package}
#' \item{integral}{use the integral equation method from Smith, 2000 to obtain the first passage time distribution}
#' \item{chisq}{simulate the first passage time distribution using the Euler-Maruyama method, and compare simulated and observed distrubtions using the chisquare statistic}
#' }
#' 
#' @usage model$set_objective(objective)
#' 
#' @param objective string; the objective function to be used; either "rtdists", "integral", or "chisq"
#' 
#' @return modifies the field \code{obj}
#' 
set_ddm_objective <- function(objective=NULL) {
  
  if (is.null(objective)) {
    if(("aprime" %in% self$par_names) | ("kappa" %in% self$par_names) | ("tc" %in% self$par_names)) {
      self$obj = private$integral_obj
    } else {
      self$obj = private$rtdists_obj
    }
  } else if (objective == "rtdists") {
    self$obj = private$rtdists_obj
  } else if (objective == "integral") {
    self$obj = private$integral_obj
  } else if (objective == "chisq") {
    self$obj = private$chisq_obj
  } else {
    if(("aprime" %in% self$par_names) | ("kappa" %in% self$par_names) | ("tc" %in% self$par_names)) {
      message("specified objective not supported, using integral method with collapsing bounds.")
      self$obj = private$integral_obj
    } else {
      message("specified objective not supported, using rtdists.")
      self$obj = private$rtdists_obj
    }
  }
  
  invisible(self)
  
}


init_diffusion_model = function(dat, model_name="ddm", include=NULL, depends_on=NULL, as_function=NULL, start_values=NULL, fixed_pars=NULL, max_time=10, extra_condition=NULL, bounds=NULL, objective=NULL, ...){
  
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
  sort_var = c(depends_on, extra_condition, as_function_vars, "correctSide", "response")
  setorderv(self$data, sort_var)
  simulate_conditions = c("correctSide", unique(c(depends_on, extra_condition, as_function_vars)))
  self$data = self$data[rt<max_time]
  q_list = get_rt_quantiles(self$data, conditions = simulate_conditions, ...)
  self$data_q = q_list[[1]]
  private$p_q = q_list[[2]]
  private$as_function = lapply(as_function, function(x) ifelse(class(x) == "function", x, x[[1]]))
  par_transform = self$data[, get(simulate_conditions), simulate_conditions]
  par_transform[, V1:=NULL]
  
  
  # set default parameter values
  all_pars = c("v", "a", "t0", "z", "dc", "sv", "st0", "sz", "aprime", "kappa", "tc")
  values = c(1, 1.5, .3, .5, 0, 0, 0, 0, 0.25, 1, 0.25)
  lower = c(-10, 0, 0, .2, -10, 0, 0, 0, 0, 0, 0)
  upper = c(10, 10, 1, .8, 10, 1, .25, .2, 1, 10, 2)
  
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
  self$set_objective(objective)
  
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
    if (!bounds %in% c("fixed", "hyperbolic", "weibull")) {
      warning("specified bounds not supported, using fixed bounds.")
      private$bounds = 0L
    } else {
      private$bounds = as.numeric(factor(bounds, levels=c("fixed", "hyperbolic", "weibull"))) - 1
    }
    
  }
  
}


#' predict diffusion model  (for internal use)
#' 
#' Predict behavior with given diffusion model parameters..
#' This function is only intended for use with a diffusion model object,
#' and should not be called directly outside of the diffusion model class.
#' 
#' @usage model$predict(pars=NULL, n=10000, ...)
#'
#' @param pars numeric vector; vector of parameters. If NULL, uses model$solution$pars.
#' @param n integer; number of decisions to simulate for each condition. If the number of conditions is equal to the length of the data, e.g. if using as_function with a continuous predictor, ignores \code{n} and simulates one decision per condition
#' @param method string; "euler" for euler-maruyama simulation or "integral" for integral equation method
#' @param ... additional arguments passed to \code{sim_ddm}
#'
#' @return data.table with simulation conditions, decision (upper or lower boundary) and response time
#' 
#' @keywords internal
#' 
predict_diffusion_model = function(pars=NULL, n=10000, method="euler", ...){
  
  if(is.null(pars)) pars = self$solution$pars
  private$set_params(pars)
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, (1:(length(private$sim_cond))) := NULL]
  
  all_sim = data.table()
  
  if (pars_only_mat[, .N] < self$data[, .N]) {
    
    for(i in 1:pars_only_mat[, .N]) {
      
      if (method == "euler") {
        this_sim = do.call(sim_ddm, c(list(n=n),
                                      as.list(pars_only_mat[i]),
                                      as.list(private$fixed),
                                      max_time=private$max_time,
                                      bounds=private$bounds,
                                      ...))
      } else if (method == "integral") {
        
        this_fpt = do.call(ddm_integral_fpt, c(as.list(pars_only_mat[i]),
                                               as.list(private$fixed),
                                               max_time=private$max_time,
                                               bounds=private$bounds,
                                               ...))
        this_sim = fpt_to_sim(this_fpt, n=n)
        
      } else {
        stop("method not supported!")
      }
      
      all_sim = rbind(all_sim, data.table(private$par_matrix[i, 1:length(private$sim_cond)], this_sim))
      
    } 
    
  } else {
    
    if (n == 10000) n = 1
    
    all_sim = data.table()
    
    for (i in 1:n) {
      
      if (method == "euler") {
        this_sim = do.call(sim_ddm_vec, c(as.list(pars_only_mat),
                                          as.list(private$fixed),
                                          max_time=private$max_time,
                                          bounds=private$bounds,
                                          ...))
      } else {
        stop("method not supported!")
      }
      
      all_sim = rbind(all_sim, cbind(private$par_matrix[, 1:length(private$sim_cond)], this_sim))
      
    }
    
  }
  
  all_sim
  
}


#' Drift Diffusion Model R6 Class
#' 
#' @description
#' 
#' R6 Class that defines a diffusion model to be applied to a set of behavioral data.
#' This class implements all of the following models:
#' \itemize{
#' \item{pure DDM}
#' \item{extended DDM}
#' \item{collapsing bounds using hyperbolic ratio (Voskuilen et al., 2016) or weibull (Hawkins et al., 2015) function}
#' }
#' 
#' @details 
#' 
#' For details regarding other methods, see:
#' \itemize{
#' \item{dm$set_objective: \code{\link{set_ddm_objective}}}
#' \item{dm$fit: \code{\link{fit_diffusion_model}}}
#' \item{dm$predict: \code{\link{predict_diffusion_model}}}
#' }
#' 
#' @usage dm <- diffusion_model$new(dat, model_name="ddm", include=NULL, depends_on=NULL, as_function=NULL, start_values=NULL, fixed_pars=NULL, max_time=10, extra_condition=NULL, bounds=NULL, objective=NULL, ...)
#' @usage dm$set_objective(objective=NULL)
#' @usage dm$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
#' @usage dm$predict(pars=NULL, n=10000, ...)
#' 
#' @param dat data table; contains at least 2 columns: rt - response time for trial, response - upper or lower boundary (1 or 0)
#' @param model_name string; name to identify model, default = "ddm"
#' @param include character; vector of parameters to include in model. drift rate v, boundary a, and non-decision time t0 are included by default always. Can specify ddm parameters starting point z, drift rate variability sv, non decision time variability st0, starting point variability sz, and wiener diffusion noise s. Also can include collapsing bound parameters: degree of collapse aprime (weibull only), slope of collapse kappa, time constant of collapse tc.
#' @param depends_on named character; if a parameter value depends on a task condition, for example drift rate depends on task difficulty, specify here: c("v" = "task difficulty"). Can specify multiple dependent parameters as a vector
#' @param as_function list; list specifying parameters whose value is a function of other parameters
#' @param start_values named numeric; to change default starting parameter values, specify as named vector. E.g. c("v"=2, "a"=3) to set drift rate to 2 and boundary to 3
#' @param fixed_pars named numeric; to fix a parameter at a specified value, use a named vector as with start_values
#' @param max_time numeric; max time to simulate a decision. Lower max time keeps computation times lower, but too low will compromise accuracy
#' @param extra_condition character; vector of task condition names. Will calculate first passage times for each condition. Recommended only when comparing a model without depends_on with a model that contains a depends_on parameter.
#' @param bounds string: either "fixed" for fixed bounds, or "weibull" or "hyperbolic" for collapsing bounds according to weibull or hyperbolic ratio functions
#' @param objective character: "rtdists" to use the rtdists package (pure and extended ddm only, will not work with collapsing bounds), "integral" to use the integral method from Voskuilen et al., 2016, or "chisquare" to use the difference in chisq from actual vs. simulated response times
#' @param ... additional arguments passed to the objective function
#'
#' @field name arbitrary name for model
#' @field data data.frame or data.table with at least three columns: correctSide, which boundary is the correct answer; response, the response of the subject on that trial; and rt, the response time on the trial
#' @field data_q response times quantile (only used if objective = "chisq")
#' @field par_names the parameter names for the model
#' @field par_corresponding which underlying diffusion model parameter does this parameter correspond to (important if as_function is used)
#' @field par_values the current parameter values
#' @field start_values starting parameter values (used for optimization)
#' @field obj the objective function used for optimization
#' @field solution information about model fit (returned from the package modelfitr). This attribute is NULL upon initiation, and is created by running the \code{fit} method
#'
#' @return \code{base_diffusion_model} object
#'
#' @export
#' 
diffusion_model = R6::R6Class("diffusion_model",
                              inherit=base_diffusion_model,
                              public=list(
                                data_q=NULL,
                                initialize=init_diffusion_model,
                                set_objective=set_ddm_objective,
                                predict=predict_diffusion_model
                              ), private=list(
                                p_q=NULL,
                                max_time=NULL,
                                as_function=NULL,
                                bounds=NULL,
                                set_params=set_dm_parameters,
                                check_par_constraints=check_ddm_constraints,
                                rtdists_obj = ddm_rtdists_nll,
                                integral_obj = ddm_integral_nll,
                                chisq_obj = ddm_sim_x2
                              ),
                              lock_objects = FALSE)