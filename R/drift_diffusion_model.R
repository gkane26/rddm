#####################################
# Helper functions for diffusion model object

check_ddm_constraints <- function(){
  par_matrix_names = names(private$par_matrix)
  checks = 0
  if ("a" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, a <= 0]) # a
  if ("t0" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, t0 < 0]) # t0
  if ("z" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (z <= 0) | (z >= 1)]) # z
  if ("sv" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (sv < 0) | (sv >= 10)]) # sz
  if ("sz" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (sz < 0) | (sz >= 1)]) # sz
  if ("st0" %in% par_matrix_names)
    checks = checks + sum(private$par_matrix[, (st0 < 0) | (st0 >= 1)]) # st0
  if ("aprime"  %in% par_matrix_names)
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
#' \item{rtdists}{use the fast-dm method implemented in the rtdists package}
#' \item{integral}{use the integral equation method from Smith, 2000 to obtain the first passage time distribution}
#' \item{chisq}{simulate the first passage time distribution using the Euler-Maruyama method, and compare simulated and observed distrubtions using the chisquare statistic}
#' }
#' 
#' @usage model$set_objective(objective)
#' 
#' @param objective string; the objective function to be used; either "rtdists", "integral", "chisq", or "qmpe"
#' 
#' @return modifies the field \code{obj}
#' 
set_ddm_objective <- function(objective=NULL) {
  
  if (is.null(objective)) {
    if (("uslope" %in% self$par_names) | ("udelay" %in% self$par_names) | ("umag" %in% self$par_names)) {
      self$obj = private$chisq_obj
    } else if(("aprime" %in% self$par_names) | ("kappa" %in% self$par_names) | ("tc" %in% self$par_names)) {
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
  } else if (objective == "qmpe") {
    self$obj = private$qmpe_obj
  } else {
    if (("uslope" %in% self$par_names) | ("udelay" %in% self$par_names) | ("umag" %in% self$par_names)) {
      message("specified objective not supported, using chisq method with urgency signal.")
      self$obj = private$chisq_obj
    } else if(("aprime" %in% self$par_names) | ("kappa" %in% self$par_names) | ("tc" %in% self$par_names)) {
      message("specified objective not supported, using integral method with collapsing bounds.")
      self$obj = private$integral_obj
    } else {
      message("specified objective not supported, using rtdists.")
      self$obj = private$rtdists_obj
    }
  }
  
  invisible(self)
  
}


init_diffusion_model = function(dat,
                                model_name="ddm",
                                include=NULL,
                                depends_on=NULL,
                                as_function=NULL,
                                start_values=NULL,
                                fixed_pars=NULL,
                                extra_condition=NULL,
                                bounds=NULL,
                                urgency=NULL,
                                objective=NULL,
                                max_time=10,
                                verbose=TRUE,
                                ...){
  
  # set default parameter values
  all_pars = c("v", "a", "t0",
               "z", "dc",
               "sv", "st0", "sz",
               "aprime", "kappa", "tc",
               "uslope", "udelay", "umag")
  values = c(1, 1.5, .3,
             .5, 0,
             0, 0, 0,
             0.25, 1, 0.25,
             0, 0, 0)
  lower = c(-100, 0, 0,
            .1, -100,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0)
  upper = c(100, 25, 5,
            .9, 100,
            10, 1, 1,
            1, 10, 2,
            10, 10, 10)
  default_pars = c("v", "a", "t0")
  start_if_include = c(sv=0.1, sz=0.1, st0=0.1, uslope=1, udelay=1, umag=1)
  
  dat$response2 = ifelse(dat$response == 0, "lower", "upper")
  
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
                   verbose=verbose,
                   ...)
  
  self$set_objective(objective)
  
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
#' @param method string; "rtdists" for rtdists package (rdiffusion), "euler" for euler-maruyama simulation, or "integral" for integral equation method
#' @param ... additional arguments passed to \code{sim_ddm}
#'
#' @return data.table with simulation conditions, decision (upper or lower boundary) and response time
#' 
#' @keywords internal
#' 
predict_diffusion_model = function(pars=NULL,
                                   n=NULL,
                                   method="rtdists",
                                   ...){
  
  if(is.null(pars)) pars = self$par_values
  private$set_params(pars)
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, (1:(length(private$sim_cond))) := NULL]
  
  if (method == "rtdists") {
    
    if (private$bounds > 0) {
      warning("can't use rtdists method with collapsing bound. Using integral method instead.")
      method = "integral"
    }
    
    legal_pars = c("v", "a", "t0", "z", "sv", "sz", "st0", "s")
    fixed = private$fixed[names(private$fixed) %in% legal_pars]
    
    if ("z" %in% names(pars_only_mat)) {
      pars_only_mat[, z := z*a]
    }
    if ("sz" %in% names(pars_only_mat)) {
      pars_only_mat[, sz := sz * a]
    }
    if ("st0" %in% names(pars_only_mat)) {
      pars_only_mat[, st0 := st0 * t0]
      pars_only_mat[, t0 := t0 - st0 / 2]
    }
    if("dc" %in% names(pars_only_mat)) {
      pars_only_mat[, v := v + dc]
      pars_only_mat[, dc := NULL]
    }
    
  } else {
    
    fixed = private$fixed
    
  }
  
  all_sim = data.table()
  
  if (pars_only_mat[, .N] < self$data[, .N]) {
    
    for(i in 1:pars_only_mat[, .N]) {
      
      if (is.null(n)) n = 10000
      
      if (method == "rtdists") {
        
        this_sim = setDT(do.call(rtdists::rdiffusion, c(list(n=n),
                                                        as.list(pars_only_mat[i]),
                                                        fixed,
                                                        maxt=private$max_time)))
        this_sim[, response := ifelse(response == "upper", 1, 0)]
        this_sim[, rt := round(rt, 3)]
        setcolorder(this_sim, c("response", "rt"))
        
      } else if (method == "euler") {
        
        this_sim = do.call(sim_ddm, c(list(n=n),
                                      as.list(pars_only_mat[i]),
                                      fixed,
                                      max_time=private$max_time,
                                      bounds=private$bounds,
                                      urgency=private$urgency,
                                      max_time=private$max_time,
                                      ...))$behavior
        
      } else if (method == "integral") {
        
        this_fpt = do.call(ddm_integral_fpt, c(as.list(pars_only_mat[i]),
                                               fixed,
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
    
    if (is.null(n)) n = 1
    
    all_sim = data.table()
    
    for (i in 1:n) {
      
      if (method == "rtdists") {
        
        this_sim = setDT(do.call(rtdists::rdiffusion, c(list(n=pars_only_mat[, .N]),
                                                        as.list(pars_only_mat),
                                                        fixed,
                                                        maxt=private$max_time)))
        this_sim[, response := ifelse(response == "upper", 1, 0)]
        this_sim[, rt := round(rt, 3)]
        setcolorder(this_sim, c("response", "rt"))
        
        
      } else if (method == "euler") {
        
        this_sim = do.call(sim_ddm_vec, c(as.list(pars_only_mat),
                                          fixed,
                                          max_time=private$max_time,
                                          bounds=private$bounds,
                                          urgency=private$urgency,
                                          ...))$behavior
        
      } else {
        
        stop("method not supported!")
        
      }
      
      all_sim = rbind(all_sim, cbind(private$par_matrix[, 1:length(private$sim_cond)], this_sim))
      
    }
    
  }
  
  all_sim
  
}

#' simulate diffusion model  (for internal use)
#' 
#' simulate DDM with given diffusion model parameters..
#' This function is only intended for use with a diffusion model object,
#' and should not be called directly outside of the diffusion model class.
#' Please use \code{sim_ddm} as a standalone function.
#' 
#' @usage model$simulate(n, par_values, par_names=NULL, ...)
#'
#' @param n integer; number of decisions to simulate for each condition. If the number of conditions is equal to the length of the data, e.g. if using as_function with a continuous predictor, ignores \code{n} and simulates one decision per condition
#' @param pars numeric vector; vector of parameters
#' @param ... additional arguments passed to \code{sim_ddm}
#'
#' @return data.table with simulation conditions, decision (upper or lower boundary) and response time
#' 
#' @keywords internal
#' 
simulate_diffusion_model = function(n=10000, pars=NULL, ...) {
  
  if(is.null(pars)) {
    if (is.null(self$par_values)) {
      pars = self$start_values
    } else {
      pars = self$par_values
    }
  }
  
  private$set_params(pars)
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, -(1:(length(private$sim_cond)))]
  
  all_sim = data.table()
  
  if (private$par_transform[, .N] < self$data[, .N]) {
    
    for (i in 1:private$par_transform[, .N]) {
      
      # get conditions
      this_sim = private$par_transform[i, 1:length(private$sim_cond)]
      
      # simulate trials
      par_list = as.list(pars_only_mat[i])
      this_sim = data.table(this_sim,
                            do.call(sim_ddm, c(n=n,
                                               par_list,
                                               private$fixed,
                                               bounds=private$bounds,
                                               urgency=private$urgency,
                                               max_time=private$max_time,
                                               ...))$behavior)
      
      all_sim = rbind(all_sim, this_sim)
      
    }
    
  } else {
    
    par_list = as.list(pars_only_mat)
    for (i in 1:n) {
      
      this_sim = data.table(private$par_transform[, 1:length(private$sim_cond)],
                            do.call(sim_ddm_vec, c(par_list,
                                                   private$fixed,
                                                   bounds=private$bounds,
                                                   urgency=private$urgency,
                                                   max_time=private$max_time,
                                                   ...))$behavior)
      
      all_sim = rbind(all_sim, this_sim)
      
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
#' \item{dm$simulate: \code{\link{simulate_diffusion_model}}}
#' }
#' 
#' @usage dm <- diffusion_model$new(dat, model_name="ddm", include=NULL, depends_on=NULL, as_function=NULL, start_values=NULL, fixed_pars=NULL, max_time=10, extra_condition=NULL, bounds=NULL, objective=NULL, ...)
#' @usage dm$set_objective(objective=NULL)
#' @usage dm$fit(use_bounds=TRUE, transform_pars=FALSE, ...)
#' @usage dm$predict(pars=NULL, n=10000, ...)
#' @usage dm$simulate(n, par_values)
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
#' @param objective character: "rtdists" to use the rtdists package (pure and extended ddm only, will not work with collapsing bounds), "integral" to use the integral method from Voskuilen et al., 2016, "chisquare" to use the difference in chisq from actual vs. simulated response times, or "qmpe" to use the quantile maximum product estimation method (comparing simulated vs. actual RT distribution) from Brown & Heathcote
#' @param verbose logical: If TRUE, print messages related to model initation. Default = TRUE
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
                                initialize=init_diffusion_model,
                                set_objective=set_ddm_objective,
                                predict=predict_diffusion_model,
                                simulate=simulate_diffusion_model
                              ), private=list(
                                check_par_constraints=check_ddm_constraints,
                                rtdists_obj = ddm_rtdists_nll,
                                integral_obj = ddm_integral_nll,
                                chisq_obj = ddm_sim_x2,
                                qmpe_obj = ddm_sim_qmpe
                              ),
                              lock_objects = FALSE)