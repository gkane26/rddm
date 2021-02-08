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
      start_vals = private$logistic_transform(start_vals)
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
      self$solution$pars = private$logistic_untransform(self$solution$pars)
    }
    
  }
  
  invisible(self)
  
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
                                      par_names=NULL,
                                      par_corresponding=NULL,
                                      par_values=NULL,
                                      start_values=NULL,
                                      obj=NULL,
                                      solution=NULL,
                                      initialize=function(data, model_name="base_dm"){
                                        if(!("subject" %in% names(data)))
                                          data$subject = NA
                                        if (!("correctSide" %in% names(data)))
                                          data$correctSide = 1
                                        self$data=data.table::copy(data.table::setDT(data))
                                        self$name=model_name
                                      },
                                      
                                      fit=fit_diffusion_model
                                    ),
                                    
                                    private=list(
                                      lower=NULL,
                                      upper=NULL,
                                      fixed=NULL,
                                      par_transform=NULL,
                                      par_matrix=NULL,
                                      sim_cond=NULL,
                                      use_weibull_bound=NULL,
                                      logistic_untransform=function(pars) (private$upper-private$lower) / (1 + exp(-pars)) + private$lower,
                                      logistic_transform=function(pars) -log((private$upper-private$lower)/(pars-private$lower) - 1)
                                    ))
