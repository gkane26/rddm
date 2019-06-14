#' Fit diffusion model (for internal use)
#' 
#' @usage diffusion_model$fit()
#'
#' @param method string; optimization method to use, see optimx method argument, default = NULL -- will use Nelder-Mead if unbounded or nlminb if bounds specified.
#' @param save_dir string; path to directory to save results. will save results only if specified. default = NA
#' @param base_file string; base name of file. if NA, will create name a default name: "fit_<modelName>_<subject>.Rdata"
#' @param check_prior_fit logical; if True, will look for and load previous fit in save_dir/base_file
#' @param refit logical; if True, if prior fit is found, will refit model starting at prior point
#' @param transform_pars logical; if True, will use logistic transformation on parameters to ensure they don't violate bounds even if using an unbounded optimization method
#' @param use_bounds logical; if True, perform bounded optimization
#' @param n_cons integer; will return when n_cons number of consecutive optimization runs without finding a better solution
#' @param tol numeric; if likelihood of current fit is within tol of last fit, considered same result
#'
#' @return definition of base diffusion model object
#' 
#' @keywords internal
#' 
#' @import optimx
#'
#' @export
fit_diffusion_model <- function(method=NULL, save_dir=NA, base_file=NA, check_prior_fit=F, refit=F,
                                transform_pars=F, use_bounds=F, n_cons=0, max_tries=25, tol=1e-5, ...){
  cat("\n")
  
  if(is.null(method))
    if(use_bounds)
      method="nlminb"
    else
      method="Nelder-Mead"
    
    ### check for prior fit of model
    ### if model exists, return if refit=F or if convergence criteria met (convergence code = 0 and kkt conditions met)
    if(check_prior_fit){
      try({
        if(!is.na(save_dir)){
          if(is.na(base_file))
            base_file = paste0("fit", "_", self$name, "_", self$data[1,subject])
          file_name = paste0(base_file, ".Rdata")
          full_path = file.path(save_dir, file_name)
          if(file.exists(full_path)){
            saved_dm = readRDS(file=full_path)
            if(!is.null(saved_dm$fit_info)){
              if(!refit){
                if(!is.null(saved_dm$fit_cons)) cons = saved_dm$fit_cons else cons = 0
                if((!is.na(saved_dm$fit_info$kkt1)) & (!is.na(saved_dm$fit_info$kkt2))){
                  if(((saved_dm$fit_info$convcode == 0) & (saved_dm$fit_info$kkt1) & (saved_dm$fit_info$kkt2)) | (cons >= n_cons)){
                    cat("Model =", self$name, "// Subject =", self$data[1,subject], "already fit, good solution found! Model saved to file =", full_path, "\n\n")
                    return()
                  }
                }else{
                  if(cons >= n_cons){
                    cat("Model =", self$name, "// Subject =", self$data[1,subject], "already fit, good solution found! Model saved to file =", full_path, "\n\n")
                    return()
                  }
                }
              }
              self$start_values = saved_dm$fit_par
              old_val = saved_dm$fit_obj
            }
          }
        }
      })
    }
    
    ### create file name to save model fit (if save_dir is specified)
    save=F
    if(!is.na(save_dir)){
      if(is.na(base_file))
        base_file = paste0("fit", "_", self$name, "_", self$data[1,subject])
      file_name = paste0(base_file, ".Rdata")
      full_path = file.path(save_dir, file_name)
      save=T
    }
    
    ### set up fit model
    
    cat("Fitting Model =", self$name, "// Subject =", self$data[1,subject], "\n\n")
    
    if(use_bounds){
      self$start_values[self$start_values < private$lower] = private$lower[self$start_values < private$lower]
      self$start_values[self$start_values > private$upper] = private$upper[self$start_values > private$upper]
    }
    
    convergence = F
    it=0
    if(!exists("cons")) cons=0
    start = as.numeric(self$start_values)
    if(transform_pars)
      start = logistic_transform(start)
    if(!exists("old_val")) old_val = Inf
    use_method = 1
    
    # fit model, check for converegence
    fit_start = as.numeric(Sys.time())
    while(!convergence & it<max_tries){
      
      # model fit
      if(use_bounds)
        fit = optimx(start, self$obj, method=method[use_method], transform_pars=transform_pars, ...)
      else
        fit = optimx(start, self$obj, method=method[use_method], transform_pars=transform_pars, ...)
      
      if(fit$value < old_val){
        self$fit_info = fit
        self$fit_obj = fit$value
        if(transform_pars)
          self$fit_par = logistic_untransform(as.numeric(fit[1,1:length(self$par_names)]))
        else
          self$fit_par = as.numeric(fit[1,1:length(self$par_names)])
        names(self$fit_par) = self$par_names
        self$fit_aic = 2*length(self$par_names) + 2*self$fit_obj
        self$fit_bic = log(self$data[,.N])*length(self$par_names) + 2*self$fit_obj
      }
      
      if(old_val-fit$value > tol){
        cons = 0
        old_val = fit$value
      }else{
        cons = cons + 1
      }
      self$fit_cons = cons
      
      # save result
      if(save) saveRDS(self, file=full_path)
      
      #check convergence
      if((!is.na(fit$kkt1)) & (!is.na(fit$kkt2)))
        convergence = ((fit$convcode == 0) & (fit$kkt1) & (fit$kkt2)) | (cons >= n_cons)
      else
        convergence = (cons >= n_cons)
      
      if(!convergence){
        start = rnorm(length(self$fit_par), self$fit_par, abs(as.numeric(self$fit_par))*.05)
        start[start>private$upper] = private$upper[start>private$upper]
        start[start<private$lower] = private$lower[start<private$lower]
        if(transform_pars)
          start = logistic_transform(start)
        
        # if(old_val-fit$value > tol){
        #   start = as.numeric(fit[1,1:length(start)])
        #   old_val = fit$value
        # }else{
        #   start = rnorm(length(diffusion_model$fit_par), diffusion_model$fit_par, abs(as.numeric(diffusion_model$fit_par))*.05)
        #   start[start>private$upper] = private$upper[start>private$upper]
        #   start[start<private$lower] = private$lower[start<private$lower]
        #   if(transform_pars)
        #     start = logistic_transform(start, diffusion_model$lower, diffusion_model$upper)
        # }
      }
      
      it = it + 1
      cat("model =", self$name, "/ iteration =", it, "/ cons =", cons, "/ convergence =", convergence, "/ obj =", fit$value, "/ fit time =", as.numeric(Sys.time())-fit_start, "\n")
      use_method = ifelse(use_method+1 <= length(method), use_method+1, 1)
      
    }
    
    #fit_time = as.numeric(Sys.time()) - fit_start
    #cat("model =", self$name, "done! // par =", self$fit_par, "// obj =", self$fit_obj, " // time =", fit_time, "\n")
}

get_aic <- function(pars=NULL, dat=NULL, ...){
  if(is.null(pars)) pars=self$fit_par
  if(is.null(dat)) dat=self$data
  obj = self$obj(pars, dat, ...)
  return(2*length(pars) + 2*obj)
}

get_bic <- function(pars=NULL, dat=NULL, ...){
  if(is.null(pars)) pars=self$fit_par
  if(is.null(dat)) dat=self$data
  obj = self$obj(pars, dat, ...)
  return(log(dat[,.N])*length(pars) + 2*obj)
}


#' Base diffusion model object
#'
#' @usage base_diffusion_model$new(dat, model_name)
#'
#' @param dat data table; contains at least 2 columns: rt - response time for trial, response - upper or lower boundary (1 or 0)
#' @param model_name string (optional); name to identify model
#'
#' @return fills in the fields: fit_info, fit_obj, fit_par, fit_cons, fit_aic, fit_bic in diffusion_model object
#'
#' @import data.table
#' 
#' @export
base_diffusion_model <- R6::R6Class("base_diffusion_model",
                                    public=list(
                                      name=NULL,
                                      data=NULL,
                                      par_names=NULL,
                                      par_values=NULL,
                                      start_values=NULL,
                                      obj=NULL,
                                      fit_info=NULL,
                                      fit_obj=NULL,
                                      fit_par=NULL,
                                      fit_cons=NULL,
                                      fit_aic=NULL,
                                      fit_bic=NULL,
                                      initialize=function(dat, model_name="base_dm"){
                                        if(!("subject" %in% names(dat)))
                                          dat$subject = NA
                                        if (!("correctSide" %in% names(dat)))
                                          dat$correctSide = 1
                                        self$data=data.table::setDT(dat)
                                        self$name=model_name
                                      },
                                      fit=fit_diffusion_model
                                    ), private=list(
                                      lower=NULL,
                                      upper=NULL,
                                      fixed=NULL,
                                      par_transform=NULL,
                                      par_matrix=NULL,
                                      sim_cond=NULL,
                                      use_weibull_bound=NULL,
                                      get_aic=get_aic,
                                      get_bic=get_bic,
                                      logistic_untransform=function(pars) (private$upper-private$lower) / (1 + exp(-pars)) + private$lower,
                                      logistic_transform=function(pars) -log((private$upper-private$lower)/(pars-private$lower) - 1)
                                    ))
