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
    private$par_matrix[[names(private$par_matrix)[[i]]]] = self$par_values[private$par_matrix[[names(private$par_matrix)[i]]]]
  }
  if("v" %in% names(private$par_matrix))
    private$par_matrix[correctSide==0, v := -v]
}

get_rt_quantiles = function(dat, qs=seq(.1, .9, .2), rt_var="rt", conditions = c("correctSide")){
  p_q = diff(c(0,qs,1))
  rt_qs = dat[, c("n_response" = .N, "p_response" = numeric(1), as.list(quantile(get(rt_var), qs))), c("response", conditions)]
  rt_qs[, p_response := n_response/sum(n_response), conditions]
  list(rt_qs, p_q)
}

init_diffusion_model = function(dat, model_name="ddm", include=NULL, depends_on=NULL, start_values=NULL, fixed_pars=NULL, max_time=10, extra_condition=NULL, use_weibull_bound=T, objective=NULL, ...){
  
  super$initialize(dat, model_name)
  
  # get task conditions and rt quantiles
  sort_var = c(depends_on, extra_condition, "correctSide", "response")
  setorderv(self$data, sort_var)
  simulate_conditions = c("correctSide", unique(c(depends_on, extra_condition)))
  self$data = self$data[rt<max_time]
  q_list = get_rt_quantiles(self$data, conditions = simulate_conditions, ...)
  self$data_q = q_list[[1]]
  private$p_q = q_list[[2]]
  par_transform = self$data[, get(simulate_conditions), simulate_conditions]
  par_transform[, V1:=NULL]
  
  # set default parameter values
  all_pars = c("v", "a", "t0", "z", "sv", "st0", "sz", "a_prime", "kappa", "tc")
  values = c(1, 1.5, .3, .5, 0, 0, 0, 0, 0, .25)
  lower = c(-10, 1, .25, .2, 0, 0, 0, 0, 0, 0)
  upper = c(10, 10, 1, .8, 1, .25, .2, 1, 5, 2)
  
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
    self$data = self$data[order(get(depends_on))]
  par_names = character()
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
  
  ddm_integral_nll = function(pars, min_p=1e-10, transform_pars=F, debug_lik=F, ...){
    
    #check params
    if(transform_pars){
      pars = private$logistic_untransform(pars)
    }else if(sum(pars<private$lower, pars>private$upper) > 0){
      return(1e10)
    }
    private$set_params(pars)
    
    ### get list of rt densities
    pars_only_mat = copy(private$par_matrix)
    pars_only_mat = pars_only_mat[, (1:(length(private$sim_cond))) := NULL]
    density_list = list()
    for(i in 1:pars_only_mat[,.N]){
      this_density = get_first_passage_density(c(as.list(pars_only_mat[i]), private$fixed), max_time=private$max_time, use_weibull_bound=private$use_weibull_bound, ...)
      density_list[[i]] = c(list(density=this_density), unlist(private$par_matrix[i, private$sim_cond, with=F]))
    }
    # density_list = foreach::foreach(i=1:pars_only_mat[,.N]) %do% {
    #   this_density = get_first_passage_density(c(as.list(pars_only_mat[i]), private$fixed), max_time=private$max_time, use_weibull_bound=private$use_weibull_bound, ...)
    #   c(list(density=this_density), unlist(private$par_matrix[i, private$sim_cond, with=F]))
    # }
    
    # get likelihood for each rt
    nll_dat = get_rt_liks(copy(self$data), density_list, min_p=min_p)
    
    # get negative log likelihood
    nll = -sum(log(nll_dat[,p_rt]))
    if(is.infinite(nll)) browser()
    if(is.na(nll)) browser()
    if(debug_lik) cat(nll, "\n")
    return(nll)
  }
  ddm_rtdists_nll = function(pars, min_p=1e-10, transform_pars=F, debug_lik=F){
    
    #check params
    if(transform_pars){
      pars = private$logistic_untransform(pars)
    }else if(sum(pars<private$lower, pars>private$upper) > 0){
      return(1e10)
    }
    private$set_params(pars)
    
    this_par_matrix = copy(private$par_matrix)
    this_par_matrix[, z := z*a]
    
    
    # loop through conditions to get likelihood
    nll = 0
    for(i in 1:this_par_matrix[, .N]){
      sub_dat = copy(self$data)
      for(j in 1:length(private$sim_cond)){
        sub_dat = sub_dat[get(private$sim_cond[j]) == this_par_matrix[i, get(private$sim_cond[j])]]
      }
      sub_p = do.call(rtdists::ddiffusion, c(list(rt=sub_dat[,rt]), list(response=ifelse(sub_dat[,response]==0, "lower", "upper")), as.list(this_par_matrix[i,(1+length(private$sim_cond)):length(this_par_matrix)])))
      sub_p[sub_p < min_p] = min_p
      nll = nll - sum(log(sub_p))
    }
    
    if(is.infinite(nll)) browser()
    if(is.na(nll)) browser()
    if(debug_lik) cat(nll, "\n")
    return(nll)
  }
  ddm_sim_x2 = function(pars, n_sim=10000, transform_pars=F, debug_lik=F, ...){
    
    #check params
    if(transform_pars){
      pars = private$logistic_untransform(pars)
    }else if(sum(pars<private$lower, pars>private$upper) > 0){
      return(1e10)
    }
    private$set_params(pars)
    
    pars_only_mat = copy(private$par_matrix)
    pars_only_mat = pars_only_mat[, (1:(length(private$sim_cond))) := NULL]
    
    # loop through conditions to get likelihood
    chisq = 0
    rt_q_idx = ((length(self$data_q)-length(private$p_q)+2):length(self$data_q))
    for(i in 1:pars_only_mat[, .N]){
      sub_q = copy(self$data_q)
      for(j in 1:length(private$sim_cond)){
        sub_q = sub_q[get(private$sim_cond[j]) == private$par_matrix[i, get(private$sim_cond[j])]]
      }
      this_sim = setDT(do.call(sim_ddm, c(n=n_sim, as.list(pars_only_mat[i], private$fixed), max_time=private$max_time, use_weibull_bound=private$use_weibull_bound, ...)))
      
      # chisq for overtime trials
      if(length(sub_q[response==-1, n_response])==0){
        chisq_ot = 0
      }else{
        exp = sub_q[response==-1, p_response] * n
        obs = this_sim[response==-1, .N]
        chisq_ot = (obs-exp)^2/exp
      }
      
      # get cumulative response times from simulated data
      sim_fpt = sim_to_fpt(this_sim)
      cum_fpt = apply(sim_fpt, 2, cumsum)
      total_response = apply(sim_fpt, 2, sum)
      
      # chisq for down boundary decisions
      n_down = sub_q[response==0, n_response]
      if(length(n_down)==0) n_down = 0
      if(n_down < 5){
        chisq_down = 0
      }else{
        exp = sub_q[response==0, p_response] * n_sim * private$p_q
        t_idx = round(as.numeric(sub_q[response==0, rt_q_idx, with=F])/.01)
        sim_hist_down = diff(c(0, cum_fpt[t_idx, 2], total_response[2]))
        chisq_down = sum((sim_hist_down - exp)^2/exp)
      }
      
      # chisq for up boundary decisions
      n_up = sub_q[response==1, n_response]
      if(length(n_up)==0) n_up = 0
      if(n_up < 5){
        chisq_up = 0
      }else{
        exp = sub_q[response==1, p_response] * n_sim* private$p_q
        t_idx = round(as.numeric(sub_q[response==1, rt_q_idx, with=F])/.01)
        sim_hist_up = diff(c(0, cum_fpt[t_idx, 1], total_response[1]))
        chisq_up = sum((sim_hist_up - exp)^2/exp)
      }
      
      chisq = chisq + chisq_ot + chisq_down + chisq_up
      
    }
    
    if(is.na(chisq)) chisq = 1e10
    if(debug_lik) cat(chisq, "\n")
    chisq
  }

  if(("a_prime" %in% par_names) | ("kappa" %in% par_names) | ("tc" %in% par_names)){
    if(is.null(objective)){
      message("Objective not specified, using integral method with collapsing bounds")
      self$obj = ddm_integral_nll
    }else if(objective == "chisq"){
      self$obj = ddm_sim_x2
    }else{
      if(objective != "integral")
        warning("WARNING :: Specified objective function is not compatible with collapsing bounds. Using integral method.")
      self$obj = ddm_integral_nll  
    }
  }else{
    if(is.null(objective)){
      message("Objective not specified, using rtdists method")
      self$obj = ddm_rtdists_nll
    }else if(objective == "chisq"){
      self$obj = ddm_sim_x2
    }else if(objective == "integral"){
      self$obj = ddm_integral_nll
    }else{
      if(objective != "rtdists")
        warning("WARNING :: Specified objective function is not compatible with collapsing bounds. Using rtdists method.")
      self$obj = ddm_rtdists_nll  
    }
  }
  
  self$par_names=par_names
  self$start_values=par_values
  private$max_time=max_time
  private$lower=par_lower
  private$upper=par_upper
  private$fixed=fixed_pars
  private$par_transform=par_transform
  private$sim_cond=simulate_conditions
  private$use_weibull_bound=use_weibull_bound
}


#' Drift diffusion model object
#'
#' @usage diffusion_model$new(dat, model_name)
#'
#' @param dat data table; contains at least 2 columns: rt - response time for trial, response - upper or lower boundary (1 or 0)
#' @param model_name string; name to identify model, default = "ddm"
#' @param include character; vector of parameters to include in model. drift rate v, boundary a, and non-decision time t0 are included by default always. Can specify ddm parameters starting point z, drift rate variability sv, non decision time variability st0, starting point variability sz, and wiener diffusion noise s. Also can include collapsing bound parameters: degree of collapse a_prime (weibull only), slope of collapse kappa, time constant of collapse tc.
#' @param depends_on named character; if a parameter value depends on a task condition, for example drift rate depends on task difficulty, specify here: c("v" = "task difficulty"). Can specify multiple dependent parameters as a vector
#' @param start_values named numeric; to change default starting parameter values, specify as named vector. E.g. c("v"=2, "a"=3) to set drift rate to 2 and boundary to 3
#' @param fixed_pars named numeric; to fix a parameter at a specified value, use a named vector as with start_values
#' @param max_time numeric; max time to simulate a decision. Lower max time keeps computation times lower, but too low will compromise accuracy
#' @param extra_condition character; vector of task condition names. Will calculate first passage times for each condition. Recommended only when comparing a model without depends_on with a model that contains a depends_on parameter.
#' @param use_weibull_bound logical: if T, use weibull function for collapsing bounds. Default = F
#' @param objective character: "rtdists" to use the rtdists package (pure and extended ddm only, will not work with collapsing bounds), "integral" to use the integral method from Voskuilen et al., 2016, or "chisquare" to use the difference in chisq from actual vs. simulated response times
#'
#' @return definition of base diffusion model object
#'#' 
#' @export
diffusion_model = R6::R6Class("diffusion_model",
                          inherit=base_diffusion_model,
                          public=list(
                            data_q=NULL,
                            initialize=init_diffusion_model
                          ), private=list(
                            p_q=NULL,
                            max_time=NULL,
                            set_params=set_dm_parameters
                          ))
