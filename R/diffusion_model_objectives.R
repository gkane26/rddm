#' @noRd
#' @importFrom foreach %dopar%
ddm_rtdists_nll = function(pars, dat=NULL, min_p=1e-10, transform_pars=F, debug_lik=F){
  
  if(is.null(dat)) dat = self$data
  
  #check params
  if(transform_pars){
    pars = private$logistic_untransform(pars)
  }else if(sum(pars<private$lower, pars>private$upper) > 0){
    return(1e10)
  }
  private$set_params(pars)
  
  this_par_matrix = copy(private$par_matrix)
  if ("z" %in% names(this_par_matrix)){
    this_par_matrix[, z := z*a]
  }
  
  if("dc" %in% names(this_par_matrix)){
    this_par_matrix[, v := v + dc]
    this_par_matrix[, dc := NULL]
  }
  
  # loop through conditions to get likelihood
  p_response = c()
  if (this_par_matrix[, .N] < dat[, .N]) {
    nll = 0
    for(i in 1:this_par_matrix[, .N]){
      sub_dat = copy(dat)
      for(j in 1:length(private$sim_cond)){
        sub_dat = sub_dat[get(private$sim_cond[j]) == this_par_matrix[i, get(private$sim_cond[j])]]
      }
      sub_p = do.call(rtdists::ddiffusion, c(list(rt=sub_dat[,rt]), list(response=ifelse(sub_dat[,response]==0, "lower", "upper")), as.list(this_par_matrix[i,(1+length(private$sim_cond)):length(this_par_matrix)])))
      p_response = c(p_response, sub_p)
    }
  } else {
    p_response = do.call(rtdists::ddiffusion, c(list(rt=dat[, rt]), list(response=dat[, ifelse(response == 0, "lower", "upper")]), as.list(this_par_matrix[, (1+length(private$sim_cond)):length(this_par_matrix)])))
  }
  
  p_response[p_response < min_p] = min_p
  nll = -sum(log(p_response))
  
  if(is.infinite(nll)) browser()
  if(is.na(nll)) browser()
  if(debug_lik) cat(nll, "\n")
  return(nll)
}

#' @noRd
#' @importFrom foreach %dopar%
ddm_integral_nll = function(pars, dat=NULL, min_p=1e-10, transform_pars=F, debug_lik=F, ...){
  
  if(is.null(dat)) dat = self$data
  
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
  
  # density_list = list()
  # for(i in 1:pars_only_mat[,.N]){
  #   this_density = get_first_passage_density(c(as.list(pars_only_mat[i]), private$fixed), max_time=private$max_time, use_weibull_bound=private$use_weibull_bound, ...)
  #   density_list[[i]] = c(list(density=this_density), unlist(private$par_matrix[i, private$sim_cond, with=F]))
  # }
  
  if (!foreach::getDoParRegistered()) foreach::registerDoSEQ()

  fixed = private$fixed
  max_time = private$max_time
  use_weibull_bound=private$use_weibull_bound
  par_matrix = private$par_matrix
  sim_cond = private$sim_cond
  density_list = foreach::foreach(i=1:pars_only_mat[,.N]) %dopar% {
    this_density = get_first_passage_density(c(as.list(pars_only_mat[i]), fixed), max_time=max_time, use_weibull_bound=use_weibull_bound, ...)
    c(list(density=this_density), unlist(par_matrix[i, sim_cond, with=F]))
  }
  
  # get likelihood for each rt
  nll_dat = get_rt_liks(copy(dat), density_list, min_p=min_p)
  
  # get negative log likelihood
  nll = -sum(log(nll_dat[,p_rt]))
  if(is.infinite(nll)) browser()
  if(is.na(nll)) browser()
  if(debug_lik) cat(nll, "\n")
  return(nll)
}

#' @noRd
#' @importFrom foreach %dopar%
ddm_sim_x2 = function(pars, dat=NULL, n_sim=10000, transform_pars=F, debug_lik=F, qs=seq(.1, .9, .2), rt_var="rt", conditions=c("correctSide"), ...){
  
  if(is.null(dat))
    dat_q = self$data_q
  else
    dat_q = get_rt_quantiles(dat, qs=qs, rt_var=rt_var, conditions=conditions)
  
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
  rt_q_idx = ((length(dat_q)-length(private$p_q)+2):length(dat_q))
  for(i in 1:pars_only_mat[, .N]){
    sub_q = copy(dat_q)
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
