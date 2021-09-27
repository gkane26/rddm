#' @noRd
#' @importFrom foreach %dopar%
ddm_rtdists_nll = function(pars,
                           dat=NULL,
                           min_p=1e-10,
                           transform_pars=F,
                           check_constraints=F,
                           debug=F){
  
  ### check constraints
  
  checks = private$objective_checks(pars,
                                    transform_pars=transform_pars,
                                    check_constraints=check_constraints,
                                    debug=debug)
  pars = checks[[1]]
  pass = checks[[2]]
  
  if (!is.na(pass) & !pass) {
    nll = 1e10
    if(debug) cat(nll, "\n")
    return(nll)
  }
  
  if(is.null(dat)) {
    dat = self$data
  }
  
  this_par_matrix = copy(private$par_matrix)
  if ("z" %in% names(this_par_matrix)){
    this_par_matrix[, z := z*a]
  }
  
  if("dc" %in% names(this_par_matrix)){
    this_par_matrix[, v := v + dc]
    this_par_matrix[, dc := NULL]
  }
  
  # loop through conditions to get likelihood
  
  legal_pars = c("v", "a", "t0", "z", "sv", "sz", "st0", "s")
  fixed = private$fixed[names(private$fixed) %in% legal_pars]
  
  p_response = c()
  if (this_par_matrix[, .N] < dat[, .N]) {
    nll = 0
    for(i in 1:this_par_matrix[, .N]){
      sub_dat = copy(dat)
      for(j in 1:length(private$sim_cond)){
        sub_dat = sub_dat[get(private$sim_cond[j]) == this_par_matrix[i, get(private$sim_cond[j])]]
      }
      sub_p = do.call(rtdists::ddiffusion, c(list(rt=sub_dat[,rt]), list(response=ifelse(sub_dat[,response]==0, "lower", "upper")), as.list(this_par_matrix[i,(1+length(private$sim_cond)):length(this_par_matrix)]), fixed))
      p_response = c(p_response, sub_p)
    }
  } else {
    p_response = do.call(rtdists::ddiffusion, c(list(rt=dat[, rt]), list(response=dat[, ifelse(response == 0, "lower", "upper")]), as.list(this_par_matrix[, (1+length(private$sim_cond)):length(this_par_matrix)]), fixed))
  }
  
  p_response[p_response < min_p] = min_p
  nll = -sum(log(p_response))
  
  if(is.na(nll)) nll = 1e10
  if(debug) cat(nll, "\n")
  nll
  
}

#' @noRd
#' @importFrom foreach %dopar%
ddm_integral_nll = function(pars,
                            dat=NULL,
                            min_p=1e-10,
                            transform_pars=F,
                            check_constraints=F,
                            debug=F,
                            ...){
  
  
  ### check constraints
  
  checks = private$objective_checks(pars,
                                    transform_pars=transform_pars,
                                    check_constraints=check_constraints,
                                    debug=debug)
  pars = checks[[1]]
  pass = checks[[2]]
  
  if (!is.na(pass) & !pass) {
    nll = 1e10
    if(debug) cat(nll, "\n")
    return(nll)
  }
  
  if(is.null(dat)) dat = self$data
  
  ### get list of rt densities
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, (1:(length(private$sim_cond))) := NULL]
  
  if (!foreach::getDoParRegistered()) foreach::registerDoSEQ()
  
  legal_pars = c("v", "a", "t0", "z", "dc", "sv", "sz", "st0", "aprime", "kappa", "tc", "s")
  fixed = private$fixed[names(private$fixed) %in% legal_pars]
  
  max_time = private$max_time
  bounds=private$bounds
  par_matrix = private$par_matrix
  sim_cond = private$sim_cond
  density_list = foreach::foreach(i=1:pars_only_mat[,.N]) %dopar% {
    this_density = get_first_passage_density(c(as.list(pars_only_mat[i]), fixed), max_time=max_time, bounds=bounds, ...)
    c(list(density=this_density), unlist(par_matrix[i, sim_cond, with=F]))
  }
  
  ### get likelihood for each rt
  
  nll_dat = get_rt_liks(copy(dat), density_list, min_p=min_p)
  
  ### get negative log likelihood
  
  nll_dat[p_rt < min_p, p_rt := min_p]
  nll = -sum(log(nll_dat[, p_rt]))
  
  if(is.na(nll)) nll = 1e10
  if(debug) cat(nll, "\n")
  nll
  
}


#' @noRd
#' @importFrom foreach %dopar%
ddm_sim_x2 = function(pars,
                      data_q=NULL,
                      n_sim=10000,
                      transform_pars=F,
                      check_constraints=F,
                      debug=F,
                      seed=-1L,
                      ...){
  
  if (seed > 0) {
    set.seed(seed)
  }
  
  ### check constraints
  
  checks = private$objective_checks(pars,
                                    transform_pars=transform_pars,
                                    check_constraints=check_constraints,
                                    debug=debug)
  pars = checks[[1]]
  pass = checks[[2]]
  
  if (!is.na(pass) & !pass) {
    chisq = 1e10
    if(debug) cat(chisq, "\n")
    return(chisq)
  }
  
  if (is.null(data_q)) {
    data_q = copy(self$data_q)
  }
  
  ### loop through conditions to get chisquare
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, -(1:(length(private$sim_cond)))]
  rt_q_cols = (length(data_q)-length(private$p_q)+2):length(data_q)
  
  chisq = 0
  
  if (private$par_transform[, .N] < self$data[, .N]) {
    
    for (i in 1:private$par_transform[, .N]) {
      
      # simulate trials
      par_list = as.list(pars_only_mat[i])
      this_sim = setDT(do.call(sim_ddm, c(n=n_sim,
                                          par_list,
                                          private$fixed,
                                          bounds=private$bounds,
                                          urgency=private$urgency,
                                          max_time=private$max_time,
                                          seed=seed,
                                          ...))$behavior)
      
      # get rt quantile matrix
      sub_q = copy(data_q)
      for(j in 1:length(private$sim_cond)) {
        sub_q = sub_q[get(private$sim_cond[j]) == private$par_matrix[i, get(private$sim_cond[j])]]
      }
      rt_q_mat = as.matrix(sub_q[, rt_q_cols, .(response), with=F])
      n_rt = sub_q[, n_response, .(response)][, n_response]
      sim_rts = list(this_sim[response == 0, rt],
                     this_sim[response == 1, rt],
                     this_sim[is.na(response), rt])
      
      chisq = chisq + quantile_chisquare(sim_rts, rt_q_mat, private$p_q, n_rt)
      
    }
    
  } else {
    
    par_list = as.list(pars_only_mat)
    this_sim = data.table()
    for (i in 1:n_sim) {
      sub_sim = setDT(do.call(sim_ddm_vec, c(par_list,
                                             private$fixed,
                                             bounds=private$bounds,
                                             urgency=private$urgency,
                                             max_time=private$max_time,
                                             seed=seed,
                                             ...))$behavior)
      sub_sim[, correctSide := private$par_transform[, correctSide]]
      this_sim = rbind(this_sim, sub_sim)
    }
    
    for (c in unique(data_q[, correctSide])) {
      
      sub_q = data_q[correctSide == c]
      sub_sim = this_sim[correctSide == c]
      rt_q_mat = as.matrix(sub_q[, rt_q_cols, .(response), with=F])
      n_rt = sub_q[, n_response, .(response)][, n_response]
      sim_rts = list(sub_sim[response == 0, rt],
                     sub_sim[response == 1, rt],
                     sub_sim[is.na(response), rt])
      
      chisq = chisq + quantile_chisquare(sim_rts, rt_q_mat, private$p_q, n_rt)
      
    }
    
  }
  
  if(is.na(chisq)) chisq = 1e10
  if(debug) cat(chisq, "\n")
  chisq
  
}

#' @noRd
#' @importFrom foreach %dopar%
ddm_sim_qmpe = function(pars,
                        data_q=NULL,
                        n_sim=10000,
                        min_p=1e-10,
                        transform_pars=F,
                        check_constraints=F,
                        debug=F,
                        seed=-1L,
                        ...){
  
  ### check constraints
  
  if (seed > 0) {
    set.seed(seed)
  }
  
  checks = private$objective_checks(pars,
                                    transform_pars=transform_pars,
                                    check_constraints=check_constraints,
                                    debug=debug)
  pars = checks[[1]]
  pass = checks[[2]]
  
  if (!is.na(pass) & !pass) {
    qmpe_nll = 1e10
    if(debug) cat(qmpe_nll, "\n")
    return(qmpe_nll)
  }
  
  if (is.null(data_q)) {
    data_q = copy(self$data_q)
  }
  
  ### loop through conditions to get qmpe
  
  pars_only_mat = copy(private$par_matrix)
  pars_only_mat = pars_only_mat[, -(1:(length(private$sim_cond)))]
  rt_q_cols = (length(data_q)-length(private$p_q)+2):length(data_q)
  
  qmpe_nll = 0
  
  if (private$par_transform[, .N] < self$data[, .N]) {
    
    for (i in 1:private$par_transform[, .N]) {
      
      # simulate trials
      par_list = as.list(pars_only_mat[i])
      this_sim = setDT(do.call(sim_ddm, c(n=n_sim,
                                          par_list,
                                          private$fixed,
                                          bounds=private$bounds,
                                          urgency=private$urgency,
                                          max_time=private$max_time,
                                          seed=seed,
                                          ...))$behavior)
      
      # get rt quantile matrix
      sub_q = copy(data_q)
      for(j in 1:length(private$sim_cond)) {
        sub_q = sub_q[get(private$sim_cond[j]) == private$par_matrix[i, get(private$sim_cond[j])]]
      }
      rt_q_mat = as.matrix(sub_q[, rt_q_cols, .(response), with=F])
      n_rt = sub_q[, n_response, .(response)][, n_response]
      sim_rts = list(this_sim[response == 0, rt],
                     this_sim[response == 1, rt],
                     this_sim[is.na(response), rt])
      
      qmpe_nll = qmpe_nll + qmpe(sim_rts, rt_q_mat, private$p_q, n_rt, min_p=min_p)
      
    }
    
  } else {
    
    par_list = as.list(pars_only_mat)
    this_sim = data.table()
    for (i in 1:n_sim) {
      sub_sim = setDT(do.call(sim_ddm_vec, c(par_list,
                                             private$fixed,
                                             bounds=private$bounds,
                                             urgency=private$urgency,
                                             max_time=private$max_time,
                                             seed=seed,
                                             ...))$behavior)
      sub_sim[, correctSide := private$par_transform[, correctSide]]
      this_sim = rbind(this_sim, sub_sim)
    }
    
    for (c in unique(data_q[, correctSide])) {
      
      sub_q = data_q[correctSide == c]
      sub_sim = this_sim[correctSide == c]
      rt_q_mat = as.matrix(sub_q[, rt_q_cols, .(response), with=F])
      n_rt = sub_q[, n_response, .(response)][, n_response]
      sim_rts = list(sub_sim[response == 0, rt],
                     sub_sim[response == 1, rt],
                     sub_sim[is.na(response), rt])
      
      qmpe_nll = qmpe_nll + qmpe(sim_rts, rt_q_mat, private$p_q, n_rt, min_p=min_p)
      
    }
    
  }
  
  if(is.na(qmpe_nll)) qmpe_nll = 1e10
  if(debug) cat(qmpe_nll, "\n")
  qmpe_nll
  
}
