#' chisquare objective function
#' 
#' @param rts list of simulated response times for different conditions
#' @param rt_qs matrix of response time quantiles for each condition (condition x rt quantile value)
#' @param p_q vector of probabilities corresponding to response time quantiles
#' @param n_rt vector with the number of empirical response times for each condition
#' 
#' @return numeric scalar, chisquare value
#' 
#' @export
quantile_chisquare = function(rts, rt_qs, p_q, n_rt){
  
  n_sim = sum(sapply(rts, length))
  n_emp = sum(n_rt)

  chisq = 0
  for(i in 1:nrow(rt_qs)) {
    
    if (any(is.na(rt_qs[i,]))) {
      expected = (n_rt[i] / n_emp) * n_sim
      observed = length(rts[[i]])
    } else {
      expected = p_q * (n_rt[i] / n_emp) * n_sim
      observed = suppressWarnings(hist(rts[[i]], c(0, rt_qs[i,], Inf), plot=F)$counts)
    }
    
    expected_denominator = sapply(expected, function(x) max(1, x))
    chisq = chisq + sum((observed - expected)^2 / expected_denominator)
    
  }
  
  if(is.na(chisq)) chisq = 1e10
  chisq
  
}

#' qmpe objective function
#' 
#' @param rts list of simulated response times for different conditions
#' @param rt_qs matrix of response time quantiles for each condition (condition x rt quantile value)
#' @param p_q vector of probabilities corresponding to response time quantiles
#' @param n_rt vector with the number of empirical response times for each condition
#' @param min_p numeric; minimum probability in bin
#' 
#' @return numeric scalar, qmpe value
#' 
#' @export
qmpe = function(rts, rt_qs, p_q, n_rt, min_p=1e-10){
  
  n_sim = sum(sapply(rts, length))
  n_emp = sum(n_rt)
  
  qmpe = 0
  for(i in 1:nrow(rt_qs)) {
    
    if (any(is.na(rt_qs[i,]))) {
      n_in_bin = n_rt[i]
      sim_in_bin = length(rts[[i]])
    } else {
      n_in_bin = p_q * n_rt[i]
      sim_in_bin = suppressWarnings(hist(rts[[i]], c(0, rt_qs[i,], Inf), plot=F)$counts)
    }
    
    p_in_bin = sim_in_bin / n_sim
    p_in_bin[p_in_bin < min_p] = min_p
    qmpe = qmpe - sum(n_in_bin * log(p_in_bin))
    
  }
  
  if(is.na(qmpe)) qmpe = 1e10
  qmpe
  
}
