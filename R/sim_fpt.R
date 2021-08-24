#' Get simulated responses and rts from first passage time density
#' 
#' @param fpt matrix; 2 X time points matrix with first passage time for upper and lower boundary crossings
#' @param dt numeric; size of time bins used to calculate first passage time density, default = .01
#' @param n integer; number of trials to simulate 
#'
#' @return data.table with two columns: response and rt
#'
#' @import data.table
#'  
#' @export
fpt_to_sim = function(fpt, dt=.01, n=10000) {
  
  tseq = seq(dt, nrow(fpt)*dt, dt)
  n_fpt = round(fpt * n)
  
  responses = c()
  rts = c()
  for(i in 1:nrow(fpt)) {
    responses = c(responses, rep(1, n_fpt[i, 1]), rep(0, n_fpt[i, 2]))
    rts = c(rts, rep(tseq[i], sum(n_fpt[i,])))
  }

  data.frame(response=responses, rt=rts)
  
}

#' Get first passage time density from simulated diffusion model data
#' 
#' @param sim_data numeric; drift rate
#' @param sim_dt numeric; time step of simulation, default = .001
#' @param fpt_dt numeric; size of time bins to calculate first passage time, default = .01
#' @param max_time numeriic; max time to calculate probability of first passage times
#'
#' @return matrix with two columns: first column is p(upper), second is p(lower)
#' 
#' @export
sim_to_fpt = function(sim_data, sim_dt=.001, fpt_dt=.01, max_time=10){
  bins = seq(0, max_time+fpt_dt, fpt_dt)
  upper = hist(as.numeric(sim_data[sim_data$response==1, "rt"]), breaks=bins, plot=F)$counts
  lower = hist(as.numeric(sim_data[sim_data$response==0, "rt"]), breaks=bins, plot=F)$counts
  cbind(upper, lower) / nrow(sim_data)
}

