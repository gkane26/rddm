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
  upper = hist(sim_data[sim_data$response==1, rt], breaks=bins, plot=F)$counts
  lower = hist(as.numeric(sim_data[sim_data$response==0, rt]), breaks=bins, plot=F)$counts
  cbind(upper, lower)
}
