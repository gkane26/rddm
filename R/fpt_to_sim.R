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
