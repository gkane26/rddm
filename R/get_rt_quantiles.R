#' Get response time quantiles
#' 
#' @param dat data.table; choice and response time data
#' @param qs numeric; vector of quantiles
#' @param rt_var character; name of rt variable in data table
#' @param conditions character vector; vector of columns indicating the conditions in which to calcluate rt quantiles
#'
#' @return list with data.table withresponse time quantiles and vector of rt quantiles
#' 
#' @export
get_rt_quantiles = function(dat, qs=seq(.1, .9, .2), rt_var="rt", conditions=NULL){
  p_q = diff(c(0, qs, 1))
  rt_qs = dat[, c("n_response" = .N, "p_response" = numeric(1), as.list(quantile(get(rt_var), qs))), c("response", conditions)]
  rt_qs[, p_response := n_response/sum(n_response), conditions]
  list(rt_qs, p_q)
}