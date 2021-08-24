#' Get response time quantiles
#' 
#' @param dat data.table; choice and response time data
#' @param qs numeric; vector of quantiles
#' @param rt_var character; name of rt variable in data table
#' @param conditions character vector; vector of columns indicating the conditions in which to calcluate rt quantiles
#'
#' @return list with data.table with response time quantiles and vector of rt quantiles
#' 
#' @export
get_rt_quantiles = function(dat, qs=seq(.1, .9, .2), rt_var="rt", conditions=NULL){
  
  dat = setDT(dat)
  p_q = diff(c(0, qs, 1))
  
  rt_qs = dat[, data.table(response = c(0, 1, NA),
                           n_response = c(sum((response == 0) & !(is.na(response))),
                                          sum((response == 1) & !(is.na(response))),
                                          sum(is.na(response))),
                           p_response = c(0, 0, 0),
                           rbind(quantile(get(rt_var)[(response == 0) & !(is.na(response))], qs),
                                 quantile(get(rt_var)[(response == 1) & !(is.na(response))], qs),
                                 rep(NA, length(qs)))),
              c(conditions)]
  
  list(rt_qs, p_q)
  
}