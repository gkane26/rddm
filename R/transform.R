#' Parameter Trasnform functions
#' 
#' @description Parameter Transform Functions
#' 
#' @param x numeric; vector of parameter values
#' @param lower numeric; parameter lower bound
#' @param upper numeric; parameter upper bound
#' 
#' @return transformed parameter vector
#' 
#' @name transform
#' @rdname transform
#'  
#' @export
logistic = function(x, lower, upper) {
  (upper - lower) / (1 + exp(-x)) + lower
}

#' @rdname transform
#' @export
inv_logistic = function(x, lower, upper) {
  -log((upper-lower)/(x-lower) - 1)
}


