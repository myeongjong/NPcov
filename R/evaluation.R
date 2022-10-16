###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description :
###
###########################################################################

###########################################################################
###
###########################################################################

#' Computing SSE of an estimator obtained by the IDEA algorithm
#'
#' @param storage data storage
#' @param x Input (e.g., distance)
#' @param y Output (e.g., correlation value)
#' @param h Bandwidth
#' @param kernel Type of kernel: epan, gaussian, or uniform
#' @param shape Shape of estimator: general or monotone
#'
#' @return SSE
#' @export
#'
#' @examples
#' # No example
eval_sse <- function(storage, x, y, h, kernel = "epan", shape = "general")
{
  obj_value   <- rep(0, ncol(storage))
  for(j in 1:ncol(storage)){
    for(k in 1:length(x)){
      obj_value[j] <- obj_value[j] + ( y[k] - ksreflect(r = x[k], data = storage[,j], h = h, kernel = kernel, shape = shape) )^2
    }
  }

  return( obj_value )
}
