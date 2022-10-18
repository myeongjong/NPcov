###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Iterated Density Estimation Evolutionary Algorithm (IDEA)
###
###########################################################################

###########################################################################
###    Building blocks
###########################################################################

.initialization <- function(n, m, dist_init = "exponential")
{
  if(dist_init == "exponential"){

    storage     <- as.data.frame(matrix(rexp(n * m), m, n))

  } else {

    storage     <- as.data.frame(matrix(rexp(n * m), m, n))

  } ; colnames(storage) <- paste0("Z", 1:n)

  return( storage )
}

.selection <- function(storage, tau, objective)
{
  storage     <- storage[,order(objective)]
  objective   <- objective[order(objective)]

  storage[,(as.integer(tau * ncol(storage)) + 1):ncol(storage)] <- NA

  return( storage )
}

.replacement <- function(storage, tau, h, kernel = "epan", method = "reflection")
{
  n           <- ncol(storage)
  m           <- nrow(storage)
  num_sample  <- as.integer(tau * n)
  data        <- as.numeric(na.omit(unlist(storage)))

  storage           <- cbind(storage[,1:num_sample], matrix(rkrnl((n - num_sample) * m, data, h, kernel, method), m, n-num_sample))
  colnames(storage) <- paste0("Z", 1:n)

  return( storage )
}

.termination <- function(storage_old, storage, h, kernel = "epan", method = "reflection") {

  data_old <- na.omit(unlist(storage_old))
  data_new <- na.omit(unlist(storage))

  abs(sum( log( krnl(data_old, h, data_old, kernel, method) / krnl(data_new, h, data_old, kernel, method) ) )) / length(data_old)

}

###########################################################################
###    IDEA for regression
###########################################################################

#' Iterated Density Estimation Evolutionary Algorithm (IDEA) for regression
#'
#' @param x x
#' @param y y
#' @param h Bandwidth
#' @param n n
#' @param m m
#' @param tau tau
#' @param kernel Type of kernel: epan, gaussian, or uniform
#' @param shape Shape of estimator: general or monotone
#' @param method kernel density estimation/sampling method: reflection or transformation
#' @param dist_init Distribution for initialization
#' @param min_iter Minimum number of iteration
#' @param max_iter Maximum number of iteration
#' @param tol Error tolerance
#' @param check Iteration tolerance
#'
#' @return List
#' @export
#'
#' @examples
#' m           <- 20
#' n           <- 100
#' tau         <- 0.2
#' h           <- 0.1
#' kernel      <- "gaussian"
#' shape       <- "general"
#' method      <- "reflection"
#' dist_init   <- "exponential"
#' min_iter    <- 20
#' max_iter    <- 100
#' wave_model  <- function(h, c_0, c_1, a) ifelse(h == 0, 0, c_0 + c_1 * (1-a*sin(h/a)/h))
#'
#' x           <- seq(0.1, 10, by = 0.1)
#' y           <- wave_model(x, 0, 1, 1)
#' y           <- 1 - y
#' y_err       <- y + rnorm(length(y), 0, 0.2)
#'
#' r           <- seq(0.05, 10, by = 0.05)
#' ptm         <- proc.time()
#' storage     <- IDEA(x, y_err, h, n, m, tau, kernel, shape, method, dist_init, min_iter, max_iter)
#' ptm         <- proc.time() - ptm
#' data        <- na.omit( unlist(storage[[2]]) )
#'
#' z           <- rep(NA, length(r))
#' for(i in 1:length(r)) z[i] <- ksreflect(r = r[i], data = data, h = h, kernel = kernel, shape = shape)
#'
#' plot(x, y_err, col = "gray")
#' lines(r, z, col = "red", lwd = 2)
#' lines(x, y, col = "black", lwd = 2)
#'
#' plot(na.omit(storage$history_obj_min), type="l", ylim = c(4, 15), lty = "dotted", xlab = "Iteration", ylab = "Objective value")
#' lines(na.omit(storage$history_obj_max_sel), col = "red")
#' lines(na.omit(storage$history_obj_mean), col = "blue")
#' lines(na.omit(storage$history_obj_max), lty = "dotted")
#'
#' ptm
IDEA <- function(x, y, h, n, m, tau, kernel = "epan", shape = "general", method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 400, tol = 1e-3, check = 5){

  history_obj_min         <- rep(NA, max_iter)
  history_obj_max         <- rep(NA, max_iter)
  history_obj_max_sel     <- rep(NA, max_iter)
  history_obj_mean        <- rep(NA, max_iter)
  history_KL              <- rep(NA, max_iter)

  storage                 <- .initialization(n, m, dist_init)
  obj_value               <- eval_sse_C(as.matrix(storage), x, y, h, kernel, shape) # eval_sse(storage, x, y, h, kernel, shape)

  history_obj_min[1]      <- min(obj_value)
  history_obj_mean[1]     <- mean(obj_value)
  history_obj_max[1]      <- max(obj_value)

  storage                 <- .selection(storage, tau, obj_value)
  history_obj_max_sel[1]  <- obj_value[order(obj_value)][(as.integer(tau * ncol(storage)))]

  for(i in 1:max_iter) {

    storage_old             <- storage
    storage                 <- .replacement(storage, tau, h, kernel, method)
    obj_value               <- eval_sse_C(as.matrix(storage), x, y, h, kernel, shape) # eval_sse_C(storage, x, y, h, kernel, shape)

    history_obj_min[i+1]    <- min(obj_value)
    history_obj_mean[i+1]   <- mean(obj_value)
    history_obj_max[i+1]    <- max(obj_value)

    storage                 <- .selection(storage, tau, obj_value)
    history_obj_max_sel[i+1]  <- obj_value[order(obj_value)][(as.integer(tau * ncol(storage)))]

    history_KL[i+1]         <- .termination(storage_old, storage, h, kernel, method)
    if(max(tail(na.omit(history_KL), n = check)) <= tol & i >= min_iter) break
  }

  output                  <- list(iteration = i, storage = storage[, 1:(as.integer(tau * ncol(storage)))], history_obj_min = history_obj_min[1:(i+1)],
                                  history_obj_max_sel = history_obj_max_sel[1:(i+1)],
                                  history_obj_mean = history_obj_mean[1:(i+1)],
                                  history_obj_max = history_obj_max[1:(i+1)],
                                  history_KL = history_KL[1:(i+1)])
  return(output)
}
