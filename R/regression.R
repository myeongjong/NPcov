###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : variable selection with cross validation
###
###########################################################################

###########################################################################
###
###########################################################################

#' Positive defnite nonparametric regression (PDNR)
#'
#' @param x x
#' @param y y
#' @param h Bandwidth
#' @param n n
#' @param m m
#' @param tau tau
#' @param eval Evaluation points
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
#' \dontrun{
#' x = seq(0.1, 10, by = 0.1)
#' y = ftn_wave(x, 0, 1, 1)
#' y = 1 - y
#' y_err = y + rnorm(length(y), 0, 0.2)
#'
#' h = 0.08
#' n = 50
#' m = 10
#' tau = 0.2
#' eval = NULL
#' kernel = "uniform"
#' shape = "general"
#' method = "reflection"
#' dist_init = "exponential"
#' min_iter = 20
#' max_iter = 100
#' tol = 1e-3
#' check = 5
#'
#' ptm <- proc.time()
#' result <- pdnr(x = x, y = y_err, h = h, n = n, m = m, tau = tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
#' ptm <- proc.time() - ptm # ~30s
#'
#' plot(x, y_err)
#' lines(x, y, type="l", col="red")
#' lines(result$eval, result$yhat)
#'
#' plot(na.omit(result$history_obj_max), type="l", ylim = c(0, 15))
#' lines(na.omit(result$history_obj_mean))
#' lines(na.omit(result$history_obj_max_sel))
#' lines(na.omit(result$history_obj_min))
#' }
pdnr <- function(x, y, h, n, m, tau, eval = NULL, kernel = "epan", shape = "general", method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 400, tol = 1e-3, check = 5)
{
  ### checkargs
  if(is.null(eval)) eval <- x

  ### IDEA
  output    <- IDEA(x = x, y = y, h = h, n = n, m = m, tau = tau, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
  pseudo    <- as.numeric(na.omit(unlist(output$storage)))

  ### ksreflect
  z         <- rep(NA, length(eval))
  for(i in 1:length(z)) z[i] <- ksreflect_C(eval[i], pseudo, h, kernel, shape)

  ### return
  result    <- list(yhat = z, eval = eval, pseudo = pseudo, iteration = output$iteration,
                    history_obj_min = output$history_obj_min,
                    history_obj_max_sel = output$history_obj_max_sel,
                    history_obj_mean = output$history_obj_mean,
                    history_obj_max = output$history_obj_max,
                    history_KL = output$history_KL)
  return( result )
}

#' Computing MS(P)E of PDNR model
#'
#' @param x_train x_train
#' @param y_train y_train
#' @param x_test x_test
#' @param y_test y_test
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
#' @return MS(P)E
#' @export
#'
#' @examples
#' # No example
mse_pdnr <- function(x_train, y_train, x_test, y_test, h, n, m, tau, kernel = "epan", shape = "general", method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 400, tol = 1e-3, check = 5)
{
  ### IDEA
  output    <- IDEA(x = x_train, y = y_train, h = h, n = n, m = m, tau = tau, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
  pseudo    <- as.numeric(na.omit(unlist(output$storage)))

  ### ksreflect
  yhat      <- rep(NA, length(x_test))
  for(i in 1:length(x_test)) yhat[i] <- ksreflect_C(x_test[i], pseudo, h, kernel, shape)

  return( mean((yhat - y_test)^2) )
}

#' Variable selection of PDNR using cross validation
#'
#' @param x x
#' @param y y
#' @param hs Possible values of h
#' @param ns Possible values of n
#' @param ms Possible values of m
#' @param taus Possible values of tau
#' @param expand Logical: Use of expand.grid or not for inputs?
#' @param eval Evaluation points
#' @param k Number of folds
#' @param kernel Type of kernel: epan, gaussian, or uniform
#' @param shape Shape of estimator: general or monotone
#' @param method kernel density estimation/sampling method: reflection or transformation
#' @param dist_init Distribution for initialization
#' @param min_iter Minimum number of iteration
#' @param max_iter Maximum number of iteration
#' @param tol Error tolerance
#' @param check Iteration tolerance
#' @param ncores Number of cores
#'
#' @return List
#' @export
#'
#' @examples
#' \dontrun{
#' x = seq(0.1, 10, by = 0.1)
#' y = ftn_wave(x, 0, 1, 1)
#' y = 1 - y
#' y_err = y + rnorm(length(y), 0, 0.2)
#' hs = c(0.1, 0.2)
#' ns = c(100, 200)
#' ms = c(10, 20)
#' taus = c(0.1, 0.2)
#' expand = TRUE
#' eval = NULL
#' k = 5
#' kernel = "epan"
#' shape = "general"
#' method = "reflection"
#' dist_init = "exponential"
#' min_iter = 20
#' max_iter = 400
#' tol = 1e-3
#' check = 5
#' ncores = NULL
#'
#' ptm <- proc.time()
#' output <- cv_pdnr(x = x, y = y_err, hs = hs, ns = ns, ms = ms, taus = taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
#' ptm <- proc.time() - ptm # ~ 1.5 hr
#' }
cv_pdnr <- function(x, y, hs, ns, ms, taus, expand = TRUE, eval = NULL, k = 5, kernel = "epan", shape = "general", method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 400, tol = 1e-3, check = 5, ncores = NULL)
{
  ### checkargs
  if(is.null(eval)) eval <- x

  if( expand == TRUE) {

    inputs            <- expand.grid(hs, ns, ms, taus)

  } else {

    inputs            <- data.frame(hs, ns, ms, taus)
  }

  colnames(inputs)  <- c("h", "n", "m", "tau")

  `%do%`            <- foreach::`%do%`
  `%dopar%`         <- foreach::`%dopar%`

  ### k-fold
  idx               <- cut(sample(1:length(x)), breaks = k, labels = FALSE)

  ### ncores
  if(is.null(ncores)) {
    no_cores          <- parallel::detectCores() - 2
  } else {
    no_cores          <- ncores
  }

  ### parallel CV
  output            <- list()
  cl                <- parallel::makeCluster(no_cores)

  doParallel::registerDoParallel(cl)
  output            <- foreach::foreach(h = inputs$h, n = inputs$n, m = inputs$m, tau = inputs$tau, .packages = c("NPcov")) %dopar% {

    result            <- rep(NA, k)
    for(i in 1:k) result[i] <- mse_pdnr(x_train = x[idx != i], y_train = y[idx != i], x_test = x[idx == i], y_test = y[idx == i], h = h, n = n, m = m, tau = tau, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)

    result
  }

  parallel::stopCluster(cl)

  ### performance
  perfmat           <- matrix(NA, nrow = nrow(inputs), ncol = k)
  for(i in 1:length(output)) perfmat[i, ] <- output[[i]]
  perfmat           <- as.data.frame(perfmat)
  colnames(perfmat) <- paste0("MSE (fold_", 1:k, ")")
  perfmat$`MSE (overall)`   <- rowMeans(perfmat)
  result            <- cbind(inputs, perfmat)

  ### best model
  idx.best          <- which.min(perfmat$`MSE (overall)`)
  inputs.best       <- inputs[idx.best, ]
  model.best        <- pdnr(x = x, y = y, h = inputs.best$h, n = inputs.best$n, m = inputs.best$m, tau = inputs.best$tau, eval = NULL, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)

  ### return
  return( list(input.best = inputs.best, model.best = model.best, cvresult = result) )
}
