###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Simulation - regression with white noise errors
###
###########################################################################

library(NPcov)
library(np)

rm(list = ls())
set.seed(11072022)

###########################################################################
### simulation function
###########################################################################

# c(2, 3, 4, 5, 6)
# seq(0.1, 0.2, by = 0.01)

simReg_without_cv <- function(nsim = 1, n_train = 200, n_test = 200, covftn = function(x) ftn_wave(x, 0, 1, 0.5), m, h, shape = "general")
{
  output        <- matrix(NA, nrow = nsim, ncol = 3)
  history_data  <- list()
  history_fit   <- list()
  for(i in 1:nsim) {

    x             <- sort(runif(n_train, min = 0.1, max = 10))
    y             <- covftn(x)
    y             <- 1 - y
    y_err         <- y + rnorm(length(y), 0, 0.2)
    x_test        <- sort(runif(n_test, min = 0.1, max = 10))
    y_test        <- covftn(x_test)
    y_test        <- 1 - y_test

    history_data[[i]] <- list(x = x, y = y, y_err = y_err, x_test = x_test, y_test = y_test)

    # inputs        <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
    # inputs$taus   <- rep(0.1, nrow(inputs))
    # inputs$ns     <- inputs$ms * 10
    # result_cv_us  <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = FALSE, eval = NULL, k = 5, kernel = "gaussian", shape = shape, method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 200, tol = 1e-3, check = 5, ncores = NULL)
    fit_us        <- pdnr(x = x, y = y_err, h = h, n = m * 10, m = m, tau = 0.1, eval = x_test, kernel = "gaussian", shape = shape, method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 200, tol = 1e-3, check = 5)

    result_cv_lc  <- npregbw(formula = y_err ~ x, regtype = "lc", bwmethod = "cv.ls", ckertype = "gaussian", folds = 5)
    fit_lc        <- npreg(result_cv_lc, newdata = data.frame(x = x_test))

    result_cv_ll  <- npregbw(formula = y_err ~ x, regtype = "ll", bwmethod = "cv.ls", ckertype = "gaussian", folds = 5)
    fit_ll        <- npreg(result_cv_ll, newdata = data.frame(x = x_test))

    history_fit[[i]] <- list(result_cv_lc = result_cv_lc, fit_lc = fit_lc, result_cv_ll = result_cv_ll, fit_ll = fit_ll, fit_us = fit_us)

    # plot(x, y_err, col = "grey80", ylim = c(-0.5, 1))
    # lines(x, y, lwd = 2)
    # lines(as.numeric(t(fit_ll$eval)), fit_ll$mean, col="chartreuse2", lwd = 2)
    # lines(as.numeric(t(fit_lc$eval)), fit_lc$mean, col="dodgerblue2", lwd = 2)
    # lines(fit_us$eval, fit_us$yhat, col = "firebrick2", lwd = 2)

    output[i, ]   <- c(sqrt(mean((fit_lc$mean - y_test)^2)), sqrt(mean((fit_ll$mean - y_test)^2)), sqrt(mean((fit_us$yhat - y_test)^2)))
  }

  output        <- as.data.frame(output)
  colnames(output) <- c("LC", "LL", "US")
  rownames(output) <- paste0("sim ", 1:nsim)
  return( list(history_data = history_data, history_fit = history_fit, output = output) )
}

simReg_with_cv <- function(nsim = 1, n_train = 200, n_test = 200, covftn = function(x) ftn_wave(x, 0, 1, 0.5), ms = c(4, 5, 6, 7, 8), hs = seq(0.15, 0.25, by = 0.01), shape = "general")
{
  output        <- matrix(NA, nrow = nsim, ncol = 3)
  history_data  <- list()
  history_fit   <- list()
  for(i in 1:nsim) {

    x             <- sort(runif(n_train, min = 0.1, max = 10))
    y             <- covftn(x)
    y             <- 1 - y
    y_err         <- y + rnorm(length(y), 0, 0.2)
    x_test        <- sort(runif(n_test, min = 0.1, max = 10))
    y_test        <- covftn(x_test)
    y_test        <- 1 - y_test

    history_data[[i]] <- list(x = x, y = y, y_err = y_err, x_test = x_test, y_test = y_test)

    inputs        <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
    inputs$taus   <- rep(0.1, nrow(inputs))
    inputs$ns     <- inputs$ms * 10
    result_cv_us  <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = FALSE, eval = NULL, k = 5, kernel = "gaussian", shape = shape, method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 200, tol = 1e-3, check = 5, ncores = NULL)
    fit_us        <- pdnr(x = x, y = y_err, h = result_cv_us$input.best$h, n = result_cv_us$input.best$n, m = result_cv_us$input.best$m, tau = result_cv_us$input.best$tau, eval = x_test, kernel = "gaussian", shape = shape, method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 200, tol = 1e-3, check = 5)

    result_cv_lc  <- npregbw(formula = y_err ~ x, regtype = "lc", bwmethod = "cv.ls", ckertype = "gaussian", folds = 5)
    fit_lc        <- npreg(result_cv_lc, newdata = data.frame(x = x_test))

    result_cv_ll  <- npregbw(formula = y_err ~ x, regtype = "ll", bwmethod = "cv.ls", ckertype = "gaussian", folds = 5)
    fit_ll        <- npreg(result_cv_ll, newdata = data.frame(x = x_test))

    history_fit[[i]] <- list(result_cv_lc = result_cv_lc, fit_lc = fit_lc, result_cv_ll = result_cv_ll, fit_ll = fit_ll, result_cv_us = result_cv_us, fit_us = fit_us)

    # plot(x, y_err, col = "grey80", ylim = c(-0.5, 1))
    # lines(x, y, lwd = 2)
    # lines(as.numeric(t(fit_ll$eval)), fit_ll$mean, col="chartreuse2", lwd = 2)
    # lines(as.numeric(t(fit_lc$eval)), fit_lc$mean, col="dodgerblue2", lwd = 2)
    # lines(fit_us$eval, fit_us$yhat, col = "firebrick2", lwd = 2)

    output[i, ]   <- c(sqrt(mean((fit_lc$mean - y_test)^2)), sqrt(mean((fit_ll$mean - y_test)^2)), sqrt(mean((fit_us$yhat - y_test)^2)))
  }

  output        <- as.data.frame(output)
  colnames(output) <- c("LC", "LL", "US")
  rownames(output) <- paste0("sim ", 1:nsim)
  return( list(history_data = history_data, history_fit = history_fit, output = output) )
}

###########################################################################
###
###########################################################################

ptm1 <- proc.time()
result_wave_complex <- simReg_without_cv(nsim = 200, n_train = 200, n_test = 200, covftn = function(x) ftn_wave(x, 0, 1, 0.5), m = 4, h = 0.16, shape = "general")
# result_wave_complex <- simReg_with_cv(nsim = 50, n_train = 200, n_test = 200, covftn = function(x) ftn_wave(x, 0, 1, 0.5), ms = c(2, 3, 4, 5, 6), hs = seq(0.1, 0.2, by = 0.01), shape = "general")
ptm1 <- proc.time() - ptm1

ptm1 # 2389.09
colMeans(result_wave_complex$output)

ptm2 <- proc.time()
result_spherical_complex <- simReg_without_cv(nsim = 200, n_train = 200, n_test = 200, covftn = function(x) ftn_spherical(x, 0, 0.8, 2), m = 6, h = 0.01, shape = "monotone")
# result_spherical_complex <- simReg_with_cv(nsim = 50, n_train = 200, n_test = 200, covftn = function(x) ftn_spherical(x, 0, 0.8, 2), ms = c(6, 7, 8, 9, 10), hs = seq(0.01, 0.1, by = 0.01), shape = "monotone")
ptm2 <- proc.time() - ptm2

ptm2 # 2619.30
colMeans(result_spherical_complex$output)

save(ptm1, ptm2, result_wave_complex, result_spherical_complex, file = "out/simout_regression_with_white_noise_errors.RData")
