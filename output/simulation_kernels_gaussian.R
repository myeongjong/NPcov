###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Gaussian kernel
###
###########################################################################

library(NPcov)

###########################################################################
### simulation 1: simple wave
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_wave(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(6, 7, 8, 9, 10)
hs          <- seq(0.1, 0.2, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
shape       <- "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 1648.21

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 22.90

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "rdata/simout_gaussian_wave_simple.RData")

###########################################################################
### simulation 2: complex wave
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_wave(x, 0, 1, 0.5)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(4, 5, 6, 7, 8)
hs          <- seq(0.15, 0.25, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
shape       <- "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 1195.69

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 43.05

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "rdata/simout_gaussian_wave_complex.RData")

###########################################################################
### simulation 3: simple spherical case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(3, 4, 5, 6, 7)
hs          <- seq(0.15, 0.25, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
shape       <- "monotone" # "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 480.74

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 7.87

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "rdata/simout_gaussian_sphe_simple.RData")

###########################################################################
### simulation 4: complex spherical case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 0.8, 2)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(6, 7, 8, 9, 10)
hs          <- seq(0.01, 0.1, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
shape       <- "monotone" # "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 889.16

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 13.16

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "rdata/simout_gaussian_sphe_complex.RData")

###########################################################################
### simulation 5: simple exponential case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_exponential(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(5, 6, 7, 8, 9)
hs          <- seq(0.2, 0.3, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
shape       <- "monotone" # "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 2070.68

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 20.25

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "rdata/simout_gaussian_expo_simple.RData")

###########################################################################
### simulation 6: complex exponential case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_exponential(x, 0, 1, 4)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(2, 3, 4, 5, 6)
hs          <- seq(0.05, 0.15, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
shape       <- "monotone" # "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = x, y = y_err, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 4372.51

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best #  21.61

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "rdata/simout_gaussian_expo_complex.RData")
