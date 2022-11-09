###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Epanechnikov kernel
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

ms          <- c(3, 4, 5, 6, 7)
hs          <- seq(0.25, 0.35, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "epan"
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
ptm_cv #

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best #

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "out/simout_epan_wave_simple.RData")

###########################################################################
### simulation 2: complex wave
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_wave(x, 0, 1, 0.5)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(16, 17, 18, 19, 20)
hs          <- seq(0.25, 0.35, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "epan"
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
ptm_cv # 29785.15

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 438.60

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "out/simout_epan_wave_complex.RData")

###########################################################################
### simulation 3: simple spherical case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(16, 17, 18, 19, 20)
hs          <- seq(0.01, 0.1, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "epan"
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
ptm_cv # 16483.33

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 466.32

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "out/simout_epan_sphe_simple.RData")

###########################################################################
### simulation 4: complex spherical case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 0.8, 2)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(16, 17, 18, 19, 20)
hs          <- seq(0.01, 0.1, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "epan"
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
ptm_cv #

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best #

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "out/simout_epan_sphe_complex.RData")

###########################################################################
### simulation 5: simple exponential case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_exponential(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(3, 4, 5, 6, 7)
hs          <- seq(0.25, 0.35, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "epan"
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
ptm_cv # 330.87

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 8.91

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "out/simout_epan_expo_simple.RData")

###########################################################################
### simulation 6: complex exponential case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_exponential(x, 0, 1, 4)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

ms          <- c(3, 4, 5, 6, 7)
hs          <- seq(0.25, 0.35, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "epan"
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
ptm_cv # 324.45

ptm_best    <- proc.time()
bestout     <- pdnr(x = x, y = y_err, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = shape, method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 4.47

plot(x, y_err, xlab = "x", ylab = "y")
abline(h = 0, lwd = 2, col = "gray", lty = "dashed")
lines(x, y, lwd = 2, col = "black")
lines(bestout$eval, bestout$yhat, lwd = 2, col = "red")

save(cvout, ptm_cv, bestout, ptm_best, file = "out/simout_epan_expo_complex.RData")
