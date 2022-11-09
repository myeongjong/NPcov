###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Real data application
###
###########################################################################

library(geoR)
library(NPcov)

rm(list = ls())
set.seed(11072022) # set.seed(11032022)

###########################################################################
### data
###########################################################################

data("SIC", package = "geoR")
rm(sic.some, sic.367)

# pdf("img/sic100.pdf")
# points(sic.100, borders = sic.borders, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5, xlab = "", ylab = "")
# dev.off()
#
# pdf("img/sicall.pdf")
# points(sic.all, borders = sic.borders, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5, xlab = "", ylab = "")
# dev.off()

extract_sic <- function(target, xmax = 8)
{
  dmat        <- dist(sic.100$coords)
  dvec        <- as.vector(dmat)
  dmat        <- as.matrix(dmat)

  vals        <- rep(NA, length(dvec))
  for(i in 1:length(vals)) {

    row         <- which(dmat == dvec[i], arr.ind = TRUE)[1,1]
    col         <- which(dmat == dvec[i], arr.ind = TRUE)[1,2]
    average     <- mean(target$data)
    vals[i]     <- (target$data[row] - average) * (target$data[col] - average)
  }

  vhat        <- var(target$data) * (length(target$data)-1) / length(target$data)
  y_err       <- vals[order(dvec)] / ( var(target$data) * (length(target$data)-1) / length(target$data) )
  x           <- dvec[order(dvec)]
  xmax0       <- max(x)
  x           <- x / xmax0 * xmax

  return( list(r = x, xmax0 = xmax0, xmax = xmax, chat = y_err, vhat = vhat) )
}

output_sic100 <- extract_sic(target = sic.100)
output_sicall <- extract_sic(target = sic.all)

# par(mfrow = c(1, 2))
# plot(output_sic100$r, output_sic100$chat, col = "gray", xlab = "r", ylab = "corr")
# plot(output_sicall$r, output_sicall$chat, col = "gray")
# par(mfrow = c(1, 1))

###########################################################################
### output_sic100
###########################################################################

expand      <- FALSE
eval        <- NULL
k           <- 5
kernel      <- "gaussian"
# shape       <- "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 200
tol         <- 1e-3
check       <- 5
ncores      <- NULL

#####

ms          <- c(7, 8, 9, 10)
hs          <- seq(0.15, 0.2, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

ptm_cv      <- proc.time()
cvout       <- cv_pdnr(x = output_sic100$r, y = output_sic100$chat, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = "general", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv      <- proc.time() - ptm_cv
ptm_cv # 31888.37

ptm_best    <- proc.time()
bestout     <- pdnr(x = output_sic100$r, y = output_sic100$chat, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = "general", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best    <- proc.time() - ptm_best
ptm_best # 2781.91

#####

ms          <- c(2, 3, 4, 5)
hs          <- seq(0.01, 0.1, by = 0.01)
inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus <- rep(0.1, nrow(inputs))
inputs$ns   <- inputs$ms * 10

ptm_cv_mon  <- proc.time()
cvout_mon   <- cv_pdnr(x = output_sic100$r, y = output_sic100$chat, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = "monotone", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv_mon  <- proc.time() - ptm_cv_mon
ptm_cv_mon # 242.58

ptm_best_mon<- proc.time()
bestout_mon <- pdnr(x = output_sic100$r, y = output_sic100$chat, h = cvout_mon$input.best$h, n = cvout_mon$input.best$n, m = cvout_mon$input.best$m, tau = cvout_mon$input.best$tau, eval = eval, kernel = kernel, shape = "monotone", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
ptm_best_mon<- proc.time() - ptm_best_mon
ptm_best_mon # 5.03

#####

plot(output_sic100$r, output_sic100$chat * output_sic100$vhat, xlab = "r", ylab = "c", col = "gray", ylim = c(-1.25, 1.25) * output_sic100$vhat)
abline(h = 0, lwd = 2, col = "black", lty = "dashed")
lines(bestout$eval, bestout$yhat * output_sic100$vhat, lwd = 2, col = "red")
lines(bestout_mon$eval, bestout_mon$yhat * output_sic100$vhat, lwd = 2, col = "blue")

save(sic.100, sic.borders, output_sic100, cvout, ptm_cv, bestout, ptm_best, cvout_mon, ptm_cv_mon, bestout_mon, ptm_best_mon, file = "out/appout_sic100_gaussian.RData")

###########################################################################
### output_sicall
###########################################################################

# ms          <- c(5, 6, 7, 8)
# hs          <- seq(0.05, 0.15, by = 0.05)
# inputs      <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
# inputs$taus <- rep(0.1, nrow(inputs))
# inputs$ns   <- inputs$ms * 10
#
# expand      <- FALSE
# eval        <- NULL
# k           <- 5
# kernel      <- "gaussian"
# # shape       <- "general"
# method      <- "reflection"
# dist_init   <- "exponential"
# min_iter    <- 20
# max_iter    <- 200
# tol         <- 1e-3
# check       <- 5
# ncores      <- NULL
#
# #####
#
# ptm_cv      <- proc.time()
# cvout       <- cv_pdnr(x = output_sicall$r, y = output_sicall$chat, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = "general", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
# ptm_cv      <- proc.time() - ptm_cv
# ptm_cv #
#
# ptm_best    <- proc.time()
# bestout     <- pdnr(x = output_sicall$r, y = output_sicall$chat, h = cvout$input.best$h, n = cvout$input.best$n, m = cvout$input.best$m, tau = cvout$input.best$tau, eval = eval, kernel = kernel, shape = "general", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
# ptm_best    <- proc.time() - ptm_best
# ptm_best #
#
# #####
#
# ptm_cv_mon  <- proc.time()
# cvout_mon   <- cv_pdnr(x = output_sicall$r, y = output_sicall$chat, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = "monotone", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
# ptm_cv_mon  <- proc.time() - ptm_cv_mon
# ptm_cv_mon #
#
# ptm_best_mon<- proc.time()
# bestout_mon <- pdnr(x = output_sicall$r, y = output_sicall$chat, h = cvout_mon$input.best$h, n = cvout_mon$input.best$n, m = cvout_mon$input.best$m, tau = cvout_mon$input.best$tau, eval = eval, kernel = kernel, shape = "monotone", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)
# ptm_best_mon<- proc.time() - ptm_best_mon
# ptm_best_mon #
#
# #####
#
# plot(output_sicall$r, output_sicall$chat * output_sicall$vhat, xlab = "r", ylab = "c", col = "gray", ylim = c(-1.25, 1.25) * output_sicall$vhat)
# abline(h = 0, lwd = 2, col = "black", lty = "dashed")
# lines(bestout$eval, bestout$yhat * output_sicall$vhat, lwd = 2, col = "red")
# lines(bestout_mon$eval, bestout_mon$yhat * output_sicall$vhat, lwd = 2, col = "blue")
#
# save(sic.all, sic.borders, output_sicall, cvout, ptm_cv, bestout, ptm_best, cvout_mon, ptm_cv_mon, bestout_mon, ptm_best_mon, file = "out/appout_sicall_gaussian.RData")
