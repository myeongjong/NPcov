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
  dmat          <- dist(sic.100$coords)
  dvec          <- as.vector(dmat)
  dmat          <- as.matrix(dmat)

  vals          <- rep(NA, length(dvec))
  for(i in 1:length(vals)) {

    row           <- which(dmat == dvec[i], arr.ind = TRUE)[1,1]
    col           <- which(dmat == dvec[i], arr.ind = TRUE)[1,2]
    average       <- mean(target$data)
    vals[i]       <- (target$data[row] - average) * (target$data[col] - average)
  }

  vhat          <- var(target$data) * (length(target$data)-1) / length(target$data)
  y_err         <- vals[order(dvec)] # / ( var(target$data) * (length(target$data)-1) / length(target$data) )
  x             <- dvec[order(dvec)]
  xmax0         <- max(x)
  x             <- x / xmax0 * xmax

  return( list(r = x, xmax0 = xmax0, xmax = xmax, chat = y_err, vhat = vhat) )
}

output_sic100 <- extract_sic(target = sic.100)
output_sicall <- extract_sic(target = sic.all)

# par(mfrow = c(1, 2))
# plot(output_sic100$r, output_sic100$chat, col = "gray", xlab = "r", ylab = "corr")
# plot(output_sicall$r, output_sicall$chat, col = "gray")
# par(mfrow = c(1, 1))

###########################################################################
### setting
###########################################################################

expand        <- FALSE
eval          <- NULL
k             <- 5
kernel        <- "gaussian"
# shape         <- "general"
method        <- "reflection"
dist_init     <- "exponential"
min_iter      <- 20
max_iter      <- 200
tol           <- 1e-3
check         <- 5
ncores        <- NULL

###########################################################################
### isotropic estimator
###########################################################################

ms            <- c(8, 9, 10)
hs            <- seq(0.15, 0.2, by = 0.01)
inputs        <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus   <- rep(0.1, nrow(inputs))
inputs$ns     <- inputs$ms * 10

ptm_cv_iso    <- proc.time()
cvout_iso     <- cv_pdnr(x = output_sic100$r, y = output_sic100$chat / output_sic100$vhat, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = "general", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv_iso    <- proc.time() - ptm_cv_iso
ptm_cv_iso # 24750.28

ptm_best_iso  <- proc.time()
vhat          <- output_sic100$vhat
for(i in 1:100) {

  message("Step ", i)

  vhat_old      <- vhat
  bestout_iso   <- pdnr(x = output_sic100$r, y = output_sic100$chat / vhat, h = cvout_iso$input.best$h, n = cvout_iso$input.best$n, m = cvout_iso$input.best$m, tau = cvout_iso$input.best$tau, eval = eval, kernel = kernel, shape = "general", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)

  vgrid         <- seq(from = max(0, vhat - 1500), to = vhat + 1500, length.out = 1000)
  msevec        <- c()
  for(j in 1:length(vgrid)) msevec[j] <- mean((output_sic100$chat - vgrid[j] * bestout_iso$yhat)^2)

  vhat          <- vgrid[which.min(msevec)]

  if(abs(vhat_old - vhat) < 0.005 * output_sic100$vhat) break
}
ptm_best_iso  <- proc.time() - ptm_best_iso
i ; ptm_best_iso # 28305.59

vhat_iso      <- vhat

###########################################################################
### monotone estimator
###########################################################################

ms            <- c(2, 3, 4, 5)
hs            <- seq(0.01, 0.1, by = 0.01)
inputs        <- expand.grid(ms, hs) ; colnames(inputs) <- c("ms", "hs")
inputs$taus   <- rep(0.1, nrow(inputs))
inputs$ns     <- inputs$ms * 10

ptm_cv_mon    <- proc.time()
cvout_mon     <- cv_pdnr(x = output_sic100$r, y = output_sic100$chat / output_sic100$vhat, hs = inputs$hs, ns = inputs$ns, ms = inputs$ms, taus = inputs$taus, expand = expand, eval = eval, k = k, kernel = kernel, shape = "monotone", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check, ncores = ncores)
ptm_cv_mon    <- proc.time() - ptm_cv_mon
ptm_cv_mon # 268.34

ptm_best_mon  <- proc.time()
vhat          <- output_sic100$vhat
for(i in 1:100) {

  message("Step ", i)

  vhat_old      <- vhat
  bestout_mon   <- pdnr(x = output_sic100$r, y = output_sic100$chat / vhat, h = cvout_mon$input.best$h, n = cvout_mon$input.best$n, m = cvout_mon$input.best$m, tau = cvout_mon$input.best$tau, eval = eval, kernel = kernel, shape = "monotone", method = method, dist_init = dist_init, min_iter = min_iter, max_iter = max_iter, tol = tol, check = check)

  vgrid         <- seq(from = max(0, vhat - 1500), to = vhat + 1500, length.out = 1000)
  msevec        <- c()
  for(j in 1:length(vgrid)) msevec[j] <- mean((output_sic100$chat - vgrid[j] * bestout_mon$yhat)^2)

  vhat          <- vgrid[which.min(msevec)]

  if(abs(vhat_old - vhat) < 0.005 * output_sic100$vhat) break
}
ptm_best_mon  <- proc.time() - ptm_best_mon
i ; ptm_best_mon # 33.19

vhat_mon      <- vhat

###########################################################################
### save
###########################################################################

plot(output_sic100$r, output_sic100$chat, xlab = "r", ylab = "c", col = "gray", ylim = c(-1.25, 1.25) * output_sic100$vhat)
abline(h = 0, lwd = 2, col = "black", lty = "dashed")
lines(bestout_iso$eval, bestout_iso$yhat * vhat_iso, lwd = 2, col = "red")
lines(bestout_mon$eval, bestout_mon$yhat * vhat_mon, lwd = 2, col = "blue")

save(sic.100, sic.borders, output_sic100, cvout_iso, ptm_cv_iso, bestout_iso, vhat_iso, ptm_best_iso, cvout_mon, ptm_cv_mon, bestout_mon, vhat_mon, ptm_best_mon, file = "out/appout_sic100_gaussian.RData")

