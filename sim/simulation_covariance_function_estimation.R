###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Simulation - covariance function estimation
###
###########################################################################

library(NPcov)
library(geoR)

rm(list = ls())

###########################################################################
### simulation function
###########################################################################

extract_covpts <- function(coords, data)
{
  dmat        <- dist(coords)
  dvec        <- as.vector(dmat)
  dmat        <- as.matrix(dmat)

  vals        <- rep(NA, length(dvec))
  for(i in 1:length(vals)) {

    row         <- which(dmat == dvec[i], arr.ind = TRUE)[1,1]
    col         <- which(dmat == dvec[i], arr.ind = TRUE)[1,2]
    average     <- mean(data)
    vals[i]     <- (data[row] - average) * (data[col] - average)
  }

  vhat        <- var(data) * (length(data)-1) / length(data)
  y_err       <- vals[order(dvec)] # / ( var(data) * (length(data)-1) / length(data) )
  x           <- dvec[order(dvec)]

  return( list(r = x, chat = y_err, vhat = vhat) )
}

###########################################################################
###
###########################################################################

set.seed(11092022)

simgrf_swave  <- grf(200, cov.pars = c(1, 1), xlims = c(0, 10/sqrt(2)), ylims = c(0, 10/sqrt(2)), cov.model = "wave", nug = 0, nsim = 1, messages = FALSE)
output_swave  <- extract_covpts(simgrf_swave$coords, simgrf_swave$data)
x             <- output_swave$r
y_err         <- output_swave$chat
vhat          <- output_swave$vhat
y             <- 1 - ftn_wave(x, 0, 1, 1)
input_swave   <- list(x = x, y = y, y_err = y_err)

###

ptm1          <- proc.time()
for(i in 1:100) {

  vhat_old      <- vhat
  fit_us_swave  <- pdnr(x = x, y = y_err / vhat, h = 0.2, n = 50, m = 5, tau = 0.1, eval = NULL, kernel = "gaussian", shape = "general", method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 200, tol = 1e-3, check = 5)

  vgrid         <- seq(0.5, 1.5, by = 0.005)
  msevec        <- c()
  for(j in 1:length(vgrid)) msevec[j] <- mean((y_err - vgrid[j] * fit_us_swave$yhat)^2)

  vhat          <- vgrid[which.min(msevec)]

  if(abs(vhat_old - vhat) < 0.01) break
}
ptm1          <- proc.time() - ptm1
i ; ptm1 # 2954.09

vhat_swave    <- vhat

###

plot(x, y_err, col = "grey80", ylim = c(-0.5, 1))
lines(x, y, lwd = 2)
lines(fit_us_swave$eval, fit_us_swave$yhat * vhat_swave, col = "firebrick2", lwd = 2)

# save(simgrf_swave, output_swave, fit_us_swave, vhat_swave, input_swave, file = "out/simout_covariance_function_estimation.RData")

###########################################################################
###
###########################################################################

set.seed(11102022)

simgrf_sexpo  <- grf(200, cov.pars = c(1, 1), xlims = c(0, 10/sqrt(2)), ylims = c(0, 10/sqrt(2)), cov.model = "exponential", nug = 0, nsim = 1, messages = FALSE)
output_sexpo  <- extract_covpts(simgrf_sexpo$coords, simgrf_sexpo$data)
x             <- output_sexpo$r
y_err         <- output_sexpo$chat
vhat          <- output_sexpo$vhat
y             <- 1 - ftn_exponential(x, 0, 1, 1)
input_sexpo   <- list(x = x, y = y, y_err = y_err)

###

ptm2          <- proc.time()
for(i in 1:100) {

  vhat_old      <- vhat
  fit_us_sexpo  <- pdnr(x = x, y = y_err, h = 0.1, n = 100, m = 10, tau = 0.1, eval = NULL, kernel = "gaussian", shape = "monotone", method = "reflection", dist_init = "exponential", min_iter = 20, max_iter = 200, tol = 1e-3, check = 5)

  vgrid         <- seq(0.5, 1.5, by = 0.005)
  msevec        <- c()
  for(j in 1:length(vgrid)) msevec[j] <- mean((y_err - vgrid[j] * fit_us_sexpo$yhat)^2)

  vhat          <- vgrid[which.min(msevec)]

  if(abs(vhat_old - vhat) < 0.01) break
}
ptm2          <- proc.time() - ptm2
i ; ptm2 # 481.84

vhat_sexpo    <- vhat

###

plot(x, y_err, col = "grey80", ylim = c(-0.5, 1))
lines(x, y, lwd = 2)
lines(fit_us_sexpo$eval, fit_us_sexpo$yhat * vhat_sexpo, col = "firebrick2", lwd = 2)

save(simgrf_swave, output_swave, fit_us_swave, vhat_swave, input_swave, simgrf_sexpo, output_sexpo, fit_us_sexpo, vhat_sexpo, input_sexpo, file = "out/simout_covariance_function_estimation.RData")
