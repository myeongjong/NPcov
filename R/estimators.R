###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Functions for implementing kernel-type estimators based on the reflection method
###
###########################################################################

###########################################################################
###
###########################################################################

#' Kernel-type estimator based on the reflection method
#'
#' @param r Distance
#' @param data Pseudo data
#' @param h Bandwidth
#' @param kernel Type of kernel: epan, gaussian, or uniform
#' @param shape Shape of estimator: general or monotone
#'
#' @return Numeric value of the estimator
#' @export
#'
#' @examples
#' r <- seq(0.05, 5, by = 0.05)
#' data <- rexp(10)
#' v1 <- c() ; v2 <- c() ; v3 <- c()
#' v4 <- c() ; v5 <- c() ; v6 <- c()
#'
#' for(i in 1:length(r)) v1[i] <- ksreflect(r = r[i], data = data, h = 1, kernel = "epan", shape = "general")
#' for(i in 1:length(r)) v2[i] <- ksreflect(r = r[i], data = data, h = 1, kernel = "gaussian", shape = "general")
#' for(i in 1:length(r)) v3[i] <- ksreflect(r = r[i], data = data, h = 1, kernel = "uniform", shape = "general")
#' for(i in 1:length(r)) v4[i] <- ksreflect(r = r[i], data = data, h = 1, kernel = "epan", shape = "monotone")
#' for(i in 1:length(r)) v5[i] <- ksreflect(r = r[i], data = data, h = 1, kernel = "gaussian", shape = "monotone")
#' for(i in 1:length(r)) v6[i] <- ksreflect(r = r[i], data = data, h = 1, kernel = "uniform", shape = "monotone")
#'
#' par(mfrow = c(1, 2))
#' plot(r, v1, type="l", col = "red", ylim = c(-0.5, 2))
#' lines(r, v2, col="blue")
#' lines(r, v3, col="green")
#'
#' plot(r, v4, type="l", col = "red", ylim = c(-0.5, 2))
#' lines(r, v5, col="blue")
#' lines(r, v6, col="green")
#' par(mfrow = c(1, 1))
ksreflect <- function(r, data, h, kernel = "epan", shape = "general")
{
  if( shape %in% c("general") ) {

    output <- .estimator_2d_reflection(r = r, data = data, h = h, kernel = kernel)

  } else if( shape %in% c("monotone") ) {

    output <- .estimator_infd_reflection(r = r, data = data, h = h, kernel = kernel)

  } else {

    stop("The input type must be either general or monotone.")
  }

  return( output )
}

###########################################################################
###
###########################################################################

.est_2d_epan_refl <- function(r, data, h)
{
  value   <- 0
  for(j in 1:length(data)) value <- value + (1-(data[j]/h)^2) * (.lmda1(r*(data[j]+h)) - .lmda1(r*(data[j]-h))) - ((r^2*(data[j]+h)^2 * .besselJ_modified(r*(data[j]+h), nu = 1) - .lmda0(r*(data[j]+h))) - (r^2*(data[j]-h)^2 * .besselJ_modified(r*(data[j]-h), nu = 1) - .lmda0(r*(data[j]-h)))) / (r * h)^2 + (2 * data[j] / r / h^2) * (r*(data[j]+h) * .besselJ_modified(r*(data[j]+h), nu = 1) - r*(data[j]-h) * .besselJ_modified(r*(data[j]-h), nu = 1))

  value   <- value * 3 / 4 / length(data) / h / r
  return( value )
}

.est_2d_gaussian_refl <- function(r, data, h)
{
  value   <- 0
  for(j in 1:length(data)) value <- value + .besselJ_generalized(r, h, data[j])

  value   <- value * exp(-0.25 * h^2 * r^2) / pi / length(data)
  return( value )
}

.est_2d_uniform_refl <- function(r, data, h)
{
  value   <- 0
  for(j in 1:length(data)) value <- value + .lmda1(r*(data[j]+h)) - .lmda1(r*(data[j]-h))

  value   <- value / 2 / length(data) / h / r
  return( value )
}

###########################################################################
###
###########################################################################

.est_infd_epan_refl <- function(r, data, h)
{
  value   <- 0
  for(j in 1:length(data)) value <- value + (1-(data[j]/h)^2) * (sqrt(pi)/r) * (pnorm(sqrt(2) * r * (data[j]+h)) - pnorm(sqrt(2) * r * (data[j]-h))) + (data[j] / h^2 / r^2) * (exp(-(r * (data[j]-h))^2) - exp(-(r * (data[j]+h))^2)) - (1/2/r^2/h^2) * (((data[j]-h)*exp(-(r * (data[j]-h))^2) - (data[j]+h)*exp(-(r * (data[j]+h))^2)) + (sqrt(pi)/r) * (pnorm(sqrt(2) * r * (data[j]+h)) - pnorm(sqrt(2) * r * (data[j]-h))))

  value   <- value * 3 / 4 / length(data) / h
  return( value )
}

.est_infd_gaussian_refl <- function(r, data, h)
{
  value   <- 0
  A       <- r^2 + 1/(2 * h^2)
  for(j in 1:length(data)) value <- value + exp((data[j]/(h^2))^2/4/A - data[j]^2/2/h^2)

  value   <- value / (length(data) * h) / (sqrt(2 * A)) ##############WARNING
  return( value )
}

.est_infd_uniform_refl <- function(r, data, h)
{
  value   <- 0
  for(j in 1:length(data)) value <- value + pnorm(sqrt(2)* r*(data[j]+h)) - pnorm(sqrt(2)* r*(data[j]-h))

  value   <- value * sqrt(pi) / 2 / length(data) / h / r
  return( value )
}

###########################################################################
###
###########################################################################

.estimator_2d_reflection <- function(r, data, h, kernel = "epan")
{
  if(kernel == "epan"){

    .est_2d_epan_refl(r, data, h)

  } else if(kernel == "gaussian"){

    .est_2d_gaussian_refl(r, data, h)

  } else if(kernel == "uniform"){

    .est_2d_uniform_refl(r, data, h)

  } else {

    return( NA )
  }
}

.estimator_infd_reflection <- function(r, data, h, kernel = "epan")
{
  if(kernel == "epan"){

    .est_infd_epan_refl(r, data, h)

  } else if(kernel == "gaussian"){

    .est_infd_gaussian_refl(r, data, h)

  } else if(kernel == "uniform"){

    .est_infd_uniform_refl(r, data, h)

  } else {

    return( NA )
  }
}

###########################################################################
###
###########################################################################


