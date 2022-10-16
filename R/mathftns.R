###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Math functions for implementing the kernel-type estimators using the reflection method
###
###########################################################################

.besselJ_modified <- function(x, nu)
{
  if(x >= 0){

    return( besselJ(x, nu) )

  } else {

    return( (-1)^nu * besselJ(-x, nu) )
  }
}

.struveH_modified <- function(x, nu)
{
  if(x >= 0){

    return( RandomFieldsUtils::struveH(x, nu) )

  } else {

    return( (-1)^(nu + 1) * RandomFieldsUtils::struveH(-x, nu) )
  }
}

.lmda0 <- function(x)
{
  output      <- (pi * x / 2) * (.besselJ_modified(x, nu = 1) * .struveH_modified(x, nu = 0) - .besselJ_modified(x, nu = 0) * .struveH_modified(x, nu = 1))

  return( output )
}

.lmda1 <- function(x)
{
  output      <- x * .besselJ_modified(x, nu = 0) + .lmda0(x)

  return( output )
}

.besselJ_generalized <- function(r, h, t)
{
  output      <- integrate(function(x) exp(0.25 * h^2 * r^2 * cos(2 * x)) * cos(r * t * sin(x)), lower = 0, upper = pi, stop.on.error = FALSE)

  return( output$value )
}
