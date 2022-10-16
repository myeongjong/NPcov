###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description :
###
###    Reference : Cressie, N. (2015). Statistics for spatial data. John Wiley & Sons.
###
###########################################################################

###########################################################################
###    Some isotropic variogram models (p61)
###########################################################################

#' Linear isotropic semivariogram model
#'
#' @param h h
#' @param c c
#' @param b b
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' x = seq(0, 15, by = 0.001)
#' y = ftn_linear(x, 0, 1)
#' plot(x, y, type="l")
ftn_linear <- function(h, c, b)
{
  ifelse(h == 0, 0, c + b * h)
}

#' Spherical isotropic semivariogram model
#'
#' @param h h
#' @param c_0 c0
#' @param c_1 c1
#' @param a a
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' x = seq(0, 15, by = 0.001)
#' y = ftn_spherical(x, 0, 1, 1)
#' plot(x, y, type="l")
ftn_spherical <- function(h, c_0, c_1, a)
{
  ifelse(h == 0, 0, ifelse(h <= a, c_0 + c_1 * (1.5*(h/a) - 0.5*(h/a)^3), c_0 + c_1))
}

#' Exponential isotropic semivariogram model
#'
#' @param h h
#' @param c_0 c0
#' @param c_1 c1
#' @param a a
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' x = seq(0, 15, by = 0.001)
#' y = ftn_exponential(x, 0, 1, 1)
#' plot(x, y, type="l")
ftn_exponential <- function(h, c_0, c_1, a)
{
  ifelse(h == 0, 0, c_0 + c_1 * (1-exp(-h/a)))
}

#' Rational quadratic semivariogram model
#'
#' @param h h
#' @param c_0 c0
#' @param c_1 c1
#' @param a a
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' x = seq(0, 15, by = 0.001)
#' y = ftn_rational_quadratic(x, 0, 1, 1)
#' plot(x, y, type="l")
ftn_rational_quadratic <- function(h, c_0, c_1, a)
{
  ifelse(h == 0, 0, c_0 + c_1 * h^2 / (1 + (h^2)/a))
}

#' Wave isotropic semivariogram model
#'
#' @param h h
#' @param c_0 c0
#' @param c_1 c1
#' @param a a
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' x = seq(0, 15, by = 0.001)
#' y = ftn_wave(x, 0, 1, 1)
#' plot(x, y, type="l")
ftn_wave <- function(h, c_0, c_1, a)
{
  ifelse(h == 0, 0, c_0 + c_1 * (1-a*sin(h/a)/h))
}

#' Power isotropic semivariogram model
#'
#' @param h h
#' @param c c
#' @param b b
#' @param lambda lambda
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' x = seq(0, 15, by = 0.001)
#' y = ftn_power(x, 0, 1, 1)
#' plot(x, y, type="l")
ftn_power <- function(h, c, b, lambda)
{
  ifelse(h == 0, 0, c + b * h^lambda)
}
