###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Drawing random sample from kernel density estimation
###
###    Reference : Silverman, B. W. (1986). Density estimation for statistics and data analysis (Vol. 26). CRC press.
###
###########################################################################

###########################################################################
###
###########################################################################

#' Drawing random sample from kernel density estimation
#'
#' @param n Desired size of random sample
#' @param data pseudo data
#' @param h Bandwidth
#' @param kernel Type of kernel: epan, gaussian, or uniform
#' @param method Method of sampling: transformation or reflection
#'
#' @return Vector of the random sample
#' @export
#'
#' @examples
#' n = 100
#' data = c(0.5, 10)
#' h = 5
#'
#' par(mfrow = c(1, 3))
#' hist(rkrl(n, data, h, "epan", "reflection"), main = "Epanechnikov")
#' hist(rkrl(n, data, h, "gaussian", "reflection"), main = "Gaussian")
#' hist(rkrl(n, data, h, "uniform", "reflection"), main = "Uniform")
#' par(mfrow = c(1, 1))
rkrl <- function(n, data, h, kernel = "epan", method = "reflection")
{
  if(kernel == "epan") {

    sam <- .generation_epanechnikov(n, data, h, method)

  } else if(kernel == "gaussian") {

    sam <- .generation_gaussian(n, data, h, method)

  } else if(kernel == "uniform") {

    sam <- .generation_uniform(n, data, h, method)

  } else {

    sam <- NULL
  }

  return( sam )
}

###########################################################################
###    Building blocks
###########################################################################

.transformation <- function(data, adj_positive) log(data + adj_positive)

.epan_dist <- function(x) ifelse(x <= -1, 0, ifelse(x >=1, 1, 0.25 * (2 + 3 * x - x^3)))

###########################################################################
###
###########################################################################

.generation_uniform <- function(n, data, h, method = "reflection")
{
  if(method == "transformation") {

    x       <- .transformation(data, 0)
    nsam    <- length(x)
    x_ind   <- sample(1:nsam, size = n, replace = TRUE)
    sam     <- x[x_ind] + 0.5 * h * runif(n = n, min = -1, max = 1)

    return( exp(sam) )

  } else {

    x       <- data
    nsam    <- length(x)
    x_ind   <- sample(1:nsam, size = n, replace = TRUE)
    sam     <- x[x_ind] + 0.5 * h * runif(n = n, min = -1, max = 1)

    return( abs(sam) )
  }
}

.generation_gaussian <- function(n, data, h, method = "reflection")
{
  if(method == "transformation") {

    x       <- .transformation(data, 0)
    nsam    <- length(x)
    x_ind   <- sample(1:nsam, size = n, replace = TRUE)
    sam     <- x[x_ind] + h * rnorm(n = n, mean = 0, sd = 1)

    return( exp(sam) )

  } else {

    x       <- data
    nsam    <- length(x)
    x_ind   <- sample(1:nsam, size = n, replace = TRUE)
    sam     <- x[x_ind] + h * rnorm(n = n, mean = 0, sd = 1)

    return( abs(sam) )
  }
}

.generation_epanechnikov <- function(n, data, h, method = "reflection")
{
  if(method == "transformation") {

    x       <- .transformation(data, 0)
    eps     <- runif(n = n, min = 0, max = 1)

    nsam    <- length(x)
    x_ind   <- sample(1:nsam, size = n, replace = TRUE)
    sam     <- x[x_ind] + h * unlist(lapply(eps, function(x) uniroot(f = function(y) .epan_dist(y)-x, interval = c(-1, 1))$root))

    return( exp(sam) )

  } else {

    x       <- data
    eps     <- runif(n = n, min = 0, max = 1)

    nsam    <- length(x)
    x_ind   <- sample(1:nsam, size = n, replace = TRUE)
    sam     <- x[x_ind] + h * unlist(lapply(eps, function(x) uniroot(f = function(y) .epan_dist(y)-x, interval = c(-1, 1))$root))

    return( abs(sam) )
  }
}
