###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Estimating kernel densities
###
###    Reference : Silverman, B. W. (1986). Density estimation for statistics and data analysis (Vol. 26). CRC press.
###
###########################################################################

###########################################################################
###
###########################################################################

#' Kernel density estimation using pseudo data
#'
#' @param x x
#' @param h Bandwidth
#' @param data pseudo data
#' @param kernel Type of kernel: epan, gaussian, or uniform
#' @param method Method of sampling: transformation or reflection
#'
#' @return Values of kernel density estimation
#' @export
#'
#' @examples
#' data <- rexp(1000)
#' x <- seq(min(data), 5, by = 0.05)
#'
#' par(mfrow = c(1, 2))
#' y <- krnl(x, 0.01, data, "epan", "reflection")
#' z <- krnl(x, 0.01, data, "epan", "transformation")
#' plot(x,y, type="l", ylim=c(0,2), col="blue")
#' lines(x,z, col="red")
#' lines(x, dexp(x))
#'
#' y <- krnl(x, 1, data, "epan", "reflection")
#' z <- krnl(x, 1, data, "epan", "transformation")
#' plot(x,y, type="l", ylim=c(0,2), col="blue")
#' lines(x,z, col="red")
#' lines(x, dexp(x))
#' par(mfrow = c(1, 1))
krnl <- function(x, h, data, kernel = "epan", method = "reflection")
{
  if(method == "transformation"){

    fhat      <- .kernEst_transformation(x, h, data, kernel)

  } else if(method == "reflection"){

    fhat      <- .kernEst_reflection(x, h, data, kernel)

  } else {

    fhat      <- NA
  }

  return( fhat )
}

###########################################################################
###    Building blocks
###########################################################################

.reflection <- function(data) list(data = data, pseudo = -data)

###########################################################################
###
###########################################################################

.kernel_uniform <- function(x, h, data)
{
  unlist(unname( lapply(x, function(x) colSums( dplyr::select(dplyr::mutate(data.frame(pt=data), value = ifelse(x-h <= pt & pt <= x+h, 0.5, 0)), value) )) ))
}

.kernel_gaussian <- function(x, h, data)
{
  unlist(unname( lapply(x, function(x) colSums( dplyr::select(dplyr::mutate(data.frame(pt=data), value = dnorm((x-pt)/h)), value) )) ))
}

.kernel_epanechnikov <- function(x, h, data)
{
  unlist(unname( lapply(x, function(x) colSums( dplyr::select(dplyr::mutate(data.frame(pt=data), value = ifelse(x-h <= pt & pt <= x+h, 0.75 * (1- ((x-pt)/h)^2), 0)), value) )) ))
}

###########################################################################
###
###########################################################################

.kernEst_transformation <- function(x, h, data, kernel = "epan", adj_positive = abs(min(data)))
{
  adj_data  <- .transformation(data, adj_positive)
  adj_x     <- .transformation(x, adj_positive)

  if(kernel == "epan"){

    fhat      <- .kernel_epanechnikov(adj_x, h, adj_data) / (length(data) * h)

  } else if(kernel == "gaussian"){

    fhat      <- .kernel_gaussian(adj_x, h, adj_data) / (length(data) * h)

  } else if(kernel == "uniform"){

    fhat      <- .kernel_uniform(adj_x, h, adj_data) / (length(data) * h)

  } else {

    fhat      <- NA
  }

  return( fhat / (exp(adj_x)) )
}

.kernEst_reflection <- function(x, h, data, kernel = "epan")
{
  adj_data  <- unlist(.reflection(data))

  if(kernel == "epan"){

    fhat      <- .kernel_epanechnikov(x, h, adj_data) / (length(data) * h)

  } else if(kernel == "gaussian"){

    fhat      <- .kernel_gaussian(x, h, adj_data) / (length(data) * h)

  } else if(kernel == "uniform"){

    fhat      <- .kernel_uniform(x, h, adj_data) / (length(data) * h)

  } else {

    fhat <- NA
  }

  return( fhat )
}
