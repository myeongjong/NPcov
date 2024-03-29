% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idea.R
\name{IDEA}
\alias{IDEA}
\title{Iterated Density Estimation Evolutionary Algorithm (IDEA) for regression}
\usage{
IDEA(
  x,
  y,
  h,
  n,
  m,
  tau,
  kernel = "epan",
  shape = "general",
  method = "reflection",
  dist_init = "exponential",
  min_iter = 20,
  max_iter = 400,
  tol = 0.001,
  check = 5
)
}
\arguments{
\item{x}{x}

\item{y}{y}

\item{h}{Bandwidth}

\item{n}{n}

\item{m}{m}

\item{tau}{tau}

\item{kernel}{Type of kernel: epan, gaussian, or uniform}

\item{shape}{Shape of estimator: general or monotone}

\item{method}{kernel density estimation/sampling method: reflection or transformation}

\item{dist_init}{Distribution for initialization}

\item{min_iter}{Minimum number of iteration}

\item{max_iter}{Maximum number of iteration}

\item{tol}{Error tolerance}

\item{check}{Iteration tolerance}
}
\value{
List
}
\description{
Iterated Density Estimation Evolutionary Algorithm (IDEA) for regression
}
\examples{
m           <- 20
n           <- 100
tau         <- 0.2
h           <- 0.1
kernel      <- "gaussian"
shape       <- "general"
method      <- "reflection"
dist_init   <- "exponential"
min_iter    <- 20
max_iter    <- 100
wave_model  <- function(h, c_0, c_1, a) ifelse(h == 0, 0, c_0 + c_1 * (1-a*sin(h/a)/h))

x           <- seq(0.1, 10, by = 0.1)
y           <- wave_model(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

r           <- seq(0.05, 10, by = 0.05)
ptm         <- proc.time()
storage     <- IDEA(x, y_err, h, n, m, tau, kernel, shape, method, dist_init, min_iter, max_iter)
ptm         <- proc.time() - ptm
data        <- na.omit( unlist(storage[[2]]) )

z           <- rep(NA, length(r))
for(i in 1:length(r)) z[i] <- ksreflect(r = r[i], data = data, h = h, kernel = kernel, shape = shape)

plot(x, y_err, col = "gray")
lines(r, z, col = "red", lwd = 2)
lines(x, y, col = "black", lwd = 2)

plot(na.omit(storage$history_obj_min), type="l", ylim = c(4, 15), lty = "dotted", xlab = "Iteration", ylab = "Objective value")
lines(na.omit(storage$history_obj_max_sel), col = "red")
lines(na.omit(storage$history_obj_mean), col = "blue")
lines(na.omit(storage$history_obj_max), lty = "dotted")

ptm
}
