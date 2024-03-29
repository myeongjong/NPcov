% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{evalEach_sse_C}
\alias{evalEach_sse_C}
\title{Computing SSE of an estimator obtained by the IDEA algorithm for each pseudo data (C++ version)}
\usage{
evalEach_sse_C(data, x, y, h, kernel, shape)
}
\arguments{
\item{data}{Pseudo data}

\item{x}{Input (e.g., distance)}

\item{y}{Output (e.g., correlation value)}

\item{h}{Bandwidth}

\item{kernel}{Type of kernel: epan, gaussian, or uniform}

\item{shape}{Shape of estimator: general or monotone}
}
\value{
SSE
}
\description{
Computing SSE of an estimator obtained by the IDEA algorithm for each pseudo data (C++ version)
}
