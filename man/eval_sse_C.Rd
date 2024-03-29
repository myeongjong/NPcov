% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{eval_sse_C}
\alias{eval_sse_C}
\title{Computing SSE of an estimator obtained by the IDEA algorithm (C++ version)}
\usage{
eval_sse_C(storage, x, y, h, kernel, shape)
}
\arguments{
\item{storage}{Data storage}

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
Computing SSE of an estimator obtained by the IDEA algorithm (C++ version)
}
