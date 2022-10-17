#ifndef RCPPLIST_H
#define RCPPLIST_H

#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <Rcpp.h>
#include <RcppNumerical.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

double struve_intern(double x, double nu, double factor_Sign, bool expscaled);

double struveH(const double x, const double nu);

double struveL(const double x, const double nu, const bool expScaled);

double struveH_modified_C(const double x, const double nu);

double besselJ(double x, double nu);

double besselJ_modified_C(const double x, const double nu);

Rcpp::List integrate_list(const double& r, const double& h, const double& t);

double besselJ_generalized_C(const double& r, const double& h, const double& t);

double lmda0_C(const double x);

double lmda1_C(const double x);

#endif
